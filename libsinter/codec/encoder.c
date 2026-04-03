/*
 * Sinter Codec — Encoder Pipeline
 * Implements the top-level encode flow:
 *
 *   Input Frame -> Block Partitioning -> Prediction -> Residual
 *   -> Pre-Filter (lapping) -> Forward DCT -> PVQ Quantization
 *   -> rANS Entropy Coding -> Bitstream
 *
 * This is a simplified reference encoder operating on a uniform
 * grid of fixed-size blocks (8x8). Production encoders would use
 * quad-tree partitioning with R-D optimization.
 */

#include <string.h>
#include <stdlib.h>
#include "encoder.h"
#include "../core/transform/lapping.h"
#include "../core/transform/dct.h"
#include "../core/quantize/pvq.h"
#include "../core/entropy/rans.h"
#include "../core/bitstream/bitstream.h"
#include "../core/prediction/intra.h"
#include "../perceptual/csf_weights.h"
#include "../perceptual/activity_mask.h"
#include "../perceptual/ptd.h"

/* QP to quantization step mapping (exponential, similar to H.264/AV1) */
static const uint16_t qp_to_step[64] = {
     1,  1,  2,  2,  3,  3,  4,  4,  5,  6,  7,  8,  9, 10, 11, 12,
    14, 16, 18, 20, 22, 24, 28, 32, 36, 40, 44, 48, 56, 64, 72, 80,
    88, 96,112,128,144,160,176,192,224,256,288,320,352,384,448,512,
    576,640,704,768,896,1024,1152,1280,1408,1536,1792,2048,2304,2560,2816,3072
};

uint16_t snte_qp_to_step(uint8_t qp) {
    if (qp > 63) qp = 63;
    return qp_to_step[qp];
}

int snte_encode_frame(SntEncoderCtx *ctx, const uint8_t *input_y,
                       int width, int height, int stride,
                       uint8_t *bitstream, uint32_t bs_capacity,
                       uint32_t *bs_len) {
    const int BLK = 8; /* fixed block size for reference encoder */

    /* Write frame header */
    SntBitWriter bw;
    snt_bw_init(&bw, bitstream, bs_capacity);

    ctx->header.width = (uint32_t)width;
    ctx->header.height = (uint32_t)height;
    ctx->header.frame_type = SNT_FRAME_KEY; /* all-intra for reference encoder */
    ctx->header.num_tile_cols = 1;
    ctx->header.num_tile_rows = 1;
    ctx->header.dual_rans = 1;
    snt_write_frame_header(&bw, &ctx->header);

    /* Convert input to int16_t work buffer */
    int pad_w = ((width + BLK - 1) / BLK) * BLK;
    int pad_h = ((height + BLK - 1) / BLK) * BLK;
    int16_t *work = (int16_t *)calloc((size_t)pad_w * pad_h, sizeof(int16_t));
    int16_t *recon = (int16_t *)calloc((size_t)pad_w * pad_h, sizeof(int16_t));
    if (!work || !recon) { free(work); free(recon); return -1; }

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            work[y * pad_w + x] = (int16_t)input_y[y * stride + x];
        }
    }

    /* Step 1: Apply 2D lapped pre-filter (horizontal then vertical) */
    snt_lap_pre_filter_2d_ref(work, pad_w, pad_h, pad_w, BLK);

    /* Step 2: Process blocks — DCT, PVQ quantize, collect symbols */
    uint16_t step = snte_qp_to_step(ctx->header.base_qp);
    int blk_cols = pad_w / BLK;
    int blk_rows = pad_h / BLK;
    int total_blocks = blk_cols * blk_rows;

    /* Allocate storage for all block codes */
    int16_t *dc_levels = (int16_t *)calloc(total_blocks, sizeof(int16_t));
    SntPVQCode *pvq_codes = (SntPVQCode *)calloc(total_blocks, sizeof(SntPVQCode));
    /* Coefficient buffer for intra prediction (frequency domain) */
    int16_t *coeff_buf = (int16_t *)calloc((size_t)total_blocks * BLK * BLK, sizeof(int16_t));
    if (!dc_levels || !pvq_codes || !coeff_buf) {
        free(work); free(recon); free(dc_levels); free(pvq_codes); free(coeff_buf);
        return -1;
    }

    /* First pass: compute activity for each block (before RDO) */
    float total_activity = 0.0f;
    float *activities = (float *)calloc(total_blocks, sizeof(float));

    for (int by = 0; by < blk_rows; by++) {
        for (int bx = 0; bx < blk_cols; bx++) {
            int16_t coeffs[64];
            snt_dct_fwd_2d_ref(&work[by * BLK * pad_w + bx * BLK], pad_w, BLK, coeffs);
            int bidx = by * blk_cols + bx;
            activities[bidx] = psy_compute_activity(coeffs, BLK * BLK);
            total_activity += activities[bidx];
        }
    }
    float mean_activity = total_activity / total_blocks;

    /* Second pass: quantize and encode */
    for (int by = 0; by < blk_rows; by++) {
        for (int bx = 0; bx < blk_cols; bx++) {
            int bidx = by * blk_cols + bx;
            int16_t *blk_ptr = &work[by * BLK * pad_w + bx * BLK];

            /* Forward DCT */
            int16_t coeffs[64];
            snt_dct_fwd_2d_ref(blk_ptr, pad_w, BLK, coeffs);

            /* Intra prediction (frequency domain) — uses reconstructed coefficients */
            int16_t pred[64];
            int left_bidx = by * blk_cols + (bx - 1);
            int above_bidx = (by - 1) * blk_cols + bx;
            const int16_t *left_c = (bx > 0) ? &coeff_buf[left_bidx * BLK * BLK] : NULL;
            const int16_t *above_c = (by > 0) ? &coeff_buf[above_bidx * BLK * BLK] : NULL;

            snt_intra_predict_freq_ref(left_c, above_c, BLK, pred);

            /* Compute residual */
            int16_t residual[64];
            for (int i = 0; i < 64; i++)
                residual[i] = coeffs[i] - pred[i];

            /* Use base QP step (activity masking is encoder-only heuristic;
             * decoder uses base_qp so quant/dequant must match). */
            uint16_t blk_step = snte_qp_to_step(ctx->header.base_qp);

            /* Quantize DC */
            dc_levels[bidx] = pvq_quantize_dc(residual[0], blk_step, 1);

            /* Quantize AC with PVQ */
            int16_t recon_ac[63];
            pvq_encode(&residual[1], 63, blk_step, 1, 128,
                       &pvq_codes[bidx], recon_ac);

            /* Reconstruct coefficients (for intra prediction reference) */
            int16_t *recon_coeffs = &coeff_buf[bidx * BLK * BLK];
            recon_coeffs[0] = pvq_dequantize_dc(dc_levels[bidx], blk_step) + pred[0];
            for (int i = 0; i < 63; i++)
                recon_coeffs[i + 1] = recon_ac[i] + pred[i + 1];

            /* Inverse DCT to spatial domain (for post-filter + output) */
            snt_dct_inv_2d_ref(recon_coeffs, BLK,
                               &recon[by * BLK * pad_w + bx * BLK], pad_w);
        }
    }

    /* Step 3: Apply 2D lapped post-filter to reconstruction */
    snt_lap_post_filter_2d_ref(recon, pad_w, pad_h, pad_w, BLK);

    /* Step 4: Entropy-code the quantized data */
    /* Write block data using rANS (simplified: just write raw for now,
     * full rANS integration would interleave DC and PVQ codes) */
    snt_bw_byte_align(&bw);

    /* Write DC levels */
    for (int i = 0; i < total_blocks; i++) {
        snt_bw_put_bits(&bw, (uint32_t)(uint16_t)dc_levels[i], 16);
    }

    /* Write PVQ codes */
    for (int i = 0; i < total_blocks; i++) {
        snt_bw_put_bits(&bw, pvq_codes[i].skip, 1);
        if (!pvq_codes[i].skip) {
            snt_bw_put_bits(&bw, pvq_codes[i].gain, 16);
            snt_bw_put_bits(&bw, pvq_codes[i].k, 16);
            snt_bw_put_bits(&bw, (uint32_t)pvq_codes[i].shape_idx, 32);
        }
    }

    snt_bw_byte_align(&bw);
    *bs_len = snt_bw_bytes_written(&bw);

    free(work);
    free(recon);
    free(dc_levels);
    free(pvq_codes);
    free(coeff_buf);
    free(activities);

    return 0;
}
