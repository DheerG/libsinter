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
#include "../core/quantize/qp_table.h"
#include "../perceptual/csf_weights.h"
#include "../perceptual/activity_mask.h"
#include "../perceptual/ptd.h"

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
    uint16_t step = snt_qp_to_step(ctx->header.base_qp);
    int blk_cols = pad_w / BLK;
    int blk_rows = pad_h / BLK;
    int total_blocks = blk_cols * blk_rows;

    /* Allocate storage for all block codes */
    int16_t *dc_levels = (int16_t *)calloc(total_blocks, sizeof(int16_t));
    SntPVQCode *pvq_codes = (SntPVQCode *)calloc(total_blocks, sizeof(SntPVQCode));
    int8_t *qp_deltas = (int8_t *)calloc(total_blocks, sizeof(int8_t));
    /* Coefficient buffer for intra prediction (frequency domain) */
    int16_t *coeff_buf = (int16_t *)calloc((size_t)total_blocks * BLK * BLK, sizeof(int16_t));
    if (!dc_levels || !pvq_codes || !coeff_buf || !qp_deltas) {
        free(work); free(recon); free(dc_levels); free(pvq_codes);
        free(coeff_buf); free(qp_deltas);
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

            /* Activity-based QP offset — signaled in bitstream so decoder matches */
            int8_t qp_delta = psy_activity_qp_offset(activities[bidx],
                                                       mean_activity,
                                                       ctx->psy_cfg.alpha);
            qp_deltas[bidx] = qp_delta;
            int eff_qp = ctx->header.base_qp + qp_delta;
            if (eff_qp < 0) eff_qp = 0;
            if (eff_qp > 63) eff_qp = 63;
            uint16_t blk_step = snt_qp_to_step((uint8_t)eff_qp);

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

    /* Step 4: Entropy-code the quantized data using rANS.
     *
     * Encoding order (forward, for CDF synchronization):
     *   For each block: qp_delta, dc_level, skip, [gain, k, shape_idx if !skip]
     *
     * rANS encodes in REVERSE but the CDF state must be computed in
     * FORWARD order. We do two passes:
     *   Pass 1 (forward): record CDF states and symbol ranges
     *   Pass 2 (reverse): rANS encode using saved ranges */
    snt_bw_byte_align(&bw);

    /* Count total symbols for rANS */
    int nsyms = 0;
    for (int i = 0; i < total_blocks; i++) {
        nsyms += 3; /* qp_delta_sign, qp_delta_abs, dc_sign */
        nsyms += 1; /* dc_abs */
        nsyms += 1; /* skip flag */
        if (!pvq_codes[i].skip) {
            nsyms += 3; /* gain, k, shape_idx placeholder */
        }
    }

    /* For this reference implementation, use a hybrid approach:
     * rANS for frequently-occurring symbols (skip flag, small values)
     * and raw bits for large values (shape_idx).
     * This demonstrates the rANS integration while keeping the
     * implementation tractable. */

    /* Initialize CDFs */
    RansCDF cdf_skip;    /* binary: skip or not */
    RansCDF cdf_gain;    /* PVQ gain low byte (0-255) */
    rans_cdf_init_uniform(&cdf_skip, 2, RANS_ADAPT_RATE_MODE);
    rans_cdf_init_uniform(&cdf_gain, 256, RANS_ADAPT_RATE_COEFF);

    /* Pass 1 (forward): record ranges for rANS reverse encoding */
    typedef struct { uint16_t cdf_lo; uint16_t freq; } SymRange;
    int max_ranges = total_blocks * 4;
    SymRange *ranges = (SymRange *)calloc(max_ranges, sizeof(SymRange));
    int range_count = 0;

    /* Temporary CDFs for forward pass */
    RansCDF fwd_skip, fwd_gain;
    rans_cdf_init_uniform(&fwd_skip, 2, RANS_ADAPT_RATE_MODE);
    rans_cdf_init_uniform(&fwd_gain, 256, RANS_ADAPT_RATE_COEFF);

    for (int i = 0; i < total_blocks; i++) {
        /* Skip flag — binary, benefits most from rANS (highly skewed at high QP) */
        uint16_t skip_sym = pvq_codes[i].skip ? 1 : 0;
        rans_cdf_symbol_range(&fwd_skip, skip_sym, &ranges[range_count].cdf_lo,
                              &ranges[range_count].freq);
        rans_cdf_adapt(&fwd_skip, skip_sym);
        range_count++;

        /* PVQ gain: code low byte through rANS (high byte raw if needed) */
        if (!pvq_codes[i].skip) {
            uint16_t gain_lo = (uint16_t)(pvq_codes[i].gain & 0xFF);
            rans_cdf_symbol_range(&fwd_gain, gain_lo, &ranges[range_count].cdf_lo,
                                  &ranges[range_count].freq);
            rans_cdf_adapt(&fwd_gain, gain_lo);
            range_count++;
        }
    }

    /* Pass 2 (reverse): rANS encode */
    uint16_t *rev_buf = (uint16_t *)calloc(65536, sizeof(uint16_t));
    RansEncoder rans_enc;
    rans_enc_init(&rans_enc, rev_buf, 65536, ctx->header.dual_rans);

    for (int r = range_count - 1; r >= 0; r--) {
        rans_enc_put(&rans_enc, ranges[r].cdf_lo, ranges[r].freq, (uint32_t)r);
    }
    rans_enc_flush(&rans_enc);

    /* Write rANS-coded data to bitstream */
    uint32_t rans_len;
    const uint16_t *rans_data = rans_enc_data(&rans_enc, &rans_len);

    /* Write rANS payload length first (so decoder knows how much to read) */
    snt_bw_put_bits(&bw, rans_len, 16);

    /* Write rANS payload as bytes */
    for (uint32_t i = 0; i < rans_len; i++) {
        snt_bw_put_bits(&bw, rans_data[i] & 0xFF, 8);
        snt_bw_put_bits(&bw, rans_data[i] >> 8, 8);
    }

    /* Write non-rANS data: QP deltas, DC levels, gain high byte, K, shape */
    for (int i = 0; i < total_blocks; i++) {
        snt_bw_put_bits(&bw, (uint32_t)(qp_deltas[i] + 8), 4);
    }
    /* DC levels as raw 16-bit (full range, not truncated) */
    for (int i = 0; i < total_blocks; i++) {
        snt_bw_put_bits(&bw, (uint32_t)(uint16_t)dc_levels[i], 16);
    }
    /* PVQ supplementary data */
    for (int i = 0; i < total_blocks; i++) {
        if (!pvq_codes[i].skip) {
            snt_bw_put_bits(&bw, (pvq_codes[i].gain >> 8) & 0xFF, 8); /* gain high byte */
            snt_bw_put_bits(&bw, pvq_codes[i].k, 16);
            snt_bw_put_bits(&bw, (uint32_t)(pvq_codes[i].shape_idx >> 32), 32);
            snt_bw_put_bits(&bw, (uint32_t)(pvq_codes[i].shape_idx & 0xFFFFFFFF), 32);
        }
    }

    snt_bw_byte_align(&bw);
    *bs_len = snt_bw_bytes_written(&bw);

    free(rev_buf);
    free(ranges);

    free(work);
    free(recon);
    free(dc_levels);
    free(pvq_codes);
    free(qp_deltas);
    free(coeff_buf);
    free(activities);

    return 0;
}
