/*
 * Sinter Codec — Simple Encoder (quality-validated pipeline)
 * Fixed 8x8 blocks, raw bit packing, no quad-tree overhead.
 * This is the quality baseline — matches the stage-by-stage diagnostic.
 */

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "encoder.h"
#include "../core/transform/lapping.h"
#include "../core/transform/dct.h"
#include "../core/quantize/pvq.h"
#include "../core/quantize/qp_table.h"
#include "../core/bitstream/bitstream.h"
#include "../util/fixed_point.h"
#include "../perceptual/activity_mask.h"

#define BLK 8

int snte_encode_frame(SntEncoderCtx *ctx, const uint8_t *input_y,
                       const uint8_t *input_u, const uint8_t *input_v,
                       int width, int height, int stride_y, int stride_uv,
                       uint8_t *bitstream, uint32_t bs_capacity,
                       uint32_t *bs_len) {
    SntBitWriter bw;
    snt_bw_init(&bw, bitstream, bs_capacity);

    ctx->header.width = (uint32_t)width;
    ctx->header.height = (uint32_t)height;
    ctx->header.frame_type = SNT_FRAME_KEY;
    snt_write_frame_header(&bw, &ctx->header);

    int pad_w = ((width + BLK - 1) / BLK) * BLK;
    int pad_h = ((height + BLK - 1) / BLK) * BLK;
    int blk_cols = pad_w / BLK;
    int blk_rows = pad_h / BLK;
    int total_blocks = blk_cols * blk_rows;

    int16_t *work = (int16_t *)calloc((size_t)pad_w * pad_h, sizeof(int16_t));
    int16_t *recon = (int16_t *)calloc((size_t)pad_w * pad_h, sizeof(int16_t));
    if (!work || !recon) { free(work); free(recon); return -1; }

    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
            work[y * pad_w + x] = (int16_t)input_y[y * stride_y + x];

    /* Pre-filter */
    snt_lap_pre_filter_2d_ref(work, pad_w, pad_h, pad_w, BLK);

    uint16_t step = snt_qp_to_step(ctx->header.base_qp);

    /* Write block count */
    snt_bw_byte_align(&bw);
    snt_bw_put_bits(&bw, (uint32_t)total_blocks, 16);

    /* Encode each block: DCT → PVQ → raw bits to bitstream */
    for (int by = 0; by < blk_rows; by++) {
        for (int bx = 0; bx < blk_cols; bx++) {
            int16_t coeffs[64];
            snt_dct_fwd_2d_ref(&work[by * BLK * pad_w + bx * BLK], pad_w, BLK, coeffs);

            /* Quantize DC */
            int16_t dc_level = pvq_quantize_dc(coeffs[0], step, 1);
            snt_bw_put_bits(&bw, (uint32_t)(uint16_t)dc_level, 16);

            /* Quantize AC via band-partitioned PVQ.
             * Split 63 AC coefficients into 4 bands of ~15-16 dims each.
             * This avoids combinatorial overflow and gives much better
             * shape resolution (K=40-60 per band vs K=14 monolithic). */
            #define NUM_BANDS 4
            static const int band_bounds[NUM_BANDS + 1] = {0, 15, 31, 47, 63};
            int16_t recon_ac[63];
            memset(recon_ac, 0, sizeof(recon_ac));

            for (int b = 0; b < NUM_BANDS; b++) {
                int start = band_bounds[b];
                int n = band_bounds[b + 1] - start;
                SntPVQCode code;
                int16_t band_recon[16];
                pvq_encode(&coeffs[1 + start], n, step, 1, 128, &code, band_recon);

                snt_bw_put_bits(&bw, code.skip, 1);
                if (!code.skip) {
                    snt_bw_put_bits(&bw, code.gain & 0xFFFF, 16);
                    snt_bw_put_bits(&bw, code.k, 16);
                    int shape_bits = pvq_shape_bits(n, code.k);
                    if (shape_bits < 1) shape_bits = 1;
                    if (shape_bits > 64) shape_bits = 64;
                    snt_bw_put_bits(&bw, (uint32_t)shape_bits, 8);
                    if (shape_bits > 32) {
                        snt_bw_put_bits(&bw, (uint32_t)(code.shape_idx >> 32), shape_bits - 32);
                        snt_bw_put_bits(&bw, (uint32_t)(code.shape_idx & 0xFFFFFFFF), 32);
                    } else {
                        snt_bw_put_bits(&bw, (uint32_t)code.shape_idx, shape_bits);
                    }
                }
                for (int i = 0; i < n; i++)
                    recon_ac[start + i] = band_recon[i];
            }

            /* Reconstruct for reference */
            int16_t recon_coeffs[64];
            recon_coeffs[0] = pvq_dequantize_dc(dc_level, step);
            for (int i = 0; i < 63; i++)
                recon_coeffs[i + 1] = recon_ac[i];

            snt_dct_inv_2d_ref(recon_coeffs, BLK,
                               &recon[by * BLK * pad_w + bx * BLK], pad_w);
        }
    }

    /* Post-filter */
    snt_lap_post_filter_2d_ref(recon, pad_w, pad_h, pad_w, BLK);

    /* Chroma (simple DC quantization) */
    snt_bw_byte_align(&bw);
    if (input_u && input_v) {
        snt_bw_put_bits(&bw, 1, 1);
        int cw = width / 2, ch = height / 2;
        int cqp = ctx->header.base_qp + ctx->header.chroma_qp_offset;
        if (cqp < 0) cqp = 0; if (cqp > 63) cqp = 63;
        uint16_t cstep = snt_qp_to_step((uint8_t)cqp);
        for (int y = 0; y < ch; y += 4)
            for (int x = 0; x < cw; x += 4) {
                int32_t sum = 0; int cnt = 0;
                for (int dy = 0; dy < 4 && y+dy < ch; dy++)
                    for (int dx = 0; dx < 4 && x+dx < cw; dx++)
                        { sum += input_u[(y+dy)*stride_uv+x+dx]; cnt++; }
                snt_bw_put_bits(&bw, (uint32_t)(uint16_t)pvq_quantize_dc((int16_t)(cnt?sum/cnt-128:0), cstep, 1), 16);
            }
        for (int y = 0; y < ch; y += 4)
            for (int x = 0; x < cw; x += 4) {
                int32_t sum = 0; int cnt = 0;
                for (int dy = 0; dy < 4 && y+dy < ch; dy++)
                    for (int dx = 0; dx < 4 && x+dx < cw; dx++)
                        { sum += input_v[(y+dy)*stride_uv+x+dx]; cnt++; }
                snt_bw_put_bits(&bw, (uint32_t)(uint16_t)pvq_quantize_dc((int16_t)(cnt?sum/cnt-128:0), cstep, 1), 16);
            }
    } else {
        snt_bw_put_bits(&bw, 0, 1);
    }

    snt_bw_byte_align(&bw);
    *bs_len = snt_bw_bytes_written(&bw);

    free(work);
    free(recon);
    return 0;
}
