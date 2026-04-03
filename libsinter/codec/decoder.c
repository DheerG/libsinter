/*
 * Sinter Codec — Simple Decoder (matches encoder_simple.c)
 */

#include <string.h>
#include <stdlib.h>
#include "decoder.h"
#include "../util/fixed_point.h"
#include "../core/transform/lapping.h"
#include "../core/transform/dct.h"
#include "../core/quantize/pvq.h"
#include "../core/quantize/qp_table.h"
#include "../core/bitstream/bitstream.h"
#include "../core/grain/grain_synth.h"
#include "../perceptual/csf_weights.h"

#define BLK 8

int sntd_decode_frame(SntDecoderCtx *ctx, const uint8_t *bitstream,
                       uint32_t bs_len, uint8_t *output_y,
                       uint8_t *output_u, uint8_t *output_v,
                       int *out_width, int *out_height,
                       int stride_y, int stride_uv) {
    SntBitReader br;
    snt_br_init(&br, bitstream, bs_len);
    snt_read_frame_header(&br, &ctx->header);

    int width = (int)ctx->header.width;
    int height = (int)ctx->header.height;
    *out_width = width;
    *out_height = height;

    int pad_w = ((width + BLK - 1) / BLK) * BLK;
    int pad_h = ((height + BLK - 1) / BLK) * BLK;
    int blk_cols = pad_w / BLK;
    int blk_rows = pad_h / BLK;

    uint16_t step = snt_qp_to_step(ctx->header.base_qp);

    snt_br_byte_align(&br);
    uint32_t total_blocks = snt_br_get_bits(&br, 16);
    (void)total_blocks; /* we know it's blk_cols * blk_rows */

    int16_t *recon = (int16_t *)calloc((size_t)pad_w * pad_h, sizeof(int16_t));
    int16_t *dc_recon_above = (int16_t *)calloc(blk_cols, sizeof(int16_t));
    int16_t dc_recon_left = 0;
    #define NUM_PRED_AC 4
    static const int pred_ac_idx[NUM_PRED_AC] = {1, 8, 2, 9};
    int16_t *ac_recon_above = (int16_t *)calloc((size_t)blk_cols * NUM_PRED_AC, sizeof(int16_t));
    int16_t ac_recon_left[NUM_PRED_AC] = {0};

    /* Decode each block */
    for (int by = 0; by < blk_rows; by++) {
        dc_recon_left = 0;
        for (int bx = 0; bx < blk_cols; bx++) {
            int16_t recon_coeffs[64];

            /* DC prediction (mirrors encoder exactly) */
            int16_t dc_pred = 0;
            if (bx > 0 && by > 0)
                dc_pred = (dc_recon_left + dc_recon_above[bx] + 1) >> 1;
            else if (bx > 0)
                dc_pred = dc_recon_left;
            else if (by > 0)
                dc_pred = dc_recon_above[bx];

            /* DC: read residual, add prediction */
            int16_t dc_level = (int16_t)(uint16_t)snt_br_get_bits(&br, 16);
            int16_t dc_recon = pvq_dequantize_dc(dc_level, step) + dc_pred;
            recon_coeffs[0] = dc_recon;
            dc_recon_left = dc_recon;
            dc_recon_above[bx] = dc_recon;

            /* AC via band-partitioned PVQ (4 bands of ~15 dims) */
            #define NUM_BANDS 16
            static const int band_bounds[NUM_BANDS + 1] = {0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,59,63};

            for (int b = 0; b < NUM_BANDS; b++) {
                int start = band_bounds[b];
                int n = band_bounds[b + 1] - start;

                uint8_t quant_mode = (uint8_t)snt_br_get_bits(&br, 1);

                if (quant_mode == 1) {
                    /* Scalar mode with CSF-weighted step (must match encoder) */
                    float csf_sum = 0;
                    for (int i = 0; i < n; i++)
                        csf_sum += psy_csf_weight(8, 1 + start + i);
                    uint8_t band_csf_d = (uint8_t)(csf_sum / n * 255.0f);
                    if (band_csf_d < 1) band_csf_d = 1;
                    uint16_t eff_step = step;
                    if (band_csf_d < 128) {
                        eff_step = (uint16_t)((uint32_t)step * 256 / band_csf_d);
                        if (eff_step > 3072) eff_step = 3072;
                    }
                    for (int i = 0; i < n; i++) {
                        int16_t level = (int16_t)(uint16_t)snt_br_get_bits(&br, 16);
                        recon_coeffs[1 + start + i] = pvq_dequantize_dc(level, eff_step);
                    }
                } else {
                    /* PVQ mode */
                    uint8_t skip = (uint8_t)snt_br_get_bits(&br, 1);
                    if (!skip) {
                        SntPVQCode code;
                        code.skip = 0;
                        code.gain = snt_br_get_bits(&br, 16);
                        code.k = (uint16_t)snt_br_get_bits(&br, 16);
                        int shape_bits = pvq_shape_bits(n, code.k);
                        if (shape_bits > 64) shape_bits = 64;
                        if (shape_bits > 32) {
                            uint64_t hi = snt_br_get_bits(&br, shape_bits - 32);
                            uint64_t lo = snt_br_get_bits(&br, 32);
                            code.shape_idx = (hi << 32) | lo;
                        } else {
                            code.shape_idx = snt_br_get_bits(&br, shape_bits);
                        }
                        int16_t band_recon[16];
                        pvq_decode(&code, n, step, band_recon);
                        for (int i = 0; i < n; i++)
                            recon_coeffs[1 + start + i] = band_recon[i];
                    } else {
                        for (int i = 0; i < n; i++)
                            recon_coeffs[1 + start + i] = 0;
                    }
                }
            }

            /* Add back AC prediction (mirrors encoder) */
            for (int p = 0; p < NUM_PRED_AC; p++) {
                int idx = pred_ac_idx[p];
                int16_t pred_val = 0;
                if (idx < 8 && bx > 0)
                    pred_val = ac_recon_left[p];
                else if (idx >= 8 && by > 0)
                    pred_val = ac_recon_above[bx * NUM_PRED_AC + p];
                recon_coeffs[idx] += pred_val;
                ac_recon_left[p] = recon_coeffs[idx];
                ac_recon_above[bx * NUM_PRED_AC + p] = recon_coeffs[idx];
            }

            snt_dct_inv_2d_ref(recon_coeffs, BLK,
                               &recon[by * BLK * pad_w + bx * BLK], pad_w);
        }
        memset(ac_recon_left, 0, sizeof(ac_recon_left));
    }

    /* Post-filter */
    snt_lap_post_filter_2d_ref(recon, pad_w, pad_h, pad_w, BLK);

    /* Copy to output */
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
            output_y[y * stride_y + x] = snt_clip_pixel((int32_t)recon[y * pad_w + x]);

    /* Decode chroma */
    snt_br_byte_align(&br);
    if (!snt_br_eof(&br)) {
        uint32_t chroma_present = snt_br_get_bits(&br, 1);
        if (chroma_present && output_u && output_v) {
            int cw = width / 2, ch = height / 2;
            int cqp = ctx->header.base_qp + ctx->header.chroma_qp_offset;
            if (cqp < 0) cqp = 0; if (cqp > 63) cqp = 63;
            uint16_t cstep = snt_qp_to_step((uint8_t)cqp);

            for (int y = 0; y < ch; y += 4)
                for (int x = 0; x < cw; x += 4) {
                    int16_t level = (int16_t)(uint16_t)snt_br_get_bits(&br, 16);
                    int16_t dc = pvq_dequantize_dc(level, cstep) + 128;
                    for (int dy = 0; dy < 4 && y+dy < ch; dy++)
                        for (int dx = 0; dx < 4 && x+dx < cw; dx++)
                            output_u[(y+dy)*stride_uv+x+dx] = snt_clip_pixel(dc);
                }
            for (int y = 0; y < ch; y += 4)
                for (int x = 0; x < cw; x += 4) {
                    int16_t level = (int16_t)(uint16_t)snt_br_get_bits(&br, 16);
                    int16_t dc = pvq_dequantize_dc(level, cstep) + 128;
                    for (int dy = 0; dy < 4 && y+dy < ch; dy++)
                        for (int dx = 0; dx < 4 && x+dx < cw; dx++)
                            output_v[(y+dy)*stride_uv+x+dx] = snt_clip_pixel(dc);
                }
        }
    }

    /* Grain synthesis */
    if (ctx->header.grain.enabled && output_u && output_v) {
        snt_grain_apply(&ctx->header.grain, output_y, output_u, output_v,
                        width, height, stride_y, stride_uv);
    }

    free(recon);
    free(dc_recon_above);
    free(ac_recon_above);
    return 0;
}
