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
#include "../perceptual/csf_weights.h"

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

    /* DC + low-AC prediction buffers */
    int16_t *dc_recon_above = (int16_t *)calloc(blk_cols, sizeof(int16_t));
    int16_t dc_recon_left = 0;
    /* Predict lowest 4 AC coefficients from neighbors:
     * AC[1] (horiz freq-1) from left, AC[8] (vert freq-1) from above,
     * AC[2] (horiz freq-2) from left, AC[9] (diag) from above-left avg */
    #define NUM_PRED_AC 4
    static const int pred_ac_idx[NUM_PRED_AC] = {1, 8, 2, 9};
    int16_t *ac_recon_above = (int16_t *)calloc((size_t)blk_cols * NUM_PRED_AC, sizeof(int16_t));
    int16_t ac_recon_left[NUM_PRED_AC] = {0};

    /* Encode each block: DCT → DC predict → PVQ → raw bits to bitstream */
    for (int by = 0; by < blk_rows; by++) {
        dc_recon_left = 0;
        for (int bx = 0; bx < blk_cols; bx++) {
            int16_t coeffs[64];
            snt_dct_fwd_2d_ref(&work[by * BLK * pad_w + bx * BLK], pad_w, BLK, coeffs);

            /* DC prediction from reconstructed neighbors (coefficient domain).
             * Average of left and above if both available. */
            int16_t dc_pred = 0;
            if (bx > 0 && by > 0)
                dc_pred = (dc_recon_left + dc_recon_above[bx] + 1) >> 1;
            else if (bx > 0)
                dc_pred = dc_recon_left;
            else if (by > 0)
                dc_pred = dc_recon_above[bx];

            /* Quantize DC residual (not absolute DC) */
            int16_t dc_residual = coeffs[0] - dc_pred;
            int16_t dc_level = pvq_quantize_dc(dc_residual, step, 1);
            snt_bw_put_sexpgolomb(&bw, dc_level); /* VLC: ~2-4 bits for small residuals */

            /* Reconstruct DC for prediction of future blocks */
            int16_t dc_recon = pvq_dequantize_dc(dc_level, step) + dc_pred;
            dc_recon_left = dc_recon;
            dc_recon_above[bx] = dc_recon;

            /* Low-AC prediction: predict lowest AC coefficients from neighbors.
             * AC[1] predicted from left (horizontal continuity),
             * AC[8] predicted from above (vertical continuity),
             * AC[2] from left, AC[9] from above. */
            int16_t ac_pred[64]; /* prediction for each coefficient */
            memset(ac_pred, 0, sizeof(ac_pred));
            for (int p = 0; p < NUM_PRED_AC; p++) {
                int idx = pred_ac_idx[p];
                int16_t pred_val = 0;
                /* AC[1], AC[2] predict from left; AC[8], AC[9] predict from above */
                if (idx < 8 && bx > 0) {
                    pred_val = ac_recon_left[p]; /* horizontal prediction */
                } else if (idx >= 8 && by > 0) {
                    pred_val = ac_recon_above[bx * NUM_PRED_AC + p]; /* vertical prediction */
                }
                ac_pred[idx] = pred_val;
                coeffs[idx] -= pred_val; /* subtract prediction from coefficient */
            }

            #define NUM_BANDS 16
            static const int band_bounds[NUM_BANDS + 1] = {0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,59,63};
            int16_t recon_ac[63];
            memset(recon_ac, 0, sizeof(recon_ac));

            for (int b = 0; b < NUM_BANDS; b++) {
                int start = band_bounds[b];
                int n = band_bounds[b + 1] - start;

                /* Compute per-band CSF weight (Jens-Rainer's fix):
                 * High-frequency bands get lower weight -> fewer K pulses -> cleaner skip.
                 * Low-frequency bands get higher weight -> more K pulses -> better fidelity. */
                float csf_sum = 0.0f;
                for (int i = 0; i < n; i++)
                    csf_sum += psy_csf_weight(BLK, 1 + start + i);
                uint8_t band_csf = (uint8_t)(csf_sum / n * 255.0f);
                if (band_csf < 1) band_csf = 1;

                /* Hybrid PVQ/Scalar mode decision per band.
                 * If energy is concentrated in 1-2 coefficients (sparse),
                 * scalar quantization is better (no shape noise).
                 * If energy is spread (texture), PVQ is better (preserves energy). */
                int16_t *band_coeffs = &coeffs[1 + start];
                int32_t max_sq = 0, sum_sq = 0;
                for (int i = 0; i < n; i++) {
                    int32_t sq = (int32_t)band_coeffs[i] * band_coeffs[i];
                    sum_sq += sq;
                    if (sq > max_sq) max_sq = sq;
                }
                /* Energy concentration: if top coeff holds >60% of energy, use scalar */
                int use_scalar = (sum_sq > 0 && max_sq * 10 > sum_sq * 6);

                int16_t band_recon[16];

                if (use_scalar) {
                    /* Scalar mode with CSF-weighted step: wider step for HF bands
                     * zeros out perceptually invisible coefficients, reducing noise.
                     * This is Sinter's equivalent of H.264's quantization matrix. */
                    uint16_t eff_step = step;
                    if (band_csf < 128) {
                        eff_step = (uint16_t)((uint32_t)step * 256 / band_csf);
                        if (eff_step > 3072) eff_step = 3072;
                    }
                    snt_bw_put_bits(&bw, 1, 1); /* mode=1: scalar */
                    for (int i = 0; i < n; i++) {
                        int16_t level = pvq_quantize_dc(band_coeffs[i], eff_step, 1);
                        snt_bw_put_sexpgolomb(&bw, level); /* VLC: ~1-5 bits */
                        band_recon[i] = pvq_dequantize_dc(level, eff_step);
                    }
                } else {
                    /* PVQ mode */
                    snt_bw_put_bits(&bw, 0, 1); /* mode=0: PVQ */
                    SntPVQCode code;
                    pvq_encode(band_coeffs, n, step, 1, band_csf, &code, band_recon);

                    snt_bw_put_bits(&bw, code.skip, 1);
                    if (!code.skip) {
                        snt_bw_put_expgolomb(&bw, code.gain);
                        snt_bw_put_expgolomb(&bw, code.k);
                        int shape_bits = pvq_shape_bits(n, code.k);
                        if (shape_bits < 1) shape_bits = 1;
                        if (shape_bits > 64) shape_bits = 64;
                        if (shape_bits > 32) {
                            snt_bw_put_bits(&bw, (uint32_t)(code.shape_idx >> 32), shape_bits - 32);
                            snt_bw_put_bits(&bw, (uint32_t)(code.shape_idx & 0xFFFFFFFF), 32);
                        } else {
                            snt_bw_put_bits(&bw, (uint32_t)code.shape_idx, shape_bits);
                        }
                    }
                }
                for (int i = 0; i < n; i++)
                    recon_ac[start + i] = band_recon[i];
            }

            /* Reconstruct for reference */
            int16_t recon_coeffs[64];
            recon_coeffs[0] = dc_recon;
            for (int i = 0; i < 63; i++)
                recon_coeffs[i + 1] = recon_ac[i];

            /* Add back AC prediction to get full reconstructed coefficients */
            for (int p = 0; p < NUM_PRED_AC; p++) {
                int idx = pred_ac_idx[p];
                recon_coeffs[idx] += ac_pred[idx];
            }

            /* Store reconstructed AC for future prediction */
            for (int p = 0; p < NUM_PRED_AC; p++) {
                int idx = pred_ac_idx[p];
                ac_recon_left[p] = recon_coeffs[idx];
                ac_recon_above[bx * NUM_PRED_AC + p] = recon_coeffs[idx];
            }

            snt_dct_inv_2d_ref(recon_coeffs, BLK,
                               &recon[by * BLK * pad_w + bx * BLK], pad_w);
        }
        /* Reset left AC prediction at row start */
        memset(ac_recon_left, 0, sizeof(ac_recon_left));
    }

    /* Post-filter */
    snt_lap_post_filter_2d_ref(recon, pad_w, pad_h, pad_w, BLK);

    /* Chroma: 4x4 DCT + scalar quantization (team consensus: Option B).
     * No PVQ, no lapping — chroma doesn't need texture preservation
     * or deblocking at 4:2:0 half-resolution. */
    snt_bw_byte_align(&bw);
    if (input_u && input_v) {
        snt_bw_put_bits(&bw, 1, 1);
        int cw = width / 2, ch = height / 2;
        int cblk = 4; /* 4x4 chroma blocks */
        int cpad_w = ((cw + cblk - 1) / cblk) * cblk;
        int cpad_h = ((ch + cblk - 1) / cblk) * cblk;
        int cblk_cols = cpad_w / cblk;
        int cqp = ctx->header.base_qp + ctx->header.chroma_qp_offset;
        if (cqp < 0) cqp = 0; if (cqp > 63) cqp = 63;
        uint16_t cstep = snt_qp_to_step((uint8_t)cqp);

        /* Encode each chroma plane: U then V */
        const uint8_t *chroma_planes[2] = {input_u, input_v};
        for (int plane = 0; plane < 2; plane++) {
            const uint8_t *cp = chroma_planes[plane];
            int16_t c_dc_left = 0;
            int16_t *c_dc_above = (int16_t *)calloc(cblk_cols, sizeof(int16_t));

            for (int by = 0; by < cpad_h; by += cblk) {
                c_dc_left = 0;
                for (int bx = 0; bx < cpad_w; bx += cblk) {
                    /* Read 4x4 block, subtract 128 to center */
                    int16_t cblock[16];
                    for (int dy = 0; dy < 4; dy++)
                        for (int dx = 0; dx < 4; dx++) {
                            int py = by + dy, px = bx + dx;
                            cblock[dy*4+dx] = (py < ch && px < cw) ?
                                (int16_t)cp[py*stride_uv+px] - 128 : 0;
                        }

                    /* Forward 4x4 DCT */
                    int16_t ccoeffs[16];
                    snt_dct_fwd_2d_ref(cblock, 4, 4, ccoeffs);

                    /* DC prediction */
                    int cbx = bx / cblk;
                    int16_t c_dc_pred = 0;
                    if (bx > 0 && by > 0)
                        c_dc_pred = (c_dc_left + c_dc_above[cbx] + 1) >> 1;
                    else if (bx > 0) c_dc_pred = c_dc_left;
                    else if (by > 0) c_dc_pred = c_dc_above[cbx];

                    /* Quantize DC residual */
                    int16_t dc_res = ccoeffs[0] - c_dc_pred;
                    int16_t dc_level = pvq_quantize_dc(dc_res, cstep, 1);
                    snt_bw_put_sexpgolomb(&bw, dc_level);

                    /* Reconstruct DC for prediction */
                    int16_t dc_recon = pvq_dequantize_dc(dc_level, cstep) + c_dc_pred;
                    c_dc_left = dc_recon;
                    c_dc_above[cbx] = dc_recon;

                    /* Scalar quantize all 15 AC coefficients */
                    for (int i = 1; i < 16; i++) {
                        int16_t ac_level = pvq_quantize_dc(ccoeffs[i], cstep, 1);
                        snt_bw_put_sexpgolomb(&bw, ac_level);
                    }
                }
            }
            free(c_dc_above);
        }
    } else {
        snt_bw_put_bits(&bw, 0, 1);
    }

    snt_bw_byte_align(&bw);
    *bs_len = snt_bw_bytes_written(&bw);

    free(work);
    free(recon);
    free(dc_recon_above);
    free(ac_recon_above);
    return 0;
}
