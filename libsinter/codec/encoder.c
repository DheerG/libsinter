/*
 * Sinter Codec — Encoder with rANS entropy coding + inter prediction
 * Two-pass: forward record + reverse rANS encode.
 * Inter spec: team consensus (Jim/Timothy/Gary/Monty/Jens-Rainer).
 *
 * Key frames: lapped transform pipeline (pre/post filter + DCT)
 * Inter frames: no lapping, MC prediction + residual coding
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
#include "../core/entropy/rans.h"
#include "../core/prediction/motion_est.h"
#include "../util/fixed_point.h"
#include "../perceptual/activity_mask.h"
#include "../perceptual/csf_weights.h"

#define BLK 8

/* ── Zigzag signed-to-unsigned mapping ────────────────────────────── */
static inline uint16_t s2u(int16_t v) {
    uint32_t u = (v >= 0) ? (uint32_t)(2 * v) : (uint32_t)(-2 * v - 1);
    return (u > 254) ? 255 : (uint16_t)u;
}

/* ── Symbol range recording ───────────────────────────────────────── */
typedef struct { uint16_t cdf_lo; uint16_t freq; } SymRange;

/*
 * Symbol emission order (normative, must match decoder exactly):
 *
 * FOR each luma block (raster: by, bx):
 *   IF inter frame:
 *     REC(cdf_inter, inter_flag)    // 0=intra, 1=inter
 *     IF inter_flag:
 *       REC(cdf_iskip, skip_flag)   // 0=coded, 1=skip
 *       IF skip: (no more symbols, copy ref block)
 *       IF coded:
 *         REC(cdf_mvx, s2u(mvd_x))
 *         REC(cdf_mvy, s2u(mvd_y))
 *   REC(cdf_dc, s2u(dc_level))     // DC of residual (inter) or pixels (intra)
 *   // Low-AC prediction subtracted before band processing (key frames only)
 *   FOR band b=0..15:
 *     REC(cdf_mode, mode)
 *     IF scalar: FOR i=0..n-1: REC(cdf_ac, s2u(level))
 *     IF PVQ:    REC(cdf_skip, skip)
 *                IF !skip: REC(cdf_gain, min(gain,255)), REC(cdf_k, min(k,255))
 * FOR each chroma plane (U, V):
 *   FOR each 4x4 block:
 *     REC(cdf_dc, s2u(dc_level))
 *     FOR i=1..15: REC(cdf_ac, s2u(level))
 */

int snte_encode_frame(SntEncoderCtx *ctx, const uint8_t *input_y,
                       const uint8_t *input_u, const uint8_t *input_v,
                       int width, int height, int stride_y, int stride_uv,
                       uint8_t *bitstream, uint32_t bs_capacity,
                       uint32_t *bs_len) {
    SntBitWriter bw;
    snt_bw_init(&bw, bitstream, bs_capacity);

    int pad_w = ((width + BLK - 1) / BLK) * BLK;
    int pad_h = ((height + BLK - 1) / BLK) * BLK;
    int blk_cols = pad_w / BLK;
    int blk_rows = pad_h / BLK;
    int total_blocks = blk_cols * blk_rows;

    /* Determine frame type */
    int is_inter = (ctx->num_refs > 0 && ctx->ref_y != NULL);
    ctx->header.width = (uint32_t)width;
    ctx->header.height = (uint32_t)height;
    ctx->header.frame_type = is_inter ? SNT_FRAME_INTER : SNT_FRAME_KEY;
    snt_write_frame_header(&bw, &ctx->header);

    uint16_t step = snt_qp_to_step(ctx->header.base_qp);

    /* Load pixels into work buffer */
    int16_t *work = (int16_t *)calloc((size_t)pad_w * pad_h, sizeof(int16_t));
    int16_t *recon = (int16_t *)calloc((size_t)pad_w * pad_h, sizeof(int16_t));
    if (!work || !recon) { free(work); free(recon); return -1; }

    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
            work[y * pad_w + x] = (int16_t)input_y[y * stride_y + x];

    /* KEY ONLY: lapping pre-filter */
    if (!is_inter)
        snt_lap_pre_filter_2d_ref(work, pad_w, pad_h, pad_w, BLK);

    /* ── SECTION A: Fixed header ──────────────────────────────────── */
    snt_bw_byte_align(&bw);
    snt_bw_put_bits(&bw, (uint32_t)total_blocks, 16);
    snt_bw_put_bits(&bw, (input_u && input_v) ? 1u : 0u, 1);
    snt_bw_byte_align(&bw);

    /* ── Allocate rANS working buffers ────────────────────────────── */
    int cw = width / 2, ch = height / 2;
    int cblk = 4;
    int cpad_w = ((cw + cblk - 1) / cblk) * cblk;
    int cpad_h = ((ch + cblk - 1) / cblk) * cblk;
    int chroma_blocks = (input_u && input_v) ? 2 * (cpad_w / cblk) * (cpad_h / cblk) : 0;

    uint32_t max_ranges = (uint32_t)(total_blocks * 85 + chroma_blocks * 16 + 256);
    SymRange *ranges = (SymRange *)malloc(max_ranges * sizeof(SymRange));
    uint32_t rc = 0;

    uint32_t ov_cap = (uint32_t)(total_blocks * 128 + 4096);
    uint8_t *ov_buf = (uint8_t *)malloc(ov_cap);
    SntBitWriter ov;
    snt_bw_init(&ov, ov_buf, ov_cap);

    /* ── Init 10 CDFs (6 existing + 4 inter) ─────────────────────── */
    RansCDF cdf_dc, cdf_ac, cdf_mode, cdf_skip, cdf_gain, cdf_k;
    RansCDF cdf_inter, cdf_iskip, cdf_mvx, cdf_mvy;
    rans_cdf_init_uniform(&cdf_dc,   256, RANS_ADAPT_RATE_COEFF);
    rans_cdf_init_uniform(&cdf_ac,   256, RANS_ADAPT_RATE_COEFF);
    rans_cdf_init_uniform(&cdf_mode,   2, RANS_ADAPT_RATE_MODE);
    rans_cdf_init_uniform(&cdf_skip,   2, RANS_ADAPT_RATE_MODE);
    rans_cdf_init_uniform(&cdf_gain, 256, RANS_ADAPT_RATE_COEFF);
    rans_cdf_init_uniform(&cdf_k,    256, RANS_ADAPT_RATE_COEFF);
    rans_cdf_init_uniform(&cdf_inter,  2, RANS_ADAPT_RATE_MODE);
    rans_cdf_init_uniform(&cdf_iskip,  2, RANS_ADAPT_RATE_MODE);
    rans_cdf_init_uniform(&cdf_mvx, 256, RANS_ADAPT_RATE_COEFF);
    rans_cdf_init_uniform(&cdf_mvy, 256, RANS_ADAPT_RATE_COEFF);

    #define REC(cdf, sym) do { \
        rans_cdf_symbol_range(&(cdf), (sym), &ranges[rc].cdf_lo, &ranges[rc].freq); \
        rans_cdf_adapt(&(cdf), (sym)); \
        rc++; \
    } while(0)

    /* ── PHASE 1: Forward pass — record symbols + overflow ─────── */
    int16_t *dc_recon_above = (int16_t *)calloc(blk_cols, sizeof(int16_t));
    int16_t dc_recon_left = 0;
    #define NUM_PRED_AC 4
    static const int pred_ac_idx[NUM_PRED_AC] = {1, 8, 2, 9};
    int16_t *ac_recon_above = (int16_t *)calloc((size_t)blk_cols * NUM_PRED_AC, sizeof(int16_t));
    int16_t ac_recon_left[NUM_PRED_AC] = {0};

    /* MV field for spatial prediction (inter only) */
    SntMV *mv_field = (SntMV *)calloc(blk_cols * blk_rows, sizeof(SntMV));

    #define NUM_BANDS 16
    static const int band_bounds[NUM_BANDS + 1] = {0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,59,63};

    for (int by = 0; by < blk_rows; by++) {
        dc_recon_left = 0;
        for (int bx = 0; bx < blk_cols; bx++) {
            int blk_idx = by * blk_cols + bx;
            int px0 = bx * BLK, py0 = by * BLK;
            int block_is_inter = 0;
            int block_is_skip = 0;
            SntMV best_mv = {0, 0};

            /* ── Inter mode decision ──────────────────────────── */
            if (is_inter) {
                /* Compute SAD for skip (MV=0).
                 * Note: input accessed at stride width (not stride_y) to match
                 * ref buffer stride. Correct when stride_y == width (CLI). */
                uint32_t skip_sad = 0;
                for (int dy = 0; dy < BLK; dy++)
                    for (int dx = 0; dx < BLK; dx++) {
                        int py = py0 + dy, px_v = px0 + dx;
                        if (py < height && px_v < width) {
                            int diff = (int)input_y[py * width + px_v]
                                     - (int)ctx->ref_y[py * width + px_v];
                            skip_sad += (uint32_t)(diff < 0 ? -diff : diff);
                        }
                    }

                /* Motion search: cur uses stride_y, ref uses width.
                 * Both must use the same stride for snte_motion_search.
                 * We pass width — correct when stride_y == width (CLI). */
                SntMV search_mv = snte_motion_search(
                    input_y, ctx->ref_y, width,
                    px0, py0, BLK, 16, width, height);
                /* Motion search returns 1/8-pel; extract integer-pel */
                int mv_ipx = search_mv.x / 8;
                int mv_ipy = search_mv.y / 8;

                /* SAD at best MV */
                uint32_t inter_sad = 0;
                for (int dy = 0; dy < BLK; dy++)
                    for (int dx = 0; dx < BLK; dx++) {
                        int py = py0 + dy, px_v = px0 + dx;
                        int ry = py + mv_ipy, rx = px_v + mv_ipx;
                        if (py < height && px_v < width &&
                            ry >= 0 && ry < height && rx >= 0 && rx < width) {
                            int diff = (int)input_y[py * width + px_v]
                                     - (int)ctx->ref_y[ry * width + rx];
                            inter_sad += (uint32_t)(diff < 0 ? -diff : diff);
                        }
                    }

                /* Rough intra SAD (vs flat DC) */
                int32_t mean = 0;
                for (int dy = 0; dy < BLK; dy++)
                    for (int dx = 0; dx < BLK; dx++) {
                        int py = py0 + dy, px_v = px0 + dx;
                        if (py < height && px_v < width)
                            mean += input_y[py * width + px_v];
                    }
                mean /= (BLK * BLK);
                uint32_t intra_sad = 0;
                for (int dy = 0; dy < BLK; dy++)
                    for (int dx = 0; dx < BLK; dx++) {
                        int py = py0 + dy, px_v = px0 + dx;
                        if (py < height && px_v < width) {
                            int diff = (int)input_y[py * width + px_v] - (int)mean;
                            intra_sad += (uint32_t)(diff < 0 ? -diff : diff);
                        }
                    }

                /* Decision: skip < inter-coded < intra, with MV cost bias */
                /* Skip threshold: if SAD < step * BLK * 2, the quantized residual
                 * would be mostly zeros — skip avoids coding overhead.
                 * (Jens-Rainer: step*BLK was too conservative, step*BLK*2 catches
                 * more near-static blocks without visible artifacts.) */
                uint32_t skip_threshold = (uint32_t)step * BLK * 2;
                uint32_t mv_cost = (uint32_t)(abs(mv_ipx) + abs(mv_ipy)) * 4;

                if (skip_sad <= skip_threshold) {
                    block_is_inter = 1;
                    block_is_skip = 1;
                    best_mv.x = 0; best_mv.y = 0;
                } else if (inter_sad + mv_cost < intra_sad) {
                    block_is_inter = 1;
                    block_is_skip = 0;
                    best_mv.x = (int16_t)mv_ipx;
                    best_mv.y = (int16_t)mv_ipy;
                } else {
                    block_is_inter = 0;
                }

                /* REC: inter_flag */
                REC(cdf_inter, block_is_inter ? 1 : 0);

                if (block_is_inter) {
                    /* REC: skip_flag */
                    REC(cdf_iskip, block_is_skip ? 1 : 0);

                    if (block_is_skip) {
                        /* Copy reference block to recon */
                        for (int dy = 0; dy < BLK; dy++)
                            for (int dx = 0; dx < BLK; dx++) {
                                int py = py0 + dy, px_v = px0 + dx;
                                if (py < height && px_v < width)
                                    recon[py * pad_w + px_v] =
                                        (int16_t)ctx->ref_y[py * width + px_v];
                            }
                        mv_field[blk_idx] = best_mv;
                        /* Reset prediction state for skip blocks.
                         * Skip blocks don't go through DCT, so there's no
                         * reconstructed DC/AC to store. Next block predicts
                         * from 0 instead of actual recon. (Gary: document this.) */
                        dc_recon_left = 0;
                        dc_recon_above[bx] = 0;
                        for (int p = 0; p < NUM_PRED_AC; p++) {
                            ac_recon_left[p] = 0;
                            ac_recon_above[bx * NUM_PRED_AC + p] = 0;
                        }
                        continue; /* skip to next block */
                    }

                    /* MV prediction from spatial neighbors */
                    const SntMV *left_mv = (bx > 0) ? &mv_field[blk_idx - 1] : NULL;
                    const SntMV *above_mv = (by > 0) ? &mv_field[blk_idx - blk_cols] : NULL;
                    const SntMV *above_right_mv = (by > 0 && bx + 1 < blk_cols) ?
                        &mv_field[blk_idx - blk_cols + 1] : NULL;
                    SntMV mv_pred = snt_mv_predict(left_mv, above_mv, above_right_mv);

                    /* MV delta (integer-pel, stored directly, not 1/8-pel) */
                    int16_t mvd_x = best_mv.x - mv_pred.x;
                    int16_t mvd_y = best_mv.y - mv_pred.y;

                    /* REC: MV deltas */
                    uint16_t mvx_sym = s2u(mvd_x);
                    REC(cdf_mvx, mvx_sym);
                    if (mvx_sym == 255)
                        snt_bw_put_sexpgolomb(&ov, mvd_x);

                    uint16_t mvy_sym = s2u(mvd_y);
                    REC(cdf_mvy, mvy_sym);
                    if (mvy_sym == 255)
                        snt_bw_put_sexpgolomb(&ov, mvd_y);

                    mv_field[blk_idx] = best_mv;
                }
            }

            /* ── DCT input: residual (inter) or pixels (intra) ──── */
            int16_t dct_input[64];
            if (block_is_inter && !block_is_skip) {
                /* Inter coded: DCT of residual = cur - MC prediction */
                for (int dy = 0; dy < BLK; dy++)
                    for (int dx = 0; dx < BLK; dx++) {
                        int py = py0 + dy, px_v = px0 + dx;
                        int ry = py + best_mv.y, rx = px_v + best_mv.x;
                        int16_t cur_val = (py < height && px_v < width) ?
                            work[py * pad_w + px_v] : 0;
                        int16_t ref_val = (ry >= 0 && ry < height && rx >= 0 && rx < width) ?
                            (int16_t)ctx->ref_y[ry * width + rx] : 0;
                        dct_input[dy * BLK + dx] = cur_val - ref_val;
                    }
            } else {
                /* Intra (key or inter-frame intra): use work buffer directly */
                for (int dy = 0; dy < BLK; dy++)
                    for (int dx = 0; dx < BLK; dx++)
                        dct_input[dy * BLK + dx] = work[(py0 + dy) * pad_w + (px0 + dx)];
            }

            int16_t coeffs[64];
            snt_dct_fwd_2d_ref(dct_input, BLK, BLK, coeffs);

            /* DC prediction */
            int16_t dc_pred = 0;
            if (bx > 0 && by > 0)
                dc_pred = (dc_recon_left + dc_recon_above[bx] + 1) >> 1;
            else if (bx > 0)
                dc_pred = dc_recon_left;
            else if (by > 0)
                dc_pred = dc_recon_above[bx];

            int16_t dc_residual = coeffs[0] - dc_pred;
            int16_t dc_level = pvq_quantize_dc(dc_residual, step, 1);

            uint16_t dc_sym = s2u(dc_level);
            REC(cdf_dc, dc_sym);
            if (dc_sym == 255)
                snt_bw_put_sexpgolomb(&ov, dc_level);

            int16_t dc_recon = pvq_dequantize_dc(dc_level, step) + dc_pred;
            dc_recon_left = dc_recon;
            dc_recon_above[bx] = dc_recon;

            /* Low-AC prediction (key frames only — inter residuals don't benefit) */
            int16_t ac_pred[64];
            memset(ac_pred, 0, sizeof(ac_pred));
            if (!is_inter) {
                for (int p = 0; p < NUM_PRED_AC; p++) {
                    int idx = pred_ac_idx[p];
                    int16_t pred_val = 0;
                    if (idx < 8 && bx > 0)
                        pred_val = ac_recon_left[p];
                    else if (idx >= 8 && by > 0)
                        pred_val = ac_recon_above[bx * NUM_PRED_AC + p];
                    ac_pred[idx] = pred_val;
                    coeffs[idx] -= pred_val;
                }
            }

            /* Band processing (same for all block types) */
            int16_t recon_ac[63];
            memset(recon_ac, 0, sizeof(recon_ac));

            for (int b = 0; b < NUM_BANDS; b++) {
                int start = band_bounds[b];
                int n = band_bounds[b + 1] - start;

                float csf_sum = 0.0f;
                for (int i = 0; i < n; i++)
                    csf_sum += psy_csf_weight(BLK, 1 + start + i);
                uint8_t band_csf = (uint8_t)(csf_sum / n * 255.0f);
                if (band_csf < 1) band_csf = 1;

                int16_t *band_coeffs = &coeffs[1 + start];
                int32_t max_sq = 0, sum_sq = 0;
                for (int i = 0; i < n; i++) {
                    int32_t sq = (int32_t)band_coeffs[i] * band_coeffs[i];
                    sum_sq += sq;
                    if (sq > max_sq) max_sq = sq;
                }
                int use_scalar = (sum_sq > 0 && max_sq * 10 > sum_sq * 6);

                uint16_t mode_sym = use_scalar ? 1 : 0;
                REC(cdf_mode, mode_sym);

                int16_t band_recon[16];

                if (use_scalar) {
                    uint16_t eff_step = step;
                    if (band_csf < 128) {
                        eff_step = (uint16_t)((uint32_t)step * 256 / band_csf);
                        if (eff_step > 3072) eff_step = 3072;
                    }
                    for (int i = 0; i < n; i++) {
                        int16_t level = pvq_quantize_dc(band_coeffs[i], eff_step, 1);
                        uint16_t ac_sym = s2u(level);
                        REC(cdf_ac, ac_sym);
                        if (ac_sym == 255)
                            snt_bw_put_sexpgolomb(&ov, level);
                        band_recon[i] = pvq_dequantize_dc(level, eff_step);
                    }
                } else {
                    SntPVQCode code;
                    pvq_encode(band_coeffs, n, step, 1, band_csf, &code, band_recon);

                    uint16_t skip_sym = code.skip ? 1 : 0;
                    REC(cdf_skip, skip_sym);

                    if (!code.skip) {
                        uint16_t gain_sym = (code.gain > 255) ? 255 : (uint16_t)code.gain;
                        REC(cdf_gain, gain_sym);
                        if (gain_sym == 255)
                            snt_bw_put_expgolomb(&ov, code.gain);

                        uint16_t k_sym = (code.k > 255) ? 255 : (uint16_t)code.k;
                        REC(cdf_k, k_sym);
                        if (k_sym == 255)
                            snt_bw_put_expgolomb(&ov, (uint32_t)code.k);

                        int shape_bits = pvq_shape_bits(n, code.k);
                        if (shape_bits < 1) shape_bits = 1;
                        if (shape_bits > 64) shape_bits = 64;
                        if (shape_bits > 32) {
                            snt_bw_put_bits(&ov, (uint32_t)(code.shape_idx >> 32), shape_bits - 32);
                            snt_bw_put_bits(&ov, (uint32_t)(code.shape_idx & 0xFFFFFFFF), 32);
                        } else {
                            snt_bw_put_bits(&ov, (uint32_t)code.shape_idx, shape_bits);
                        }
                    }
                }
                for (int i = 0; i < n; i++)
                    recon_ac[start + i] = band_recon[i];
            }

            /* Reconstruct coefficients */
            int16_t recon_coeffs[64];
            recon_coeffs[0] = dc_recon;
            for (int i = 0; i < 63; i++)
                recon_coeffs[i + 1] = recon_ac[i];

            /* Add back AC prediction (key frames only) */
            if (!is_inter) {
                for (int p = 0; p < NUM_PRED_AC; p++) {
                    int idx = pred_ac_idx[p];
                    recon_coeffs[idx] += ac_pred[idx];
                }
            }

            /* Store reconstructed AC for future prediction (key frames only —
             * inter frames don't use AC prediction, so this avoids dead writes) */
            if (!is_inter) {
                for (int p = 0; p < NUM_PRED_AC; p++) {
                    int idx = pred_ac_idx[p];
                    ac_recon_left[p] = recon_coeffs[idx];
                    ac_recon_above[bx * NUM_PRED_AC + p] = recon_coeffs[idx];
                }
            }

            /* IDCT → recon buffer */
            int16_t recon_blk[64];
            snt_dct_inv_2d_ref(recon_coeffs, BLK, recon_blk, BLK);

            /* Write to recon: add MC prediction for inter coded blocks */
            for (int dy = 0; dy < BLK; dy++)
                for (int dx = 0; dx < BLK; dx++) {
                    int py = py0 + dy, px_v = px0 + dx;
                    int16_t val = recon_blk[dy * BLK + dx];
                    if (block_is_inter && !block_is_skip) {
                        int ry = py + best_mv.y, rx = px_v + best_mv.x;
                        int16_t ref_val = (ry >= 0 && ry < height && rx >= 0 && rx < width) ?
                            (int16_t)ctx->ref_y[ry * width + rx] : 0;
                        val += ref_val;
                    }
                    recon[py * pad_w + px_v] = val;
                }
        }
        memset(ac_recon_left, 0, sizeof(ac_recon_left));
    }

    /* KEY ONLY: lapping post-filter */
    if (!is_inter)
        snt_lap_post_filter_2d_ref(recon, pad_w, pad_h, pad_w, BLK);

    /* Chroma: 4x4 DCT + scalar quantization (same for key and inter) */
    if (input_u && input_v) {
        int cblk_cols = cpad_w / cblk;
        int cqp = ctx->header.base_qp + ctx->header.chroma_qp_offset;
        if (cqp < 0) cqp = 0; if (cqp > 63) cqp = 63;
        uint16_t cstep = snt_qp_to_step((uint8_t)cqp);

        const uint8_t *chroma_planes[2] = {input_u, input_v};
        for (int plane = 0; plane < 2; plane++) {
            const uint8_t *cp = chroma_planes[plane];
            int16_t c_dc_left = 0;
            int16_t *c_dc_above = (int16_t *)calloc(cblk_cols, sizeof(int16_t));

            for (int cby = 0; cby < cpad_h; cby += cblk) {
                c_dc_left = 0;
                for (int cbx = 0; cbx < cpad_w; cbx += cblk) {
                    int16_t cblock[16];
                    for (int dy = 0; dy < 4; dy++)
                        for (int dx = 0; dx < 4; dx++) {
                            int py = cby + dy, px = cbx + dx;
                            cblock[dy*4+dx] = (py < ch && px < cw) ?
                                (int16_t)cp[py*stride_uv+px] - 128 : 0;
                        }

                    int16_t ccoeffs[16];
                    snt_dct_fwd_2d_ref(cblock, 4, 4, ccoeffs);

                    int ci = cbx / cblk;
                    int16_t c_dc_pred = 0;
                    if (cbx > 0 && cby > 0)
                        c_dc_pred = (c_dc_left + c_dc_above[ci] + 1) >> 1;
                    else if (cbx > 0) c_dc_pred = c_dc_left;
                    else if (cby > 0) c_dc_pred = c_dc_above[ci];

                    int16_t dc_res = ccoeffs[0] - c_dc_pred;
                    int16_t cdc_level = pvq_quantize_dc(dc_res, cstep, 1);

                    uint16_t cdc_sym = s2u(cdc_level);
                    REC(cdf_dc, cdc_sym);
                    if (cdc_sym == 255)
                        snt_bw_put_sexpgolomb(&ov, cdc_level);

                    int16_t cdc_recon = pvq_dequantize_dc(cdc_level, cstep) + c_dc_pred;
                    c_dc_left = cdc_recon;
                    c_dc_above[ci] = cdc_recon;

                    for (int i = 1; i < 16; i++) {
                        int16_t ac_level = pvq_quantize_dc(ccoeffs[i], cstep, 1);
                        uint16_t ac_sym = s2u(ac_level);
                        REC(cdf_ac, ac_sym);
                        if (ac_sym == 255)
                            snt_bw_put_sexpgolomb(&ov, ac_level);
                    }
                }
            }
            free(c_dc_above);
        }
    }

    /* ── PHASE 2: Reverse rANS encode ──────────────────────────────── */
    uint32_t rev_cap = rc * 2 + 64;
    uint16_t *rev_buf = (uint16_t *)malloc(rev_cap * sizeof(uint16_t));
    RansEncoder renc;
    rans_enc_init(&renc, rev_buf, rev_cap, 1);

    for (int32_t i = (int32_t)rc - 1; i >= 0; i--)
        rans_enc_put(&renc, ranges[i].cdf_lo, ranges[i].freq, (uint32_t)i);
    rans_enc_flush(&renc);

    /* ── PHASE 3: Write to output bitstream ────────────────────────── */
    uint32_t rans_words;
    const uint16_t *rans_data = rans_enc_data(&renc, &rans_words);
    uint32_t rans_bytes = rans_words * 2;

    snt_bw_put_bits(&bw, rans_bytes, 32);
    for (uint32_t i = 0; i < rans_words; i++) {
        snt_bw_put_byte(&bw, (uint8_t)(rans_data[i] & 0xFF));
        snt_bw_put_byte(&bw, (uint8_t)(rans_data[i] >> 8));
    }

    snt_bw_byte_align(&ov);
    uint32_t ov_len = snt_bw_bytes_written(&ov);
    for (uint32_t i = 0; i < ov_len; i++)
        snt_bw_put_byte(&bw, ov_buf[i]);

    snt_bw_byte_align(&bw);
    *bs_len = snt_bw_bytes_written(&bw);

    /* ── Store reconstructed frame as reference ────────────────────── */
    if (!ctx->ref_y || ctx->ref_w != width || ctx->ref_h != height) {
        free(ctx->ref_y);
        ctx->ref_y = (uint8_t *)malloc((size_t)width * height);
        ctx->ref_w = width;
        ctx->ref_h = height;
    }
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
            ctx->ref_y[y * width + x] = snt_clip_pixel((int32_t)recon[y * pad_w + x]);
    ctx->num_refs = 1;

    /* Cleanup */
    free(work);
    free(recon);
    free(dc_recon_above);
    free(ac_recon_above);
    free(mv_field);
    free(ranges);
    free(ov_buf);
    free(rev_buf);
    return 0;
}
