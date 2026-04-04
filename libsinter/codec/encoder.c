/*
 * Sinter Codec — Encoder with rANS entropy coding
 * Two-pass: forward record + reverse rANS encode.
 * Spec v3: team consensus from roundtable (Jim/Timothy/Gary/Monty/Jens-Rainer).
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
#include "../util/fixed_point.h"
#include "../perceptual/activity_mask.h"
#include "../perceptual/csf_weights.h"

#define BLK 8

/* ── Zigzag signed-to-unsigned mapping (Timothy) ──────────────────── */
static inline uint16_t s2u(int16_t v) {
    uint32_t u = (v >= 0) ? (uint32_t)(2 * v) : (uint32_t)(-2 * v - 1);
    return (u > 254) ? 255 : (uint16_t)u;
}

/* ── Symbol range recording (Monty's REC macro) ──────────────────── */
typedef struct { uint16_t cdf_lo; uint16_t freq; } SymRange;

/*
 * Symbol emission order (normative, must match decoder exactly):
 *
 * FOR each luma block (raster: by, bx):
 *   REC(cdf_dc, s2u(dc_level))
 *   // Note: low-AC prediction (coeffs 1,8,2,9) is subtracted from coeffs
 *   // BEFORE band processing. Predicted residuals flow through their bands'
 *   // scalar/PVQ paths — they are NOT separate rANS symbols.
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

    /* ── SECTION A: Fixed header (Jim) ─────────────────────────────── */
    /* Written ONCE to bw. Nothing else touches bw until Phase 3. */
    snt_bw_byte_align(&bw);
    snt_bw_put_bits(&bw, (uint32_t)total_blocks, 16);
    snt_bw_put_bits(&bw, (input_u && input_v) ? 1u : 0u, 1); /* chroma_present */
    snt_bw_byte_align(&bw);

    /* ── Allocate rANS working buffers ────────────────────────────── */
    /* Chroma block count */
    int cw = width / 2, ch = height / 2;
    int cblk = 4;
    int cpad_w = ((cw + cblk - 1) / cblk) * cblk;
    int cpad_h = ((ch + cblk - 1) / cblk) * cblk;
    int chroma_blocks = (input_u && input_v) ? 2 * (cpad_w / cblk) * (cpad_h / cblk) : 0;

    uint32_t max_ranges = (uint32_t)(total_blocks * 80 + chroma_blocks * 16 + 256);
    SymRange *ranges = (SymRange *)malloc(max_ranges * sizeof(SymRange));
    uint32_t rc = 0;

    /* Overflow buffer (separate from bw!) */
    uint32_t ov_cap = (uint32_t)(total_blocks * 128 + 4096);
    uint8_t *ov_buf = (uint8_t *)malloc(ov_cap);
    SntBitWriter ov;
    snt_bw_init(&ov, ov_buf, ov_cap);

    /* ── Init 6 CDFs (Gary/Timothy) ──────────────────────────────── */
    RansCDF cdf_dc, cdf_ac, cdf_mode, cdf_skip, cdf_gain, cdf_k;
    rans_cdf_init_uniform(&cdf_dc,   256, RANS_ADAPT_RATE_COEFF);
    rans_cdf_init_uniform(&cdf_ac,   256, RANS_ADAPT_RATE_COEFF);
    rans_cdf_init_uniform(&cdf_mode,   2, RANS_ADAPT_RATE_MODE);
    rans_cdf_init_uniform(&cdf_skip,   2, RANS_ADAPT_RATE_MODE);
    rans_cdf_init_uniform(&cdf_gain, 256, RANS_ADAPT_RATE_COEFF);
    rans_cdf_init_uniform(&cdf_k,    256, RANS_ADAPT_RATE_COEFF);

    /* REC macro: record symbol range + adapt CDF (Monty) */
    #define REC(cdf, sym) do { \
        rans_cdf_symbol_range(&(cdf), (sym), &ranges[rc].cdf_lo, &ranges[rc].freq); \
        rans_cdf_adapt(&(cdf), (sym)); \
        rc++; \
    } while(0)

    /* ── PHASE 1: Forward pass — record symbols + overflow ─────── */
    /* DC + low-AC prediction buffers */
    int16_t *dc_recon_above = (int16_t *)calloc(blk_cols, sizeof(int16_t));
    int16_t dc_recon_left = 0;
    #define NUM_PRED_AC 4
    static const int pred_ac_idx[NUM_PRED_AC] = {1, 8, 2, 9};
    int16_t *ac_recon_above = (int16_t *)calloc((size_t)blk_cols * NUM_PRED_AC, sizeof(int16_t));
    int16_t ac_recon_left[NUM_PRED_AC] = {0};

    for (int by = 0; by < blk_rows; by++) {
        dc_recon_left = 0;
        for (int bx = 0; bx < blk_cols; bx++) {
            int16_t coeffs[64];
            snt_dct_fwd_2d_ref(&work[by * BLK * pad_w + bx * BLK], pad_w, BLK, coeffs);

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

            /* REC: DC residual → cdf_dc (zigzag mapped) */
            uint16_t dc_sym = s2u(dc_level);
            REC(cdf_dc, dc_sym);
            if (dc_sym == 255)
                snt_bw_put_sexpgolomb(&ov, dc_level); /* escape to overflow */

            /* Reconstruct DC */
            int16_t dc_recon = pvq_dequantize_dc(dc_level, step) + dc_pred;
            dc_recon_left = dc_recon;
            dc_recon_above[bx] = dc_recon;

            /* Low-AC prediction */
            int16_t ac_pred[64];
            memset(ac_pred, 0, sizeof(ac_pred));
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

            /* Band processing */
            #define NUM_BANDS 16
            static const int band_bounds[NUM_BANDS + 1] = {0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,59,63};
            int16_t recon_ac[63];
            memset(recon_ac, 0, sizeof(recon_ac));

            for (int b = 0; b < NUM_BANDS; b++) {
                int start = band_bounds[b];
                int n = band_bounds[b + 1] - start;

                /* CSF weight */
                float csf_sum = 0.0f;
                for (int i = 0; i < n; i++)
                    csf_sum += psy_csf_weight(BLK, 1 + start + i);
                uint8_t band_csf = (uint8_t)(csf_sum / n * 255.0f);
                if (band_csf < 1) band_csf = 1;

                /* Mode decision: scalar if energy concentrated */
                int16_t *band_coeffs = &coeffs[1 + start];
                int32_t max_sq = 0, sum_sq = 0;
                for (int i = 0; i < n; i++) {
                    int32_t sq = (int32_t)band_coeffs[i] * band_coeffs[i];
                    sum_sq += sq;
                    if (sq > max_sq) max_sq = sq;
                }
                int use_scalar = (sum_sq > 0 && max_sq * 10 > sum_sq * 6);

                /* REC: mode → cdf_mode */
                uint16_t mode_sym = use_scalar ? 1 : 0;
                REC(cdf_mode, mode_sym);

                int16_t band_recon[16];

                if (use_scalar) {
                    /* Scalar mode with CSF-weighted step */
                    uint16_t eff_step = step;
                    if (band_csf < 128) {
                        eff_step = (uint16_t)((uint32_t)step * 256 / band_csf);
                        if (eff_step > 3072) eff_step = 3072;
                    }
                    for (int i = 0; i < n; i++) {
                        int16_t level = pvq_quantize_dc(band_coeffs[i], eff_step, 1);
                        /* REC: scalar level → cdf_ac (zigzag) */
                        uint16_t ac_sym = s2u(level);
                        REC(cdf_ac, ac_sym);
                        if (ac_sym == 255)
                            snt_bw_put_sexpgolomb(&ov, level);
                        band_recon[i] = pvq_dequantize_dc(level, eff_step);
                    }
                } else {
                    /* PVQ mode */
                    SntPVQCode code;
                    pvq_encode(band_coeffs, n, step, 1, band_csf, &code, band_recon);

                    /* REC: skip → cdf_skip */
                    uint16_t skip_sym = code.skip ? 1 : 0;
                    REC(cdf_skip, skip_sym);

                    if (!code.skip) {
                        /* REC: gain → cdf_gain */
                        uint16_t gain_sym = (code.gain > 255) ? 255 : (uint16_t)code.gain;
                        REC(cdf_gain, gain_sym);
                        if (gain_sym == 255)
                            snt_bw_put_expgolomb(&ov, code.gain);

                        /* REC: k → cdf_k */
                        uint16_t k_sym = (code.k > 255) ? 255 : (uint16_t)code.k;
                        REC(cdf_k, k_sym);
                        if (k_sym == 255)
                            snt_bw_put_expgolomb(&ov, (uint32_t)code.k);

                        /* Shape index → raw bits to overflow */
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

            /* Reconstruct for reference */
            int16_t recon_coeffs[64];
            recon_coeffs[0] = dc_recon;
            for (int i = 0; i < 63; i++)
                recon_coeffs[i + 1] = recon_ac[i];

            /* Add back AC prediction */
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
        memset(ac_recon_left, 0, sizeof(ac_recon_left));
    }

    /* Post-filter */
    snt_lap_post_filter_2d_ref(recon, pad_w, pad_h, pad_w, BLK);

    /* Chroma: 4x4 DCT + scalar quantization */
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

                    /* DC prediction */
                    int ci = cbx / cblk;
                    int16_t c_dc_pred = 0;
                    if (cbx > 0 && cby > 0)
                        c_dc_pred = (c_dc_left + c_dc_above[ci] + 1) >> 1;
                    else if (cbx > 0) c_dc_pred = c_dc_left;
                    else if (cby > 0) c_dc_pred = c_dc_above[ci];

                    int16_t dc_res = ccoeffs[0] - c_dc_pred;
                    int16_t cdc_level = pvq_quantize_dc(dc_res, cstep, 1);

                    /* REC: chroma DC → cdf_dc (shared with luma) */
                    uint16_t cdc_sym = s2u(cdc_level);
                    REC(cdf_dc, cdc_sym);
                    if (cdc_sym == 255)
                        snt_bw_put_sexpgolomb(&ov, cdc_level);

                    int16_t cdc_recon = pvq_dequantize_dc(cdc_level, cstep) + c_dc_pred;
                    c_dc_left = cdc_recon;
                    c_dc_above[ci] = cdc_recon;

                    /* REC: chroma AC → cdf_ac (shared with luma) */
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
    rans_enc_init(&renc, rev_buf, rev_cap, 1); /* dual interleaved */

    for (int32_t i = (int32_t)rc - 1; i >= 0; i--)
        rans_enc_put(&renc, ranges[i].cdf_lo, ranges[i].freq, (uint32_t)i);
    rans_enc_flush(&renc);

    /* ── PHASE 3: Write to output bitstream ────────────────────────── */
    /* Section B: rANS payload */
    uint32_t rans_words;
    const uint16_t *rans_data = rans_enc_data(&renc, &rans_words);
    uint32_t rans_bytes = rans_words * 2;

    snt_bw_put_bits(&bw, rans_bytes, 32);
    for (uint32_t i = 0; i < rans_words; i++) {
        snt_bw_put_byte(&bw, (uint8_t)(rans_data[i] & 0xFF));
        snt_bw_put_byte(&bw, (uint8_t)(rans_data[i] >> 8));
    }

    /* Section C: overflow */
    snt_bw_byte_align(&ov);
    uint32_t ov_len = snt_bw_bytes_written(&ov);
    for (uint32_t i = 0; i < ov_len; i++)
        snt_bw_put_byte(&bw, ov_buf[i]);

    snt_bw_byte_align(&bw);
    *bs_len = snt_bw_bytes_written(&bw);

    /* Cleanup */
    free(work);
    free(recon);
    free(dc_recon_above);
    free(ac_recon_above);
    free(ranges);
    free(ov_buf);
    free(rev_buf);
    return 0;
}
