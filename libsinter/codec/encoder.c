/*
 * Sinter Codec — Encoder (rANS adaptive entropy coding)
 * Fixed 8x8 blocks, two-pass rANS encoding for luma + chroma.
 * This replaces the exp-Golomb VLC bitstream with rANS.
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

/* ── Signed-to-unsigned zigzag mapping for rANS symbols ────────────── */
/*
 * Maps signed coefficient values to unsigned symbols [0..255].
 * Symbol 255 is the escape code: the true value is written in the
 * raw overflow section that follows the rANS payload.
 *
 *  val ->  0 => sym 0
 *  val ->  1 => sym 2,  val -> -1 => sym 1
 *  val ->  2 => sym 4,  val -> -2 => sym 3
 *  ...
 *  val -> 127 => sym 254, val -> -128 => sym 255 (escape)
 */
static uint16_t signed_to_sym(int16_t val) {
    uint32_t uval = (val >= 0) ? (uint32_t)(2 * val) : (uint32_t)(-2 * val - 1);
    return (uval > 255) ? 255 : (uint16_t)uval;
}

static int escaped(int16_t val) {
    uint32_t uval = (val >= 0) ? (uint32_t)(2 * val) : (uint32_t)(-2 * val - 1);
    return uval > 255;
}

/* ── Two-pass rANS symbol recording ────────────────────────────────── */
typedef struct { uint16_t cdf_lo; uint16_t freq; } SymRange;

#define RECORD_SYM(cdf_ptr, sym) do { \
    rans_cdf_symbol_range(cdf_ptr, sym, &ranges[range_count].cdf_lo, \
                          &ranges[range_count].freq); \
    rans_cdf_adapt(cdf_ptr, sym); \
    range_count++; \
} while(0)

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

    /* ================================================================
     * SYMBOL ORDER (normative — must match decoder exactly):
     *
     * LUMA (per block, blocks in raster order):
     *   dc_sym          (cdf_dc, 256 symbols)   zigzag mapped
     *   4 x ac_pred_sym (cdf_low_ac, 256)        zigzag mapped
     *   For each of 16 bands:
     *     mode_sym       (cdf_mode, 2)            0=PVQ, 1=scalar
     *     If PVQ mode:
     *       skip_sym     (cdf_skip, 2)            0=active, 1=skip
     *       If !skip: gain_lo (cdf_gain, 256), k_lo (cdf_k, 256)
     *     If scalar mode:
     *       n x ac_level_sym (cdf_scalar_ac, 256) zigzag mapped
     *
     * CHROMA (per plane, U then V; per 4x4 block in raster order):
     *   dc_sym           (cdf_chroma_dc, 256)     zigzag mapped
     *   15 x ac_sym      (cdf_chroma_ac, 256)     zigzag mapped
     *
     * RAW OVERFLOW (after rANS payload, in forward order):
     *   For each luma block:
     *     dc: if escaped -> sign(1 bit) + exp-golomb(|overflow|)
     *     4 x ac_pred: if escaped -> sign(1) + exp-golomb(|overflow|)
     *     For each of 16 bands:
     *       If PVQ & !skip: gain_hi(8 bits if gain>255),
     *                        k_hi(8 bits if k>255),
     *                        shape_idx(pvq_shape_bits raw)
     *       If scalar: per coeff if escaped -> sign(1) + exp-golomb
     *   For each chroma block (U then V):
     *     dc + 15 AC: if escaped -> sign(1) + exp-golomb(|overflow|)
     * ================================================================ */

    /* DC + low-AC prediction buffers */
    int16_t *dc_recon_above = (int16_t *)calloc(blk_cols, sizeof(int16_t));
    int16_t dc_recon_left = 0;
    #define NUM_PRED_AC 4
    static const int pred_ac_idx[NUM_PRED_AC] = {1, 8, 2, 9};
    int16_t *ac_recon_above = (int16_t *)calloc((size_t)blk_cols * NUM_PRED_AC, sizeof(int16_t));
    int16_t ac_recon_left[NUM_PRED_AC] = {0};

    /* Allocate symbol range buffer for two-pass rANS.
     * Worst case per luma block: 1 DC + 4 AC_pred + 16*(1 mode + 1 skip + 2 + 4 scalar) = ~133
     * Per chroma block: 16 symbols. Budget generously. */
    uint32_t max_luma_syms = (uint32_t)total_blocks * 140;
    /* Chroma: two planes, each up to (cpad_w/4)*(cpad_h/4) blocks * 16 symbols */
    int cw = width / 2, ch = height / 2;
    int cblk = 4;
    int cpad_w = ((cw + cblk - 1) / cblk) * cblk;
    int cpad_h = ((ch + cblk - 1) / cblk) * cblk;
    uint32_t max_chroma_syms = (uint32_t)(cpad_w / cblk) * (cpad_h / cblk) * 16 * 2;
    uint32_t max_syms = max_luma_syms + max_chroma_syms;

    SymRange *ranges = (SymRange *)malloc(max_syms * sizeof(SymRange));
    if (!ranges) { free(work); free(recon); free(dc_recon_above); free(ac_recon_above); return -1; }
    uint32_t range_count = 0;

    /* Overflow bitwriter: raw bits for escaped values, shape indices, etc. */
    uint32_t overflow_cap = bs_capacity / 2;
    uint8_t *overflow_buf = (uint8_t *)calloc(overflow_cap, 1);
    if (!overflow_buf) { free(ranges); free(work); free(recon); free(dc_recon_above); free(ac_recon_above); return -1; }
    SntBitWriter ov;
    snt_bw_init(&ov, overflow_buf, overflow_cap);

    /* Initialize CDFs */
    RansCDF cdf_dc, cdf_low_ac, cdf_mode, cdf_skip, cdf_gain, cdf_k, cdf_scalar_ac;
    RansCDF cdf_chroma_dc, cdf_chroma_ac;
    rans_cdf_init_uniform(&cdf_dc, 256, RANS_ADAPT_RATE_COEFF);
    rans_cdf_init_uniform(&cdf_low_ac, 256, RANS_ADAPT_RATE_COEFF);
    rans_cdf_init_uniform(&cdf_mode, 2, RANS_ADAPT_RATE_MODE);
    rans_cdf_init_uniform(&cdf_skip, 2, RANS_ADAPT_RATE_MODE);
    rans_cdf_init_uniform(&cdf_gain, 256, RANS_ADAPT_RATE_COEFF);
    rans_cdf_init_uniform(&cdf_k, 256, RANS_ADAPT_RATE_COEFF);
    rans_cdf_init_uniform(&cdf_scalar_ac, 256, RANS_ADAPT_RATE_COEFF);
    rans_cdf_init_uniform(&cdf_chroma_dc, 256, RANS_ADAPT_RATE_COEFF);
    rans_cdf_init_uniform(&cdf_chroma_ac, 256, RANS_ADAPT_RATE_COEFF);

    /* ════════════════════════════════════════════════════════════════
     *  PASS 1 (forward): Process blocks, record symbol ranges,
     *                     write overflow to ov bitwriter.
     * ════════════════════════════════════════════════════════════════ */

    #define NUM_BANDS 16
    static const int band_bounds[NUM_BANDS + 1] = {0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,59,63};

    for (int by = 0; by < blk_rows; by++) {
        dc_recon_left = 0;
        for (int bx = 0; bx < blk_cols; bx++) {
            int16_t coeffs[64];
            snt_dct_fwd_2d_ref(&work[by * BLK * pad_w + bx * BLK], pad_w, BLK, coeffs);

            /* DC prediction from reconstructed neighbors */
            int16_t dc_pred = 0;
            if (bx > 0 && by > 0)
                dc_pred = (dc_recon_left + dc_recon_above[bx] + 1) >> 1;
            else if (bx > 0)
                dc_pred = dc_recon_left;
            else if (by > 0)
                dc_pred = dc_recon_above[bx];

            /* Quantize DC residual */
            int16_t dc_residual = coeffs[0] - dc_pred;
            int16_t dc_level = pvq_quantize_dc(dc_residual, step, 1);

            /* Record DC symbol via rANS */
            uint16_t dc_sym = signed_to_sym(dc_level);
            RECORD_SYM(&cdf_dc, dc_sym);

            /* Write DC overflow if escaped */
            if (escaped(dc_level)) {
                uint32_t uval = (dc_level >= 0) ? (uint32_t)(2 * dc_level) : (uint32_t)(-2 * dc_level - 1);
                uint32_t overflow_val = uval - 255;
                snt_bw_put_bits(&ov, dc_level < 0 ? 1 : 0, 1);
                snt_bw_put_expgolomb(&ov, overflow_val);
            }

            /* Reconstruct DC for prediction */
            int16_t dc_recon = pvq_dequantize_dc(dc_level, step) + dc_pred;
            dc_recon_left = dc_recon;
            dc_recon_above[bx] = dc_recon;

            /* Low-AC prediction */
            int16_t ac_pred[64];
            memset(ac_pred, 0, sizeof(ac_pred));
            for (int p = 0; p < NUM_PRED_AC; p++) {
                int idx = pred_ac_idx[p];
                int16_t pred_val = 0;
                if (idx < 8 && bx > 0) {
                    pred_val = ac_recon_left[p];
                } else if (idx >= 8 && by > 0) {
                    pred_val = ac_recon_above[bx * NUM_PRED_AC + p];
                }
                ac_pred[idx] = pred_val;
                coeffs[idx] -= pred_val;
            }

            /* Record 4 low-AC prediction symbols */
            for (int p = 0; p < NUM_PRED_AC; p++) {
                int idx = pred_ac_idx[p];
                int16_t ac_level = pvq_quantize_dc(coeffs[idx], step, 1);
                uint16_t ac_sym = signed_to_sym(ac_level);
                RECORD_SYM(&cdf_low_ac, ac_sym);

                if (escaped(ac_level)) {
                    uint32_t uval = (ac_level >= 0) ? (uint32_t)(2 * ac_level) : (uint32_t)(-2 * ac_level - 1);
                    uint32_t overflow_val = uval - 255;
                    snt_bw_put_bits(&ov, ac_level < 0 ? 1 : 0, 1);
                    snt_bw_put_expgolomb(&ov, overflow_val);
                }
                /* Store quantized level back for band processing */
                coeffs[idx] = ac_level;
            }

            /* Band processing */
            int16_t recon_ac[63];
            memset(recon_ac, 0, sizeof(recon_ac));

            for (int b = 0; b < NUM_BANDS; b++) {
                int start = band_bounds[b];
                int n = band_bounds[b + 1] - start;

                /* Per-band CSF weight */
                float csf_sum = 0.0f;
                for (int i = 0; i < n; i++)
                    csf_sum += psy_csf_weight(BLK, 1 + start + i);
                uint8_t band_csf = (uint8_t)(csf_sum / n * 255.0f);
                if (band_csf < 1) band_csf = 1;

                /* Mode decision: scalar vs PVQ */
                int16_t *band_coeffs = &coeffs[1 + start];
                int32_t max_sq = 0, sum_sq = 0;
                for (int i = 0; i < n; i++) {
                    int32_t sq = (int32_t)band_coeffs[i] * band_coeffs[i];
                    sum_sq += sq;
                    if (sq > max_sq) max_sq = sq;
                }
                int use_scalar = (sum_sq > 0 && max_sq * 10 > sum_sq * 6);

                int16_t band_recon[16];

                /* Record mode symbol */
                uint16_t mode_sym = use_scalar ? 1 : 0;
                RECORD_SYM(&cdf_mode, mode_sym);

                if (use_scalar) {
                    /* Scalar mode with CSF-weighted step */
                    uint16_t eff_step = step;
                    if (band_csf < 128) {
                        eff_step = (uint16_t)((uint32_t)step * 256 / band_csf);
                        if (eff_step > 3072) eff_step = 3072;
                    }
                    for (int i = 0; i < n; i++) {
                        int16_t level = pvq_quantize_dc(band_coeffs[i], eff_step, 1);
                        uint16_t sym = signed_to_sym(level);
                        RECORD_SYM(&cdf_scalar_ac, sym);

                        if (escaped(level)) {
                            uint32_t uval = (level >= 0) ? (uint32_t)(2 * level) : (uint32_t)(-2 * level - 1);
                            uint32_t overflow_val = uval - 255;
                            snt_bw_put_bits(&ov, level < 0 ? 1 : 0, 1);
                            snt_bw_put_expgolomb(&ov, overflow_val);
                        }
                        band_recon[i] = pvq_dequantize_dc(level, eff_step);
                    }
                } else {
                    /* PVQ mode */
                    SntPVQCode code;
                    pvq_encode(band_coeffs, n, step, 1, band_csf, &code, band_recon);

                    /* Record skip symbol */
                    uint16_t skip_sym = code.skip;
                    RECORD_SYM(&cdf_skip, skip_sym);

                    if (!code.skip) {
                        /* Record gain_lo and k_lo (clamp to 255, overflow in raw section) */
                        uint16_t gain_lo = (code.gain > 255) ? 255 : (uint16_t)code.gain;
                        uint16_t k_lo = (code.k > 255) ? 255 : (uint16_t)code.k;
                        RECORD_SYM(&cdf_gain, gain_lo);
                        RECORD_SYM(&cdf_k, k_lo);

                        /* Raw overflow: gain_hi, k_hi, shape_idx */
                        if (code.gain > 255)
                            snt_bw_put_bits(&ov, (uint32_t)(code.gain - 255), 8);
                        if (code.k > 255)
                            snt_bw_put_bits(&ov, (uint32_t)(code.k - 255), 8);

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

            /* Reconstruct for reference (same as original) */
            int16_t recon_coeffs[64];
            recon_coeffs[0] = dc_recon;
            for (int i = 0; i < 63; i++)
                recon_coeffs[i + 1] = recon_ac[i];

            /* Add back AC prediction for low-AC coefficients.
             * Note: for predicted coefficients, recon_ac stored the quantized
             * level (not dequantized), so we need to dequantize first. */
            for (int p = 0; p < NUM_PRED_AC; p++) {
                int idx = pred_ac_idx[p];
                recon_coeffs[idx] = pvq_dequantize_dc(coeffs[idx], step) + ac_pred[idx];
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

    /* Chroma: 4x4 DCT + scalar quantization, rANS coded */
    int has_chroma = (input_u && input_v);
    if (has_chroma) {
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
                    int cbxi = cbx / cblk;
                    int16_t c_dc_pred = 0;
                    if (cbx > 0 && cby > 0)
                        c_dc_pred = (c_dc_left + c_dc_above[cbxi] + 1) >> 1;
                    else if (cbx > 0) c_dc_pred = c_dc_left;
                    else if (cby > 0) c_dc_pred = c_dc_above[cbxi];

                    int16_t dc_res = ccoeffs[0] - c_dc_pred;
                    int16_t dc_lev = pvq_quantize_dc(dc_res, cstep, 1);

                    /* Record chroma DC */
                    uint16_t dc_s = signed_to_sym(dc_lev);
                    RECORD_SYM(&cdf_chroma_dc, dc_s);

                    if (escaped(dc_lev)) {
                        uint32_t uval = (dc_lev >= 0) ? (uint32_t)(2 * dc_lev) : (uint32_t)(-2 * dc_lev - 1);
                        uint32_t overflow_val = uval - 255;
                        snt_bw_put_bits(&ov, dc_lev < 0 ? 1 : 0, 1);
                        snt_bw_put_expgolomb(&ov, overflow_val);
                    }

                    /* Reconstruct DC for prediction */
                    int16_t c_dc_recon = pvq_dequantize_dc(dc_lev, cstep) + c_dc_pred;
                    c_dc_left = c_dc_recon;
                    c_dc_above[cbxi] = c_dc_recon;

                    /* Record 15 chroma AC coefficients */
                    for (int i = 1; i < 16; i++) {
                        int16_t ac_lev = pvq_quantize_dc(ccoeffs[i], cstep, 1);
                        uint16_t ac_s = signed_to_sym(ac_lev);
                        RECORD_SYM(&cdf_chroma_ac, ac_s);

                        if (escaped(ac_lev)) {
                            uint32_t uval = (ac_lev >= 0) ? (uint32_t)(2 * ac_lev) : (uint32_t)(-2 * ac_lev - 1);
                            uint32_t overflow_val = uval - 255;
                            snt_bw_put_bits(&ov, ac_lev < 0 ? 1 : 0, 1);
                            snt_bw_put_expgolomb(&ov, overflow_val);
                        }
                    }
                }
            }
            free(c_dc_above);
        }
    }

    /* ════════════════════════════════════════════════════════════════
     *  PASS 2 (reverse): Feed all recorded ranges into rANS encoder.
     * ════════════════════════════════════════════════════════════════ */

    uint32_t rev_cap = range_count * 2 + 64; /* generous: each symbol may emit <=2 words */
    uint16_t *rev_buf = (uint16_t *)malloc(rev_cap * sizeof(uint16_t));
    if (!rev_buf) {
        free(ranges); free(overflow_buf); free(work); free(recon);
        free(dc_recon_above); free(ac_recon_above);
        return -1;
    }

    RansEncoder rans_enc;
    rans_enc_init(&rans_enc, rev_buf, rev_cap, ctx->header.dual_rans);

    /* Encode in reverse order (rANS is LIFO) */
    for (int32_t i = (int32_t)range_count - 1; i >= 0; i--) {
        rans_enc_put(&rans_enc, ranges[i].cdf_lo, ranges[i].freq, (uint32_t)i);
    }
    rans_enc_flush(&rans_enc);

    /* Get encoded rANS payload */
    uint32_t rans_words;
    const uint16_t *rans_data = rans_enc_data(&rans_enc, &rans_words);
    uint32_t rans_bytes = rans_words * 2;

    /* ════════════════════════════════════════════════════════════════
     *  WRITE to output bitstream:
     *    [block_count: 16 bits]
     *    [chroma_present: 1 bit]
     *    [rans_payload_length: 16 bits (in bytes)]
     *    [rans_payload: rans_bytes raw bytes]
     *    [overflow section: variable]
     * ════════════════════════════════════════════════════════════════ */

    snt_bw_byte_align(&bw);
    snt_bw_put_bits(&bw, (uint32_t)total_blocks, 16);
    snt_bw_put_bits(&bw, has_chroma ? 1 : 0, 1);
    snt_bw_byte_align(&bw);
    snt_bw_put_bits(&bw, rans_bytes, 16);

    /* Write rANS payload bytes (little-endian 16-bit words) */
    for (uint32_t i = 0; i < rans_words; i++) {
        snt_bw_put_byte(&bw, (uint8_t)(rans_data[i] & 0xFF));
        snt_bw_put_byte(&bw, (uint8_t)(rans_data[i] >> 8));
    }

    /* Write overflow section */
    snt_bw_byte_align(&ov);
    uint32_t overflow_len = snt_bw_bytes_written(&ov);
    for (uint32_t i = 0; i < overflow_len; i++) {
        snt_bw_put_byte(&bw, overflow_buf[i]);
    }

    snt_bw_byte_align(&bw);
    *bs_len = snt_bw_bytes_written(&bw);

    free(rev_buf);
    free(ranges);
    free(overflow_buf);
    free(work);
    free(recon);
    free(dc_recon_above);
    free(ac_recon_above);
    return 0;
}
