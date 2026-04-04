/*
 * Sinter Codec — Decoder with rANS entropy coding
 * Mirrors encoder symbol order exactly (spec v3).
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
#include "../core/entropy/rans.h"
#include "../core/grain/grain_synth.h"
#include "../perceptual/csf_weights.h"

#define BLK 8

/* ── Zigzag unsigned-to-signed mapping (Timothy) ──────────────────── */
static inline int16_t u2s(uint16_t s) {
    return (s % 2 == 0) ? (int16_t)(s / 2) : -(int16_t)((s + 1) / 2);
}

/*
 * Symbol consumption order (must match encoder exactly):
 *
 * FOR each luma block (raster: by, bx):
 *   DEC(cdf_dc) → dc_level
 *   FOR band b=0..15:
 *     DEC(cdf_mode) → mode
 *     IF scalar: FOR i=0..n-1: DEC(cdf_ac) → level
 *     IF PVQ:    DEC(cdf_skip) → skip
 *                IF !skip: DEC(cdf_gain) → gain, DEC(cdf_k) → k
 * FOR each chroma plane (U, V):
 *   FOR each 4x4 block:
 *     DEC(cdf_dc) → dc_level
 *     FOR i=1..15: DEC(cdf_ac) → level
 */

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

    /* ── SECTION A: Read fixed header ──────────────────────────────── */
    snt_br_byte_align(&br);
    uint32_t total_blocks = snt_br_get_bits(&br, 16);
    (void)total_blocks;
    uint32_t chroma_present = snt_br_get_bits(&br, 1);
    snt_br_byte_align(&br);

    /* ── SECTION B: Read rANS payload ──────────────────────────────── */
    uint32_t rans_bytes = snt_br_get_bits(&br, 32);
    uint8_t *rans_buf = (uint8_t *)malloc(rans_bytes + 4);
    for (uint32_t i = 0; i < rans_bytes; i++)
        rans_buf[i] = snt_br_get_byte(&br);

    /* ── SECTION C: Overflow reader starts at current position ─────── */
    uint32_t ov_start = br.byte_pos;
    SntBitReader ovr;
    snt_br_init(&ovr, bitstream + ov_start, bs_len - ov_start);

    /* ── Init rANS decoder ─────────────────────────────────────────── */
    RansDecoder rdec;
    rans_dec_init(&rdec, rans_buf, rans_bytes, 1); /* dual interleaved */

    /* ── Init 6 CDFs (must match encoder exactly) ──────────────────── */
    RansCDF cdf_dc, cdf_ac, cdf_mode, cdf_skip, cdf_gain, cdf_k;
    rans_cdf_init_uniform(&cdf_dc,   256, RANS_ADAPT_RATE_COEFF);
    rans_cdf_init_uniform(&cdf_ac,   256, RANS_ADAPT_RATE_COEFF);
    rans_cdf_init_uniform(&cdf_mode,   2, RANS_ADAPT_RATE_MODE);
    rans_cdf_init_uniform(&cdf_skip,   2, RANS_ADAPT_RATE_MODE);
    rans_cdf_init_uniform(&cdf_gain, 256, RANS_ADAPT_RATE_COEFF);
    rans_cdf_init_uniform(&cdf_k,    256, RANS_ADAPT_RATE_COEFF);

    /* DEC helper: decode symbol + adapt CDF (Monty — mirrors REC exactly) */
    uint16_t _dec_sym; /* shared temp for DEC macro */
    #define DEC(cdf) (_dec_sym = rans_dec_symbol(&rdec, &(cdf)), \
        rans_cdf_adapt(&(cdf), _dec_sym), _dec_sym)

    /* ── Decode luma blocks ────────────────────────────────────────── */
    int16_t *recon = (int16_t *)calloc((size_t)pad_w * pad_h, sizeof(int16_t));
    int16_t *dc_recon_above = (int16_t *)calloc(blk_cols, sizeof(int16_t));
    int16_t dc_recon_left = 0;
    #define NUM_PRED_AC 4
    static const int pred_ac_idx[NUM_PRED_AC] = {1, 8, 2, 9};
    int16_t *ac_recon_above = (int16_t *)calloc((size_t)blk_cols * NUM_PRED_AC, sizeof(int16_t));
    int16_t ac_recon_left[NUM_PRED_AC] = {0};

    for (int by = 0; by < blk_rows; by++) {
        dc_recon_left = 0;
        for (int bx = 0; bx < blk_cols; bx++) {
            int16_t recon_coeffs[64];

            /* DC prediction */
            int16_t dc_pred = 0;
            if (bx > 0 && by > 0)
                dc_pred = (dc_recon_left + dc_recon_above[bx] + 1) >> 1;
            else if (bx > 0)
                dc_pred = dc_recon_left;
            else if (by > 0)
                dc_pred = dc_recon_above[bx];

            /* DEC: DC residual ← cdf_dc */
            uint16_t dc_sym = DEC(cdf_dc);
            int16_t dc_level;
            if (dc_sym == 255)
                dc_level = (int16_t)snt_br_get_sexpgolomb(&ovr);
            else
                dc_level = u2s(dc_sym);

            int16_t dc_recon = pvq_dequantize_dc(dc_level, step) + dc_pred;
            recon_coeffs[0] = dc_recon;
            dc_recon_left = dc_recon;
            dc_recon_above[bx] = dc_recon;

            /* Band processing */
            #define NUM_BANDS 16
            static const int band_bounds[NUM_BANDS + 1] = {0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,59,63};

            for (int b = 0; b < NUM_BANDS; b++) {
                int start = band_bounds[b];
                int n = band_bounds[b + 1] - start;

                /* DEC: mode ← cdf_mode */
                uint16_t mode_sym = DEC(cdf_mode);

                if (mode_sym == 1) {
                    /* Scalar mode — CSF-weighted step (must match encoder) */
                    float csf_sum = 0;
                    for (int i = 0; i < n; i++)
                        csf_sum += psy_csf_weight(8, 1 + start + i);
                    uint8_t band_csf = (uint8_t)(csf_sum / n * 255.0f);
                    if (band_csf < 1) band_csf = 1;
                    uint16_t eff_step = step;
                    if (band_csf < 128) {
                        eff_step = (uint16_t)((uint32_t)step * 256 / band_csf);
                        if (eff_step > 3072) eff_step = 3072;
                    }
                    for (int i = 0; i < n; i++) {
                        /* DEC: scalar level ← cdf_ac */
                        uint16_t ac_sym = DEC(cdf_ac);
                        int16_t level;
                        if (ac_sym == 255)
                            level = (int16_t)snt_br_get_sexpgolomb(&ovr);
                        else
                            level = u2s(ac_sym);
                        recon_coeffs[1 + start + i] = pvq_dequantize_dc(level, eff_step);
                    }
                } else {
                    /* PVQ mode */
                    /* DEC: skip ← cdf_skip */
                    uint16_t skip_sym = DEC(cdf_skip);

                    if (!skip_sym) {
                        SntPVQCode code;
                        code.skip = 0;

                        /* DEC: gain ← cdf_gain */
                        uint16_t gain_sym = DEC(cdf_gain);
                        if (gain_sym == 255)
                            code.gain = snt_br_get_expgolomb(&ovr);
                        else
                            code.gain = (uint32_t)gain_sym;

                        /* DEC: k ← cdf_k */
                        uint16_t k_sym = DEC(cdf_k);
                        if (k_sym == 255)
                            code.k = (uint16_t)snt_br_get_expgolomb(&ovr);
                        else
                            code.k = k_sym;

                        /* Shape index from overflow (raw bits) */
                        int shape_bits = pvq_shape_bits(n, code.k);
                        if (shape_bits < 1) shape_bits = 1;
                        if (shape_bits > 64) shape_bits = 64;
                        if (shape_bits > 32) {
                            uint64_t hi = snt_br_get_bits(&ovr, shape_bits - 32);
                            uint64_t lo = snt_br_get_bits(&ovr, 32);
                            code.shape_idx = (hi << 32) | lo;
                        } else {
                            code.shape_idx = snt_br_get_bits(&ovr, shape_bits);
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
    if (chroma_present && output_u && output_v) {
        int cw = width / 2, ch_h = height / 2;
        int cblk = 4;
        int cpad_w = ((cw + cblk - 1) / cblk) * cblk;
        int cpad_h = ((ch_h + cblk - 1) / cblk) * cblk;
        int cblk_cols = cpad_w / cblk;
        int cqp = ctx->header.base_qp + ctx->header.chroma_qp_offset;
        if (cqp < 0) cqp = 0; if (cqp > 63) cqp = 63;
        uint16_t cstep = snt_qp_to_step((uint8_t)cqp);

        uint8_t *chroma_out[2] = {output_u, output_v};
        for (int plane = 0; plane < 2; plane++) {
            int16_t c_dc_left = 0;
            int16_t *c_dc_above = (int16_t *)calloc(cblk_cols, sizeof(int16_t));

            for (int cby = 0; cby < cpad_h; cby += cblk) {
                c_dc_left = 0;
                for (int cbx = 0; cbx < cpad_w; cbx += cblk) {
                    int ci = cbx / cblk;

                    /* DC prediction */
                    int16_t c_dc_pred = 0;
                    if (cbx > 0 && cby > 0)
                        c_dc_pred = (c_dc_left + c_dc_above[ci] + 1) >> 1;
                    else if (cbx > 0) c_dc_pred = c_dc_left;
                    else if (cby > 0) c_dc_pred = c_dc_above[ci];

                    /* DEC: chroma DC ← cdf_dc */
                    uint16_t cdc_sym = DEC(cdf_dc);
                    int16_t cdc_level;
                    if (cdc_sym == 255)
                        cdc_level = (int16_t)snt_br_get_sexpgolomb(&ovr);
                    else
                        cdc_level = u2s(cdc_sym);

                    int16_t cdc_recon = pvq_dequantize_dc(cdc_level, cstep) + c_dc_pred;
                    c_dc_left = cdc_recon;
                    c_dc_above[ci] = cdc_recon;

                    /* DEC: chroma AC ← cdf_ac */
                    int16_t ccoeffs[16];
                    ccoeffs[0] = cdc_recon;
                    for (int i = 1; i < 16; i++) {
                        uint16_t ac_sym = DEC(cdf_ac);
                        int16_t ac_level;
                        if (ac_sym == 255)
                            ac_level = (int16_t)snt_br_get_sexpgolomb(&ovr);
                        else
                            ac_level = u2s(ac_sym);
                        ccoeffs[i] = pvq_dequantize_dc(ac_level, cstep);
                    }

                    /* Inverse 4x4 DCT */
                    int16_t cblock[16];
                    snt_dct_inv_2d_ref(ccoeffs, 4, cblock, 4);

                    for (int dy = 0; dy < 4; dy++)
                        for (int dx = 0; dx < 4; dx++) {
                            int py = cby + dy, px = cbx + dx;
                            if (py < ch_h && px < cw)
                                chroma_out[plane][py*stride_uv+px] =
                                    snt_clip_pixel((int32_t)cblock[dy*4+dx] + 128);
                        }
                }
            }
            free(c_dc_above);
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
    free(rans_buf);
    return 0;
}
