/*
 * Sinter Codec — rANS Decoder (matches rANS encoder)
 * Two-pass approach: decode all rANS symbols first, then read raw overflow bits.
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
#include "../core/entropy/rans.h"

#define BLK 8
#define NUM_BANDS 16
#define NUM_PRED_AC 4

/* Inverse of encoder's signed_to_sym: unsigned sym back to signed value.
 * sym=0 -> 0, sym=1 -> -1, sym=2 -> +1, sym=3 -> -2, sym=4 -> +2, ...
 * sym=255 is the escape code: actual value comes from raw overflow. */
static int16_t sym_to_signed(uint16_t sym) {
    return (sym % 2 == 0) ? (int16_t)(sym / 2) : -(int16_t)((sym + 1) / 2);
}

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

    static const int band_bounds[NUM_BANDS + 1] = {
        0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,59,63
    };
    static const int pred_ac_idx[NUM_PRED_AC] = {1, 8, 2, 9};

    /* ── Read bitstream header (must match encoder write order) ────── */
    snt_br_byte_align(&br);
    uint32_t total_blocks = snt_br_get_bits(&br, 16);
    uint32_t chroma_present = snt_br_get_bits(&br, 1);
    snt_br_byte_align(&br);
    uint32_t rans_payload_len = snt_br_get_bits(&br, 16); /* rANS payload in bytes */

    /* rANS payload starts at current byte position */
    uint32_t rans_start_byte = br.byte_pos;
    const uint8_t *rans_buf = bitstream + rans_start_byte;

    /* ── Initialize CDFs identically to encoder ────────────────────── */
    RansCDF cdf_dc, cdf_low_ac, cdf_mode, cdf_skip, cdf_gain, cdf_k, cdf_scalar_ac;
    rans_cdf_init_uniform(&cdf_dc,        256, RANS_ADAPT_RATE_COEFF);
    rans_cdf_init_uniform(&cdf_low_ac,    256, RANS_ADAPT_RATE_COEFF);
    rans_cdf_init_uniform(&cdf_mode,        2, RANS_ADAPT_RATE_MODE);
    rans_cdf_init_uniform(&cdf_skip,        2, RANS_ADAPT_RATE_MODE);
    rans_cdf_init_uniform(&cdf_gain,      256, RANS_ADAPT_RATE_COEFF);
    rans_cdf_init_uniform(&cdf_k,         256, RANS_ADAPT_RATE_COEFF);
    rans_cdf_init_uniform(&cdf_scalar_ac, 256, RANS_ADAPT_RATE_COEFF);

    RansCDF cdf_chroma_dc, cdf_chroma_ac;
    rans_cdf_init_uniform(&cdf_chroma_dc, 256, RANS_ADAPT_RATE_COEFF);
    rans_cdf_init_uniform(&cdf_chroma_ac, 256, RANS_ADAPT_RATE_COEFF);

    /* ── Initialize rANS decoder from payload ──────────────────────── */
    RansDecoder rdec;
    rans_dec_init(&rdec, rans_buf, rans_payload_len, ctx->header.dual_rans);

    /* ── Allocate symbol storage ──────────────────────────────────── */
    /* rans_syms stores DC, ac_pred, and scalar AC symbols (in order).
     * PVQ symbols (mode, skip, gain, k) are stored in separate arrays
     * indexed by [blk_idx * NUM_BANDS + band]. */
    uint32_t max_rans_syms = total_blocks * 128 + 4096;
    uint16_t *rans_syms = (uint16_t *)malloc(max_rans_syms * sizeof(uint16_t));
    uint32_t rans_sym_count = 0;

    uint8_t  *block_band_modes = (uint8_t *)malloc((size_t)total_blocks * NUM_BANDS);
    uint8_t  *pvq_skip_flags   = (uint8_t *)malloc((size_t)total_blocks * NUM_BANDS);
    uint16_t *pvq_gain_syms    = (uint16_t *)malloc((size_t)total_blocks * NUM_BANDS * sizeof(uint16_t));
    uint16_t *pvq_k_syms       = (uint16_t *)malloc((size_t)total_blocks * NUM_BANDS * sizeof(uint16_t));

    /* ── PASS 1: Decode ALL rANS symbols in forward order ─────────── */
    for (uint32_t blk_idx = 0; blk_idx < total_blocks; blk_idx++) {
        /* dc_sym (cdf_dc, 256 symbols) */
        uint16_t dc_sym = rans_dec_symbol(&rdec, &cdf_dc);
        rans_cdf_adapt(&cdf_dc, dc_sym);
        rans_syms[rans_sym_count++] = dc_sym;

        /* 4 x ac_pred_sym (cdf_low_ac, 256 symbols) */
        for (int p = 0; p < NUM_PRED_AC; p++) {
            uint16_t ac_sym = rans_dec_symbol(&rdec, &cdf_low_ac);
            rans_cdf_adapt(&cdf_low_ac, ac_sym);
            rans_syms[rans_sym_count++] = ac_sym;
        }

        /* 16 bands */
        for (int b = 0; b < NUM_BANDS; b++) {
            int n = band_bounds[b + 1] - band_bounds[b];

            /* mode_sym (cdf_mode, 2 symbols): 0=PVQ, 1=scalar */
            uint16_t mode_sym = rans_dec_symbol(&rdec, &cdf_mode);
            rans_cdf_adapt(&cdf_mode, mode_sym);
            block_band_modes[blk_idx * NUM_BANDS + b] = (uint8_t)mode_sym;

            if (mode_sym == 0) {
                /* PVQ mode */
                uint16_t skip_sym = rans_dec_symbol(&rdec, &cdf_skip);
                rans_cdf_adapt(&cdf_skip, skip_sym);
                pvq_skip_flags[blk_idx * NUM_BANDS + b] = (uint8_t)skip_sym;

                if (!skip_sym) {
                    uint16_t gain_lo = rans_dec_symbol(&rdec, &cdf_gain);
                    rans_cdf_adapt(&cdf_gain, gain_lo);
                    pvq_gain_syms[blk_idx * NUM_BANDS + b] = gain_lo;

                    uint16_t k_lo = rans_dec_symbol(&rdec, &cdf_k);
                    rans_cdf_adapt(&cdf_k, k_lo);
                    pvq_k_syms[blk_idx * NUM_BANDS + b] = k_lo;
                }
            } else {
                /* Scalar mode: n x ac_level_sym (cdf_scalar_ac, 256 symbols) */
                for (int i = 0; i < n; i++) {
                    uint16_t ac_sym = rans_dec_symbol(&rdec, &cdf_scalar_ac);
                    rans_cdf_adapt(&cdf_scalar_ac, ac_sym);
                    rans_syms[rans_sym_count++] = ac_sym;
                }
            }
        }
    }

    /* Chroma rANS symbols */
    int cw = width / 2, ch = height / 2;
    int cblk = 4;
    int cpad_w = ((cw + cblk - 1) / cblk) * cblk;
    int cpad_h = ((ch + cblk - 1) / cblk) * cblk;
    int cblk_cols = cpad_w / cblk;
    int chroma_blocks_per_plane = (cpad_w / cblk) * (cpad_h / cblk);

    uint32_t max_chroma_syms = (uint32_t)chroma_blocks_per_plane * 2 * 16 + 256;
    uint16_t *chroma_syms = (uint16_t *)malloc(max_chroma_syms * sizeof(uint16_t));
    uint32_t chroma_sym_count = 0;

    for (int plane = 0; plane < 2; plane++) {
        for (int bi = 0; bi < chroma_blocks_per_plane; bi++) {
            uint16_t dc_sym = rans_dec_symbol(&rdec, &cdf_chroma_dc);
            rans_cdf_adapt(&cdf_chroma_dc, dc_sym);
            chroma_syms[chroma_sym_count++] = dc_sym;

            for (int i = 0; i < 15; i++) {
                uint16_t ac_sym = rans_dec_symbol(&rdec, &cdf_chroma_ac);
                rans_cdf_adapt(&cdf_chroma_ac, ac_sym);
                chroma_syms[chroma_sym_count++] = ac_sym;
            }
        }
    }

    /* ── PASS 2: Read raw overflow section ────────────────────────── */
    /* Advance bitreader past rANS payload to the raw overflow section */
    br.byte_pos = rans_start_byte + rans_payload_len;
    br.bit_pos = 0;

    /* overflow_vals stores: escaped DC/AC/scalar values (sexpgolomb),
     * AND PVQ gain/k actual values (always stored for non-skip PVQ bands).
     * All in the order the encoder wrote them. */
    uint32_t max_overflow = rans_sym_count + chroma_sym_count +
                            total_blocks * NUM_BANDS * 2 + 256;
    int32_t  *overflow_vals  = (int32_t *)malloc(max_overflow * sizeof(int32_t));
    uint64_t *shape_indices  = (uint64_t *)malloc((size_t)total_blocks * NUM_BANDS * sizeof(uint64_t));
    uint32_t ov_count = 0;

    /* Walk rans_syms in the same order as pass 1 to identify escapes */
    uint32_t si = 0; /* index into rans_syms */

    for (uint32_t blk_idx = 0; blk_idx < total_blocks; blk_idx++) {
        /* dc overflow */
        if (rans_syms[si] == 255)
            overflow_vals[ov_count++] = snt_br_get_sexpgolomb(&br);
        si++;

        /* 4 x ac_pred overflow */
        for (int p = 0; p < NUM_PRED_AC; p++) {
            if (rans_syms[si] == 255)
                overflow_vals[ov_count++] = snt_br_get_sexpgolomb(&br);
            si++;
        }

        /* 16 bands */
        for (int b = 0; b < NUM_BANDS; b++) {
            int n = band_bounds[b + 1] - band_bounds[b];
            uint8_t mode = block_band_modes[blk_idx * NUM_BANDS + b];

            if (mode == 0) {
                /* PVQ mode */
                uint8_t skip = pvq_skip_flags[blk_idx * NUM_BANDS + b];
                if (!skip) {
                    uint16_t gain_sym = pvq_gain_syms[blk_idx * NUM_BANDS + b];
                    uint16_t k_sym    = pvq_k_syms[blk_idx * NUM_BANDS + b];

                    /* Resolve gain: use sym directly, or read overflow */
                    uint32_t gain_val = (gain_sym == 255)
                        ? snt_br_get_expgolomb(&br) : (uint32_t)gain_sym;
                    overflow_vals[ov_count++] = (int32_t)gain_val;

                    /* Resolve k: use sym directly, or read overflow */
                    uint16_t k_val = (k_sym == 255)
                        ? (uint16_t)snt_br_get_expgolomb(&br) : k_sym;
                    overflow_vals[ov_count++] = (int32_t)k_val;

                    /* Read shape index (raw bits) */
                    int shape_bits = pvq_shape_bits(n, k_val);
                    if (shape_bits > 64) shape_bits = 64;
                    uint64_t shape_idx;
                    if (shape_bits > 32) {
                        uint64_t hi = snt_br_get_bits(&br, shape_bits - 32);
                        uint64_t lo = snt_br_get_bits(&br, 32);
                        shape_idx = (hi << 32) | lo;
                    } else if (shape_bits > 0) {
                        shape_idx = snt_br_get_bits(&br, shape_bits);
                    } else {
                        shape_idx = 0;
                    }
                    shape_indices[blk_idx * NUM_BANDS + b] = shape_idx;
                }
            } else {
                /* Scalar mode: read overflows for escaped symbols */
                for (int i = 0; i < n; i++) {
                    if (rans_syms[si + i] == 255)
                        overflow_vals[ov_count++] = snt_br_get_sexpgolomb(&br);
                }
                si += n;
            }
        }
    }

    /* Chroma overflows */
    for (uint32_t ci = 0; ci < chroma_sym_count; ci++) {
        if (chroma_syms[ci] == 255)
            overflow_vals[ov_count++] = snt_br_get_sexpgolomb(&br);
    }

    /* ── PASS 3: Reconstruct luma blocks ──────────────────────────── */
    int16_t *recon = (int16_t *)calloc((size_t)pad_w * pad_h, sizeof(int16_t));
    int16_t *dc_recon_above = (int16_t *)calloc(blk_cols, sizeof(int16_t));
    int16_t dc_recon_left = 0;
    int16_t *ac_recon_above = (int16_t *)calloc((size_t)blk_cols * NUM_PRED_AC, sizeof(int16_t));
    int16_t ac_recon_left[NUM_PRED_AC] = {0};

    /* Re-walk rans_syms and overflow_vals in the same order as pass 2 */
    si = 0;
    uint32_t ov_idx = 0;

    for (int by = 0; by < blk_rows; by++) {
        dc_recon_left = 0;
        for (int bx = 0; bx < blk_cols; bx++) {
            uint32_t blk_idx = (uint32_t)(by * blk_cols + bx);
            int16_t recon_coeffs[64];

            /* DC prediction (mirrors encoder exactly) */
            int16_t dc_pred = 0;
            if (bx > 0 && by > 0)
                dc_pred = (dc_recon_left + dc_recon_above[bx] + 1) >> 1;
            else if (bx > 0)
                dc_pred = dc_recon_left;
            else if (by > 0)
                dc_pred = dc_recon_above[bx];

            /* DC: sym -> signed level -> dequant + prediction */
            uint16_t dc_sym = rans_syms[si++];
            int16_t dc_level = (dc_sym == 255)
                ? (int16_t)overflow_vals[ov_idx++] : sym_to_signed(dc_sym);
            int16_t dc_recon = pvq_dequantize_dc(dc_level, step) + dc_pred;
            recon_coeffs[0] = dc_recon;
            dc_recon_left = dc_recon;
            dc_recon_above[bx] = dc_recon;

            /* 4 AC prediction residual levels */
            int16_t ac_pred_levels[NUM_PRED_AC];
            for (int p = 0; p < NUM_PRED_AC; p++) {
                uint16_t ac_sym = rans_syms[si++];
                ac_pred_levels[p] = (ac_sym == 255)
                    ? (int16_t)overflow_vals[ov_idx++] : sym_to_signed(ac_sym);
            }

            /* Decode band data */
            for (int b = 0; b < NUM_BANDS; b++) {
                int start = band_bounds[b];
                int n = band_bounds[b + 1] - start;
                uint8_t mode = block_band_modes[blk_idx * NUM_BANDS + b];

                if (mode == 1) {
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
                        uint16_t ac_sym = rans_syms[si++];
                        int16_t level = (ac_sym == 255)
                            ? (int16_t)overflow_vals[ov_idx++] : sym_to_signed(ac_sym);
                        recon_coeffs[1 + start + i] = pvq_dequantize_dc(level, eff_step);
                    }
                } else {
                    /* PVQ mode */
                    uint8_t skip = pvq_skip_flags[blk_idx * NUM_BANDS + b];
                    if (!skip) {
                        SntPVQCode code;
                        code.skip = 0;
                        code.gain = (uint32_t)overflow_vals[ov_idx++];
                        code.k    = (uint16_t)overflow_vals[ov_idx++];
                        code.shape_idx = shape_indices[blk_idx * NUM_BANDS + b];

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
                recon_coeffs[idx] = pvq_dequantize_dc(ac_pred_levels[p], step) + pred_val;
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

    /* Copy luma to output */
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
            output_y[y * stride_y + x] = snt_clip_pixel((int32_t)recon[y * pad_w + x]);

    /* ── Decode chroma using pre-decoded rANS symbols ─────────────── */
    if (output_u && output_v) {
        int cqp = ctx->header.base_qp + ctx->header.chroma_qp_offset;
        if (cqp < 0) cqp = 0; if (cqp > 63) cqp = 63;
        uint16_t cstep = snt_qp_to_step((uint8_t)cqp);

        uint8_t *chroma_out[2] = {output_u, output_v};
        uint32_t c_si = 0;

        for (int plane = 0; plane < 2; plane++) {
            int16_t c_dc_left = 0;
            int16_t *c_dc_above = (int16_t *)calloc(cblk_cols, sizeof(int16_t));

            for (int cby = 0; cby < cpad_h; cby += cblk) {
                c_dc_left = 0;
                for (int cbx = 0; cbx < cpad_w; cbx += cblk) {
                    int cbxi = cbx / cblk;

                    /* DC prediction */
                    int16_t c_dc_pred = 0;
                    if (cbx > 0 && cby > 0)
                        c_dc_pred = (c_dc_left + c_dc_above[cbxi] + 1) >> 1;
                    else if (cbx > 0) c_dc_pred = c_dc_left;
                    else if (cby > 0) c_dc_pred = c_dc_above[cbxi];

                    /* DC from rANS sym */
                    uint16_t dc_sym = chroma_syms[c_si++];
                    int16_t dc_level = (dc_sym == 255)
                        ? (int16_t)overflow_vals[ov_idx++] : sym_to_signed(dc_sym);
                    int16_t c_dc_recon = pvq_dequantize_dc(dc_level, cstep) + c_dc_pred;
                    c_dc_left = c_dc_recon;
                    c_dc_above[cbxi] = c_dc_recon;

                    /* 15 AC coefficients from rANS syms */
                    int16_t ccoeffs[16];
                    ccoeffs[0] = c_dc_recon;
                    for (int i = 1; i < 16; i++) {
                        uint16_t ac_sym = chroma_syms[c_si++];
                        int16_t ac_level = (ac_sym == 255)
                            ? (int16_t)overflow_vals[ov_idx++] : sym_to_signed(ac_sym);
                        ccoeffs[i] = pvq_dequantize_dc(ac_level, cstep);
                    }

                    /* Inverse 4x4 DCT */
                    int16_t cblock[16];
                    snt_dct_inv_2d_ref(ccoeffs, 4, cblock, 4);

                    /* Write to output (+128 to uncenter) */
                    for (int dy = 0; dy < 4; dy++)
                        for (int dx = 0; dx < 4; dx++) {
                            int py = cby + dy, px = cbx + dx;
                            if (py < ch && px < cw)
                                chroma_out[plane][py * stride_uv + px] =
                                    snt_clip_pixel((int32_t)cblock[dy * 4 + dx] + 128);
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

    /* Cleanup */
    free(recon);
    free(dc_recon_above);
    free(ac_recon_above);
    free(rans_syms);
    free(chroma_syms);
    free(overflow_vals);
    free(shape_indices);
    free(block_band_modes);
    free(pvq_gain_syms);
    free(pvq_k_syms);
    free(pvq_skip_flags);
    return 0;
}
