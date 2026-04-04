/*
 * Sinter Codec — Decoder with rANS entropy coding + inter prediction
 * Mirrors encoder symbol order exactly.
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
#include "../core/prediction/motion_est.h"
#include "../core/grain/grain_synth.h"
#include "../perceptual/csf_weights.h"

#define BLK 8

/* ── Zigzag unsigned-to-signed mapping ────────────────────────────── */
static inline int16_t u2s(uint16_t s) {
    return (s % 2 == 0) ? (int16_t)(s / 2) : -(int16_t)((s + 1) / 2);
}

/*
 * Symbol consumption order (must match encoder exactly):
 *
 * FOR each luma block (raster: by, bx):
 *   IF inter frame:
 *     DEC(cdf_inter) → inter_flag
 *     IF inter_flag:
 *       DEC(cdf_iskip) → skip_flag
 *       IF skip: copy ref block, continue
 *       IF coded: DEC(cdf_mvx), DEC(cdf_mvy) → MV delta
 *   DEC(cdf_dc) → dc_level (residual for inter, absolute for intra)
 *   FOR band b=0..15: same as before
 * FOR each chroma plane: same as before
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

    int is_inter = (ctx->header.frame_type == SNT_FRAME_INTER);
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

    /* ── SECTION C: Overflow reader ────────────────────────────────── */
    uint32_t ov_start = br.byte_pos;
    SntBitReader ovr;
    snt_br_init(&ovr, bitstream + ov_start, bs_len - ov_start);

    /* ── Init rANS decoder ─────────────────────────────────────────── */
    RansDecoder rdec;
    rans_dec_init(&rdec, rans_buf, rans_bytes, 1);

    /* ── Init 10 CDFs (must match encoder exactly) ─────────────────── */
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

    /* Cross-frame CDF carry: restore from previous frame for inter */
    if (is_inter && ctx->cdf_state_valid && ctx->cdf_state) {
        RansCDF *saved = (RansCDF *)ctx->cdf_state;
        cdf_dc = saved[0]; cdf_ac = saved[1]; cdf_mode = saved[2];
        cdf_skip = saved[3]; cdf_gain = saved[4]; cdf_k = saved[5];
        cdf_inter = saved[6]; cdf_iskip = saved[7]; cdf_mvx = saved[8]; cdf_mvy = saved[9];
    }

    uint16_t _dec_sym;
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

    /* MV field and block flags for chroma inter derivation */
    SntMV *mv_field = (SntMV *)calloc(blk_cols * blk_rows, sizeof(SntMV));
    int *blk_is_inter = (int *)calloc(blk_cols * blk_rows, sizeof(int));
    int *blk_is_skip = (int *)calloc(blk_cols * blk_rows, sizeof(int));

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

            /* ── Inter mode decode ────────────────────────────── */
            if (is_inter) {
                uint16_t inter_sym = DEC(cdf_inter);
                block_is_inter = (inter_sym == 1);
                blk_is_inter[blk_idx] = block_is_inter;

                if (block_is_inter) {
                    uint16_t skip_sym = DEC(cdf_iskip);
                    block_is_skip = (skip_sym == 1);
                    blk_is_skip[blk_idx] = block_is_skip;

                    if (block_is_skip) {
                        /* Copy reference block to recon */
                        for (int dy = 0; dy < BLK; dy++)
                            for (int dx = 0; dx < BLK; dx++) {
                                int py = py0 + dy, px_v = px0 + dx;
                                if (py < height && px_v < width && ctx->ref_y)
                                    recon[py * pad_w + px_v] =
                                        (int16_t)ctx->ref_y[py * width + px_v];
                            }
                        mv_field[blk_idx] = best_mv;
                        /* Reset prediction state for skip blocks —
                         * no DCT path, so no reconstructed DC/AC to store. */
                        dc_recon_left = 0;
                        dc_recon_above[bx] = 0;
                        for (int p = 0; p < NUM_PRED_AC; p++) {
                            ac_recon_left[p] = 0;
                            ac_recon_above[bx * NUM_PRED_AC + p] = 0;
                        }
                        continue;
                    }

                    /* MV prediction from spatial neighbors */
                    const SntMV *left_mv = (bx > 0) ? &mv_field[blk_idx - 1] : NULL;
                    const SntMV *above_mv = (by > 0) ? &mv_field[blk_idx - blk_cols] : NULL;
                    const SntMV *above_right_mv = (by > 0 && bx + 1 < blk_cols) ?
                        &mv_field[blk_idx - blk_cols + 1] : NULL;
                    SntMV mv_pred = snt_mv_predict(left_mv, above_mv, above_right_mv);

                    /* Decode MV deltas */
                    uint16_t mvx_sym = DEC(cdf_mvx);
                    int16_t mvd_x;
                    if (mvx_sym == 255)
                        mvd_x = (int16_t)snt_br_get_sexpgolomb(&ovr);
                    else
                        mvd_x = u2s(mvx_sym);

                    uint16_t mvy_sym = DEC(cdf_mvy);
                    int16_t mvd_y;
                    if (mvy_sym == 255)
                        mvd_y = (int16_t)snt_br_get_sexpgolomb(&ovr);
                    else
                        mvd_y = u2s(mvy_sym);

                    best_mv.x = mvd_x + mv_pred.x;
                    best_mv.y = mvd_y + mv_pred.y;
                    mv_field[blk_idx] = best_mv;
                }
            }

            /* ── Decode DC + bands (same for intra and inter coded) ── */
            int16_t recon_coeffs[64];

            /* DC prediction */
            int16_t dc_pred = 0;
            if (bx > 0 && by > 0)
                dc_pred = (dc_recon_left + dc_recon_above[bx] + 1) >> 1;
            else if (bx > 0)
                dc_pred = dc_recon_left;
            else if (by > 0)
                dc_pred = dc_recon_above[bx];

            uint16_t dc_sym = DEC(cdf_dc);
            int16_t dc_level;
            if (dc_sym == 255)
                dc_level = (int16_t)snt_br_get_sexpgolomb(&ovr);
            else
                dc_level = u2s(dc_sym);

            int16_t dc_recon_val = pvq_dequantize_dc(dc_level, step) + dc_pred;
            recon_coeffs[0] = dc_recon_val;
            dc_recon_left = dc_recon_val;
            dc_recon_above[bx] = dc_recon_val;

            /* Band processing */
            for (int b = 0; b < NUM_BANDS; b++) {
                int start = band_bounds[b];
                int n = band_bounds[b + 1] - start;

                uint16_t mode_sym = DEC(cdf_mode);

                if (mode_sym == 1) {
                    /* Scalar mode */
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
                    uint16_t skip_sym = DEC(cdf_skip);

                    if (!skip_sym) {
                        SntPVQCode code;
                        code.skip = 0;

                        uint16_t gain_sym = DEC(cdf_gain);
                        if (gain_sym == 255)
                            code.gain = snt_br_get_expgolomb(&ovr);
                        else
                            code.gain = (uint32_t)gain_sym;

                        uint16_t k_sym = DEC(cdf_k);
                        if (k_sym == 255)
                            code.k = (uint16_t)snt_br_get_expgolomb(&ovr);
                        else
                            code.k = k_sym;

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

            /* Add back AC prediction (key frames only) */
            if (!is_inter) {
                for (int p = 0; p < NUM_PRED_AC; p++) {
                    int idx = pred_ac_idx[p];
                    int16_t pred_val = 0;
                    if (idx < 8 && bx > 0)
                        pred_val = ac_recon_left[p];
                    else if (idx >= 8 && by > 0)
                        pred_val = ac_recon_above[bx * NUM_PRED_AC + p];
                    recon_coeffs[idx] += pred_val;
                }
            }

            /* Store reconstructed AC for future prediction (key frames only) */
            if (!is_inter) {
                for (int p = 0; p < NUM_PRED_AC; p++) {
                    int idx = pred_ac_idx[p];
                    ac_recon_left[p] = recon_coeffs[idx];
                    ac_recon_above[bx * NUM_PRED_AC + p] = recon_coeffs[idx];
                }
            }

            /* IDCT */
            int16_t recon_blk[64];
            snt_dct_inv_2d_ref(recon_coeffs, BLK, recon_blk, BLK);

            /* Write to recon: add MC prediction for inter coded blocks */
            for (int dy = 0; dy < BLK; dy++)
                for (int dx = 0; dx < BLK; dx++) {
                    int py = py0 + dy, px_v = px0 + dx;
                    int16_t val = recon_blk[dy * BLK + dx];
                    if (block_is_inter && !block_is_skip && ctx->ref_y) {
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

    /* Copy to output */
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
            output_y[y * stride_y + x] = snt_clip_pixel((int32_t)recon[y * pad_w + x]);

    /* Decode chroma with inter prediction support */
    if (chroma_present && output_u && output_v) {
        int cw = width / 2, ch_h = height / 2;
        int cblk = 4;
        int cpad_w = ((cw + cblk - 1) / cblk) * cblk;
        int cpad_h = ((ch_h + cblk - 1) / cblk) * cblk;
        int cblk_cols = cpad_w / cblk;
        /* Reduced chroma QP offset for inter frames (matches encoder) */
        int cqp_off = is_inter ? 1 : ctx->header.chroma_qp_offset;
        int cqp = ctx->header.base_qp + cqp_off;
        if (cqp < 0) cqp = 0; if (cqp > 63) cqp = 63;
        uint16_t cstep = snt_qp_to_step((uint8_t)cqp);

        uint8_t *chroma_out[2] = {output_u, output_v};
        const uint8_t *chroma_ref[2] = {ctx->ref_u, ctx->ref_v};
        for (int plane = 0; plane < 2; plane++) {
            const uint8_t *cr = chroma_ref[plane];
            int16_t c_dc_left = 0;
            int16_t *c_dc_above = (int16_t *)calloc(cblk_cols, sizeof(int16_t));

            for (int cby = 0; cby < cpad_h; cby += cblk) {
                c_dc_left = 0;
                for (int cbx = 0; cbx < cpad_w; cbx += cblk) {
                    int ci = cbx / cblk;
                    int cj = cby / cblk;

                    /* Check co-located luma block's inter/skip status */
                    int luma_blk_idx = (ci < blk_cols && cj < blk_rows) ?
                        cj * blk_cols + ci : 0;
                    int c_is_skip = is_inter && blk_is_skip[luma_blk_idx];
                    int c_is_inter = is_inter && blk_is_inter[luma_blk_idx];

                    /* Chroma skip: copy from reference */
                    if (c_is_skip && cr) {
                        SntMV luma_mv = mv_field[luma_blk_idx]; (void)luma_mv;
                        for (int dy = 0; dy < 4; dy++)
                            for (int dx = 0; dx < 4; dx++) {
                                int py = cby + dy, px = cbx + dx;
                                if (py < ch_h && px < cw)
                                    chroma_out[plane][py*stride_uv+px] = cr[py*cw+px];
                            }
                        c_dc_left = 0;
                        c_dc_above[ci] = 0;
                        continue;
                    }

                    int16_t c_dc_pred = 0;
                    if (cbx > 0 && cby > 0)
                        c_dc_pred = (c_dc_left + c_dc_above[ci] + 1) >> 1;
                    else if (cbx > 0) c_dc_pred = c_dc_left;
                    else if (cby > 0) c_dc_pred = c_dc_above[ci];

                    uint16_t cdc_sym = DEC(cdf_dc);
                    int16_t cdc_level;
                    if (cdc_sym == 255)
                        cdc_level = (int16_t)snt_br_get_sexpgolomb(&ovr);
                    else
                        cdc_level = u2s(cdc_sym);

                    int16_t cdc_recon = pvq_dequantize_dc(cdc_level, cstep) + c_dc_pred;
                    c_dc_left = cdc_recon;
                    c_dc_above[ci] = cdc_recon;

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

                    int16_t cblock[16];
                    snt_dct_inv_2d_ref(ccoeffs, 4, cblock, 4);

                    /* Output: add MC reference for inter, add 128 for intra */
                    for (int dy = 0; dy < 4; dy++)
                        for (int dx = 0; dx < 4; dx++) {
                            int py = cby + dy, px = cbx + dx;
                            if (py < ch_h && px < cw) {
                                int32_t val = (int32_t)cblock[dy*4+dx];
                                if (c_is_inter && cr) {
                                    SntMV luma_mv = mv_field[luma_blk_idx];
                                    int cmvx = luma_mv.x / 2, cmvy = luma_mv.y / 2;
                                    int ry = py + cmvy, rx = px + cmvx;
                                    int16_t ref_val = (ry >= 0 && ry < ch_h && rx >= 0 && rx < cw) ?
                                        (int16_t)cr[ry*cw+rx] : 128;
                                    val += ref_val;
                                } else {
                                    val += 128;
                                }
                                chroma_out[plane][py*stride_uv+px] = snt_clip_pixel(val);
                            }
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

    /* ── Save CDF state for cross-frame carry ─────────────────────── */
    if (!ctx->cdf_state)
        ctx->cdf_state = malloc(10 * sizeof(RansCDF));
    {
        RansCDF *saved = (RansCDF *)ctx->cdf_state;
        saved[0] = cdf_dc; saved[1] = cdf_ac; saved[2] = cdf_mode;
        saved[3] = cdf_skip; saved[4] = cdf_gain; saved[5] = cdf_k;
        saved[6] = cdf_inter; saved[7] = cdf_iskip; saved[8] = cdf_mvx; saved[9] = cdf_mvy;
    }
    ctx->cdf_state_valid = 1;

    /* ── Store reconstructed frame as reference ────────────────────── */
    int cw_dec = width / 2, ch_dec = height / 2;
    if (!ctx->ref_y || ctx->ref_w != width || ctx->ref_h != height) {
        free(ctx->ref_y); free(ctx->ref_u); free(ctx->ref_v);
        ctx->ref_y = (uint8_t *)malloc((size_t)width * height);
        ctx->ref_u = (uint8_t *)malloc((size_t)cw_dec * ch_dec);
        ctx->ref_v = (uint8_t *)malloc((size_t)cw_dec * ch_dec);
        ctx->ref_w = width;
        ctx->ref_h = height;
    }
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
            ctx->ref_y[y * width + x] = output_y[y * stride_y + x];
    if (output_u && output_v && ctx->ref_u && ctx->ref_v) {
        for (int y = 0; y < ch_dec; y++)
            for (int x = 0; x < cw_dec; x++) {
                ctx->ref_u[y * cw_dec + x] = output_u[y * stride_uv + x];
                ctx->ref_v[y * cw_dec + x] = output_v[y * stride_uv + x];
            }
    }
    ctx->num_refs = 1;

    free(recon);
    free(dc_recon_above);
    free(ac_recon_above);
    free(mv_field);
    free(blk_is_inter);
    free(blk_is_skip);
    free(rans_buf);
    return 0;
}
