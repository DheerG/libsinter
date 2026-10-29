/*
 * Sinter Codec — Decoder with rANS entropy coding + inter prediction
 * Mirrors encoder symbol order exactly.
 * No lapping post-filter: the encoder does not pre-filter (see encoder.c).
 */

#include <string.h>
#include <stdlib.h>
#include "decoder.h"
#include "cdf_set.h"
#include "scan.h"
#include "../util/fixed_point.h"
#include "../core/transform/dct.h"
#include "../core/quantize/pvq.h"
#include "../core/quantize/qp_table.h"
#include "../core/bitstream/bitstream.h"
#include "../core/entropy/rans.h"
#include "../core/prediction/motion_est.h"
#include "../core/prediction/mc.h"
#include "../core/filter/deblock.h"
#include "../core/prediction/intra_spatial.h"
#include "../core/grain/grain_synth.h"
#include "../perceptual/csf_weights.h"

#define BLK 8

/* ── Zigzag unsigned-to-signed mapping ────────────────────────────── */
static inline int16_t u2s(uint16_t s) {
    return (s % 2 == 0) ? (int16_t)(s / 2) : -(int16_t)((s + 1) / 2);
}

/* Decode one MV-delta component (escape: signed exp-Golomb overflow) */
static int16_t dec_mvd(RansDecoder *rd, RansCDF *cdf, SntBitReader *ovr) {
    uint16_t sym = rans_dec_symbol(rd, cdf);
    rans_cdf_adapt(cdf, sym);
    return sym == 255 ? (int16_t)snt_br_get_sexpgolomb(ovr) : u2s(sym);
}

/*
 * Symbol consumption order (must match encoder exactly):
 *
 * FOR each 16x16 macroblock (raster):
 *   IF inter frame AND full macroblock: DEC(mbmode); IF INTER16: DEC(mvx, mvy)
 *   FOR each 8x8 block in Z order (SKIP16/INTER16: no per-block mode symbols):
 *   IF inter frame (split macroblock):
 *     DEC(inter) → inter_flag
 *     IF inter_flag:
 *       DEC(iskip) → skip_flag
 *       IF skip: MC prediction at the predicted MV, continue
 *       IF coded: DEC(mvx), DEC(mvy) → MV delta
 *   IF intra block: DEC(imode[is_inter_frame]) → mode (SNT_IPRED_*)
 *   DEC(dc[t]) → dc_level (residual for inter, absolute for intra)
 *   FOR band b=0..15: mode[t][b], then ac[t][b/4] or skip[t][b] + gain[t]
 * FOR each chroma plane: cdc[ct], then cac[ct][snt_cac_group(i)]
 * (contexts: see cdf_set.h and the encoder's symbol-order comment)
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

    /* ── Entropy contexts: reset on key frames, carried across inter ── */
    if (!ctx->cdf_state)
        ctx->cdf_state = malloc(sizeof(SntCdfSet));
    SntCdfSet *cs = (SntCdfSet *)ctx->cdf_state;
    if (!is_inter || !ctx->cdf_state_valid)
        snt_cdf_set_init(cs);

    uint16_t _dec_sym;
    #define DEC(cdf) (_dec_sym = rans_dec_symbol(&rdec, &(cdf)), \
        rans_cdf_adapt(&(cdf), _dec_sym), _dec_sym)

    /* ── Decode luma blocks ────────────────────────────────────────── */
    int16_t *recon = (int16_t *)calloc((size_t)pad_w * pad_h, sizeof(int16_t));
    #define NUM_PRED_AC 4
    static const int pred_ac_idx[NUM_PRED_AC] = {1, 8, 2, 9};
    int total_blocks_d = blk_cols * blk_rows;
    /* Per-block state for neighbours' prediction: pixel-domain DC, low AC
     * (key frames), and whether the block is coded yet */
    int16_t *dc_pix = (int16_t *)calloc(total_blocks_d, sizeof(int16_t));
    int16_t *ac_state = (int16_t *)calloc((size_t)total_blocks_d * NUM_PRED_AC, sizeof(int16_t));
    uint8_t *coded = (uint8_t *)calloc(total_blocks_d, 1);

    /* MV field and block flags for chroma inter derivation */
    SntMV *mv_field = (SntMV *)calloc(blk_cols * blk_rows, sizeof(SntMV));
    int *blk_is_inter = (int *)calloc(blk_cols * blk_rows, sizeof(int));
    int *blk_is_skip = (int *)calloc(blk_cols * blk_rows, sizeof(int));

    #define NUM_BANDS 16
    static const int band_bounds[NUM_BANDS + 1] = {0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,59,63};

    /* Blocks are coded in 16x16 macroblocks in raster order, and the four
     * 8x8 blocks of a macroblock in Z order (top-left, top-right,
     * bottom-left, bottom-right). */
    int mb_cols = (blk_cols + 1) / 2, mb_rows = (blk_rows + 1) / 2;
    for (int mb = 0; mb < mb_cols * mb_rows; mb++) {
        /* Macroblock mode (inter frames, full macroblocks only). The 16x16
         * predictor comes from neighbours outside the macroblock. */
        int bx0 = (mb % mb_cols) * 2, by0 = (mb / mb_cols) * 2;
        int mb_mode = SNT_MB_SPLIT;
        SntMV mv_pred16 = {0, 0}, mv16 = {0, 0};
        if (is_inter && bx0 + 1 < blk_cols && by0 + 1 < blk_rows) {
            mv_pred16 = snt_mv_predict_at(mv_field, coded, blk_cols, bx0, by0);
            mb_mode = DEC(cs->mbmode);
            if (mb_mode == SNT_MB_INTER16) {
                mv16.x = dec_mvd(&rdec, &cs->mvx, &ovr) + mv_pred16.x;
                mv16.y = dec_mvd(&rdec, &cs->mvy, &ovr) + mv_pred16.y;
            }
        }
        for (int z = 0; z < 4; z++) {
            int bx = (mb % mb_cols) * 2 + (z & 1), by = (mb / mb_cols) * 2 + (z >> 1);
            if (bx >= blk_cols || by >= blk_rows) continue;
            int blk_idx = by * blk_cols + bx;
            int px0 = bx * BLK, py0 = by * BLK;
            int block_is_inter = 0;
            int block_is_skip = 0;
            SntMV best_mv = {0, 0};
            uint8_t mc_blk[BLK * BLK];  /* MC prediction (skip and inter coded blocks) */

            /* ── Inter mode decode ────────────────────────────── */
            if (is_inter) {
                if (mb_mode == SNT_MB_SPLIT) {
                    block_is_inter = DEC(cs->inter) == 1;
                    if (block_is_inter) {
                        /* MV prediction from spatial neighbors */
                        SntMV mv_pred = snt_mv_predict_at(mv_field, coded, blk_cols, bx, by);
                        block_is_skip = DEC(cs->iskip) == 1;
                        best_mv = mv_pred;  /* skip predicts at mv_pred */
                        if (!block_is_skip) {
                            best_mv.x = dec_mvd(&rdec, &cs->mvx, &ovr) + mv_pred.x;
                            best_mv.y = dec_mvd(&rdec, &cs->mvy, &ovr) + mv_pred.y;
                        }
                    }
                } else {
                    /* 16x16 skip or inter: every block takes the MB's MV */
                    block_is_inter = 1;
                    block_is_skip = (mb_mode == SNT_MB_SKIP16);
                    best_mv = block_is_skip ? mv_pred16 : mv16;
                }
                blk_is_inter[blk_idx] = block_is_inter;
                blk_is_skip[blk_idx] = block_is_skip;
                if (block_is_inter)
                    mv_field[blk_idx] = best_mv;

                if (block_is_skip) {
                    /* Copy the MC prediction to recon */
                    if (ctx->ref_y)
                        snt_mc_luma(ctx->ref_y, width, height, px0, py0, best_mv, BLK, BLK, mc_blk);
                    else
                        memset(mc_blk, 0, sizeof(mc_blk));  /* malformed: inter frame, no reference */
                    for (int dy = 0; dy < BLK; dy++)
                        for (int dx = 0; dx < BLK; dx++) {
                            int py = py0 + dy, px_v = px0 + dx;
                            if (py < height && px_v < width)
                                recon[py * pad_w + px_v] = mc_blk[dy * BLK + dx];
                        }
                    /* Neighbours predict their DC from this block's
                     * pixel-domain DC; AC prediction is key-frame only. */
                    dc_pix[blk_idx] = snt_pred_dc8(mc_blk);
                    memset(&ac_state[blk_idx * NUM_PRED_AC], 0, NUM_PRED_AC * sizeof(int16_t));
                    coded[blk_idx] = 1;
                    continue;
                }
            }

            if (block_is_inter && !block_is_skip) {
                if (ctx->ref_y)
                    snt_mc_luma(ctx->ref_y, width, height, px0, py0, best_mv, BLK, BLK, mc_blk);
                else
                    memset(mc_blk, 0, sizeof(mc_blk));  /* malformed: inter frame, no reference */
            }

            /* Intra mode; a spatial prediction goes in mc_blk */
            int imode = SNT_IPRED_TD;
            if (!block_is_inter) {
                imode = DEC(cs->imode[is_inter]);
                if (imode != SNT_IPRED_TD)
                    snt_intra_spatial_pred8(recon, pad_w, px0, py0, imode, mc_blk);
            }
            /* mc_blk holds this block's prediction (inter or spatial intra) */
            int has_pred_blk = (block_is_inter && !block_is_skip) || imode != SNT_IPRED_TD;
            int td_ac = !is_inter && imode == SNT_IPRED_TD;

            /* ── Decode DC + bands (same for intra and inter coded) ── */
            int16_t recon_coeffs[64];

            /* DC prediction: only transform-domain intra blocks predict DC
             * (from the neighbours' pixel-domain DC); inter and
             * spatial-intra residuals are zero-mean. */
            int16_t dc_pred = 0;
            if (has_pred_blk)
                dc_pred = 0;
            else if (bx > 0 && by > 0)
                dc_pred = (dc_pix[blk_idx - 1] + dc_pix[blk_idx - blk_cols] + 1) >> 1;
            else if (bx > 0)
                dc_pred = dc_pix[blk_idx - 1];
            else if (by > 0)
                dc_pred = dc_pix[blk_idx - blk_cols];

            uint16_t dc_sym = DEC(cs->dc[block_is_inter]);
            int16_t dc_level;
            if (dc_sym == 255)
                dc_level = (int16_t)snt_br_get_sexpgolomb(&ovr);
            else
                dc_level = u2s(dc_sym);

            int16_t dc_recon_val = pvq_dequantize_dc(dc_level, step) + dc_pred;
            recon_coeffs[0] = dc_recon_val;
            /* Neighbours see the pixel-domain DC (add the prediction block's DC) */
            dc_pix[blk_idx] = has_pred_blk ? dc_recon_val + snt_pred_dc8(mc_blk) : dc_recon_val;

            /* Band processing */
            for (int b = 0; b < NUM_BANDS; b++) {
                int start = band_bounds[b];
                int n = band_bounds[b + 1] - start;

                float csf_sum = 0;
                for (int i = 0; i < n; i++)
                    csf_sum += psy_csf_weight(8, snt_zigzag8[1 + start + i]);
                uint8_t band_csf = (uint8_t)(csf_sum / n * 255.0f);
                if (band_csf < 1) band_csf = 1;

                uint16_t mode_sym = DEC(cs->mode[block_is_inter][b]);

                if (mode_sym == 1) {
                    /* Scalar mode */
                    uint16_t eff_step = step;
                    if (band_csf < 128) {
                        eff_step = (uint16_t)((uint32_t)step * 256 / band_csf);
                        if (eff_step > 3072) eff_step = 3072;
                    }
                    for (int i = 0; i < n; i++) {
                        uint16_t ac_sym = DEC(cs->ac[block_is_inter][b / 4]);
                        int16_t level;
                        if (ac_sym == 255)
                            level = (int16_t)snt_br_get_sexpgolomb(&ovr);
                        else
                            level = u2s(ac_sym);
                        recon_coeffs[snt_zigzag8[1 + start + i]] = pvq_dequantize_dc(level, eff_step);
                    }
                } else {
                    /* PVQ mode */
                    uint16_t skip_sym = DEC(cs->skip[block_is_inter][b]);

                    if (!skip_sym) {
                        SntPVQCode code;
                        code.skip = 0;

                        uint16_t gain_sym = DEC(cs->gain[block_is_inter]);
                        if (gain_sym == 255)
                            code.gain = snt_br_get_expgolomb(&ovr);
                        else
                            code.gain = (uint32_t)gain_sym;

                        /* K is derived, not transmitted (same formula as encoder) */
                        code.k = pvq_compute_k(code.gain, n, band_csf);

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
                            recon_coeffs[snt_zigzag8[1 + start + i]] = band_recon[i];
                    } else {
                        for (int i = 0; i < n; i++)
                            recon_coeffs[snt_zigzag8[1 + start + i]] = 0;
                    }
                }
            }

            /* Add back AC prediction (transform-domain intra in key frames) */
            if (td_ac) {
                for (int p = 0; p < NUM_PRED_AC; p++) {
                    int idx = pred_ac_idx[p];
                    int16_t pred_val = 0;
                    if (idx < 8 && bx > 0)
                        pred_val = ac_state[(blk_idx - 1) * NUM_PRED_AC + p];
                    else if (idx >= 8 && by > 0)
                        pred_val = ac_state[(blk_idx - blk_cols) * NUM_PRED_AC + p];
                    recon_coeffs[idx] += pred_val;
                }
            }

            /* IDCT */
            int16_t recon_blk[64];
            snt_dct_inv_2d_ref(recon_coeffs, BLK, recon_blk, BLK);

            /* Write to recon: add the prediction block (inter or spatial intra) */
            for (int dy = 0; dy < BLK; dy++)
                for (int dx = 0; dx < BLK; dx++) {
                    int py = py0 + dy, px_v = px0 + dx;
                    int16_t val = recon_blk[dy * BLK + dx];
                    if (has_pred_blk)
                        val += mc_blk[dy * BLK + dx];
                    recon[py * pad_w + px_v] = val;
                    recon_blk[dy * BLK + dx] = val;
                }

            /* Key frames: store reconstructed low AC for the neighbours'
             * prediction (spatial-intra blocks: from the reconstructed block) */
            if (!is_inter) {
                const int16_t *src = recon_coeffs;
                int16_t full[64];
                if (imode != SNT_IPRED_TD) {
                    snt_dct_fwd_2d_ref(recon_blk, BLK, BLK, full);
                    src = full;
                }
                for (int p = 0; p < NUM_PRED_AC; p++)
                    ac_state[blk_idx * NUM_PRED_AC + p] = src[pred_ac_idx[p]];
            }
            coded[blk_idx] = 1;
        }
    }

    /* Copy to output */
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
            output_y[y * stride_y + x] = snt_clip_pixel((int32_t)recon[y * pad_w + x]);
    snt_deblock_luma(output_y, width, height, stride_y, step,
                     is_inter ? blk_is_skip : NULL, mv_field, blk_cols);

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

                    /* Chroma skip: MC prediction at the luma block's MV */
                    if (c_is_skip && cr) {
                        uint8_t c_skip[16];
                        snt_mc_chroma(cr, cw, ch_h, cbx, cby, mv_field[luma_blk_idx], 4, 4, c_skip);
                        for (int dy = 0; dy < 4; dy++)
                            for (int dx = 0; dx < 4; dx++) {
                                int py = cby + dy, px = cbx + dx;
                                if (py < ch_h && px < cw)
                                    chroma_out[plane][py*stride_uv+px] = c_skip[dy*4+dx];
                            }
                        c_dc_left = 0;
                        c_dc_above[ci] = 0;
                        continue;
                    }

                    int c_bt = (c_is_inter && cr) ? 1 : 0;  /* chroma context: inter-predicted? */
                    int16_t c_dc_pred = 0;
                    if (cbx > 0 && cby > 0)
                        c_dc_pred = (c_dc_left + c_dc_above[ci] + 1) >> 1;
                    else if (cbx > 0) c_dc_pred = c_dc_left;
                    else if (cby > 0) c_dc_pred = c_dc_above[ci];

                    uint16_t cdc_sym = DEC(cs->cdc[c_bt]);
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
                        uint16_t ac_sym = DEC(cs->cac[c_bt][snt_cac_group(i)]);
                        int16_t ac_level;
                        if (ac_sym == 255)
                            ac_level = (int16_t)snt_br_get_sexpgolomb(&ovr);
                        else
                            ac_level = u2s(ac_sym);
                        ccoeffs[i] = pvq_dequantize_dc(ac_level, cstep);
                    }

                    int16_t cblock[16];
                    snt_dct_inv_2d_ref(ccoeffs, 4, cblock, 4);

                    uint8_t c_mc[16];  /* chroma MC prediction: luma MV as 1/8-pel */
                    if (c_is_inter && cr)
                        snt_mc_chroma(cr, cw, ch_h, cbx, cby, mv_field[luma_blk_idx], 4, 4, c_mc);

                    /* Output: add MC reference for inter, add 128 for intra */
                    for (int dy = 0; dy < 4; dy++)
                        for (int dx = 0; dx < 4; dx++) {
                            int py = cby + dy, px = cbx + dx;
                            if (py < ch_h && px < cw) {
                                int32_t val = (int32_t)cblock[dy*4+dx];
                                if (c_is_inter && cr) {
                                    val += c_mc[dy*4+dx];
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
        snt_deblock_chroma(output_u, cw, ch_h, stride_uv, cstep, is_inter ? blk_is_skip : NULL, mv_field, blk_cols);
        snt_deblock_chroma(output_v, cw, ch_h, stride_uv, cstep, is_inter ? blk_is_skip : NULL, mv_field, blk_cols);
    }

    /* Grain synthesis */
    if (ctx->header.grain.enabled && output_u && output_v) {
        snt_grain_apply(&ctx->header.grain, output_y, output_u, output_v,
                        width, height, stride_y, stride_uv);
    }

    /* CDFs in ctx->cdf_state now carry into the next inter frame */
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
    free(dc_pix);
    free(ac_state);
    free(coded);
    free(mv_field);
    free(blk_is_inter);
    free(blk_is_skip);
    free(rans_buf);
    return 0;
}
