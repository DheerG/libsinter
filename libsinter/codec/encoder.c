/*
 * Sinter Codec — Encoder with rANS entropy coding + inter prediction
 * Two-pass: forward record + reverse rANS encode.
 * Inter spec: team consensus (Jim/Timothy/Gary/Monty/Jens-Rainer).
 *
 * Key frames: 8x8 block DCT with DC/low-AC prediction
 * Inter frames: MC prediction + residual coding
 *
 * Lapping (core/transform/lapping.c) is not applied. Its lifting network
 * has no butterfly/scaling stage, so it is not DC-preserving: a flat 128
 * image comes out with a +/-11 sawtooth at every block edge, which the DCT
 * then has to code. Disabling it saved 25-41% of intra bytes and gained
 * 2.5-3.3 dB PSNR. Re-enable only with a proper TDLT (e.g. Daala's
 * od_pre_filter); PSNR does not measure the blocking this lets through.
 */

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "encoder.h"
#include "cdf_set.h"
#include "scan.h"
#include "../core/transform/dct.h"
#include "../core/quantize/pvq.h"
#include "../core/quantize/qp_table.h"
#include "../core/bitstream/bitstream.h"
#include "../core/entropy/rans.h"
#include "../core/prediction/motion_est.h"
#include "../core/prediction/mc.h"
#include "../core/filter/deblock.h"
#include "../core/prediction/intra_spatial.h"
#include "../util/fixed_point.h"
#include "../perceptual/activity_mask.h"
#include "../perceptual/csf_weights.h"

#define BLK 8

/* ── Zigzag signed-to-unsigned mapping ────────────────────────────── */
static inline uint16_t s2u(int16_t v) {
    uint32_t u = (v >= 0) ? (uint32_t)(2 * v) : (uint32_t)(-2 * v - 1);
    return (u > 254) ? 255 : (uint16_t)u;
}

/* ── Rate estimates for encoder decisions (non-normative) ──────────── */
#define RD_LAMBDA 0.14  /* lambda = RD_LAMBDA * step^2 (orthonormal DCT, SSE) */

/* Key frames are coded this many QP steps finer than inter frames: every
 * frame until the next key frame predicts from them, so their quality pays
 * off repeatedly. The frame header carries the QP actually used. */
#define KEY_QP_OFFSET 4

/* Inter-frame intra decision: a block is coded intra only if its SAD from
 * a flat block is under 1/INTRA_SAD_BIAS of the inter SAD + MV cost. The
 * flat-block SAD makes intra look cheaper than it codes; 2 was best by
 * BD-rate. */
#define INTRA_SAD_BIAS 2

/* Ideal cost in bits of coding sym with cdf in its current state. */
static double cdf_bits(const RansCDF *cdf, uint16_t sym) {
    uint16_t lo, freq;
    rans_cdf_symbol_range(cdf, sym, &lo, &freq);
    return -log2((double)freq / RANS_PROB_TOTAL);
}

/* Length of the order-0 exp-Golomb code for v. */
static int expgolomb_len(uint32_t v) {
    int prefix = 0;
    for (uint32_t c = v + 1; c > 1; c >>= 1) prefix++;
    return 2 * prefix + 1;
}

/* Bits for a signed level: rANS symbol, plus overflow bits on escape. */
static double level_bits(const RansCDF *cdf, int16_t level) {
    uint16_t sym = s2u(level);
    double bits = cdf_bits(cdf, sym);
    if (sym == 255)
        bits += expgolomb_len(level >= 0 ? 2u * level : (uint32_t)(-2 * level - 1));
    return bits;
}

static double band_sse(const int16_t *a, const int16_t *b, int n) {
    double sse = 0;
    for (int i = 0; i < n; i++) {
        double d = (double)a[i] - b[i];
        sse += d * d;
    }
    return sse;
}

/* ── Symbol recording ──────────────────────────────────────────────── */
typedef struct { uint16_t cdf_lo; uint16_t freq; } SymRange;

/* Forward-pass symbol sink: records each symbol's rANS range, adapts its
 * CDF, and holds the overflow bit writer for escapes and PVQ shapes. */
typedef struct {
    SymRange     *ranges;
    uint32_t      rc;
    SntBitWriter *ov;
} SymSink;

static void rec(SymSink *s, RansCDF *cdf, uint16_t sym) {
    rans_cdf_symbol_range(cdf, sym, &s->ranges[s->rc].cdf_lo, &s->ranges[s->rc].freq);
    rans_cdf_adapt(cdf, sym);
    s->rc++;
}

/* ── Luma block coding ─────────────────────────────────────────────── */
#define NUM_BANDS 16
static const int band_bounds[NUM_BANDS + 1] = {0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,59,63};
#define NUM_PRED_AC 4
static const int pred_ac_idx[NUM_PRED_AC] = {1, 8, 2, 9};

/* Result of coding one 8x8 luma block's DC and bands */
typedef struct {
    int16_t recon_coeffs[64];  /* reconstructed coefficients (TD prediction included) */
    int16_t recon[64];         /* reconstructed pixels (prediction block added) */
    int16_t dc_recon;          /* reconstructed DC coefficient */
    double  bits;              /* estimated bits (estimate mode only) */
} LumaCode;

/* Code one 8x8 luma block's DC and bands.
 * dct_in: pixels (transform-domain intra) or a residual (all other modes).
 * td_pred: transform-domain prediction, used when pred is NULL: the DC,
 *   plus the low AC when td_ac. pred: prediction block added back on
 *   reconstruction (NULL for transform-domain intra).
 * t: context block type (0 intra, 1 inter).
 * With a sink, symbols are emitted and CDFs adapt as they go. Without
 * one, nothing is emitted, the CDFs are untouched, and out->bits
 * estimates the cost from their current state. */
static void code_luma(SntCdfSet *cs, SymSink *sink, uint16_t step, int t,
                      const int16_t *dct_in, const int16_t *td_pred, int td_ac,
                      const uint8_t *pred, LumaCode *out) {
    int16_t coeffs[64];
    snt_dct_fwd_2d_ref(dct_in, BLK, BLK, coeffs);
    double bits = 0;
    double lambda = RD_LAMBDA * (double)step * step;

    /* DC: transform-domain intra predicts it; other residuals are zero-mean */
    int16_t dc_pred = pred ? 0 : td_pred[0];
    int16_t dc_level = pvq_quantize_dc(coeffs[0] - dc_pred, step, 1);
    if (sink) {
        uint16_t dc_sym = s2u(dc_level);
        rec(sink, &cs->dc[t], dc_sym);
        if (dc_sym == 255)
            snt_bw_put_sexpgolomb(sink->ov, dc_level);
    } else {
        bits += level_bits(&cs->dc[t], dc_level);
    }
    out->dc_recon = pvq_dequantize_dc(dc_level, step) + dc_pred;

    /* Low-AC prediction (transform-domain intra in key frames) */
    if (td_ac)
        for (int p = 0; p < NUM_PRED_AC; p++)
            coeffs[pred_ac_idx[p]] -= td_pred[pred_ac_idx[p]];

    int16_t recon_ac[63];
    memset(recon_ac, 0, sizeof(recon_ac));
    for (int b = 0; b < NUM_BANDS; b++) {
        int start = band_bounds[b];
        int n = band_bounds[b + 1] - start;

        float csf_sum = 0.0f;
        for (int i = 0; i < n; i++)
            csf_sum += psy_csf_weight(BLK, snt_zigzag8[1 + start + i]);
        uint8_t band_csf = (uint8_t)(csf_sum / n * 255.0f);
        if (band_csf < 1) band_csf = 1;

        int16_t band_coeffs[16];  /* this band, in zigzag order */
        for (int i = 0; i < n; i++)
            band_coeffs[i] = coeffs[snt_zigzag8[1 + start + i]];

        /* Scalar candidate (CSF-weighted step on high bands) */
        uint16_t eff_step = step;
        if (band_csf < 128) {
            eff_step = (uint16_t)((uint32_t)step * 256 / band_csf);
            if (eff_step > 3072) eff_step = 3072;
        }
        int16_t s_level[16], s_recon[16];
        double s_bits = cdf_bits(&cs->mode[t][b], 1);
        for (int i = 0; i < n; i++) {
            s_level[i] = pvq_quantize_dc(band_coeffs[i], eff_step, 1);
            s_recon[i] = pvq_dequantize_dc(s_level[i], eff_step);
            s_bits += level_bits(&cs->ac[t][b / 4], s_level[i]);
        }

        /* PVQ candidate */
        SntPVQCode code;
        int16_t p_recon[16];
        pvq_encode(band_coeffs, n, step, 1, band_csf, &code, p_recon);
        uint16_t gain_sym = (code.gain > 255) ? 255 : (uint16_t)code.gain;
        int shape_bits = 0;
        double p_bits = cdf_bits(&cs->mode[t][b], 0)
                      + cdf_bits(&cs->skip[t][b], code.skip ? 1 : 0);
        if (!code.skip) {
            /* K = pvq_compute_k(gain, n, band_csf): the decoder
             * has all three inputs, so K is not transmitted. */
            shape_bits = pvq_shape_bits(n, code.k);
            if (shape_bits < 1) shape_bits = 1;
            if (shape_bits > 64) shape_bits = 64;
            p_bits += cdf_bits(&cs->gain[t], gain_sym) + shape_bits;
            if (gain_sym == 255) p_bits += expgolomb_len(code.gain);
        }

        /* Mode decision: minimum distortion + lambda * bits.
         * (The old energy-concentration heuristic picked PVQ for
         * bands where it cost far more than it saved.) */
        int use_scalar = band_sse(band_coeffs, s_recon, n) + lambda * s_bits
                      <= band_sse(band_coeffs, p_recon, n) + lambda * p_bits;
        const int16_t *band_recon = use_scalar ? s_recon : p_recon;

        if (!sink) {
            bits += use_scalar ? s_bits : p_bits;
        } else {
            rec(sink, &cs->mode[t][b], use_scalar ? 1 : 0);
            if (use_scalar) {
                for (int i = 0; i < n; i++) {
                    uint16_t ac_sym = s2u(s_level[i]);
                    rec(sink, &cs->ac[t][b / 4], ac_sym);
                    if (ac_sym == 255)
                        snt_bw_put_sexpgolomb(sink->ov, s_level[i]);
                }
            } else {
                rec(sink, &cs->skip[t][b], code.skip ? 1 : 0);
                if (!code.skip) {
                    rec(sink, &cs->gain[t], gain_sym);
                    if (gain_sym == 255)
                        snt_bw_put_expgolomb(sink->ov, code.gain);
                    if (shape_bits > 32) {
                        snt_bw_put_bits(sink->ov, (uint32_t)(code.shape_idx >> 32), shape_bits - 32);
                        snt_bw_put_bits(sink->ov, (uint32_t)(code.shape_idx & 0xFFFFFFFF), 32);
                    } else {
                        snt_bw_put_bits(sink->ov, (uint32_t)code.shape_idx, shape_bits);
                    }
                }
            }
        }
        for (int i = 0; i < n; i++)
            recon_ac[start + i] = band_recon[i];
    }

    /* Reconstruct coefficients, then pixels */
    out->recon_coeffs[0] = out->dc_recon;
    for (int i = 0; i < 63; i++)
        out->recon_coeffs[snt_zigzag8[i + 1]] = recon_ac[i];
    if (td_ac)
        for (int p = 0; p < NUM_PRED_AC; p++)
            out->recon_coeffs[pred_ac_idx[p]] += td_pred[pred_ac_idx[p]];
    snt_dct_inv_2d_ref(out->recon_coeffs, BLK, out->recon, BLK);
    if (pred)
        for (int k = 0; k < 64; k++)
            out->recon[k] += pred[k];
    out->bits = bits;
}

/*
 * Symbol emission order (normative, must match decoder exactly):
 *
 * Contexts (cdf_set.h): t = block type (0 intra, 1 inter-coded),
 * ct = chroma block type (1 if inter-predicted).
 *
 * FOR each luma block (raster: by, bx):
 *   IF inter frame:
 *     REC(inter, inter_flag)        // 0=intra, 1=inter
 *     IF inter_flag:
 *       REC(iskip, skip_flag)       // 0=coded, 1=skip
 *       IF skip: (no more symbols, MC prediction at the predicted MV)
 *       IF coded:
 *         REC(mvx, s2u(mvd_x))
 *         REC(mvy, s2u(mvd_y))
 *   IF intra block: REC(imode[is_inter_frame], mode)  // SNT_IPRED_*
 *   REC(dc[t], s2u(dc_level))       // DC of residual (inter) or pixels (intra)
 *   // Low-AC prediction subtracted before band processing (key frames only)
 *   FOR band b=0..15:
 *     REC(mode[t][b], mode)
 *     IF scalar: FOR i=0..n-1: REC(ac[t][b/4], s2u(level))
 *     IF PVQ:    REC(skip[t][b], skip)
 *                IF !skip: REC(gain[t], min(gain,255))
 *                          (K is not sent: decoder derives it with pvq_compute_k)
 * FOR each chroma plane (U, V):
 *   FOR each 4x4 block:
 *     REC(cdc[ct], s2u(dc_level))
 *     FOR i=1..15: REC(cac[ct][snt_cac_group(i)], s2u(level))
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
    SntFrameHeader fhdr = ctx->header;  /* this frame's header (QP may differ) */
    if (!is_inter)
        fhdr.base_qp = fhdr.base_qp > KEY_QP_OFFSET ? fhdr.base_qp - KEY_QP_OFFSET : 0;
    snt_write_frame_header(&bw, &fhdr);

    uint16_t step = snt_qp_to_step(fhdr.base_qp);

    /* Load pixels into work buffer */
    int16_t *work = (int16_t *)calloc((size_t)pad_w * pad_h, sizeof(int16_t));
    int16_t *recon = (int16_t *)calloc((size_t)pad_w * pad_h, sizeof(int16_t));
    if (!work || !recon) { free(work); free(recon); return -1; }

    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
            work[y * pad_w + x] = (int16_t)input_y[y * stride_y + x];

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

    uint32_t ov_cap = (uint32_t)(total_blocks * 128 + 4096);
    uint8_t *ov_buf = (uint8_t *)malloc(ov_cap);
    SntBitWriter ov;
    snt_bw_init(&ov, ov_buf, ov_cap);

    /* ── Entropy contexts: reset on key frames, carried across inter ── */
    if (!ctx->cdf_state)
        ctx->cdf_state = malloc(sizeof(SntCdfSet));
    SntCdfSet *cs = (SntCdfSet *)ctx->cdf_state;
    if (!is_inter || !ctx->cdf_state_valid)
        snt_cdf_set_init(cs);

    SymSink sink = { ranges, 0, &ov };
    #define REC(cdf, sym) rec(&sink, &(cdf), (sym))

    /* ── PHASE 1: Forward pass — record symbols + overflow ─────── */
    /* Per-block mode flags for chroma inter derivation */
    int *blk_is_inter = (int *)calloc(total_blocks, sizeof(int));
    int *blk_is_skip = (int *)calloc(total_blocks, sizeof(int));

    int16_t *dc_recon_above = (int16_t *)calloc(blk_cols, sizeof(int16_t));
    int16_t dc_recon_left = 0;
    int16_t *ac_recon_above = (int16_t *)calloc((size_t)blk_cols * NUM_PRED_AC, sizeof(int16_t));
    int16_t ac_recon_left[NUM_PRED_AC] = {0};

    /* MV field for spatial prediction (inter only) */
    SntMV *mv_field = (SntMV *)calloc(blk_cols * blk_rows, sizeof(SntMV));

    for (int by = 0; by < blk_rows; by++) {
        dc_recon_left = 0;
        for (int bx = 0; bx < blk_cols; bx++) {
            int blk_idx = by * blk_cols + bx;
            int px0 = bx * BLK, py0 = by * BLK;
            int block_is_inter = 0;
            int block_is_skip = 0;
            SntMV best_mv = {0, 0};
            uint8_t mc_blk[BLK * BLK];  /* MC prediction (skip and inter coded blocks) */

            /* ── Inter mode decision ──────────────────────────── */
            if (is_inter) {
                /* MV prediction from spatial neighbors (also used by the
                 * search, which charges each candidate for its MV delta) */
                const SntMV *left_mv = (bx > 0) ? &mv_field[blk_idx - 1] : NULL;
                const SntMV *above_mv = (by > 0) ? &mv_field[blk_idx - blk_cols] : NULL;
                const SntMV *above_right_mv = (by > 0 && bx + 1 < blk_cols) ?
                    &mv_field[blk_idx - blk_cols + 1] : NULL;
                SntMV mv_pred = snt_mv_predict(left_mv, above_mv, above_right_mv);

                /* SAD for skip, which predicts at mv_pred (so a block moving
                 * with its neighbours costs only the two flags).
                 * Note: input accessed at stride width (not stride_y) to match
                 * ref buffer stride. Correct when stride_y == width (CLI). */
                snt_mc_luma(ctx->ref_y, width, height, px0, py0, mv_pred, BLK, BLK, mc_blk);
                uint32_t skip_sad = 0;
                for (int dy = 0; dy < BLK; dy++)
                    for (int dx = 0; dx < BLK; dx++) {
                        int py = py0 + dy, px_v = px0 + dx;
                        if (py < height && px_v < width) {
                            int diff = (int)input_y[py * width + px_v]
                                     - (int)mc_blk[dy * BLK + dx];
                            skip_sad += (uint32_t)(diff < 0 ? -diff : diff);
                        }
                    }

                /* lambda for SAD-domain decisions = sqrt(mode lambda) */
                double lambda_sad = sqrt(RD_LAMBDA) * step;

                /* Motion search: cur uses stride_y, ref uses width.
                 * Both must use the same stride for snte_motion_search.
                 * We pass width — correct when stride_y == width (CLI). */
                SntMV search_mv = snte_motion_search(
                    input_y, ctx->ref_y, width,
                    px0, py0, BLK, 16, width, height, mv_pred, lambda_sad);
                /* Refine to quarter-pel; inter_sad is the SAD at the result */
                uint32_t inter_sad;
                SntMV qmv = snte_subpel_refine(input_y, ctx->ref_y, width, height,
                                               px0, py0, BLK, search_mv, mv_pred,
                                               lambda_sad, &inter_sad);

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
                double mv_cost = lambda_sad * (snte_mvd_bits(qmv.x - mv_pred.x) +
                                               snte_mvd_bits(qmv.y - mv_pred.y));

                if (skip_sad <= skip_threshold) {
                    block_is_inter = 1;
                    block_is_skip = 1;
                    best_mv = mv_pred;
                } else if (inter_sad + mv_cost < INTRA_SAD_BIAS * (double)intra_sad) {
                    block_is_inter = 1;
                    block_is_skip = 0;
                    best_mv = qmv;
                } else {
                    block_is_inter = 0;
                }

                /* Store flags for chroma inter derivation */
                blk_is_inter[blk_idx] = block_is_inter;
                blk_is_skip[blk_idx] = block_is_skip;

                /* REC: inter_flag */
                REC(cs->inter, block_is_inter ? 1 : 0);

                if (block_is_inter) {
                    /* REC: skip_flag */
                    REC(cs->iskip, block_is_skip ? 1 : 0);

                    if (block_is_skip) {
                        /* Copy the MC prediction to recon */
                        snt_mc_luma(ctx->ref_y, width, height, px0, py0, best_mv, BLK, BLK, mc_blk);
                        for (int dy = 0; dy < BLK; dy++)
                            for (int dx = 0; dx < BLK; dx++) {
                                int py = py0 + dy, px_v = px0 + dx;
                                if (py < height && px_v < width)
                                    recon[py * pad_w + px_v] = mc_blk[dy * BLK + dx];
                            }
                        mv_field[blk_idx] = best_mv;
                        /* Neighbours predict their DC from this block's
                         * pixel-domain DC; AC prediction is key-frame only. */
                        dc_recon_left = snt_pred_dc8(mc_blk);
                        dc_recon_above[bx] = dc_recon_left;
                        for (int p = 0; p < NUM_PRED_AC; p++) {
                            ac_recon_left[p] = 0;
                            ac_recon_above[bx * NUM_PRED_AC + p] = 0;
                        }
                        continue; /* skip to next block */
                    }

                    /* MV delta (quarter-pel) */
                    int16_t mvd_x = best_mv.x - mv_pred.x;
                    int16_t mvd_y = best_mv.y - mv_pred.y;

                    /* REC: MV deltas */
                    uint16_t mvx_sym = s2u(mvd_x);
                    REC(cs->mvx, mvx_sym);
                    if (mvx_sym == 255)
                        snt_bw_put_sexpgolomb(&ov, mvd_x);

                    uint16_t mvy_sym = s2u(mvd_y);
                    REC(cs->mvy, mvy_sym);
                    if (mvy_sym == 255)
                        snt_bw_put_sexpgolomb(&ov, mvd_y);

                    mv_field[blk_idx] = best_mv;
                }
            }

            /* ── DCT input: residual (inter) or pixels (intra) ──── */
            int16_t dct_input[64];
            if (block_is_inter && !block_is_skip) {
                snt_mc_luma(ctx->ref_y, width, height, px0, py0, best_mv, BLK, BLK, mc_blk);
                /* Inter coded: DCT of residual = cur - MC prediction
                 * (zero outside the frame) */
                for (int dy = 0; dy < BLK; dy++)
                    for (int dx = 0; dx < BLK; dx++) {
                        int py = py0 + dy, px_v = px0 + dx;
                        dct_input[dy * BLK + dx] = (py < height && px_v < width) ?
                            work[py * pad_w + px_v] - mc_blk[dy * BLK + dx] : 0;
                    }
            } else {
                /* Intra (key or inter-frame intra): use work buffer directly */
                for (int dy = 0; dy < BLK; dy++)
                    for (int dx = 0; dx < BLK; dx++)
                        dct_input[dy * BLK + dx] = work[(py0 + dy) * pad_w + (px0 + dx)];
            }

            int16_t coeffs[64];
            snt_dct_fwd_2d_ref(dct_input, BLK, BLK, coeffs);

            /* Transform-domain prediction for intra blocks: DC from the
             * neighbours' pixel-domain DC; low AC (key frames) from their
             * reconstructed coefficients. */
            int16_t td_pred[64];
            memset(td_pred, 0, sizeof(td_pred));
            if (!block_is_inter) {
                if (bx > 0 && by > 0)
                    td_pred[0] = (dc_recon_left + dc_recon_above[bx] + 1) >> 1;
                else if (bx > 0)
                    td_pred[0] = dc_recon_left;
                else if (by > 0)
                    td_pred[0] = dc_recon_above[bx];
                if (!is_inter)
                    for (int p = 0; p < NUM_PRED_AC; p++) {
                        int idx = pred_ac_idx[p];
                        if (idx < 8 && bx > 0)
                            td_pred[idx] = ac_recon_left[p];
                        else if (idx >= 8 && by > 0)
                            td_pred[idx] = ac_recon_above[bx * NUM_PRED_AC + p];
                    }
            }

            /* Intra mode: transform-domain prediction or a spatial
             * predictor, whichever has the lower transform-domain cost
             * (sum of |residual coefficients| + lambda * mode bits).
             * A spatial prediction goes in mc_blk, like an MC prediction;
             * luma_in becomes its residual. */
            int16_t luma_in[64];
            memcpy(luma_in, dct_input, sizeof(luma_in));
            int imode = SNT_IPRED_TD;
            if (!block_is_inter) {
                double lambda_sad = sqrt(RD_LAMBDA) * step;
                double best_cost = lambda_sad * cdf_bits(&cs->imode[is_inter], SNT_IPRED_TD);
                for (int k = 0; k < 64; k++)
                    best_cost += abs(coeffs[k] - td_pred[k]);
                for (int m = SNT_IPRED_V; m < SNT_NUM_IPRED; m++) {
                    uint8_t p[BLK * BLK];
                    int16_t res[64], rcf[64];
                    snt_intra_spatial_pred8(recon, pad_w, px0, py0, m, p);
                    for (int k = 0; k < 64; k++)
                        res[k] = dct_input[k] - p[k];
                    snt_dct_fwd_2d_ref(res, BLK, BLK, rcf);
                    double cost = lambda_sad * cdf_bits(&cs->imode[is_inter], (uint16_t)m);
                    for (int k = 0; k < 64; k++)
                        cost += abs(rcf[k]);
                    if (cost < best_cost) {
                        best_cost = cost;
                        imode = m;
                        memcpy(mc_blk, p, sizeof(p));
                        memcpy(luma_in, res, sizeof(res));
                    }
                }
                REC(cs->imode[is_inter], (uint16_t)imode);
            }
            /* mc_blk holds this block's prediction (inter or spatial intra) */
            int has_pred_blk = (block_is_inter && !block_is_skip) || imode != SNT_IPRED_TD;
            int td_ac = !is_inter && imode == SNT_IPRED_TD;

            LumaCode lc;
            code_luma(cs, &sink, step, block_is_inter, luma_in, td_pred, td_ac,
                      has_pred_blk ? mc_blk : NULL, &lc);

            /* Neighbours see the pixel-domain DC (add the prediction block's DC) */
            dc_recon_left = has_pred_blk ? lc.dc_recon + snt_pred_dc8(mc_blk) : lc.dc_recon;
            dc_recon_above[bx] = dc_recon_left;

            for (int dy = 0; dy < BLK; dy++)
                for (int dx = 0; dx < BLK; dx++)
                    recon[(py0 + dy) * pad_w + px0 + dx] = lc.recon[dy * BLK + dx];

            /* Key frames: store reconstructed low AC for the neighbours'
             * prediction. Spatial-intra blocks coded a residual, so take
             * the coefficients of the reconstructed block itself. */
            if (!is_inter) {
                const int16_t *src = lc.recon_coeffs;
                int16_t full[64];
                if (imode != SNT_IPRED_TD) {
                    snt_dct_fwd_2d_ref(lc.recon, BLK, BLK, full);
                    src = full;
                }
                for (int p = 0; p < NUM_PRED_AC; p++) {
                    ac_recon_left[p] = src[pred_ac_idx[p]];
                    ac_recon_above[bx * NUM_PRED_AC + p] = src[pred_ac_idx[p]];
                }
            }
        }
        memset(ac_recon_left, 0, sizeof(ac_recon_left));
    }

    /* Ensure a chroma reference exists to predict from (zero-filled on the first frame) */
    if (input_u && input_v) {
        if (!ctx->ref_u || ctx->ref_w != width || ctx->ref_h != height) {
            free(ctx->ref_u); free(ctx->ref_v);
            ctx->ref_u = (uint8_t *)calloc((size_t)cw * ch, 1);
            ctx->ref_v = (uint8_t *)calloc((size_t)cw * ch, 1);
        }
    }

    /* Chroma: 4x4 DCT + scalar quantization, with inter prediction.
     * Reconstruct into fresh buffers and swap them in at the end. Writing
     * into ctx->ref_u/v while predicting from them lets MVs that point
     * up/left read current-frame pixels the decoder never sees, which was
     * the real cause of the chroma "drift". */
    if (input_u && input_v) {
        uint8_t *chroma_rec[2] = {(uint8_t *)calloc((size_t)cw * ch, 1),
                                  (uint8_t *)calloc((size_t)cw * ch, 1)};
        int cblk_cols = cpad_w / cblk;
        /* Reduce chroma QP offset for inter frames to prevent quantization
         * cascade (drift). +4 for key frames (no drift), +1 for inter. */
        int cqp_off = is_inter ? 1 : ctx->header.chroma_qp_offset;
        int cqp = fhdr.base_qp + cqp_off;
        if (cqp < 0) cqp = 0; if (cqp > 63) cqp = 63;
        uint16_t cstep = snt_qp_to_step((uint8_t)cqp);

        const uint8_t *chroma_in[2] = {input_u, input_v};
        const uint8_t *chroma_ref[2] = {ctx->ref_u, ctx->ref_v};
        for (int plane = 0; plane < 2; plane++) {
            const uint8_t *cp = chroma_in[plane];
            const uint8_t *cr = chroma_ref[plane]; /* may be NULL on first frame */
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

                    /* Chroma skip (implicit from luma skip): MC prediction at
                     * the luma block's MV, which is its predicted MV */
                    if (c_is_skip && cr) {
                        uint8_t *c_ref_out = chroma_rec[plane];
                        uint8_t c_skip[16];
                        snt_mc_chroma(cr, cw, ch, cbx, cby, mv_field[luma_blk_idx], 4, 4, c_skip);
                        for (int dy = 0; dy < 4; dy++)
                            for (int dx = 0; dx < 4; dx++) {
                                int py = cby + dy, px = cbx + dx;
                                if (py < ch && px < cw && c_ref_out)
                                    c_ref_out[py*cw+px] = c_skip[dy*4+dx];
                            }
                        c_dc_left = 0;
                        c_dc_above[ci] = 0;
                        continue;
                    }

                    /* Build 4x4 block input */
                    int16_t cblock[16];
                    uint8_t c_mc[16];  /* chroma MC prediction (inter only) */
                    if (c_is_inter && cr) {
                        /* Chroma inter: the luma MV (quarter-pel) is the
                         * chroma MV in 1/8-pel units */
                        snt_mc_chroma(cr, cw, ch, cbx, cby, mv_field[luma_blk_idx], 4, 4, c_mc);
                        for (int dy = 0; dy < 4; dy++)
                            for (int dx = 0; dx < 4; dx++) {
                                int py = cby + dy, px = cbx + dx;
                                cblock[dy*4+dx] = (py < ch && px < cw) ?
                                    (int16_t)cp[py*stride_uv+px] - c_mc[dy*4+dx] : 0;
                            }
                    } else {
                        /* Intra: center by subtracting 128 */
                        for (int dy = 0; dy < 4; dy++)
                            for (int dx = 0; dx < 4; dx++) {
                                int py = cby + dy, px = cbx + dx;
                                cblock[dy*4+dx] = (py < ch && px < cw) ?
                                    (int16_t)cp[py*stride_uv+px] - 128 : 0;
                            }
                    }

                    int c_bt = (c_is_inter && cr) ? 1 : 0;  /* chroma context: inter-predicted? */
                    int16_t ccoeffs[16];
                    snt_dct_fwd_2d_ref(cblock, 4, 4, ccoeffs);

                    int16_t c_dc_pred = 0;
                    if (cbx > 0 && cby > 0)
                        c_dc_pred = (c_dc_left + c_dc_above[ci] + 1) >> 1;
                    else if (cbx > 0) c_dc_pred = c_dc_left;
                    else if (cby > 0) c_dc_pred = c_dc_above[ci];

                    int16_t dc_res = ccoeffs[0] - c_dc_pred;
                    int16_t cdc_level = pvq_quantize_dc(dc_res, cstep, 1);

                    uint16_t cdc_sym = s2u(cdc_level);
                    REC(cs->cdc[c_bt], cdc_sym);
                    if (cdc_sym == 255)
                        snt_bw_put_sexpgolomb(&ov, cdc_level);

                    int16_t cdc_recon = pvq_dequantize_dc(cdc_level, cstep) + c_dc_pred;
                    c_dc_left = cdc_recon;
                    c_dc_above[ci] = cdc_recon;

                    int16_t crecon_coeffs[16];
                    crecon_coeffs[0] = cdc_recon;
                    for (int i = 1; i < 16; i++) {
                        int16_t ac_level = pvq_quantize_dc(ccoeffs[i], cstep, 1);
                        uint16_t ac_sym = s2u(ac_level);
                        REC(cs->cac[c_bt][snt_cac_group(i)], ac_sym);
                        if (ac_sym == 255)
                            snt_bw_put_sexpgolomb(&ov, ac_level);
                        crecon_coeffs[i] = pvq_dequantize_dc(ac_level, cstep);
                    }

                    /* Reconstruct chroma for reference frame */
                    int16_t crecon_blk[16];
                    snt_dct_inv_2d_ref(crecon_coeffs, 4, crecon_blk, 4);
                    uint8_t *c_ref_out = chroma_rec[plane];
                    if (c_ref_out) {
                        for (int dy = 0; dy < 4; dy++)
                            for (int dx = 0; dx < 4; dx++) {
                                int py = cby + dy, px = cbx + dx;
                                if (py < ch && px < cw) {
                                    int32_t val = (int32_t)crecon_blk[dy*4+dx];
                                    if (c_is_inter && cr) {
                                        val += c_mc[dy*4+dx];
                                    } else {
                                        val += 128;
                                    }
                                    c_ref_out[py*cw+px] = snt_clip_pixel(val);
                                }
                            }
                    }
                }
            }
            free(c_dc_above);
        }
        free(ctx->ref_u);
        free(ctx->ref_v);
        ctx->ref_u = chroma_rec[0];
        ctx->ref_v = chroma_rec[1];
        snt_deblock_chroma(ctx->ref_u, cw, ch, cw, cstep, is_inter ? blk_is_skip : NULL, mv_field, blk_cols);
        snt_deblock_chroma(ctx->ref_v, cw, ch, cw, cstep, is_inter ? blk_is_skip : NULL, mv_field, blk_cols);
    }

    /* ── PHASE 2: Reverse rANS encode ──────────────────────────────── */
    uint32_t rev_cap = sink.rc * 2 + 64;
    uint16_t *rev_buf = (uint16_t *)malloc(rev_cap * sizeof(uint16_t));
    RansEncoder renc;
    rans_enc_init(&renc, rev_buf, rev_cap, 1);

    for (int32_t i = (int32_t)sink.rc - 1; i >= 0; i--)
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

    /* CDFs in ctx->cdf_state now carry into the next inter frame */
    ctx->cdf_state_valid = 1;

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
    snt_deblock_luma(ctx->ref_y, width, height, width, step,
                     is_inter ? blk_is_skip : NULL, mv_field, blk_cols);
    /* Chroma reference is stored during chroma encoding (reconstructed, not original) */
    ctx->num_refs = 1;

    /* Cleanup */
    free(work);
    free(recon);
    free(dc_recon_above);
    free(ac_recon_above);
    free(mv_field);
    free(blk_is_inter);
    free(blk_is_skip);
    free(ranges);
    free(ov_buf);
    free(rev_buf);
    return 0;
}
