/*
 * Sinter Codec — Encoder with range coding + inter prediction
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
#include "../core/entropy/rangecoder.h"
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

/* B-frames are coded this many QP steps coarser: nothing predicts from
 * them, so their errors don't propagate. 4 was best by BD-rate (0, 2,
 * 5 and 6 tried). */
#define B_QP_OFFSET 4

/* Ideal cost in bits of coding sym with model m in its current state. */
static double cdf_bits(const SntModel *m, uint16_t sym) {
    return -log2((double)m->freq[sym] / m->total);
}

/* Length of the order-0 exp-Golomb code for v. */
static int expgolomb_len(uint32_t v) {
    int prefix = 0;
    for (uint32_t c = v + 1; c > 1; c >>= 1) prefix++;
    return 2 * prefix + 1;
}

/* Bits for a signed level: rANS symbol, plus overflow bits on escape. */
static double level_bits(const SntModel *cdf, int16_t level) {
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

/* Squared error of a reconstruction (clipped, as it will be output)
 * against the source over the block's in-frame vw x vh region. */
static double block_sse(const int16_t *src, const int16_t *rec, int vw, int vh) {
    double sse = 0;
    for (int y = 0; y < vh; y++)
        for (int x = 0; x < vw; x++) {
            double d = (double)src[y * 8 + x] - snt_clip_pixel(rec[y * 8 + x]);
            sse += d * d;
        }
    return sse;
}

/* ── Symbol output ─────────────────────────────────────────────────── */
/* Symbol sink: the range coder (each symbol adapts its model), and the
 * overflow bit writer for escapes and PVQ shapes. */
typedef struct {
    SntRangeEnc  *rc;
    SntBitWriter *ov;
} SymSink;

static void rec(SymSink *s, SntModel *m, uint16_t sym) {
    snt_rc_encode_sym(s->rc, m, sym);
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
 * With a sink, symbols are emitted and models adapt as they go. Without
 * one, nothing is emitted, the models are untouched, and out->bits
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
    /* Band flag context: is the previous band (for band 0, the DC level)
     * nonzero? Zero bands come in runs. */
    int prev_nz = dc_level != 0;
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
        double s_bits = cdf_bits(&cs->mode[t][b][prev_nz], 1);
        int s_nz = 0;
        for (int i = 0; i < n; i++) {
            s_level[i] = pvq_quantize_dc(band_coeffs[i], eff_step, 1);
            s_recon[i] = pvq_dequantize_dc(s_level[i], eff_step);
            s_bits += level_bits(&cs->ac[t][b / 4], s_level[i]);
            s_nz |= s_level[i] != 0;
        }

        /* PVQ candidate */
        SntPVQCode code;
        int16_t p_recon[16];
        pvq_encode(band_coeffs, n, step, 1, band_csf, &code, p_recon);
        uint16_t gain_sym = (code.gain > 255) ? 255 : (uint16_t)code.gain;
        int shape_bits = 0;
        double p_bits = cdf_bits(&cs->mode[t][b][prev_nz], 0)
                      + cdf_bits(&cs->skip[t][b][prev_nz], code.skip ? 1 : 0);
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
            rec(sink, &cs->mode[t][b][prev_nz], use_scalar ? 1 : 0);
            if (use_scalar) {
                for (int i = 0; i < n; i++) {
                    uint16_t ac_sym = s2u(s_level[i]);
                    rec(sink, &cs->ac[t][b / 4], ac_sym);
                    if (ac_sym == 255)
                        snt_bw_put_sexpgolomb(sink->ov, s_level[i]);
                }
            } else {
                rec(sink, &cs->skip[t][b][prev_nz], code.skip ? 1 : 0);
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
        prev_nz = use_scalar ? s_nz : !code.skip;
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

/* ── Per-frame luma coding state ───────────────────────────────────── */
typedef struct {
    SntEncoderCtx *ctx;
    SntCdfSet     *cs;
    SymSink       *sink;
    const uint8_t *input_y;  /* inter search reads it at stride width (CLI) */
    int16_t       *work;     /* source, pad_w stride, 0 outside the frame */
    int16_t       *recon;    /* reconstruction so far, pad_w stride */
    int            width, height, pad_w, blk_cols;
    int            is_inter, is_bidi;
    const SntSplinePlane *ref_f;  /* forward (past) luma reference */
    const SntSplinePlane *ref_b;  /* backward (future) luma reference, B-frames */
    uint16_t       step;
    double         lambda;   /* RD_LAMBDA * step^2 */
    SntMV         *mv_field;  /* per block: forward MV */
    SntMV         *mvb_field; /* per block: backward MV (B-frames) */
    uint8_t       *blk_dir;   /* per block: SNT_DIR_*, 0 for intra */
    int           *blk_is_inter, *blk_is_skip;
    int16_t       *dc_pix;   /* per block: pixel-domain DC */
    int16_t       *ac_state; /* per block: reconstructed low AC (key frames) */
    uint8_t       *coded;    /* per block: coded yet? */
} EncFrame;

enum { CAND_SKIP, CAND_INTER, CAND_INTRA };

/* A block's chosen mode */
typedef struct {
    int   kind;       /* CAND_* */
    int   imode;      /* intra mode (CAND_INTRA) */
    int   dir;        /* SNT_DIR_* (CAND_INTER) */
    SntMV mvf, mvb;   /* MVs (an unused direction keeps its predictor) */
} BlockChoice;

/* Source block (0 outside the frame) and its in-frame size */
static void load_src(const EncFrame *f, int bx, int by, int16_t *src, int *vw, int *vh) {
    int px0 = bx * BLK, py0 = by * BLK;
    for (int dy = 0; dy < BLK; dy++)
        for (int dx = 0; dx < BLK; dx++)
            src[dy * BLK + dx] = f->work[(py0 + dy) * f->pad_w + px0 + dx];
    *vw = f->width - px0 < BLK ? f->width - px0 : BLK;
    *vh = f->height - py0 < BLK ? f->height - py0 : BLK;
}

/* Transform-domain prediction for intra: DC from the neighbours'
 * pixel-domain DC; low AC (key frames) from their reconstructed AC. */
static void td_prediction(const EncFrame *f, int bx, int by, int16_t *td_pred) {
    int i = by * f->blk_cols + bx, up = i - f->blk_cols;
    memset(td_pred, 0, 64 * sizeof(int16_t));
    if (bx > 0 && by > 0)
        td_pred[0] = (f->dc_pix[i - 1] + f->dc_pix[up] + 1) >> 1;
    else if (bx > 0)
        td_pred[0] = f->dc_pix[i - 1];
    else if (by > 0)
        td_pred[0] = f->dc_pix[up];
    if (!f->is_inter)
        for (int p = 0; p < NUM_PRED_AC; p++) {
            int idx = pred_ac_idx[p];
            if (idx < 8 && bx > 0)
                td_pred[idx] = f->ac_state[(i - 1) * NUM_PRED_AC + p];
            else if (idx >= 8 && by > 0)
                td_pred[idx] = f->ac_state[up * NUM_PRED_AC + p];
        }
}

/* Predicted forward and backward MVs of block (bx, by) */
static void mv_predictors(const EncFrame *f, int bx, int by, SntMV *pf, SntMV *pb) {
    *pf = snt_mv_predict_at(f->mv_field, f->coded, f->blk_cols, bx, by);
    *pb = f->is_bidi ? snt_mv_predict_at(f->mvb_field, f->coded, f->blk_cols, bx, by)
                     : (SntMV){0, 0};
}

/* Inter prediction of a block and its residual (0 outside the frame) */
static void mc_residual(const EncFrame *f, int bx, int by, int dir, SntMV mvf, SntMV mvb,
                        const int16_t *src, int vw, int vh, uint8_t *pred, int16_t *res) {
    snt_pred_luma8(f->ref_f, f->ref_b, bx * BLK, by * BLK, dir, mvf, mvb, pred);
    for (int dy = 0; dy < BLK; dy++)
        for (int dx = 0; dx < BLK; dx++)
            res[dy * BLK + dx] = (dy < vh && dx < vw) ? src[dy * BLK + dx] - pred[dy * BLK + dx] : 0;
}

/* Quarter-pel MV for the blk x blk block at (px0, py0) in ref, searched
 * around pred; lambda for SAD decisions = sqrt(lambda) */
static SntMV search_mv(const EncFrame *f, const SntSplinePlane *ref, int px0, int py0, int blk,
                       SntMV pred) {
    double lambda_sad = sqrt(RD_LAMBDA) * f->step;
    SntMV s = snte_motion_search(f->input_y, ref->pix, f->width, px0, py0, blk, 16,
                                 f->width, f->height, pred, lambda_sad);
    return snte_subpel_refine(f->input_y, ref, f->width, f->height, px0, py0, blk,
                              s, pred, lambda_sad, NULL);
}

/* Bits for an inter block's direction (B-frames) and MV deltas */
static double inter_mv_bits(const EncFrame *f, int dir, SntMV mvf, SntMV mvb,
                            SntMV pf, SntMV pb) {
    double b = f->is_bidi ? cdf_bits(&f->cs->bdir, (uint16_t)(dir - 1)) : 0;
    if (dir & SNT_DIR_FWD)
        b += level_bits(&f->cs->mvx, mvf.x - pf.x) + level_bits(&f->cs->mvy, mvf.y - pf.y);
    if (dir & SNT_DIR_BWD)
        b += level_bits(&f->cs->mvx, mvb.x - pb.x) + level_bits(&f->cs->mvy, mvb.y - pb.y);
    return b;
}

/* Choose a block's mode by trial encoding: each candidate (skip at the
 * predicted MVs; inter at the searched MVs, forward and in B-frames also
 * backward and bi-directional; intra in each mode) is coded without
 * emitting anything, and the least distortion + lambda * bits wins.
 * Returns that cost. With `provisional`, the winner's reconstruction,
 * DC, MVs and coded flag are written so later blocks of the same
 * macroblock can be decided; emission overwrites them. */
static double decide_block(EncFrame *f, int bx, int by, int provisional, BlockChoice *bc) {
    SntCdfSet *cs = f->cs;
    int blk = by * f->blk_cols + bx, px0 = bx * BLK, py0 = by * BLK;
    int16_t src[64], td_pred[64], best_rec[64];
    int vw, vh;
    load_src(f, bx, by, src, &vw, &vh);
    td_prediction(f, bx, by, td_pred);

    double best_j = 1e300;
    int16_t best_dc = 0;
    SntMV zero = {0, 0};
    bc->kind = CAND_INTRA; bc->imode = SNT_IPRED_TD; bc->dir = 0;
    bc->mvf = zero; bc->mvb = zero;

    if (f->is_inter) {
        SntMV pf, pb;
        mv_predictors(f, bx, by, &pf, &pb);
        uint8_t p[64];
        int16_t p16[64], res[64];

        /* Skip: prediction at the predicted MV(s), nothing else coded */
        int sdir = f->is_bidi ? SNT_DIR_BI : SNT_DIR_FWD;
        snt_pred_luma8(f->ref_f, f->ref_b, px0, py0, sdir, pf, pb, p);
        for (int k = 0; k < 64; k++) p16[k] = p[k];
        double j = block_sse(src, p16, vw, vh)
                 + f->lambda * (cdf_bits(&cs->inter, 1) + cdf_bits(&cs->iskip, 1));
        if (j < best_j) {
            best_j = j; bc->kind = CAND_SKIP; bc->dir = sdir; bc->mvf = pf; bc->mvb = pb;
            memcpy(best_rec, p16, sizeof(best_rec)); best_dc = snt_pred_dc8(p);
        }

        /* Inter coded at the searched MVs */
        SntMV qf = search_mv(f, f->ref_f, px0, py0, BLK, pf);
        SntMV qb = f->is_bidi ? search_mv(f, f->ref_b, px0, py0, BLK, pb) : zero;
        for (int d = SNT_DIR_FWD; d <= (f->is_bidi ? SNT_DIR_BI : SNT_DIR_FWD); d++) {
            SntMV mvf = (d & SNT_DIR_FWD) ? qf : pf, mvb = (d & SNT_DIR_BWD) ? qb : pb;
            LumaCode lc;
            mc_residual(f, bx, by, d, mvf, mvb, src, vw, vh, p, res);
            code_luma(cs, NULL, f->step, 1, res, td_pred, 0, p, &lc);
            j = block_sse(src, lc.recon, vw, vh)
              + f->lambda * (cdf_bits(&cs->inter, 1) + cdf_bits(&cs->iskip, 0)
                             + inter_mv_bits(f, d, mvf, mvb, pf, pb) + lc.bits);
            if (j < best_j) {
                best_j = j; bc->kind = CAND_INTER; bc->dir = d; bc->mvf = mvf; bc->mvb = mvb;
                memcpy(best_rec, lc.recon, sizeof(best_rec)); best_dc = lc.dc_recon + snt_pred_dc8(p);
            }
        }
    }

    for (int m = SNT_IPRED_TD; m < SNT_NUM_IPRED; m++) {
        uint8_t p[64];
        int16_t in[64];
        if (m == SNT_IPRED_TD) {
            memcpy(in, src, sizeof(in));
        } else {
            snt_intra_spatial_pred8(f->recon, f->pad_w, px0, py0, m, p);
            for (int k = 0; k < 64; k++)
                in[k] = src[k] - p[k];
        }
        LumaCode lc;
        code_luma(cs, NULL, f->step, 0, in, td_pred, !f->is_inter && m == SNT_IPRED_TD,
                  m == SNT_IPRED_TD ? NULL : p, &lc);
        double bits = cdf_bits(&cs->imode[f->is_inter], (uint16_t)m) + lc.bits
                    + (f->is_inter ? cdf_bits(&cs->inter, 0) : 0);
        double j = block_sse(src, lc.recon, vw, vh) + f->lambda * bits;
        if (j < best_j) {
            best_j = j; bc->kind = CAND_INTRA; bc->imode = m; bc->dir = 0;
            bc->mvf = zero; bc->mvb = zero;
            memcpy(best_rec, lc.recon, sizeof(best_rec));
            best_dc = m == SNT_IPRED_TD ? lc.dc_recon : lc.dc_recon + snt_pred_dc8(p);
        }
    }

    if (provisional) {
        for (int dy = 0; dy < BLK; dy++)
            for (int dx = 0; dx < BLK; dx++)
                f->recon[(py0 + dy) * f->pad_w + px0 + dx] = best_rec[dy * BLK + dx];
        f->dc_pix[blk] = best_dc;
        f->mv_field[blk] = bc->mvf;
        f->mvb_field[blk] = bc->mvb;
        f->blk_dir[blk] = (uint8_t)bc->dir;
        f->coded[blk] = 1;
    }
    return best_j;
}

/* Record a block's state after coding: reconstruction, DC, MVs, flags */
static void finish_block(EncFrame *f, int bx, int by, const int16_t *rec, int16_t dc_pix,
                         int is_skip, int dir, SntMV mvf, SntMV mvb) {
    int blk = by * f->blk_cols + bx, px0 = bx * BLK, py0 = by * BLK;
    for (int dy = 0; dy < BLK; dy++)
        for (int dx = 0; dx < BLK; dx++)
            f->recon[(py0 + dy) * f->pad_w + px0 + dx] = rec[dy * BLK + dx];
    f->dc_pix[blk] = dc_pix;
    f->mv_field[blk] = mvf;
    f->mvb_field[blk] = mvb;
    f->blk_dir[blk] = (uint8_t)dir;
    f->blk_is_inter[blk] = dir != 0;
    f->blk_is_skip[blk] = is_skip;
    f->coded[blk] = 1;
}

/* Emit an inter-coded or skipped block's residual (no mode symbols) */
static void emit_inter_body(EncFrame *f, int bx, int by, int dir, SntMV mvf, SntMV mvb,
                            int skip) {
    int16_t src[64], res[64], rec[64], td_pred[64];
    uint8_t pred[64];
    int vw, vh;
    load_src(f, bx, by, src, &vw, &vh);
    mc_residual(f, bx, by, dir, mvf, mvb, src, vw, vh, pred, res);
    int blk = by * f->blk_cols + bx;
    memset(&f->ac_state[blk * NUM_PRED_AC], 0, NUM_PRED_AC * sizeof(int16_t));
    if (skip) {
        for (int k = 0; k < 64; k++) rec[k] = pred[k];
        finish_block(f, bx, by, rec, snt_pred_dc8(pred), 1, dir, mvf, mvb);
        return;
    }
    LumaCode lc;
    memset(td_pred, 0, sizeof(td_pred));
    code_luma(f->cs, f->sink, f->step, 1, res, td_pred, 0, pred, &lc);
    finish_block(f, bx, by, lc.recon, lc.dc_recon + snt_pred_dc8(pred), 0, dir, mvf, mvb);
}

static void emit_mvd(EncFrame *f, SntMV mv, SntMV mv_pred) {
    int16_t mvd_x = mv.x - mv_pred.x, mvd_y = mv.y - mv_pred.y;
    uint16_t sx = s2u(mvd_x), sy = s2u(mvd_y);
    rec(f->sink, &f->cs->mvx, sx);
    if (sx == 255) snt_bw_put_sexpgolomb(f->sink->ov, mvd_x);
    rec(f->sink, &f->cs->mvy, sy);
    if (sy == 255) snt_bw_put_sexpgolomb(f->sink->ov, mvd_y);
}

/* Emit an inter block's direction (B-frames) and MV deltas */
static void emit_inter_mvs(EncFrame *f, int dir, SntMV mvf, SntMV mvb, SntMV pf, SntMV pb) {
    if (f->is_bidi)
        rec(f->sink, &f->cs->bdir, (uint16_t)(dir - 1));
    if (dir & SNT_DIR_FWD) emit_mvd(f, mvf, pf);
    if (dir & SNT_DIR_BWD) emit_mvd(f, mvb, pb);
}

/* Emit one 8x8 block in its chosen mode, with its mode symbols */
static void emit_block(EncFrame *f, int bx, int by, const BlockChoice *bc) {
    SntCdfSet *cs = f->cs;
    int blk = by * f->blk_cols + bx, px0 = bx * BLK, py0 = by * BLK;

    if (f->is_inter) {
        SntMV pf, pb;
        mv_predictors(f, bx, by, &pf, &pb);
        rec(f->sink, &cs->inter, bc->kind != CAND_INTRA);
        if (bc->kind != CAND_INTRA) {
            rec(f->sink, &cs->iskip, bc->kind == CAND_SKIP);
            if (bc->kind == CAND_SKIP) {
                /* Skip predicts at the predicted MV(s) */
                emit_inter_body(f, bx, by, f->is_bidi ? SNT_DIR_BI : SNT_DIR_FWD, pf, pb, 1);
            } else {
                /* An unused direction keeps its predictor, as in the decoder */
                SntMV mvf = (bc->dir & SNT_DIR_FWD) ? bc->mvf : pf;
                SntMV mvb = (bc->dir & SNT_DIR_BWD) ? bc->mvb : pb;
                emit_inter_mvs(f, bc->dir, mvf, mvb, pf, pb);
                emit_inter_body(f, bx, by, bc->dir, mvf, mvb, 0);
            }
            return;
        }
    }

    /* Intra: transform-domain or spatial prediction */
    int16_t src[64], td_pred[64], in[64];
    uint8_t pred[64];
    int vw, vh;
    load_src(f, bx, by, src, &vw, &vh);
    td_prediction(f, bx, by, td_pred);
    int imode = bc->imode;
    rec(f->sink, &cs->imode[f->is_inter], (uint16_t)imode);
    if (imode == SNT_IPRED_TD) {
        memcpy(in, src, sizeof(in));
    } else {
        snt_intra_spatial_pred8(f->recon, f->pad_w, px0, py0, imode, pred);
        for (int k = 0; k < 64; k++)
            in[k] = src[k] - pred[k];
    }
    int td_ac = !f->is_inter && imode == SNT_IPRED_TD;
    LumaCode lc;
    code_luma(cs, f->sink, f->step, 0, in, td_pred, td_ac,
              imode == SNT_IPRED_TD ? NULL : pred, &lc);
    finish_block(f, bx, by, lc.recon,
                 imode == SNT_IPRED_TD ? lc.dc_recon : lc.dc_recon + snt_pred_dc8(pred),
                 0, 0, (SntMV){0, 0}, (SntMV){0, 0});

    /* Key frames: store reconstructed low AC for the neighbours'
     * prediction. Spatial-intra blocks coded a residual, so take the
     * coefficients of the reconstructed block itself. */
    if (!f->is_inter) {
        const int16_t *srcc = lc.recon_coeffs;
        int16_t full[64];
        if (imode != SNT_IPRED_TD) {
            snt_dct_fwd_2d_ref(lc.recon, BLK, BLK, full);
            srcc = full;
        }
        for (int p = 0; p < NUM_PRED_AC; p++)
            f->ac_state[blk * NUM_PRED_AC + p] = srcc[pred_ac_idx[p]];
    }
}

/* A 16x16 inter candidate */
typedef struct {
    int   dir;
    SntMV mvf, mvb;
    double j;
} Mb16Choice;

/* 16x16 candidates for the full macroblock whose top-left block is
 * (bx0, by0): skip at the predicted MVs (pf, pb), and inter at
 * 16x16-searched MVs with each 8x8 residual coded (forward; in B-frames
 * also backward and bi-directional). Returns the skip cost and the best
 * inter candidate. */
static double eval_mb16(EncFrame *f, int bx0, int by0, SntMV pf, SntMV pb, Mb16Choice *best) {
    SntCdfSet *cs = f->cs;
    int px0 = bx0 * BLK, py0 = by0 * BLK;
    int sdir = f->is_bidi ? SNT_DIR_BI : SNT_DIR_FWD;
    SntMV zero = {0, 0};
    SntMV qf = search_mv(f, f->ref_f, px0, py0, 16, pf);
    SntMV qb = f->is_bidi ? search_mv(f, f->ref_b, px0, py0, 16, pb) : zero;

    double d_skip = 0;
    best->j = 1e300;
    for (int d = SNT_DIR_FWD; d <= (f->is_bidi ? SNT_DIR_BI : SNT_DIR_FWD); d++) {
        SntMV mvf = (d & SNT_DIR_FWD) ? qf : pf, mvb = (d & SNT_DIR_BWD) ? qb : pb;
        double d_inter = 0, bits = 0;
        for (int z = 0; z < 4; z++) {
            int bx = bx0 + (z & 1), by = by0 + (z >> 1);
            int16_t src[64], res[64], p16[64], td0[64];
            uint8_t p[64];
            int vw, vh;
            load_src(f, bx, by, src, &vw, &vh);
            if (d == SNT_DIR_FWD) {
                mc_residual(f, bx, by, sdir, pf, pb, src, vw, vh, p, res);
                for (int k = 0; k < 64; k++) p16[k] = p[k];
                d_skip += block_sse(src, p16, vw, vh);
            }
            LumaCode lc;
            memset(td0, 0, sizeof(td0));
            mc_residual(f, bx, by, d, mvf, mvb, src, vw, vh, p, res);
            code_luma(cs, NULL, f->step, 1, res, td0, 0, p, &lc);
            d_inter += block_sse(src, lc.recon, vw, vh);
            bits += lc.bits;
        }
        double j = d_inter + f->lambda * (cdf_bits(&cs->mbmode, SNT_MB_INTER16) + bits
                                          + inter_mv_bits(f, d, mvf, mvb, pf, pb));
        if (j < best->j) {
            best->j = j; best->dir = d; best->mvf = mvf; best->mvb = mvb;
        }
    }
    return d_skip + f->lambda * cdf_bits(&cs->mbmode, SNT_MB_SKIP16);
}

/* Code one 16x16 macroblock (mb_x, mb_y in macroblock units) */
static void code_macroblock(EncFrame *f, int mb_x, int mb_y, int blk_rows) {
    int bx0 = mb_x * 2, by0 = mb_y * 2;
    int full = bx0 + 1 < f->blk_cols && by0 + 1 < blk_rows;

    if (!f->is_inter || !full) {
        /* Key frames and partial macroblocks: four independent 8x8 blocks */
        for (int z = 0; z < 4; z++) {
            int bx = bx0 + (z & 1), by = by0 + (z >> 1);
            if (bx >= f->blk_cols || by >= blk_rows) continue;
            BlockChoice bc;
            decide_block(f, bx, by, 0, &bc);
            emit_block(f, bx, by, &bc);
        }
        return;
    }

    /* Inter frame, full macroblock: split into 8x8 blocks, or one 16x16
     * skip / inter block, whichever costs least. The 16x16 predictors
     * come from neighbours outside the macroblock. */
    SntMV pf, pb;
    mv_predictors(f, bx0, by0, &pf, &pb);
    BlockChoice bc[4];
    double j_split = f->lambda * cdf_bits(&f->cs->mbmode, SNT_MB_SPLIT);
    for (int z = 0; z < 4; z++)
        j_split += decide_block(f, bx0 + (z & 1), by0 + (z >> 1), 1, &bc[z]);
    Mb16Choice m16;
    double j_skip = eval_mb16(f, bx0, by0, pf, pb, &m16);

    int mode = SNT_MB_SPLIT;
    if (j_skip < j_split && j_skip <= m16.j) mode = SNT_MB_SKIP16;
    else if (m16.j < j_split) mode = SNT_MB_INTER16;

    rec(f->sink, &f->cs->mbmode, (uint16_t)mode);
    if (mode == SNT_MB_SPLIT) {
        for (int z = 0; z < 4; z++)
            emit_block(f, bx0 + (z & 1), by0 + (z >> 1), &bc[z]);
    } else if (mode == SNT_MB_SKIP16) {
        for (int z = 0; z < 4; z++)
            emit_inter_body(f, bx0 + (z & 1), by0 + (z >> 1),
                            f->is_bidi ? SNT_DIR_BI : SNT_DIR_FWD, pf, pb, 1);
    } else {
        emit_inter_mvs(f, m16.dir, m16.mvf, m16.mvb, pf, pb);
        for (int z = 0; z < 4; z++)
            emit_inter_body(f, bx0 + (z & 1), by0 + (z >> 1), m16.dir, m16.mvf, m16.mvb, 0);
    }
}

/*
 * Symbol emission order (normative, must match decoder exactly):
 *
 * Contexts (cdf_set.h): t = block type (0 intra, 1 inter-coded),
 * ct = chroma block type (1 if inter-predicted).
 *
 * FOR each 16x16 macroblock (raster):
 *   IF inter frame AND all four 8x8 blocks inside the frame:
 *     REC(mbmode, m)                // SNT_MB_SPLIT / SKIP16 / INTER16
 *     IF SKIP16:  (all four blocks: MC prediction at the MB's predicted MV)
 *     IF INTER16: [B: REC(bdir, dir - 1)], REC(mvx, mvy) per used direction,
 *                 then FOR each 8x8 block in Z order:
 *                 DC + bands as an inter-coded block (below)
 *   FOR each 8x8 block in Z order (split macroblocks only):
 *   IF inter frame:
 *     REC(inter, inter_flag)        // 0=intra, 1=inter
 *     IF inter_flag:
 *       REC(iskip, skip_flag)       // 0=coded, 1=skip
 *       IF skip: (no more symbols; prediction at the predicted MV(s),
 *                 bi-directional in B-frames)
 *       IF coded:
 *         [B-frame: REC(bdir, dir - 1)]
 *         FOR each used direction (forward, then backward):
 *           REC(mvx, s2u(mvd_x)), REC(mvy, s2u(mvd_y))
 *   IF intra block: REC(imode[is_inter_frame], mode)  // SNT_IPRED_*
 *   REC(dc[t], s2u(dc_level))       // DC of residual (inter) or pixels (intra)
 *   // Low-AC prediction subtracted before band processing (key frames only)
 *   FOR band b=0..15:
 *     REC(mode[t][b][nz], mode)       // nz: previous band (band 0: DC level) nonzero
 *     IF scalar: FOR i=0..n-1: REC(ac[t][b/4], s2u(level))
 *     IF PVQ:    REC(skip[t][b][nz], skip)
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

    /* Determine frame type: key without a reference; B-frame when the
     * caller asks and an older anchor exists; P-frame otherwise */
    int is_inter = (ctx->num_refs > 0 && ctx->ref_y != NULL);
    int is_bidi = is_inter && ctx->code_bidi && ctx->ref_old_y != NULL;
    ctx->header.width = (uint32_t)width;
    ctx->header.height = (uint32_t)height;
    ctx->header.frame_type = !is_inter ? SNT_FRAME_KEY : (is_bidi ? SNT_FRAME_BIDI : SNT_FRAME_INTER);
    SntFrameHeader fhdr = ctx->header;  /* this frame's header (QP may differ) */
    if (!is_inter)
        fhdr.base_qp = fhdr.base_qp > KEY_QP_OFFSET ? fhdr.base_qp - KEY_QP_OFFSET : 0;
    if (is_bidi)
        fhdr.base_qp = fhdr.base_qp + B_QP_OFFSET < 63 ? fhdr.base_qp + B_QP_OFFSET : 63;
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

    /* ── Range coder and overflow buffers ──────────────────────────── */
    int cw = width / 2, ch = height / 2;
    int cblk = 4;
    int cpad_w = ((cw + cblk - 1) / cblk) * cblk;
    int cpad_h = ((ch + cblk - 1) / cblk) * cblk;

    uint8_t *rc_buf = (uint8_t *)malloc(bs_capacity);
    SntRangeEnc renc;
    snt_rc_enc_init(&renc, rc_buf, bs_capacity);

    uint32_t ov_cap = (uint32_t)(total_blocks * 128 + 4096);
    uint8_t *ov_buf = (uint8_t *)malloc(ov_cap);
    SntBitWriter ov;
    snt_bw_init(&ov, ov_buf, ov_cap);

    /* ── Entropy contexts: reset on key frames, carried across inter ──
     * B-frames keep their own set, carried from B-frame to B-frame. */
    if (!ctx->cdf_state)
        ctx->cdf_state = malloc(sizeof(SntCdfSet));
    if (!ctx->cdf_state_b)
        ctx->cdf_state_b = malloc(sizeof(SntCdfSet));
    if (!is_inter)
        ctx->cdf_state_b_valid = 0;
    SntCdfSet *cs = (SntCdfSet *)(is_bidi ? ctx->cdf_state_b : ctx->cdf_state);
    if (!is_inter || !(is_bidi ? ctx->cdf_state_b_valid : ctx->cdf_state_valid))
        snt_cdf_set_init(cs);

    SymSink sink = { &renc, &ov };
    #define REC(cdf, sym) rec(&sink, &(cdf), (sym))

    /* Per-block mode flags for chroma inter derivation */
    int *blk_is_inter = (int *)calloc(total_blocks, sizeof(int));
    int *blk_is_skip = (int *)calloc(total_blocks, sizeof(int));

    /* Per-block state for neighbours' prediction: pixel-domain DC, low AC
     * (key frames), and whether the block is coded yet */
    int16_t *dc_pix = (int16_t *)calloc(total_blocks, sizeof(int16_t));
    int16_t *ac_state = (int16_t *)calloc((size_t)total_blocks * NUM_PRED_AC, sizeof(int16_t));
    uint8_t *coded = (uint8_t *)calloc(total_blocks, 1);

    /* MV fields (forward, backward) and prediction directions, for MV
     * prediction, chroma and deblocking */
    SntMV *mv_field = (SntMV *)calloc(blk_cols * blk_rows, sizeof(SntMV));
    SntMV *mvb_field = (SntMV *)calloc(blk_cols * blk_rows, sizeof(SntMV));
    uint8_t *blk_dir = (uint8_t *)calloc(blk_cols * blk_rows, 1);

    /* Luma: 16x16 macroblocks in raster order; the 8x8 blocks of a split
     * macroblock in Z order (top-left, top-right, bottom-left, bottom-right).
     * P-frames predict from the newest anchor; B-frames forward from the
     * older anchor and backward from the newest. */
    SntSplinePlane sp_new = { ctx->ref_y, ctx->ref_coef, width, height };
    SntSplinePlane sp_old = { ctx->ref_old_y, ctx->ref_old_coef, width, height };
    EncFrame ef = { ctx, cs, &sink, input_y, work, recon, width, height, pad_w, blk_cols,
                    is_inter, is_bidi,
                    is_bidi ? &sp_old : &sp_new, is_bidi ? &sp_new : NULL,
                    step, RD_LAMBDA * (double)step * step, mv_field, mvb_field, blk_dir,
                    blk_is_inter, blk_is_skip, dc_pix, ac_state, coded };
    int mb_cols = (blk_cols + 1) / 2, mb_rows = (blk_rows + 1) / 2;
    for (int mb = 0; mb < mb_cols * mb_rows; mb++)
        code_macroblock(&ef, mb % mb_cols, mb / mb_cols, blk_rows);

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
    uint8_t *chroma_rec[2] = {NULL, NULL};
    if (input_u && input_v) {
        chroma_rec[0] = (uint8_t *)calloc((size_t)cw * ch, 1);
        chroma_rec[1] = (uint8_t *)calloc((size_t)cw * ch, 1);
        int cblk_cols = cpad_w / cblk;
        /* Reduce chroma QP offset for inter frames to prevent quantization
         * cascade (drift). +4 for key frames (no drift), +1 for inter. */
        int cqp_off = is_inter ? 1 : ctx->header.chroma_qp_offset;
        int cqp = fhdr.base_qp + cqp_off;
        if (cqp < 0) cqp = 0; if (cqp > 63) cqp = 63;
        uint16_t cstep = snt_qp_to_step((uint8_t)cqp);

        const uint8_t *chroma_in[2] = {input_u, input_v};
        const uint8_t *chroma_ref[2] = {is_bidi ? ctx->ref_old_u : ctx->ref_u,
                                         is_bidi ? ctx->ref_old_v : ctx->ref_v};
        const uint8_t *chroma_ref_b[2] = {is_bidi ? ctx->ref_u : NULL, is_bidi ? ctx->ref_v : NULL};
        for (int plane = 0; plane < 2; plane++) {
            const uint8_t *cp = chroma_in[plane];
            const uint8_t *cr = chroma_ref[plane]; /* may be NULL on first frame */
            const uint8_t *crb = chroma_ref_b[plane];
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
                        snt_pred_chroma4(cr, crb, cw, ch, cbx, cby, blk_dir[luma_blk_idx],
                                         mv_field[luma_blk_idx], mvb_field[luma_blk_idx], c_skip);
                        for (int dy = 0; dy < 4; dy++)
                            for (int dx = 0; dx < 4; dx++) {
                                int py = cby + dy, px = cbx + dx;
                                if (py < ch && px < cw && c_ref_out)
                                    c_ref_out[py*cw+px] = c_skip[dy*4+dx];
                            }
                        c_dc_left = snt_pred_dc4c(c_skip);
                        c_dc_above[ci] = c_dc_left;
                        continue;
                    }

                    /* Build 4x4 block input */
                    int16_t cblock[16];
                    uint8_t c_mc[16];  /* chroma MC prediction (inter only) */
                    if (c_is_inter && cr) {
                        /* Chroma inter: the luma MV (quarter-pel) is the
                         * chroma MV in 1/8-pel units */
                        snt_pred_chroma4(cr, crb, cw, ch, cbx, cby, blk_dir[luma_blk_idx],
                                         mv_field[luma_blk_idx], mvb_field[luma_blk_idx], c_mc);
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

                    /* DC prediction: intra blocks from the neighbours'
                     * pixel-domain DC; inter residuals are zero-mean */
                    int16_t c_dc_pred = 0;
                    if (c_bt)
                        c_dc_pred = 0;
                    else if (cbx > 0 && cby > 0)
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
                    /* Neighbours see the pixel-domain DC (inter: add the prediction's) */
                    c_dc_left = c_bt ? cdc_recon + snt_pred_dc4c(c_mc) : cdc_recon;
                    c_dc_above[ci] = c_dc_left;

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
        snt_deblock_chroma(chroma_rec[0], cw, ch, cw, cstep, is_inter ? blk_is_skip : NULL, mv_field, blk_cols);
        snt_deblock_chroma(chroma_rec[1], cw, ch, cw, cstep, is_inter ? blk_is_skip : NULL, mv_field, blk_cols);
    }

    /* ── Write the range-coded bytes, then the overflow bits ───────── */
    uint32_t rc_bytes = snt_rc_enc_finish(&renc);
    snt_bw_put_bits(&bw, rc_bytes, 32);
    for (uint32_t i = 0; i < rc_bytes; i++)
        snt_bw_put_byte(&bw, rc_buf[i]);

    snt_bw_byte_align(&ov);
    uint32_t ov_len = snt_bw_bytes_written(&ov);
    for (uint32_t i = 0; i < ov_len; i++)
        snt_bw_put_byte(&bw, ov_buf[i]);

    snt_bw_byte_align(&bw);
    *bs_len = snt_bw_bytes_written(&bw);

    /* Models now carry into the next frame of the same kind */
    if (is_bidi)
        ctx->cdf_state_b_valid = 1;
    else
        ctx->cdf_state_valid = 1;

    /* ── Reconstruction: kept in last_*; anchors (key and P-frames) also
     *    become the newest reference, and the previous newest the older.
     *    B-frames are not references. ─────────────────────────────────── */
    uint8_t *rec_y = (uint8_t *)malloc((size_t)width * height);
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
            rec_y[y * width + x] = snt_clip_pixel((int32_t)recon[y * pad_w + x]);
    snt_deblock_luma(rec_y, width, height, width, step,
                     is_inter ? blk_is_skip : NULL, mv_field, blk_cols);
    if (!ctx->last_y || ctx->ref_w != width || ctx->ref_h != height) {
        free(ctx->last_y); free(ctx->last_u); free(ctx->last_v);
        ctx->last_y = (uint8_t *)malloc((size_t)width * height);
        ctx->last_u = (uint8_t *)calloc((size_t)cw * ch, 1);
        ctx->last_v = (uint8_t *)calloc((size_t)cw * ch, 1);
    }
    memcpy(ctx->last_y, rec_y, (size_t)width * height);
    if (chroma_rec[0]) {
        memcpy(ctx->last_u, chroma_rec[0], (size_t)cw * ch);
        memcpy(ctx->last_v, chroma_rec[1], (size_t)cw * ch);
    }
    if (!is_bidi) {
        free(ctx->ref_old_y);
        free(ctx->ref_old_coef);
        ctx->ref_old_y = (ctx->ref_w == width && ctx->ref_h == height) ? ctx->ref_y : NULL;
        ctx->ref_old_coef = ctx->ref_old_y ? ctx->ref_coef : NULL;
        if (!ctx->ref_old_y) { free(ctx->ref_y); free(ctx->ref_coef); }
        ctx->ref_y = rec_y;
        ctx->ref_coef = snt_spline_coefs(rec_y, width, height);
        if (chroma_rec[0]) {
            free(ctx->ref_old_u); free(ctx->ref_old_v);
            ctx->ref_old_u = ctx->ref_u;
            ctx->ref_old_v = ctx->ref_v;
            ctx->ref_u = chroma_rec[0];
            ctx->ref_v = chroma_rec[1];
        }
        ctx->ref_w = width;
        ctx->ref_h = height;
        ctx->num_refs = 1;
    } else {
        free(rec_y);
        free(chroma_rec[0]);
        free(chroma_rec[1]);
    }

    /* Cleanup */
    free(work);
    free(recon);
    free(dc_pix);
    free(ac_state);
    free(coded);
    free(mv_field);
    free(mvb_field);
    free(blk_dir);
    free(blk_is_inter);
    free(blk_is_skip);
    free(rc_buf);
    free(ov_buf);
    return 0;
}
