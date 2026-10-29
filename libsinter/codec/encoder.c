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
#ifdef SNT_EXPERIMENT_CLIPTABLES
#include "cliptables.h"
#endif
#ifdef SNT_EXPERIMENT_MIXING
#include "mixer.h"
/* Mixing state of the model set in use (one encoder per process) */
static SntMixState enc_mix[2];
static SntMixState *cur_mix;
#endif
#include "scan.h"
#include "../core/transform/dct.h"
#include "../core/quantize/pvq.h"
#include "../core/quantize/qp_table.h"
#include "../core/bitstream/bitstream.h"
#include "../core/entropy/rangecoder.h"
#include "../core/prediction/motion_est.h"
#include "../core/prediction/mc.h"
#include "../core/prediction/global_motion.h"
#include "../core/filter/deblock.h"
#include "../core/prediction/intra_acdc.h"
#include "../util/fixed_point.h"

#define BLK 8

/* ── Zigzag signed-to-unsigned mapping ────────────────────────────── */
static inline uint16_t s2u(int16_t v) {
    uint32_t u = (v >= 0) ? (uint32_t)(2 * v) : (uint32_t)(-2 * v - 1);
    return (u > 254) ? 255 : (uint16_t)u;
}

/* ── Rate estimates for encoder decisions (non-normative) ──────────── */
#ifndef RD_LAMBDA
#define RD_LAMBDA 0.14  /* lambda = RD_LAMBDA * step^2 (orthonormal DCT, SSE) */
#endif

/* Key frames are coded this many QP steps finer than inter frames: every
 * frame until the next key frame predicts from them, so their quality pays
 * off repeatedly. The frame header carries the QP actually used. */
#define KEY_QP_OFFSET 4

/* B-frames are coded this many QP steps coarser: nothing predicts from
 * them, so their errors don't propagate. 4 was best by BD-rate (0, 2,
 * 5 and 6 tried). */
#define B_QP_OFFSET 4

/* Global motion estimation counts a block as predicted well enough to
 * skip when its SSE is below GM_SKIP_SSE * lambda */
#ifndef GM_SKIP_SSE
#define GM_SKIP_SSE 30
#endif

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

/* Bits for a luma AC level at zigzag index k, coded as binary decisions
 * (levels.h) with the models' current state */
static double ac_level_bits(SntCdfSet *cs, int t, int k, int pc, int16_t v, int z_known,
                            int band, int nzprev) {
    SntLevelBin bins[2 * SNT_LEVEL_MAX_CAT + 4];
    int n = snt_level_bins(v, z_known, bins);
    double bits = 0;
    for (int i = 0; i < n; i++) {
#ifdef SNT_EXPERIMENT_MIXING
        if (bins[i].kind == 0 || bins[i].kind == 2) {
            SntMixBin mb;
            snt_mix_prepare(cur_mix, cs, bins[i].kind == 2, t, k, band, pc, nzprev, &mb);
            bits += -log2((bins[i].bit ? mb.p : 4096 - mb.p) / 4096.0);
            continue;
        }
#endif
        bits += cdf_bits(snt_level_model(cs, t, k, pc, &bins[i]), bins[i].bit);
    }
    (void)band; (void)nzprev;
    return bits;
}

/* Number of nonzero levels among the two before zigzag index k */
static inline int nz_prev2(const int16_t *lv, int k) {
    return (lv[k - 1] != 0) + (k >= 2 && lv[k - 2] != 0);
}

/* RD-optimized quantization: the candidates for a coefficient's level are
 * the rounded level l0, one step toward zero, and zero (K. Ramchandran and
 * M. Vetterli, "Rate-distortion optimal fast thresholding with complete
 * JPEG/MPEG decoder compatibility", IEEE TIP 3(5), 1994). Returns their
 * number; zero is left out when no_zero. */
static int rd_candidates(int16_t l0, int no_zero, int16_t *cand) {
    int n = 0;
    if (l0 != 0) cand[n++] = l0;
    if (l0 > 1 || l0 < -1) cand[n++] = (int16_t)(l0 > 0 ? l0 - 1 : l0 + 1);
    if (!no_zero || n == 0) cand[n++] = 0;
    return n;
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
#ifdef SNT_TRAIN_STATS
    SntModel     *base;   /* first model of the set in use */
    uint32_t    (*hist)[256];
    int           set;
#endif
} SymSink;

#ifdef SNT_TRAIN_STATS
/* Symbol counts per model set, QP class at the set's restart, and model
 * (tools/train/collect_stats.c) */
static uint32_t train_hist[2][SNT_MODEL_CLASSES][SNT_NUM_MODELS][256];
static int train_cls[2];

uint32_t snte_train_count(int set, int k, int sym) {
    uint32_t n = 0;
    for (int c = 0; c < SNT_MODEL_CLASSES; c++) n += train_hist[set][c][k][sym];
    return n;
}

void snte_train_reset(void) { memset(train_hist, 0, sizeof(train_hist)); }

void snte_train_dump(FILE *f) {
    SntCdfSet sizes;
    snt_cdf_set_init_uniform(&sizes);
    for (int set = 0; set < 2; set++)
        for (int c = 0; c < SNT_MODEL_CLASSES; c++)
            for (int k = 0; k < SNT_NUM_MODELS; k++) {
                fprintf(f, "%d %d %d %d", set, c, k, ((SntModel *)&sizes)[k].num_symbols);
                for (int i = 0; i < 256; i++) fprintf(f, " %u", train_hist[set][c][k][i]);
                fprintf(f, "\n");
            }
}
#endif

#ifdef SNT_TRAIN_STATS
/* Ideal cost in bits of every coded symbol, by model set and model; and
 * bytes outside the range coder (headers, escapes) by frame type */
static double train_bits[2][SNT_NUM_MODELS];
static double train_other[3];
void snte_bits_dump(FILE *f) {
    for (int set = 0; set < 2; set++)
        for (int k = 0; k < SNT_NUM_MODELS; k++)
            fprintf(f, "bits %d %d %.0f\n", set, k, train_bits[set][k]);
    for (int t = 0; t < 3; t++) fprintf(f, "other %d %.0f\n", t, train_other[t]);
}
#endif

static void rec(SymSink *s, SntModel *m, uint16_t sym) {
#ifdef SNT_TRAIN_STATS
    s->hist[m - s->base][sym]++;
    train_bits[s->set][m - s->base] += -log2((double)m->freq[sym] / m->total);
#endif
    snt_rc_encode_sym(s->rc, m, sym);
}

/* ── Luma block coding ─────────────────────────────────────────────── */
#define NUM_BANDS 16
static const int band_bounds[NUM_BANDS + 1] = {0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,59,63};

/* Result of coding one 8x8 luma block's DC and bands */
typedef struct {
    int16_t recon_coeffs[64];  /* reconstructed coefficients (TD prediction included) */
    int16_t recon[64];         /* reconstructed pixels (prediction block added) */
    int16_t dc_recon;          /* reconstructed DC coefficient */
    double  bits;              /* estimated bits (estimate mode only) */
} LumaCode;

/* Code one 8x8 luma block's DC and bands.
 * dct_in: pixels (intra) or a residual (inter).
 * td_pred: coefficient prediction (intra_acdc.h), used when pred is NULL:
 *   the DC, plus the first row or column when td_ac. pred: prediction
 *   block added back on reconstruction (NULL for intra).
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

    /* DC: transform-domain intra predicts it; other residuals are zero-mean.
     * Keep the best zero and the best nonzero candidate: whether the DC
     * level is zero is the first band flag's context, so it is decided
     * with the bands below. */
    int16_t dc_pred = pred ? 0 : td_pred[0];
    int16_t dc_cand[2] = {0, 0};  /* by nonzero? */
    double dc_j[2] = {1e300, 1e300};
    {
        int16_t cand[3];
        int nc = rd_candidates(pvq_quantize_dc(coeffs[0] - dc_pred, step, 1), 0, cand);
        for (int c = 0; c < nc; c++) {
            double d = (double)(coeffs[0] - dc_pred) - pvq_dequantize_dc(cand[c], step);
            double j = d * d + lambda * level_bits(&cs->dc[t], cand[c]);
            int nz = cand[c] != 0;
            if (j < dc_j[nz]) { dc_j[nz] = j; dc_cand[nz] = cand[c]; }
        }
    }

    /* AC prediction (intra VERT and HORZ modes) */
    if (td_ac)
        for (int i = 1; i < 64; i++)
            coeffs[i] -= td_pred[i];

    /* Contexts from the prediction (inter blocks; levels.h) */
    uint8_t kcls[64] = {0}, bcls[NUM_BANDS] = {0};
    if (pred && t == 1) {
        int16_t p16[64], pdct[64];
        for (int i = 0; i < 64; i++) p16[i] = pred[i];
        snt_dct_fwd_2d_ref(p16, BLK, BLK, pdct);
        snt_pred_classes(pdct, step, snt_zigzag8, band_bounds, NUM_BANDS, kcls, bcls);
    }

    /* Each band is either all zero (flag 0) or has its levels coded (flag
     * 1). First each band's levels if it is coded, chosen coefficient by
     * coefficient; zigzag order, index 1 first. */
    int16_t ac[63], s_level[63], s_recon[63];
    double j_coded[NUM_BANDS], j_zero[NUM_BANDS], lev_bits[NUM_BANDS];
    int band_nz[NUM_BANDS];
    int16_t lv[64] = {0};  /* levels so far, zigzag order */
    for (int b = 0; b < NUM_BANDS; b++) {
        int start = band_bounds[b], n = band_bounds[b + 1] - start;
        int s_nz = 0;
        lev_bits[b] = 0;
        for (int i = start; i < start + n; i++) {
            int k = 1 + i;
            ac[i] = coeffs[snt_zigzag8[k]];
            /* The last level of a coded band is known nonzero if the
             * others are zero, so it can't be zero then */
            int z_known = i == start + n - 1 && !s_nz;
            int16_t cand[3];
            int nc = rd_candidates(pvq_quantize_dc(ac[i], step, 1), z_known, cand);
            double best = 1e300, best_bits = 0;
            for (int c = 0; c < nc; c++) {
                double d = (double)ac[i] - pvq_dequantize_dc(cand[c], step);
                double cb = ac_level_bits(cs, t, k, kcls[k], cand[c], z_known, b, nz_prev2(lv, k));
                if (d * d + lambda * cb < best) {
                    best = d * d + lambda * cb; best_bits = cb; s_level[i] = cand[c];
                }
            }
            s_recon[i] = pvq_dequantize_dc(s_level[i], step);
            lev_bits[b] += best_bits;
            lv[k] = s_level[i];
            s_nz |= s_level[i] != 0;
        }
        int16_t zeros[16] = {0};
        band_nz[b] = s_nz;
        j_coded[b] = band_sse(&ac[start], &s_recon[start], n) + lambda * lev_bits[b];
        j_zero[b] = band_sse(&ac[start], zeros, n);
    }

    /* Then coded or zero for all bands at once. A band flag's model
     * depends on whether the band before it (for band 0, the DC level) is
     * nonzero, since zero bands come in runs; the least total distortion
     * + lambda * bits is a shortest path through two states per band
     * (Viterbi), as in trellis quantization (J. Wen, M. Luttrell and
     * J. Villasenor, "Trellis-based R-D optimal quantization in H.263+",
     * IEEE TIP 9(8), 2000). */
    double cost[2] = { dc_j[0], dc_j[1] };
    uint8_t from[NUM_BANDS][2];
    for (int b = 0; b < NUM_BANDS; b++) {
        double next[2] = {1e300, 1e300};
        for (int st = 0; st < 2; st++) {
            if (cost[st] >= 1e300) continue;
            const SntModel *flag = &cs->band[t][b][st][bcls[b]];
            double jz = cost[st] + j_zero[b] + lambda * cdf_bits(flag, 0);
            if (jz < next[0]) { next[0] = jz; from[b][0] = (uint8_t)st; }
            if (band_nz[b]) {
                double jc = cost[st] + j_coded[b] + lambda * cdf_bits(flag, 1);
                if (jc < next[1]) { next[1] = jc; from[b][1] = (uint8_t)st; }
            }
        }
        cost[0] = next[0]; cost[1] = next[1];
    }
    int coded[NUM_BANDS];
    int st = cost[1] < cost[0];
    for (int b = NUM_BANDS - 1; b >= 0; b--) {
        coded[b] = st;
        st = from[b][st];
    }
    int16_t dc_level = dc_cand[st];

    if (sink) {
        uint16_t dc_sym = s2u(dc_level);
        rec(sink, &cs->dc[t], dc_sym);
        if (dc_sym == 255)
            snt_bw_put_sexpgolomb(sink->ov, dc_level);
    } else {
        bits += level_bits(&cs->dc[t], dc_level);
    }
    out->dc_recon = pvq_dequantize_dc(dc_level, step) + dc_pred;

    int16_t recon_ac[63];
    memset(recon_ac, 0, sizeof(recon_ac));
    memset(lv, 0, sizeof(lv));
    lv[0] = dc_level;
    int prev_nz = dc_level != 0;
    for (int b = 0; b < NUM_BANDS; b++) {
        int start = band_bounds[b], n = band_bounds[b + 1] - start;
        SntModel *flag = &cs->band[t][b][prev_nz][bcls[b]];
        if (!sink) {
            bits += coded[b] ? cdf_bits(flag, 1) + lev_bits[b] : cdf_bits(flag, 0);
        } else {
            rec(sink, flag, (uint16_t)coded[b]);
            if (coded[b]) {
                int any = 0;
                for (int i = start; i < start + n; i++) {
                    int k = 1 + i;
                    SntLevelBin bins[2 * SNT_LEVEL_MAX_CAT + 4];
                    int nb = snt_level_bins(s_level[i], i == start + n - 1 && !any, bins);
                    for (int j = 0; j < nb; j++) {
#ifdef SNT_EXPERIMENT_MIXING
                        if (bins[j].kind == 0 || bins[j].kind == 2) {
                            SntMixBin mb;
                            snt_mix_prepare(cur_mix, cs, bins[j].kind == 2, t, k, b, kcls[k],
                                            nz_prev2(lv, k), &mb);
                            snt_rc_encode_bit(sink->rc, bins[j].bit, mb.p);
                            snt_mix_update(&mb, bins[j].bit);
                            continue;
                        }
#endif
                        rec(sink, snt_level_model(cs, t, k, kcls[k], &bins[j]), bins[j].bit);
                    }
                    any |= s_level[i] != 0;
                    lv[k] = s_level[i];
                }
            }
        }
        prev_nz = coded[b];
        if (coded[b])
            for (int i = start; i < start + n; i++)
                recon_ac[i] = s_recon[i];
    }

    /* Reconstruct coefficients, then pixels */
    out->recon_coeffs[0] = out->dc_recon;
    for (int i = 0; i < 63; i++)
        out->recon_coeffs[snt_zigzag8[i + 1]] = recon_ac[i];
    if (td_ac)
        for (int i = 1; i < 64; i++)
            out->recon_coeffs[i] += td_pred[i];
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
    SntMV          gmf, gmb;  /* global motion: skipped blocks predict with it */
    uint16_t       step;
    double         lambda;   /* RD_LAMBDA * step^2 */
    SntMV         *mv_field;  /* per block: forward MV */
    SntMV         *mvb_field; /* per block: backward MV (B-frames) */
    uint8_t       *blk_dir;   /* per block: SNT_DIR_*, 0 for intra */
    int           *blk_is_inter, *blk_is_skip;
    int16_t       *acdc;     /* per block: SNT_ACDC_STATE coefficients of the
                                reconstruction, for intra prediction */
    uint8_t       *coded;    /* per block: coded yet? */
    /* Extra MV candidates for the blocks being decided (encoder decisions):
     * the macroblock's 16x16-searched MVs for its 8x8 blocks, the 8x8
     * blocks' MVs for the 16x16 candidate */
    SntMV          xf[4], xb[4];
    int            nx;
} EncFrame;

enum { CAND_SKIP, CAND_INTER, CAND_INTRA };

/* A block's chosen mode */
typedef struct {
    int   kind;       /* CAND_* */
    int   imode;      /* intra mode (CAND_INTRA) */
    int   dir;        /* SNT_DIR_* (CAND_INTER) */
    SntMV mvf, mvb;   /* MVs (zero for an unused direction) */
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

/* Coefficient prediction for an intra block in mode m (intra_acdc.h) */
static void td_prediction(const EncFrame *f, int bx, int by, int m, int16_t *td_pred) {
    int i = by * f->blk_cols + bx;
    snt_acdc_predict(m, bx > 0 ? &f->acdc[(i - 1) * SNT_ACDC_STATE] : NULL,
                     by > 0 ? &f->acdc[(i - f->blk_cols) * SNT_ACDC_STATE] : NULL, td_pred);
}

/* Record a reconstructed block's coefficients for its neighbours' intra
 * prediction */
static void store_acdc(EncFrame *f, int blk, const int16_t *rec) {
    int16_t c[64];
    snt_dct_fwd_2d_ref(rec, BLK, BLK, c);
    snt_acdc_store(c, &f->acdc[blk * SNT_ACDC_STATE]);
}

/* Predicted forward and backward MVs of block (bx, by) */
static void mv_predictors(const EncFrame *f, int bx, int by, SntMV *pf, SntMV *pb) {
    *pf = snt_mv_predict_at(f->mv_field, f->blk_cols, bx, by);
    *pb = f->is_bidi ? snt_mv_predict_at(f->mvb_field, f->blk_cols, bx, by)
                     : (SntMV){0, 0};
}

/* Prediction of a skipped block (bx, by): as the block to its left if
 * that block is inter-predicted (MPEG-1's rule for skipped macroblocks
 * in B-pictures), otherwise with the global motion */
static void skip_motion(const EncFrame *f, int bx, int by, int *dir, SntMV *mvf, SntMV *mvb) {
    int left = by * f->blk_cols + bx - 1;
    if (bx > 0 && f->blk_dir[left]) {
        *dir = f->blk_dir[left];
        *mvf = f->mv_field[left];
        *mvb = f->mvb_field[left];
    } else {
        *dir = f->is_bidi ? SNT_DIR_BI : SNT_DIR_FWD;
        *mvf = f->gmf;
        *mvb = f->gmb;
    }
}

/* Prediction agreement class of 8x8 block (bx, by) (cdf_set.h), with its
 * skip motion and predictors */
static uint32_t agree_sad8(const EncFrame *f, int bx, int by, int sdir, SntMV smf, SntMV smb,
                           SntMV pf, SntMV pb) {
    SntMV zero = {0, 0};
    if (f->is_bidi)
        return snt_pred_sad8(f->ref_f, f->ref_b, bx * BLK, by * BLK,
                             SNT_DIR_FWD, pf, zero, SNT_DIR_BWD, zero, pb);
    return snt_pred_sad8(f->ref_f, f->ref_b, bx * BLK, by * BLK,
                         sdir, smf, smb, SNT_DIR_FWD, pf, zero);
}

static int agree_class_block(const EncFrame *f, int bx, int by) {
    int sdir;
    SntMV smf, smb, pf, pb;
    skip_motion(f, bx, by, &sdir, &smf, &smb);
    mv_predictors(f, bx, by, &pf, &pb);
    return snt_agree_class(agree_sad8(f, bx, by, sdir, smf, smb, pf, pb), 64, f->step);
}

/* The same for the 16x16 macroblock whose top-left block is (bx0, by0),
 * with the macroblock's skip motion and predictors */
static int agree_class_mb(const EncFrame *f, int bx0, int by0) {
    int sdir;
    SntMV smf, smb, pf, pb;
    skip_motion(f, bx0, by0, &sdir, &smf, &smb);
    mv_predictors(f, bx0, by0, &pf, &pb);
    uint32_t sad = 0;
    for (int z = 0; z < 4; z++)
        sad += agree_sad8(f, bx0 + (z & 1), by0 + (z >> 1), sdir, smf, smb, pf, pb);
    return snt_agree_class(sad, 256, f->step);
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

/* MV candidates for an inter block whose top-left 8x8 block is (bx, by),
 * besides the searched ones (cf[0], cb[0]): the predicted MVs (no MV
 * delta), the skip motion, the MVs of the inter blocks to the left, above
 * and above right, and f's extra candidates. Returns their number, the
 * searched ones included. w_blocks: the block's width in 8x8 blocks. */
#define MAX_MV_CAND 10
static int mv_candidates(const EncFrame *f, int bx, int by, int w_blocks,
                         SntMV qf, SntMV qb, SntMV pf, SntMV pb, SntMV *cf, SntMV *cb) {
    int n = 0, sdir;
    SntMV smf, smb;
    cf[n] = qf; cb[n++] = qb;
    cf[n] = pf; cb[n++] = pb;
    skip_motion(f, bx, by, &sdir, &smf, &smb);
    cf[n] = smf; cb[n++] = smb;
    int i = by * f->blk_cols + bx;
    if (bx > 0 && f->blk_dir[i - 1]) { cf[n] = f->mv_field[i - 1]; cb[n++] = f->mvb_field[i - 1]; }
    if (by > 0) {
        int a = i - f->blk_cols, ar = a + w_blocks;
        if (f->blk_dir[a]) { cf[n] = f->mv_field[a]; cb[n++] = f->mvb_field[a]; }
        if (bx + w_blocks < f->blk_cols && f->coded[ar] && f->blk_dir[ar]) {
            cf[n] = f->mv_field[ar]; cb[n++] = f->mvb_field[ar];
        }
    }
    for (int x = 0; x < f->nx; x++) { cf[n] = f->xf[x]; cb[n++] = f->xb[x]; }
    return n;
}

static inline int mv_eq(SntMV a, SntMV b) { return a.x == b.x && a.y == b.y; }

/* Is candidate c, in direction d, the same as an earlier one? */
static int mv_cand_dup(int d, const SntMV *cf, const SntMV *cb, int c) {
    for (int e = 0; e < c; e++)
        if ((!(d & SNT_DIR_FWD) || mv_eq(cf[e], cf[c])) && (!(d & SNT_DIR_BWD) || mv_eq(cb[e], cb[c])))
            return 1;
    return 0;
}

/* Trial order of the inter candidates: each candidate forward, then
 * backward, then bi-directional (B-frames), then bi-directional from the
 * best forward and the best backward candidate, held in slot ncand */
static void mv_cand_index(int ncand, int cand, int *d, int *c) {
    if (cand < 3 * ncand) { *d = SNT_DIR_FWD + cand / ncand; *c = cand % ncand; }
    else { *d = SNT_DIR_BI; *c = ncand; }
}

/* Keep the best forward and backward candidates (cost j) in slot ncand */
static void mv_cand_best(int d, double j, SntMV mvf, SntMV mvb, double *jdir,
                         SntMV *cf, SntMV *cb, int ncand) {
    if (d == SNT_DIR_FWD && j < jdir[0]) { jdir[0] = j; cf[ncand] = mvf; }
    if (d == SNT_DIR_BWD && j < jdir[1]) { jdir[1] = j; cb[ncand] = mvb; }
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

/* Cost (distortion + lambda * bits) of 8x8 block (bx, by) coded inter in
 * direction d at (mvf, mvb), with its mode flags (agreement class ac) and
 * MV deltas; its reconstruction in rec */
static double inter8_cost(const EncFrame *f, int bx, int by, const int16_t *src, int vw, int vh,
                          int ac, int d, SntMV mvf, SntMV mvb, SntMV pf, SntMV pb, int16_t *rec) {
    uint8_t p[64];
    int16_t res[64], td0[64] = {0};
    LumaCode lc;
    mc_residual(f, bx, by, d, mvf, mvb, src, vw, vh, p, res);
    code_luma(f->cs, NULL, f->step, 1, res, td0, 0, p, &lc);
    memcpy(rec, lc.recon, sizeof(lc.recon));
    return block_sse(src, lc.recon, vw, vh)
         + f->lambda * (cdf_bits(&f->cs->inter[ac], 1) + cdf_bits(&f->cs->iskip[ac], 0)
                        + inter_mv_bits(f, d, mvf, mvb, pf, pb) + lc.bits);
}

/* The same for the macroblock whose top-left block is (bx0, by0), coded
 * INTER16 (macroblock mode model by agreement class mac) */
static double inter16_cost(const EncFrame *f, int bx0, int by0, int mac, int d, SntMV mvf,
                           SntMV mvb, SntMV pf, SntMV pb) {
    double dist = 0, bits = 0;
    for (int z = 0; z < 4; z++) {
        int bx = bx0 + (z & 1), by = by0 + (z >> 1), vw, vh;
        int16_t src[64], res[64], td0[64] = {0};
        uint8_t p[64];
        LumaCode lc;
        load_src(f, bx, by, src, &vw, &vh);
        mc_residual(f, bx, by, d, mvf, mvb, src, vw, vh, p, res);
        code_luma(f->cs, NULL, f->step, 1, res, td0, 0, p, &lc);
        dist += block_sse(src, lc.recon, vw, vh);
        bits += lc.bits;
    }
    return dist + f->lambda * (cdf_bits(&f->cs->mbmode[mac], SNT_MB_INTER16) + bits
                               + inter_mv_bits(f, d, mvf, mvb, pf, pb));
}

/* Neighbour k (0-7) of mv at a quarter sample */
static inline SntMV mv_neighbour(SntMV mv, int k) {
    if (k >= 4) k++;
    return (SntMV){ (int16_t)(mv.x + k % 3 - 1), (int16_t)(mv.y + k / 3 - 1) };
}

/* Choose a block's mode by trial encoding: each candidate (skip with the
 * global motion; inter at the searched MVs, forward and in B-frames also
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

    double best_j = 1e300;
    SntMV zero = {0, 0};
    bc->kind = CAND_INTRA; bc->imode = SNT_IPRED_DC; bc->dir = 0;
    bc->mvf = zero; bc->mvb = zero;

    int ac = f->is_inter ? agree_class_block(f, bx, by) : 0;
    if (f->is_inter) {
        SntMV pf, pb;
        mv_predictors(f, bx, by, &pf, &pb);
        uint8_t p[64];
        int16_t p16[64];

        /* Skip: prediction with the global motion (bi-directional in
         * B-frames), nothing else coded */
        int sdir;
        SntMV smf, smb;
        skip_motion(f, bx, by, &sdir, &smf, &smb);
        snt_pred_luma8(f->ref_f, f->ref_b, px0, py0, sdir, smf, smb, p);
        for (int k = 0; k < 64; k++) p16[k] = p[k];
        double j = block_sse(src, p16, vw, vh)
                 + f->lambda * (cdf_bits(&cs->inter[ac], 1) + cdf_bits(&cs->iskip[ac], 1));
        if (j < best_j) {
            best_j = j; bc->kind = CAND_SKIP; bc->dir = sdir; bc->mvf = smf; bc->mvb = smb;
            memcpy(best_rec, p16, sizeof(best_rec));
        }

        /* Inter coded at the searched MVs and the other candidates */
        SntMV qf = search_mv(f, f->ref_f, px0, py0, BLK, pf);
        SntMV qb = f->is_bidi ? search_mv(f, f->ref_b, px0, py0, BLK, pb) : zero;
        SntMV cf[MAX_MV_CAND + 1], cb[MAX_MV_CAND + 1];
        int ncand = mv_candidates(f, bx, by, 1, qf, qb, pf, pb, cf, cb);
        double jdir[2] = {1e300, 1e300}, ij = 1e300;
        int id = 0;
        SntMV imf = zero, imb = zero;
        int16_t irec[64], trec[64];
        for (int cand = 0; cand < (f->is_bidi ? 3 * ncand + 1 : ncand); cand++) {
            int d, c;
            mv_cand_index(ncand, cand, &d, &c);
            if (mv_cand_dup(d, cf, cb, c)) continue;
            /* An unused direction's MV is zero (H.263 Annex O) */
            SntMV mvf = (d & SNT_DIR_FWD) ? cf[c] : zero;
            SntMV mvb = (d & SNT_DIR_BWD) ? cb[c] : zero;
            j = inter8_cost(f, bx, by, src, vw, vh, ac, d, mvf, mvb, pf, pb, trec);
            mv_cand_best(d, j, mvf, mvb, jdir, cf, cb, ncand);
            if (j < ij) {
                ij = j; id = d; imf = mvf; imb = mvb;
                memcpy(irec, trec, sizeof(irec));
            }
        }

        /* The best inter candidate, refined by the same cost: the 8
         * quarter-sample neighbours of each MV it uses, the other MV held
         * (for bi-directional prediction, as in S.-W. Wu and A. Gersho,
         * "Joint estimation of forward and backward motion vectors for
         * interpolative prediction of video", IEEE TIP 3(5), 1994) */
        for (int dirbit = SNT_DIR_FWD; dirbit <= SNT_DIR_BWD; dirbit <<= 1) {
            if (!(id & dirbit)) continue;
            SntMV c0 = dirbit == SNT_DIR_FWD ? imf : imb, f0 = imf, b0 = imb;
            for (int k = 0; k < 8; k++) {
                SntMV mvf = dirbit == SNT_DIR_FWD ? mv_neighbour(c0, k) : f0;
                SntMV mvb = dirbit == SNT_DIR_BWD ? mv_neighbour(c0, k) : b0;
                j = inter8_cost(f, bx, by, src, vw, vh, ac, id, mvf, mvb, pf, pb, trec);
                if (j < ij) {
                    ij = j; imf = mvf; imb = mvb;
                    memcpy(irec, trec, sizeof(irec));
                }
            }
        }
        if (ij < best_j) {
            best_j = ij; bc->kind = CAND_INTER; bc->dir = id; bc->mvf = imf; bc->mvb = imb;
            memcpy(best_rec, irec, sizeof(best_rec));
        }
    }

    for (int m = 0; m < SNT_NUM_IPRED; m++) {
        td_prediction(f, bx, by, m, td_pred);
        LumaCode lc;
        code_luma(cs, NULL, f->step, 0, src, td_pred, 1, NULL, &lc);
        double bits = cdf_bits(&cs->imode[f->is_inter], (uint16_t)m) + lc.bits
                    + (f->is_inter ? cdf_bits(&cs->inter[ac], 0) : 0);
        double j = block_sse(src, lc.recon, vw, vh) + f->lambda * bits;
        if (j < best_j) {
            best_j = j; bc->kind = CAND_INTRA; bc->imode = m; bc->dir = 0;
            bc->mvf = zero; bc->mvb = zero;
            memcpy(best_rec, lc.recon, sizeof(best_rec));
        }
    }

    if (provisional) {
        for (int dy = 0; dy < BLK; dy++)
            for (int dx = 0; dx < BLK; dx++)
                f->recon[(py0 + dy) * f->pad_w + px0 + dx] = best_rec[dy * BLK + dx];
        store_acdc(f, blk, best_rec);
        f->mv_field[blk] = bc->mvf;
        f->mvb_field[blk] = bc->mvb;
        f->blk_dir[blk] = (uint8_t)bc->dir;
        f->coded[blk] = 1;
    }
    return best_j;
}

/* Record a block's state after coding: reconstruction, coefficients for
 * intra prediction, MVs, flags */
static void finish_block(EncFrame *f, int bx, int by, const int16_t *rec,
                         int is_skip, int dir, SntMV mvf, SntMV mvb) {
    int blk = by * f->blk_cols + bx, px0 = bx * BLK, py0 = by * BLK;
    for (int dy = 0; dy < BLK; dy++)
        for (int dx = 0; dx < BLK; dx++)
            f->recon[(py0 + dy) * f->pad_w + px0 + dx] = rec[dy * BLK + dx];
    store_acdc(f, blk, rec);
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
    if (skip) {
        for (int k = 0; k < 64; k++) rec[k] = pred[k];
        finish_block(f, bx, by, rec, 1, dir, mvf, mvb);
        return;
    }
    LumaCode lc;
    memset(td_pred, 0, sizeof(td_pred));
    code_luma(f->cs, f->sink, f->step, 1, res, td_pred, 0, pred, &lc);
    finish_block(f, bx, by, lc.recon, 0, dir, mvf, mvb);
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

    if (f->is_inter) {
        SntMV pf, pb;
        mv_predictors(f, bx, by, &pf, &pb);
        int ac = agree_class_block(f, bx, by);
        rec(f->sink, &cs->inter[ac], bc->kind != CAND_INTRA);
        if (bc->kind != CAND_INTRA) {
            rec(f->sink, &cs->iskip[ac], bc->kind == CAND_SKIP);
            if (bc->kind == CAND_SKIP) {
                /* Skip predicts with the global motion */
                int sdir;
                SntMV smf, smb;
                skip_motion(f, bx, by, &sdir, &smf, &smb);
                emit_inter_body(f, bx, by, sdir, smf, smb, 1);
            } else {
                /* An unused direction's MV is zero, as in the decoder */
                SntMV zero = {0, 0};
                SntMV mvf = (bc->dir & SNT_DIR_FWD) ? bc->mvf : zero;
                SntMV mvb = (bc->dir & SNT_DIR_BWD) ? bc->mvb : zero;
                emit_inter_mvs(f, bc->dir, mvf, mvb, pf, pb);
                emit_inter_body(f, bx, by, bc->dir, mvf, mvb, 0);
            }
            return;
        }
    }

    /* Intra: coefficient prediction in the chosen mode */
    int16_t src[64], td_pred[64];
    int vw, vh;
    load_src(f, bx, by, src, &vw, &vh);
    rec(f->sink, &cs->imode[f->is_inter], (uint16_t)bc->imode);
    td_prediction(f, bx, by, bc->imode, td_pred);
    LumaCode lc;
    code_luma(cs, f->sink, f->step, 0, src, td_pred, 1, NULL, &lc);
    finish_block(f, bx, by, lc.recon, 0, 0, (SntMV){0, 0}, (SntMV){0, 0});
}

/* A 16x16 inter candidate */
typedef struct {
    int   dir;
    SntMV mvf, mvb;
    double j;
} Mb16Choice;

/* 16x16 candidates for the full macroblock whose top-left block is
 * (bx0, by0): skip (motion from skip_motion of the top-left block), and
 * inter at
 * 16x16-searched MVs with each 8x8 residual coded (forward; in B-frames
 * also backward and bi-directional). Returns the skip cost and the best
 * inter candidate. */
static double eval_mb16(EncFrame *f, int bx0, int by0, SntMV pf, SntMV pb, SntMV qf, SntMV qb,
                        int mac, Mb16Choice *best) {
    SntCdfSet *cs = f->cs;
    int sdir;
    SntMV smf, smb, zero = {0, 0};
    skip_motion(f, bx0, by0, &sdir, &smf, &smb);

    double d_skip = 0;
    for (int z = 0; z < 4; z++) {
        int bx = bx0 + (z & 1), by = by0 + (z >> 1), vw, vh;
        int16_t src[64], res[64], p16[64];
        uint8_t p[64];
        load_src(f, bx, by, src, &vw, &vh);
        mc_residual(f, bx, by, sdir, smf, smb, src, vw, vh, p, res);
        for (int k = 0; k < 64; k++) p16[k] = p[k];
        d_skip += block_sse(src, p16, vw, vh);
    }

    /* At the searched MVs and the other candidates */
    best->j = 1e300;
    SntMV cf[MAX_MV_CAND + 1], cb[MAX_MV_CAND + 1];
    int ncand = mv_candidates(f, bx0, by0, 2, qf, qb, pf, pb, cf, cb);
    double jdir[2] = {1e300, 1e300};
    for (int cand = 0; cand < (f->is_bidi ? 3 * ncand + 1 : ncand); cand++) {
        int d, c;
        mv_cand_index(ncand, cand, &d, &c);
        if (mv_cand_dup(d, cf, cb, c)) continue;
        SntMV mvf = (d & SNT_DIR_FWD) ? cf[c] : zero;
        SntMV mvb = (d & SNT_DIR_BWD) ? cb[c] : zero;
        double j = inter16_cost(f, bx0, by0, mac, d, mvf, mvb, pf, pb);
        mv_cand_best(d, j, mvf, mvb, jdir, cf, cb, ncand);
        if (j < best->j) {
            best->j = j; best->dir = d; best->mvf = mvf; best->mvb = mvb;
        }
    }

    /* Refined as for 8x8 blocks (decide_block) */
    for (int dirbit = SNT_DIR_FWD; dirbit <= SNT_DIR_BWD; dirbit <<= 1) {
        if (!(best->dir & dirbit)) continue;
        SntMV c0 = dirbit == SNT_DIR_FWD ? best->mvf : best->mvb, f0 = best->mvf, b0 = best->mvb;
        for (int k = 0; k < 8; k++) {
            SntMV mvf = dirbit == SNT_DIR_FWD ? mv_neighbour(c0, k) : f0;
            SntMV mvb = dirbit == SNT_DIR_BWD ? mv_neighbour(c0, k) : b0;
            double j = inter16_cost(f, bx0, by0, mac, best->dir, mvf, mvb, pf, pb);
            if (j < best->j) { best->j = j; best->mvf = mvf; best->mvb = mvb; }
        }
    }
    return d_skip + f->lambda * cdf_bits(&cs->mbmode[mac], SNT_MB_SKIP16);
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
    int mac = agree_class_mb(f, bx0, by0);  /* before the blocks' provisional state */
    SntMV zero = {0, 0};
    SntMV qf = search_mv(f, f->ref_f, bx0 * BLK, by0 * BLK, 16, pf);
    SntMV qb = f->is_bidi ? search_mv(f, f->ref_b, bx0 * BLK, by0 * BLK, 16, pb) : zero;
    BlockChoice bc[4];
    double j_split = f->lambda * cdf_bits(&f->cs->mbmode[mac], SNT_MB_SPLIT);
    f->xf[0] = qf; f->xb[0] = qb; f->nx = 1;
    for (int z = 0; z < 4; z++)
        j_split += decide_block(f, bx0 + (z & 1), by0 + (z >> 1), 1, &bc[z]);
    f->nx = 0;
    for (int z = 0; z < 4; z++)
        if (bc[z].kind == CAND_INTER) { f->xf[f->nx] = bc[z].mvf; f->xb[f->nx++] = bc[z].mvb; }
    Mb16Choice m16;
    double j_skip = eval_mb16(f, bx0, by0, pf, pb, qf, qb, mac, &m16);
    f->nx = 0;

    int mode = SNT_MB_SPLIT;
    if (j_skip < j_split && j_skip <= m16.j) mode = SNT_MB_SKIP16;
    else if (m16.j < j_split) mode = SNT_MB_INTER16;

    rec(f->sink, &f->cs->mbmode[mac], (uint16_t)mode);
    if (mode == SNT_MB_SPLIT) {
        for (int z = 0; z < 4; z++)
            emit_block(f, bx0 + (z & 1), by0 + (z >> 1), &bc[z]);
    } else if (mode == SNT_MB_SKIP16) {
        int sdir;
        SntMV smf, smb;
        skip_motion(f, bx0, by0, &sdir, &smf, &smb);
        for (int z = 0; z < 4; z++)
            emit_inter_body(f, bx0 + (z & 1), by0 + (z >> 1), sdir, smf, smb, 1);
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
 *     IF SKIP16:  (all four blocks: MC prediction with the global motion)
 *     IF INTER16: [B: REC(bdir, dir - 1)], REC(mvx, mvy) per used direction,
 *                 then FOR each 8x8 block in Z order:
 *                 DC + bands as an inter-coded block (below)
 *   FOR each 8x8 block in Z order (split macroblocks only):
 *   IF inter frame:
 *     REC(inter, inter_flag)        // 0=intra, 1=inter
 *     IF inter_flag:
 *       REC(iskip, skip_flag)       // 0=coded, 1=skip
 *       IF skip: (no more symbols; prediction with the global motion,
 *                 bi-directional in B-frames)
 *       IF coded:
 *         [B-frame: REC(bdir, dir - 1)]
 *         FOR each used direction (forward, then backward):
 *           REC(mvx, s2u(mvd_x)), REC(mvy, s2u(mvd_y))
 *   IF intra block: REC(imode[is_inter_frame], mode)  // SNT_IPRED_* (intra_acdc.h)
 *   REC(dc[t], s2u(dc_level))       // DC of residual (inter) or pixels (intra)
 *   // Intra: the mode's coefficient prediction is subtracted first
 *   FOR band b=0..15:
 *     REC(band[t][b][nz][pc], coded)  // nz: previous band (band 0: DC level) nonzero;
 *                                     // pc: prediction class (inter, levels.h)
 *     IF coded: FOR i=0..n-1: the level's binary decisions (levels.h);
 *               the last level's Z decision is skipped if the others are 0
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
    if (!is_inter)
        ctx->header.key_frame_number = ctx->header.frame_number;
    SntFrameHeader fhdr = ctx->header;  /* this frame's header (QP may differ) */
    if (!is_inter)
        fhdr.base_qp = fhdr.base_qp > KEY_QP_OFFSET ? fhdr.base_qp - KEY_QP_OFFSET : 0;
    if (is_bidi)
        fhdr.base_qp = fhdr.base_qp + B_QP_OFFSET < 63 ? fhdr.base_qp + B_QP_OFFSET : 63;

    uint16_t step = snt_qp_to_step(fhdr.base_qp);

    /* Load pixels into work buffer */
    int16_t *work = (int16_t *)calloc((size_t)pad_w * pad_h, sizeof(int16_t));
    int16_t *recon = (int16_t *)calloc((size_t)pad_w * pad_h, sizeof(int16_t));
    if (!work || !recon) { free(work); free(recon); return -1; }

    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
            work[y * pad_w + x] = (int16_t)input_y[y * stride_y + x];

    /* Global motion of this frame against each reference (encoder's
     * estimate); skipped blocks predict with it */
    SntSplinePlane sp_new = { ctx->ref_y, ctx->ref_coef, width, height };
    SntSplinePlane sp_old = { ctx->ref_old_y, ctx->ref_old_coef, width, height };
    SntMV gmf = {0, 0}, gmb = {0, 0};
    if (is_inter) {
        double skip_sse = GM_SKIP_SSE * RD_LAMBDA * (double)step * step;
        gmf = snte_global_motion(input_y, is_bidi ? &sp_old : &sp_new, width, height, skip_sse);
        if (is_bidi)
            gmb = snte_global_motion(input_y, &sp_new, width, height, skip_sse);
    }

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
    if (!is_inter || !(is_bidi ? ctx->cdf_state_b_valid : ctx->cdf_state_valid)) {
        snt_cdf_set_init(cs, is_bidi ? SNT_SET_B : SNT_SET_P, fhdr.base_qp);
#ifdef SNT_EXPERIMENT_MIXING
        snt_mix_init(&enc_mix[is_bidi]);
#endif
#ifdef SNT_EXPERIMENT_CLIPTABLES
        if (ctx->clip_tables)
            snt_clip_apply(cs, (const SntClipTables *)ctx->clip_tables, is_bidi);
#endif
    }
#ifdef SNT_EXPERIMENT_CLIPTABLES
    /* Key frames start with a flag: per-clip tables follow */
    if (!is_inter) {
        int send = ctx->clip_tables && ctx->clip_tables_send;
        snt_rc_encode_bit(&renc, send, 2048);
        if (send) snt_clip_encode(&renc, (const SntClipTables *)ctx->clip_tables);
        ctx->clip_tables_send = 0;
    }
#endif
#ifdef SNT_EXPERIMENT_MIXING
    cur_mix = &enc_mix[is_bidi];
#endif

#ifdef SNT_TRAIN_STATS
    int set_id = is_bidi ? SNT_SET_B : SNT_SET_P;
    if (!is_inter || !(is_bidi ? ctx->cdf_state_b_valid : ctx->cdf_state_valid))
        train_cls[set_id] = snt_model_class(fhdr.base_qp);
    SymSink sink = { &renc, &ov, (SntModel *)cs, train_hist[set_id][train_cls[set_id]], set_id };
#else
    SymSink sink = { &renc, &ov };
#endif
    #define REC(cdf, sym) rec(&sink, &(cdf), (sym))

    /* Per-block mode flags for chroma inter derivation */
    int *blk_is_inter = (int *)calloc(total_blocks, sizeof(int));
    int *blk_is_skip = (int *)calloc(total_blocks, sizeof(int));

    /* Per-block state for neighbours' prediction: coefficients of the
     * reconstruction, and whether the block is coded yet */
    int16_t *acdc = (int16_t *)calloc((size_t)total_blocks * SNT_ACDC_STATE, sizeof(int16_t));
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
    EncFrame ef = { ctx, cs, &sink, input_y, work, recon, width, height, pad_w, blk_cols,
                    is_inter, is_bidi,
                    is_bidi ? &sp_old : &sp_new, is_bidi ? &sp_new : NULL, gmf, gmb,
                    step, RD_LAMBDA * (double)step * step, mv_field, mvb_field, blk_dir,
                    blk_is_inter, blk_is_skip, acdc, coded, {{0, 0}}, {{0, 0}}, 0 };
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
                     * the luma block's MV, which is zero */
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
        snt_deblock_chroma(chroma_rec[0], cw, ch, cw, cstep, is_inter ? blk_is_skip : NULL, blk_cols);
        snt_deblock_chroma(chroma_rec[1], cw, ch, cw, cstep, is_inter ? blk_is_skip : NULL, blk_cols);
    }

    /* ── Frame: header, overflow bits, then the range-coded bytes up to
     *    the end of the frame ───────────────────────────────────────── */
    uint32_t rc_bytes = snt_rc_enc_finish(&renc);
    snt_bw_byte_align(&ov);
    uint32_t ov_len = snt_bw_bytes_written(&ov);

    snt_write_frame_header(&bw, &fhdr);
    snt_bw_put_bits(&bw, (input_u && input_v) ? 1u : 0u, 1);
    if (is_inter) {
        snt_bw_put_sexpgolomb(&bw, gmf.x);
        snt_bw_put_sexpgolomb(&bw, gmf.y);
    }
    if (is_bidi) {
        snt_bw_put_sexpgolomb(&bw, gmb.x);
        snt_bw_put_sexpgolomb(&bw, gmb.y);
    }
    snt_bw_put_expgolomb(&bw, ov_len);
    snt_bw_byte_align(&bw);
    for (uint32_t i = 0; i < ov_len; i++)
        snt_bw_put_byte(&bw, ov_buf[i]);
    for (uint32_t i = 0; i < rc_bytes; i++)
        snt_bw_put_byte(&bw, rc_buf[i]);
    *bs_len = snt_bw_bytes_written(&bw);
#ifdef SNT_TRAIN_STATS
    train_other[fhdr.frame_type] += 8.0 * (*bs_len - rc_bytes);
#endif

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
                         is_inter ? blk_is_skip : NULL, blk_cols);
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
    free(acdc);
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
