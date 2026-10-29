/*
 * Sinter Codec — Entropy-coding context set (normative)
 *
 * Every adaptive CDF used by the frame syntax lives in one struct, so the
 * encoder and decoder initialize it identically. It is kept in the codec
 * context: reset on key frames, carried across inter frames.
 *
 * Contexts separate symbols whose statistics differ: block type (0 = intra,
 * 1 = inter residual), band index, and luma vs chroma. One shared CDF has to
 * average distributions that have little in common; with separate CDFs, the
 * mode/skip flags of mostly-empty bands become nearly free.
 */
#ifndef SNT_CDF_SET_H
#define SNT_CDF_SET_H

#include "../core/entropy/rans.h"
#include "../core/prediction/intra_spatial.h"

#define SNT_NUM_BANDS   16
#define SNT_AC_GROUPS    4  /* luma AC contexts: band / 4 */
#define SNT_CAC_GROUPS   3  /* chroma AC contexts: by frequency, see snt_cac_group */

/* Macroblock modes in inter frames (full 16x16 macroblocks only) */
enum {
    SNT_MB_SPLIT   = 0,  /* four 8x8 blocks, each with its own mode */
    SNT_MB_SKIP16  = 1,  /* all four blocks predicted at the MB's predicted MV */
    SNT_MB_INTER16 = 2,  /* one MV for the MB; each 8x8 residual coded */
    SNT_NUM_MB_MODES = 3
};

typedef struct {
    RansCDF dc[2];                    /* luma DC residual, by block type */
    RansCDF ac[2][SNT_AC_GROUPS];     /* luma scalar AC level */
    RansCDF mode[2][SNT_NUM_BANDS];   /* band mode: 0 = PVQ, 1 = scalar */
    RansCDF skip[2][SNT_NUM_BANDS];   /* PVQ band skip */
    RansCDF gain[2];                  /* PVQ gain */
    RansCDF imode[2];                 /* intra mode, by frame type (0 key, 1 inter) */
    RansCDF mbmode;                   /* macroblock mode (SNT_MB_*) */
    RansCDF bdir;                     /* B-frame prediction direction - 1 (SNT_DIR_*) */
    RansCDF inter, iskip, mvx, mvy;   /* inter-frame block modes and MV deltas */
    RansCDF cdc[2];                   /* chroma DC residual, by block type */
    RansCDF cac[2][SNT_CAC_GROUPS];   /* chroma AC level */
} SntCdfSet;

/* Chroma 4x4 AC context from raster index 1..15: low (x+y <= 1),
 * mid (x+y <= 3), high. */
static inline int snt_cac_group(int i) {
    int f = (i & 3) + (i >> 2);
    return f <= 1 ? 0 : (f <= 3 ? 1 : 2);
}

static inline void snt_cdf_set_init(SntCdfSet *s) {
    for (int t = 0; t < 2; t++) {
        rans_cdf_init_uniform(&s->dc[t],   256, RANS_ADAPT_RATE_COEFF);
        rans_cdf_init_uniform(&s->gain[t], 256, RANS_ADAPT_RATE_COEFF);
        rans_cdf_init_uniform(&s->imode[t], SNT_NUM_IPRED, RANS_ADAPT_RATE_MODE);
        rans_cdf_init_uniform(&s->cdc[t],  256, RANS_ADAPT_RATE_COEFF);
        for (int g = 0; g < SNT_AC_GROUPS; g++)
            rans_cdf_init_uniform(&s->ac[t][g], 256, RANS_ADAPT_RATE_COEFF);
        for (int g = 0; g < SNT_CAC_GROUPS; g++)
            rans_cdf_init_uniform(&s->cac[t][g], 256, RANS_ADAPT_RATE_COEFF);
        for (int b = 0; b < SNT_NUM_BANDS; b++) {
            rans_cdf_init_uniform(&s->mode[t][b], 2, RANS_ADAPT_RATE_MODE);
            rans_cdf_init_uniform(&s->skip[t][b], 2, RANS_ADAPT_RATE_MODE);
        }
    }
    rans_cdf_init_uniform(&s->mbmode, SNT_NUM_MB_MODES, RANS_ADAPT_RATE_MODE);
    rans_cdf_init_uniform(&s->bdir, 3, RANS_ADAPT_RATE_MODE);
    rans_cdf_init_uniform(&s->inter,  2, RANS_ADAPT_RATE_MODE);
    rans_cdf_init_uniform(&s->iskip,  2, RANS_ADAPT_RATE_MODE);
    rans_cdf_init_uniform(&s->mvx,  256, RANS_ADAPT_RATE_COEFF);
    rans_cdf_init_uniform(&s->mvy,  256, RANS_ADAPT_RATE_COEFF);
}

#endif /* SNT_CDF_SET_H */
