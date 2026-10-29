/*
 * Sinter Codec — Entropy-coding context set (normative)
 *
 * Every adaptive model used by the frame syntax lives in one struct, so the
 * encoder and decoder initialize it identically. It is kept in the codec
 * context: reset on key frames, carried across inter frames.
 *
 * Contexts separate symbols whose statistics differ: block type (0 = intra,
 * 1 = inter residual), band index, and luma vs chroma. One shared model has
 * to average distributions that have little in common; with separate models,
 * the flags of mostly-empty bands become nearly free. Conditioning each
 * symbol's model on earlier symbols is context modeling as in Rissanen &
 * Langdon, "Universal modeling and coding" (IEEE Trans. IT, 1981).
 */
#ifndef SNT_CDF_SET_H
#define SNT_CDF_SET_H

#include "../core/entropy/rangecoder.h"
#include "../core/prediction/intra_acdc.h"

#define SNT_NUM_BANDS   16
#define SNT_AC_GROUPS    4  /* luma AC contexts: band / 4 */
#define SNT_CAC_GROUPS   3  /* chroma AC contexts: by frequency, see snt_cac_group */

/* Macroblock modes in inter frames (full 16x16 macroblocks only) */
enum {
    SNT_MB_SPLIT   = 0,  /* four 8x8 blocks, each with its own mode */
    SNT_MB_SKIP16  = 1,  /* all four blocks predicted at zero motion */
    SNT_MB_INTER16 = 2,  /* one MV for the MB; each 8x8 residual coded */
    SNT_NUM_MB_MODES = 3
};

typedef struct {
    SntModel dc[2];                    /* luma DC residual, by block type */
    SntModel ac[2][SNT_AC_GROUPS];     /* luma scalar AC level */
    SntModel band[2][SNT_NUM_BANDS][2];  /* band coded (1) or all zero (0); by previous band nonzero */
    SntModel imode[2];                 /* intra mode, by frame type (0 key, 1 inter) */
    SntModel mbmode;                   /* macroblock mode (SNT_MB_*) */
    SntModel bdir;                     /* B-frame prediction direction - 1 (SNT_DIR_*) */
    SntModel inter, iskip, mvx, mvy;   /* inter-frame block modes and MV deltas */
    SntModel cdc[2];                   /* chroma DC residual, by block type */
    SntModel cac[2][SNT_CAC_GROUPS];   /* chroma AC level */
} SntCdfSet;

/* Chroma 4x4 AC context from raster index 1..15: low (x+y <= 1),
 * mid (x+y <= 3), high. */
static inline int snt_cac_group(int i) {
    int f = (i & 3) + (i >> 2);
    return f <= 1 ? 0 : (f <= 3 ? 1 : 2);
}

/* Number of models in a set: the struct holds nothing else */
#define SNT_NUM_MODELS ((int)(sizeof(SntCdfSet) / sizeof(SntModel)))

/* Model sets: one for key and P-frames, one for B-frames */
enum { SNT_SET_P = 0, SNT_SET_B = 1 };

/* Starting tables are trained per QP class of the frame where the set
 * restarts: QP 0-11, 12-19, 20-27, 28+ */
#define SNT_MODEL_CLASSES 4
static inline int snt_model_class(int qp) {
    int c = (qp - 4) / 8;
    return c < 0 ? 0 : (c >= SNT_MODEL_CLASSES ? SNT_MODEL_CLASSES - 1 : c);
}

static inline void snt_cdf_set_init_uniform(SntCdfSet *s) {
    for (int t = 0; t < 2; t++) {
        snt_model_init(&s->dc[t], 256, SNT_RC_INIT_COEFF, SNT_RC_INC_COEFF, SNT_RC_LIMIT_COEFF);
        snt_model_init(&s->imode[t], SNT_NUM_IPRED, SNT_RC_INIT_MODE, SNT_RC_INC_MODE, SNT_RC_LIMIT_MODE);
        snt_model_init(&s->cdc[t], 256, SNT_RC_INIT_COEFF, SNT_RC_INC_COEFF, SNT_RC_LIMIT_COEFF);
        for (int g = 0; g < SNT_AC_GROUPS; g++)
            snt_model_init(&s->ac[t][g], 256, SNT_RC_INIT_COEFF, SNT_RC_INC_COEFF, SNT_RC_LIMIT_COEFF);
        for (int g = 0; g < SNT_CAC_GROUPS; g++)
            snt_model_init(&s->cac[t][g], 256, SNT_RC_INIT_COEFF, SNT_RC_INC_COEFF, SNT_RC_LIMIT_COEFF);
        for (int b = 0; b < SNT_NUM_BANDS; b++) {
            for (int z = 0; z < 2; z++)
                snt_model_init(&s->band[t][b][z], 2, SNT_RC_INIT_MODE, SNT_RC_INC_MODE, SNT_RC_LIMIT_MODE);
        }
    }
    snt_model_init(&s->mbmode, SNT_NUM_MB_MODES, SNT_RC_INIT_MODE, SNT_RC_INC_MODE, SNT_RC_LIMIT_MODE);
    snt_model_init(&s->bdir, 3, SNT_RC_INIT_MODE, SNT_RC_INC_MODE, SNT_RC_LIMIT_MODE);
    snt_model_init(&s->inter, 2, SNT_RC_INIT_MODE, SNT_RC_INC_MODE, SNT_RC_LIMIT_MODE);
    snt_model_init(&s->iskip, 2, SNT_RC_INIT_MODE, SNT_RC_INC_MODE, SNT_RC_LIMIT_MODE);
    snt_model_init(&s->mvx, 256, SNT_RC_INIT_COEFF, SNT_RC_INC_COEFF, SNT_RC_LIMIT_COEFF);
    snt_model_init(&s->mvy, 256, SNT_RC_INIT_COEFF, SNT_RC_INC_COEFF, SNT_RC_LIMIT_COEFF);
}

#include "model_tables.h"

/* Weight of a trained starting table, in observations (increments):
 * chosen by BD-rate on the tuning clips (16 to 256 were within 0.1%) */
#define SNT_TRAINED_OBS_COEFF 64
#define SNT_TRAINED_OBS_MODE  16

/* Restart a model set for a frame of the given QP: each model starts
 * from its trained table (model_tables.h) for the set and the QP's
 * class, scaled to a total worth SNT_TRAINED_OBS_* observations, with
 * every count at least 1. Models without enough training data start
 * uniform. */
static inline void snt_cdf_set_init(SntCdfSet *s, int set, int qp) {
    snt_cdf_set_init_uniform(s);
    int cls = snt_model_class(qp);
    SntModel *m = (SntModel *)s;
    for (int k = 0; k < SNT_NUM_MODELS; k++) {
        const uint16_t *p = &snt_trained_prob[set][cls][snt_trained_offset[k]];
        if (!snt_trained_valid[set][cls][k]) continue;
        uint32_t obs = m[k].num_symbols == 256 ? SNT_TRAINED_OBS_COEFF : SNT_TRAINED_OBS_MODE;
        uint32_t strength = obs * m[k].inc, total = 0;
        if (strength > m[k].limit / 2) strength = m[k].limit / 2;
        for (int i = 0; i < m[k].num_symbols; i++) {
            uint32_t f = ((uint32_t)p[i] * strength + 16384) >> 15;
            m[k].freq[i] = (uint16_t)(f < 1 ? 1 : f);
            total += m[k].freq[i];
        }
        m[k].total = total;
    }
}

#endif /* SNT_CDF_SET_H */
