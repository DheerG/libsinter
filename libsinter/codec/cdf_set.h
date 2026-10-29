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
#include "levels.h"

#define SNT_NUM_BANDS   16

/* Prediction agreement (novel tool 4, docs/NOVEL.md): before a block's
 * mode is read, the decoder forms two predictions it could make from the
 * references alone (P-frames: the skip prediction and the prediction at
 * the median MV predictor; B-frames: forward at the forward predictor and
 * backward at the backward predictor) and classes their mean absolute
 * difference against the quantizer step: below step/8, below step/2,
 * above. The class selects the mode models. */
#define SNT_AGREE_CLASSES 3
#define SNT_BDIR_DIRECT   3   /* bdir symbol: direct mode (snt_direct_mvs), no MV deltas */
#define SNT_FLAG_CONTEXTS (SNT_AGREE_CLASSES * 3)  /* mode-flag contexts, see snt_flag_context */
#define SNT_DC_CONTEXTS  10  /* luma DC contexts, see snt_dc_context */
static inline int snt_agree_class(uint32_t sad, int npix, uint16_t step) {
    uint64_t s8 = 8ull * sad, s2 = 2ull * sad, ns = (uint64_t)npix * step;
    return s8 < ns ? 0 : (s2 < ns ? 1 : 2);
}

/* Macroblock modes in inter frames (full 16x16 macroblocks only) */
enum {
    SNT_MB_SPLIT   = 0,  /* four 8x8 blocks, each with its own mode */
    SNT_MB_SKIP16  = 1,  /* all four blocks predicted at zero motion */
    SNT_MB_INTER16 = 2,  /* one MV for the MB; each 8x8 residual coded */
    SNT_NUM_MB_MODES = 3
};

typedef struct {
    /* luma DC level as binary decisions (levels.h), by block type; Z, S and
     * G also by the DC context (snt_dc_context), G also by sign */
    SntModel dz[2][SNT_DC_CONTEXTS], ds[2][SNT_DC_CONTEXTS], dg[2][SNT_DC_CONTEXTS][2];
    SntModel dx[2][SNT_LEVEL_MAX_CAT];       /* X: category */
    SntModel dm[2][SNT_LEVEL_MAX_CAT + 1];   /* M: magnitude bits, by category */
    /* luma AC levels as binary decisions (levels.h), by block type */
    SntModel acz[2][64][SNT_PRED_CLASSES];  /* Z: nonzero?, by zigzag index, prediction class */
    SntModel acs[2];                   /* S: negative? */
    SntModel acg[2][64][SNT_PRED_CLASSES];  /* G: |v| > 1?, by zigzag index, prediction class */
    SntModel acx[2][2][SNT_LEVEL_MAX_CAT];      /* X: category, by low/high index */
    SntModel acm[2][SNT_LEVEL_MAX_CAT + 1];     /* M: magnitude bits, by category */
    SntModel band[2][SNT_NUM_BANDS][2][SNT_PRED_CLASSES];  /* band coded (1) or all zero (0);
                                          by previous band nonzero, prediction class */
    SntModel imode[2];                 /* intra mode, by frame type (0 key, 1 inter) */
    SntModel mbmode[SNT_FLAG_CONTEXTS];  /* macroblock mode (SNT_MB_*), by snt_flag_context */
    SntModel bdir;                     /* B-frame prediction direction - 1 (SNT_DIR_*), or
                                          SNT_BDIR_DIRECT */
    SntModel ref1[2];                  /* P-frame inter block: from a reference other than the
                                          newest? by whether the left block is */
    SntModel ref2[2];                  /* ... the third (the key frame) rather than the second?
                                          by whether the left block uses the third */
    SntModel inter[SNT_FLAG_CONTEXTS];   /* inter-frame block: inter (1) or intra, by snt_flag_context */
    SntModel iskip[SNT_FLAG_CONTEXTS];   /* inter block: skipped (1) or coded, by snt_flag_context */
    /* MV deltas as binary decisions (levels.h), by component (x, y); Z and G
     * also by block size (8x8, 16x16) and, for y, whether x's delta is zero */
    SntModel vz[2][2][2][SNT_AGREE_CLASSES], vs[2], vg[2][2][2][SNT_AGREE_CLASSES];
    SntModel vx[2][SNT_LEVEL_MAX_CAT];       /* X: category */
    SntModel vm[2][SNT_LEVEL_MAX_CAT + 1];   /* M: magnitude bits, by category */
    /* Chroma 4x4 blocks: a coded flag, then the 16 levels in zigzag order,
     * DC first, as binary decisions (levels.h), by block type */
    SntModel cflag[2][2][3];           /* coded (1) or all zero, by luma block coded,
                                          plane (U; V after an uncoded U; V after a coded U) */
    SntModel cz[2][16], cs[2], cg[2][16];   /* Z, S, G; Z and G by zigzag index */
    SntModel cx[2][2][SNT_LEVEL_MAX_CAT];   /* X: category, by DC or AC */
    SntModel cm[2][SNT_LEVEL_MAX_CAT + 1];  /* M: magnitude bits, by category */
} SntCdfSet;

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
        for (int a = 0; a < SNT_DC_CONTEXTS; a++) {
            snt_model_init(&s->dz[t][a], 2, SNT_RC_INIT_MODE, SNT_RC_INC_MODE, SNT_RC_LIMIT_MODE);
            snt_model_init(&s->ds[t][a], 2, SNT_RC_INIT_MODE, SNT_RC_INC_MODE, SNT_RC_LIMIT_MODE);
            for (int n = 0; n < 2; n++)
                snt_model_init(&s->dg[t][a][n], 2, SNT_RC_INIT_MODE, SNT_RC_INC_MODE, SNT_RC_LIMIT_MODE);
        }
        for (int i = 0; i < SNT_LEVEL_MAX_CAT; i++)
            snt_model_init(&s->dx[t][i], 2, SNT_RC_INIT_MODE, SNT_RC_INC_MODE, SNT_RC_LIMIT_MODE);
        for (int c = 0; c <= SNT_LEVEL_MAX_CAT; c++)
            snt_model_init(&s->dm[t][c], 2, SNT_RC_INIT_MODE, SNT_RC_INC_MODE, SNT_RC_LIMIT_MODE);
        snt_model_init(&s->imode[t], SNT_NUM_IPRED, SNT_RC_INIT_MODE, SNT_RC_INC_MODE, SNT_RC_LIMIT_MODE);
        for (int k = 0; k < 64; k++) {
            for (int p = 0; p < SNT_PRED_CLASSES; p++)
                snt_model_init(&s->acz[t][k][p], 2, SNT_RC_INIT_MODE, SNT_RC_INC_MODE, SNT_RC_LIMIT_MODE);
            for (int p = 0; p < SNT_PRED_CLASSES; p++)
                snt_model_init(&s->acg[t][k][p], 2, SNT_RC_INIT_MODE, SNT_RC_INC_MODE, SNT_RC_LIMIT_MODE);
        }
        snt_model_init(&s->acs[t], 2, SNT_RC_INIT_MODE, SNT_RC_INC_MODE, SNT_RC_LIMIT_MODE);
        for (int i = 0; i < SNT_LEVEL_MAX_CAT; i++)
            for (int h = 0; h < 2; h++)
                snt_model_init(&s->acx[t][h][i], 2, SNT_RC_INIT_MODE, SNT_RC_INC_MODE, SNT_RC_LIMIT_MODE);
        for (int c = 0; c <= SNT_LEVEL_MAX_CAT; c++)
            snt_model_init(&s->acm[t][c], 2, SNT_RC_INIT_MODE, SNT_RC_INC_MODE, SNT_RC_LIMIT_MODE);
        for (int l = 0; l < 2; l++)
            for (int c = 0; c < 3; c++)
                snt_model_init(&s->cflag[t][l][c], 2, SNT_RC_INIT_MODE, SNT_RC_INC_MODE, SNT_RC_LIMIT_MODE);
        for (int k = 0; k < 16; k++) {
            snt_model_init(&s->cz[t][k], 2, SNT_RC_INIT_MODE, SNT_RC_INC_MODE, SNT_RC_LIMIT_MODE);
            snt_model_init(&s->cg[t][k], 2, SNT_RC_INIT_MODE, SNT_RC_INC_MODE, SNT_RC_LIMIT_MODE);
        }
        snt_model_init(&s->cs[t], 2, SNT_RC_INIT_MODE, SNT_RC_INC_MODE, SNT_RC_LIMIT_MODE);
        for (int i = 0; i < SNT_LEVEL_MAX_CAT; i++)
            for (int h = 0; h < 2; h++)
                snt_model_init(&s->cx[t][h][i], 2, SNT_RC_INIT_MODE, SNT_RC_INC_MODE, SNT_RC_LIMIT_MODE);
        for (int c = 0; c <= SNT_LEVEL_MAX_CAT; c++)
            snt_model_init(&s->cm[t][c], 2, SNT_RC_INIT_MODE, SNT_RC_INC_MODE, SNT_RC_LIMIT_MODE);
        for (int b = 0; b < SNT_NUM_BANDS; b++) {
            for (int z = 0; z < 2; z++)
                for (int p = 0; p < SNT_PRED_CLASSES; p++)
                    snt_model_init(&s->band[t][b][z][p], 2, SNT_RC_INIT_MODE, SNT_RC_INC_MODE, SNT_RC_LIMIT_MODE);
        }
    }
    for (int a = 0; a < SNT_FLAG_CONTEXTS; a++) {
        snt_model_init(&s->mbmode[a], SNT_NUM_MB_MODES, SNT_RC_INIT_MODE, SNT_RC_INC_MODE, SNT_RC_LIMIT_MODE);
        snt_model_init(&s->inter[a], 2, SNT_RC_INIT_MODE, SNT_RC_INC_MODE, SNT_RC_LIMIT_MODE);
        snt_model_init(&s->iskip[a], 2, SNT_RC_INIT_MODE, SNT_RC_INC_MODE, SNT_RC_LIMIT_MODE);
    }
    snt_model_init(&s->bdir, 4, SNT_RC_INIT_MODE, SNT_RC_INC_MODE, SNT_RC_LIMIT_MODE);
    for (int r = 0; r < 2; r++) {
        snt_model_init(&s->ref1[r], 2, SNT_RC_INIT_MODE, SNT_RC_INC_MODE, SNT_RC_LIMIT_MODE);
        snt_model_init(&s->ref2[r], 2, SNT_RC_INIT_MODE, SNT_RC_INC_MODE, SNT_RC_LIMIT_MODE);
    }
    for (int c = 0; c < 2; c++) {
        for (int z = 0; z < 2; z++)
            for (int x = 0; x < 2; x++)
                for (int a = 0; a < SNT_AGREE_CLASSES; a++) {
                    snt_model_init(&s->vz[c][z][x][a], 2, SNT_RC_INIT_MODE, SNT_RC_INC_MODE, SNT_RC_LIMIT_MODE);
                    snt_model_init(&s->vg[c][z][x][a], 2, SNT_RC_INIT_MODE, SNT_RC_INC_MODE, SNT_RC_LIMIT_MODE);
                }
        snt_model_init(&s->vs[c], 2, SNT_RC_INIT_MODE, SNT_RC_INC_MODE, SNT_RC_LIMIT_MODE);
        for (int i = 0; i < SNT_LEVEL_MAX_CAT; i++)
            snt_model_init(&s->vx[c][i], 2, SNT_RC_INIT_MODE, SNT_RC_INC_MODE, SNT_RC_LIMIT_MODE);
        for (int i = 0; i <= SNT_LEVEL_MAX_CAT; i++)
            snt_model_init(&s->vm[c][i], 2, SNT_RC_INIT_MODE, SNT_RC_INC_MODE, SNT_RC_LIMIT_MODE);
    }
}

#ifndef SNT_NO_TRAINED_TABLES
#include "model_tables.h"
#endif

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
#ifdef SNT_NO_TRAINED_TABLES
    (void)set; (void)qp;  /* training builds start uniform */
#else
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
#endif
}

/* Model for decision b of a luma AC level at zigzag index k with
 * prediction class pc (levels.h) */
static inline SntModel *snt_level_model(SntCdfSet *s, int t, int k, int pc, const SntLevelBin *b) {
    switch (b->kind) {
    case 0:  return &s->acz[t][k][pc];
    case 1:  return &s->acs[t];
    case 2:  return &s->acg[t][k][pc];
    case 3:  return &s->acx[t][k > SNT_LEVEL_KLOW][b->sub];
    default: return &s->acm[t][b->sub];
    }
}

/* Context of a mode flag (macroblock mode, inter, skip): the prediction
 * agreement class ac, and how many of the blocks left of and above the
 * (macro)block's top-left 8x8 block are skipped (0-2; a missing neighbour
 * counts as not skipped). A skip map is a binary image, which JBIG (ITU-T
 * T.82, 1993) codes with an adaptive model chosen by a template of
 * neighbouring pixels. */
static inline int snt_flag_context(int ac, int left_skip, int above_skip) {
    return ac * 3 + left_skip + above_skip;
}

/* Class of a DC level for the next block's DC context, as JPEG classifies
 * the previous DC difference (T.81 F.1.4.4.1, with its default bounds
 * L = 0, U = 1): zero, small positive (1, 2), small negative, large
 * positive, large negative */
static inline int snt_dc_class(int v) {
    if (v == 0) return 0;
    if (v > 0) return v <= 2 ? 1 : 3;
    return v >= -2 ? 2 : 4;
}

/* DC context of a luma block: the class of the left block's DC level,
 * and whether the above block's is nonzero (0 for a missing neighbour) */
static inline int snt_dc_context(int left, int above) {
    return snt_dc_class(left) + 5 * (above != 0);
}

/* Model for decision b of a luma DC level v in DC context da */
static inline SntModel *snt_dclevel_model(SntCdfSet *s, int t, int da, int16_t v, const SntLevelBin *b) {
    switch (b->kind) {
    case 0:  return &s->dz[t][da];
    case 1:  return &s->ds[t][da];
    case 2:  return &s->dg[t][da][v < 0];
    case 3:  return &s->dx[t][b->sub];
    default: return &s->dm[t][b->sub];
    }
}

/* Model for decision b of an MV-delta component c (0 x, 1 y) of a block of
 * size16 (0 8x8, 1 16x16); xz: for y, whether x's delta is zero */
static inline SntModel *snt_mvd_model(SntCdfSet *s, int c, int size16, int xz, int ac,
                                      const SntLevelBin *b) {
    switch (b->kind) {
    case 0:  return &s->vz[c][size16][xz][ac];
    case 1:  return &s->vs[c];
    case 2:  return &s->vg[c][size16][xz][ac];
    case 3:  return &s->vx[c][b->sub];
    default: return &s->vm[c][b->sub];
    }
}

/* Model for decision b of a chroma level at zigzag index i (levels.h) */
static inline SntModel *snt_clevel_model(SntCdfSet *s, int t, int i, const SntLevelBin *b) {
    switch (b->kind) {
    case 0:  return &s->cz[t][i];
    case 1:  return &s->cs[t];
    case 2:  return &s->cg[t][i];
    case 3:  return &s->cx[t][i > 0][b->sub];
    default: return &s->cm[t][b->sub];
    }
}

#endif /* SNT_CDF_SET_H */
