/*
 * Sinter Codec — Logistic mixing of context models (experiment)
 *
 * Phase 3, novel tool 1 (docs/NOVEL.md, section 1). Not part of the
 * format: the patent search found the same idea (adaptive mixing of
 * several models' probabilities for coefficient coding). Built only with
 * -DSNT_EXPERIMENT_MIXING, to measure what it would have saved; streams
 * from such a build don't decode with a normal build.
 *
 * The zero (Z) and "above one" (G) decisions of luma AC levels are coded
 * with a probability mixed in the logistic domain from four predictors:
 *   0  the level's usual model (cdf_set.h, adaptive counts)
 *   1  a fast estimate per block type and zigzag index (shift 4)
 *   2  a slower estimate per block type, index, prediction class and the
 *      number of nonzero levels among the two before it (shift 6)
 *   3  the same by band instead of index (shift 5)
 * plus a bias. p = squash(sum w_i stretch(p_i)); after each bit the
 * weights move by stretch(p_i) * (bit - p), the error's gradient for the
 * coding cost, and every predictor adapts. All integer arithmetic.
 */
#ifndef SNT_MIXER_H
#define SNT_MIXER_H

#include <stdint.h>
#include <string.h>
#include "cdf_set.h"
#include "mixer_tables.h"

#define SNT_MIX_INPUTS 5
/* Weight learning rate 2^-10 (2^-8 .. 2^-16 tried) */
#ifndef SNT_MIX_LR_SHIFT
#define SNT_MIX_LR_SHIFT 10
#endif

typedef struct {
    uint16_t p1[2][2][64];                         /* [kind Z/G][t][k] */
    uint16_t p2[2][2][64][SNT_PRED_CLASSES][3];    /* [kind][t][k][pc][nzprev] */
    uint16_t p3[2][2][16][SNT_PRED_CLASSES][3];    /* [kind][t][band][pc][nzprev] */
    int32_t  w[2][2][8][SNT_MIX_INPUTS];           /* [kind][t][k / 8][input] */
} SntMixState;

static inline void snt_mix_init(SntMixState *m) {
    uint16_t *p = (uint16_t *)m;
    size_t np = (sizeof(m->p1) + sizeof(m->p2) + sizeof(m->p3)) / sizeof(uint16_t);
    for (size_t i = 0; i < np; i++) p[i] = 32768;
    for (int a = 0; a < 2; a++)
        for (int t = 0; t < 2; t++)
            for (int g = 0; g < 8; g++) {
                int32_t *w = m->w[a][t][g];
                w[0] = 65536; w[1] = 0; w[2] = 0; w[3] = 0; w[4] = 0;
            }
}

/* One mixed decision: its inputs and where they live */
typedef struct {
    int16_t   st[SNT_MIX_INPUTS];
    int32_t  *w;
    SntModel *m0;
    uint16_t *s1, *s2, *s3;
    uint32_t  p;          /* mixed probability of a 1, 1/4096ths */
} SntMixBin;

static inline int snt_mix_stretch16(uint16_t p16) { return snt_stretch_tab[p16 >> 4]; }

/* kind 0: Z (model acz), 1: G (model acg) */
static inline void snt_mix_prepare(SntMixState *ms, SntCdfSet *cs, int kind, int t, int k,
                                   int band, int pc, int nzprev, SntMixBin *b) {
    b->m0 = kind == 0 ? &cs->acz[t][k][pc] : &cs->acg[t][k][pc];
    b->s1 = &ms->p1[kind][t][k];
    b->s2 = &ms->p2[kind][t][k][pc][nzprev];
    b->s3 = &ms->p3[kind][t][band][pc][nzprev];
    b->w = ms->w[kind][t][k >> 3];
    uint32_t p0 = (uint32_t)(((uint64_t)b->m0->freq[1] << 12) / b->m0->total);
    b->st[0] = snt_stretch_tab[p0 < 1 ? 1 : (p0 > 4095 ? 4095 : p0)];
    b->st[1] = (int16_t)snt_mix_stretch16(*b->s1);
    b->st[2] = (int16_t)snt_mix_stretch16(*b->s2);
    b->st[3] = (int16_t)snt_mix_stretch16(*b->s3);
    b->st[4] = 256;
    int64_t dot = 0;
    for (int i = 0; i < SNT_MIX_INPUTS; i++) dot += (int64_t)b->w[i] * b->st[i];
    int32_t x = (int32_t)(dot >= 0 ? dot >> 16 : -((-dot) >> 16));
    if (x < -2047) x = -2047;
    if (x > 2047) x = 2047;
    b->p = snt_squash_tab[x + 2048];
}

static inline void snt_mix_update(const SntMixBin *b, int bit) {
    int32_t err = ((int32_t)bit << 12) - (int32_t)b->p;
    for (int i = 0; i < SNT_MIX_INPUTS; i++) {
        int32_t d = b->st[i] * err;
        b->w[i] += d >= 0 ? d >> SNT_MIX_LR_SHIFT : -((-d) >> SNT_MIX_LR_SHIFT);
    }
    int32_t target = bit ? 65535 : 0;
    *b->s1 = (uint16_t)(*b->s1 + ((target - *b->s1) >> 4));
    *b->s2 = (uint16_t)(*b->s2 + ((target - *b->s2) >> 6));
    *b->s3 = (uint16_t)(*b->s3 + ((target - *b->s3) >> 5));
    snt_model_update(b->m0, bit);
}

#endif /* SNT_MIXER_H */
