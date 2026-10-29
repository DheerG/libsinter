/*
 * Sinter Codec — Per-clip starting tables (experiment)
 *
 * Phase 3, novel tool 3 (docs/NOVEL.md, section 3). Not part of the
 * format: sending starting probabilities in the header is VP8's design,
 * an excluded source. Built only with -DSNT_EXPERIMENT_CLIPTABLES, to
 * measure what it would have saved (tools/experiments/cliptables_eval.c).
 *
 * A first encoding pass counts each binary model's symbols over the clip.
 * The first key frame then carries, for each binary model of both model
 * sets, an update flag and, if set, a starting probability as a 5-bit
 * step (-15..15) in the logistic domain, 128/256 apart. Encoder and
 * decoder apply these starts, at the trained tables' total, whenever a
 * set restarts.
 */
#ifndef SNT_CLIPTABLES_H
#define SNT_CLIPTABLES_H

#include <stdint.h>
#include "cdf_set.h"
#include "mixer_tables.h"

#define SNT_CLIP_NONE  (-16)
#define SNT_CLIP_QSTEP 128

typedef struct {
    int8_t q[2][SNT_NUM_MODELS];  /* per set and model: SNT_CLIP_NONE or -15..15 */
} SntClipTables;

static inline uint32_t snt_clip_p1(int q) {  /* 1/4096ths */
    return snt_squash_tab[q * SNT_CLIP_QSTEP + 2048];
}

static inline void snt_clip_apply(SntCdfSet *s, const SntClipTables *ct, int set) {
    SntModel *m = (SntModel *)s;
    for (int k = 0; k < SNT_NUM_MODELS; k++) {
        if (m[k].num_symbols != 2 || ct->q[set][k] == SNT_CLIP_NONE) continue;
        uint32_t f1 = (snt_clip_p1(ct->q[set][k]) * m[k].total + 2048) >> 12;
        if (f1 < 1) f1 = 1;
        if (f1 > m[k].total - 1) f1 = m[k].total - 1;
        m[k].freq[1] = (uint16_t)f1;
        m[k].freq[0] = (uint16_t)(m[k].total - f1);
    }
}

/* Which models are binary (the same for every set) */
static inline int snt_clip_binary(int k) {
    static SntCdfSet sizes;
    static int ready;
    if (!ready) { snt_cdf_set_init_uniform(&sizes); ready = 1; }
    return ((SntModel *)&sizes)[k].num_symbols == 2;
}

static inline void snt_clip_encode(SntRangeEnc *e, const SntClipTables *ct) {
    for (int set = 0; set < 2; set++)
        for (int k = 0; k < SNT_NUM_MODELS; k++) {
            if (!snt_clip_binary(k)) continue;
            int q = ct->q[set][k];
            snt_rc_encode_bit(e, q != SNT_CLIP_NONE, 128);
            if (q != SNT_CLIP_NONE)
                for (int b = 4; b >= 0; b--)
                    snt_rc_encode_bit(e, ((q + 15) >> b) & 1, 2048);
        }
}

static inline void snt_clip_decode(SntRangeDec *d, SntClipTables *ct) {
    for (int set = 0; set < 2; set++)
        for (int k = 0; k < SNT_NUM_MODELS; k++) {
            ct->q[set][k] = SNT_CLIP_NONE;
            if (!snt_clip_binary(k) || !snt_rc_decode_bit(d, 128)) continue;
            int v = 0;
            for (int b = 4; b >= 0; b--) v = (v << 1) | snt_rc_decode_bit(d, 2048);
            ct->q[set][k] = (int8_t)(v > 30 ? 15 : v - 15);
        }
}

#endif /* SNT_CLIPTABLES_H */
