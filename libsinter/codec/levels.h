/*
 * Sinter Codec — Luma AC levels as binary decisions (normative)
 *
 * Each AC level v of a coded band is coded as a series of binary
 * decisions, each with its own adaptive model, as in the arithmetic-
 * coding mode of JPEG (ITU-T T.81 Annex F.1.4.4.2, 1992):
 *   Z   v != 0?                       model by block type and zigzag index
 *   S   v < 0?                        model by block type
 *   G   |v| > 1?                      model by block type and zigzag index
 *   X1, X2, ...  for a = |v| - 2 >= 0, unary decisions "a + 1 >= 2^i?"
 *                giving the category c (2^(c-1) <= a + 1 < 2^c)
 *                                     model by block type, low/high zigzag
 *                                     index and i
 *   M   the c - 1 bits of a + 1 below its top bit, high to low
 *                                     model by block type and c
 * Changes from the source (see PROVENANCE.md): JPEG keys Z and G to the
 * zigzag index only through its statistical bins and codes an end-of-
 * block decision; here a band's coded flag replaces end-of-block, and
 * when every earlier level of a coded band is zero the last one is known
 * to be nonzero, so its Z decision is not coded. The models are
 * Witten, Neal & Cleary's counts like every other model.
 */
#ifndef SNT_LEVELS_H
#define SNT_LEVELS_H

#include <stdint.h>

#define SNT_LEVEL_MAX_CAT 15   /* a + 1 < 2^15 */
#define SNT_LEVEL_KLOW    5    /* zigzag index <= 5: low-frequency X models */

/* The decisions of one level, in coding order; returns their count.
 * kind[]: 0 Z, 1 S, 2 G, 3 X (sub = i - 1), 4 M (sub = c). */
typedef struct {
    uint8_t kind, sub, bit;
} SntLevelBin;

static inline int snt_level_bins(int16_t v, int z_known, SntLevelBin *out) {
    int n = 0;
    if (!z_known)
        out[n++] = (SntLevelBin){ 0, 0, (uint8_t)(v != 0) };
    if (v == 0) return n;
    out[n++] = (SntLevelBin){ 1, 0, (uint8_t)(v < 0) };
    uint32_t mag = (uint32_t)(v < 0 ? -v : v);
    out[n++] = (SntLevelBin){ 2, 0, (uint8_t)(mag > 1) };
    if (mag == 1) return n;
    uint32_t a1 = mag - 1;  /* a + 1 >= 1 */
    int c = 1;
    while (c < SNT_LEVEL_MAX_CAT && a1 >= (1u << c)) {
        out[n++] = (SntLevelBin){ 3, (uint8_t)(c - 1), 1 };
        c++;
    }
    if (c < SNT_LEVEL_MAX_CAT)
        out[n++] = (SntLevelBin){ 3, (uint8_t)(c - 1), 0 };
    for (int i = c - 2; i >= 0; i--)
        out[n++] = (SntLevelBin){ 4, (uint8_t)c, (uint8_t)((a1 >> i) & 1) };
    return n;
}

/* Contexts from the prediction (inter blocks): the DCT of the motion-
 * compensated prediction block, compared with the quantizer step, per
 * zigzag index and per band. A residual is roughly the prediction's
 * misalignment times its gradient, so it tends to be nonzero where the
 * prediction has energy. Classes 0, 1, 2; intra blocks use 0. */
#define SNT_PRED_CLASSES 3

static inline void snt_pred_classes(const int16_t *pred_dct, uint16_t step, const uint8_t *zigzag,
                                    const int *band_bounds, int num_bands,
                                    uint8_t *kcls, uint8_t *bcls) {
    kcls[0] = 0;
    for (int k = 1; k < 64; k++) {
        int e = pred_dct[zigzag[k]];
        if (e < 0) e = -e;
        kcls[k] = (uint8_t)(e < step ? 0 : (e < 4 * step ? 1 : 2));
    }
    for (int b = 0; b < num_bands; b++) {
        int sum = 0;
        for (int k = band_bounds[b]; k < band_bounds[b + 1]; k++) {
            int e = pred_dct[zigzag[1 + k]];
            sum += e < 0 ? -e : e;
        }
        bcls[b] = (uint8_t)(sum < 2 * step ? 0 : (sum < 8 * step ? 1 : 2));
    }
}

#endif /* SNT_LEVELS_H */
