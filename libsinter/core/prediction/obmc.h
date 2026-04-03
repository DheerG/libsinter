/*
 * Sinter Codec — Overlapped Block Motion Compensation (OBMC)
 * Implements Section 1.5: Inter Prediction
 */
#ifndef SNT_OBMC_H
#define SNT_OBMC_H

#include <stdint.h>
#include "../sinter_types.h"

/* OBMC blending weights: 6-bit fixed-point, sum to 64 at every pixel.
 * Overlap region is 4 pixels (power of 2 for SIMD). */
#define OBMC_OVERLAP    4
#define OBMC_WEIGHT_BITS 6
#define OBMC_WEIGHT_SUM  (1 << OBMC_WEIGHT_BITS)  /* 64 */

/* ── Raised-cosine OBMC weights ──────────────────────────────────── */
/* For a 4-pixel overlap, weights for the current block:
 * w[0]=16, w[1]=32, w[2]=48, w[3]=56
 * Neighbor weight = 64 - w[i]. Sum = 64 at every position. */
extern const uint8_t obmc_weights_cur[OBMC_OVERLAP];
extern const uint8_t obmc_weights_nbr[OBMC_OVERLAP];

/* ── Motion compensation ──────────────────────────────────────────── */

/* Sub-pixel interpolation: 8-tap luma, 4-tap chroma.
 * Horizontal-first, vertical-second.
 * Coefficients sum to 128 (Q7). Rounding: add 64, shift 7. */

/* Perform motion compensation for a single block.
 * ref: reference frame luma plane. ref_stride: stride of ref.
 * mv: motion vector (1/8-pel). dst: output prediction. dst_stride.
 * blk_w, blk_h: block dimensions. */
void snt_mc_luma_ref(const uint8_t *ref, int ref_stride,
                      SntMV mv, uint8_t *dst, int dst_stride,
                      int blk_w, int blk_h);

/* OBMC blend: blend current block prediction with neighbor predictions
 * at block boundaries. Uses raised-cosine weights.
 * cur: current block's MC prediction.
 * left/above: neighbor predictions (NULL if at frame edge).
 * dst: output blended prediction.
 * All buffers have the same stride and dimensions. */
void obmc_blend_ref(const uint8_t *cur, const uint8_t *left,
                     const uint8_t *above, uint8_t *dst,
                     int stride, int blk_w, int blk_h);

/* ── Checkerboard processing order ────────────────────────────────
 *
 * Even superblocks are processed first (independent of each other),
 * then odd superblocks (which can reference even neighbors).
 *
 *   +---+---+---+---+
 *   | E | O | E | O |    E = even (pass 1)
 *   +---+---+---+---+    O = odd  (pass 2)
 *   | O | E | O | E |
 *   +---+---+---+---+
 *   | E | O | E | O |
 *   +---+---+---+---+
 *
 * is_even_block(col, row) = ((col + row) % 2 == 0)
 */
static inline int obmc_is_even_block(int col, int row) {
    return ((col + row) % 2) == 0;
}

#endif /* SNT_OBMC_H */
