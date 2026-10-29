/*
 * Sinter Codec — Sub-pixel motion compensation (normative)
 *
 * Luma MVs are in quarter-pel units. Chroma (4:2:0, half resolution) uses
 * the same MV value, which is then in 1/8-pel chroma units.
 * Reference samples outside the frame are clamped to the nearest edge.
 */
#ifndef SNT_MC_H
#define SNT_MC_H

#include <stdint.h>
#include "../sinter_types.h"
#include "spline.h"

#define SNT_MC_MAX_BLK 16  /* largest block these functions predict */

/* Predict a bw x bh luma block whose top-left is (x0, y0), displaced by
 * mv (quarter-pel), from the B-spline of ref. dst: stride bw. */
void snt_mc_luma(const SntSplinePlane *ref, int x0, int y0,
                 SntMV mv, int bw, int bh, uint8_t *dst);

/* Predict a bw x bh chroma block at (x0, y0), displaced by mv (1/8-pel).
 * ref: w x h plane, stride w. dst: stride bw. */
void snt_mc_chroma(const uint8_t *ref, int w, int h, int x0, int y0,
                   SntMV mv, int bw, int bh, uint8_t *dst);

/* Inter prediction directions: forward (past reference), backward
 * (future reference, B-frames), or the rounded average of both. */
enum { SNT_DIR_FWD = 1, SNT_DIR_BWD = 2, SNT_DIR_BI = 3 };

/* Predict the 8x8 luma block at (x0, y0) in direction dir: forward from
 * ref_f at mvf, backward from ref_b at mvb, or (a + b + 1) >> 1. */
void snt_pred_luma8(const SntSplinePlane *ref_f, const SntSplinePlane *ref_b,
                    int x0, int y0, int dir, SntMV mvf, SntMV mvb, uint8_t *pred);

/* The same for a 4x4 chroma block, with the luma MVs read as 1/8-pel */
void snt_pred_chroma4(const uint8_t *ref_f, const uint8_t *ref_b, int w, int h,
                      int x0, int y0, int dir, SntMV mvf, SntMV mvb, uint8_t *pred);

/* DC of an 8x8 prediction block on the scale of the 8x8 orthonormal DCT's
 * DC coefficient (sum / 8). Neighbouring blocks predict their DC from it. */
static inline int16_t snt_pred_dc8(const uint8_t *blk) {
    int32_t s = 0;
    for (int i = 0; i < 64; i++) s += blk[i];
    return (int16_t)((s + 4) >> 3);
}

/* DC of a 4x4 chroma prediction block relative to 128, on the scale of
 * the 4x4 orthonormal DCT's DC coefficient (sum / 4). Chroma intra blocks
 * code pixels - 128, so neighbours compare like with like. */
static inline int16_t snt_pred_dc4c(const uint8_t *blk) {
    int32_t s = -128 * 16;
    for (int i = 0; i < 16; i++) s += blk[i];
    return (int16_t)(s >= 0 ? (s + 2) / 4 : -((-s + 2) / 4));
}

#endif /* SNT_MC_H */
