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

#define SNT_MC_MAX_BLK 16  /* largest block these functions predict */

/* Predict a bw x bh luma block whose top-left is (x0, y0), displaced by
 * mv (quarter-pel). ref: w x h plane, stride w. dst: stride bw. */
void snt_mc_luma(const uint8_t *ref, int w, int h, int x0, int y0,
                 SntMV mv, int bw, int bh, uint8_t *dst);

/* Predict a bw x bh chroma block at (x0, y0), displaced by mv (1/8-pel).
 * ref: w x h plane, stride w. dst: stride bw. */
void snt_mc_chroma(const uint8_t *ref, int w, int h, int x0, int y0,
                   SntMV mv, int bw, int bh, uint8_t *dst);

/* DC of an 8x8 prediction block on the scale of the 8x8 orthonormal DCT's
 * DC coefficient (sum / 8). Neighbouring blocks predict their DC from it. */
static inline int16_t snt_pred_dc8(const uint8_t *blk) {
    int32_t s = 0;
    for (int i = 0; i < 64; i++) s += blk[i];
    return (int16_t)((s + 4) >> 3);
}

#endif /* SNT_MC_H */
