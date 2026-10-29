/*
 * Sinter Codec — Motion Estimation (Encoder-Side)
 * Implements inter prediction search for optimal motion vectors.
 */
#ifndef SNT_MOTION_EST_H
#define SNT_MOTION_EST_H

#include <stdint.h>
#include "../sinter_types.h"
#include "spline.h"

/* Full-search integer-pel motion estimation.
 * Searches [-search_range, +search_range] (plus the predicted MV) for the
 * MV minimizing SAD + lambda * bits(mv - pred), with bits approximated by
 * the exp-Golomb length of each zigzagged component. lambda = 0 gives a
 * plain minimum-SAD search.
 * cur: current frame block. ref: reference frame. stride: row stride.
 * bx, by: block position. blk_size: block dimension.
 * search_range: half-width of search window in pixels.
 * pred: predicted MV in quarter-pel units (as coded).
 * Returns: best integer-pel MV, in quarter-pel units. */
SntMV snte_motion_search(const uint8_t *cur, const uint8_t *ref,
                          int stride, int bx, int by,
                          int blk_size, int search_range,
                          int frame_w, int frame_h,
                          SntMV pred, double lambda);

/* Approximate bits to code one MV-delta component (exp-Golomb length). */
int snte_mvd_bits(int d);

/* Refine mv to quarter-pel: a half-pel step, then a quarter-pel step,
 * each trying the 8 neighbours under SAD + lambda * bits(mv - pred), with
 * predictions from snt_mc_luma. cur is a frame_w x frame_h plane with
 * stride frame_w. All MVs in quarter-pel units. Returns the refined
 * MV; *sad_out (if non-NULL) receives its SAD. */
SntMV snte_subpel_refine(const uint8_t *cur, const SntSplinePlane *ref,
                          int frame_w, int frame_h, int bx, int by,
                          int blk_size, SntMV mv, SntMV pred,
                          double lambda, uint32_t *sad_out);

/* Median of three neighbours' MVs, component by component. */
SntMV snt_mv_predict(const SntMV *left_mv, const SntMV *above_mv,
                      const SntMV *above_right_mv);

/* MV predictor of 8x8 block (bx, by), as in H.263 (section 6.1.1 and
 * Annex F, 1996): the median of three candidates, MV1 left, MV2 above
 * and MV3 above-right, except that within a 16x16 macroblock the top-left
 * block takes MV3 from the next macroblock (two blocks right) and the
 * bottom-right block takes it from the top-left block. At the left edge
 * MV1 is zero; at the top edge MV2 and MV3 are MV1; at the right edge MV3
 * is zero. Intra and skipped blocks hold zero MVs in mv_field, so they
 * count as zero candidates. A 16x16 macroblock's predictor is its
 * top-left block's. */
static inline SntMV snt_mv_predict_at(const SntMV *mv_field, int blk_cols, int bx, int by) {
    SntMV zero = {0, 0};
    int i = by * blk_cols + bx;
    int zx = bx & 1, zy = by & 1;
    SntMV mv1 = bx > 0 ? mv_field[i - 1] : zero;
    if (by == 0 && !zy)
        return mv1;  /* MV2 = MV3 = MV1: the median is MV1 */
    SntMV mv2 = mv_field[i - blk_cols];
    int x3 = zx == 0 && zy == 0 ? bx + 2 : (zx && zy ? bx - 1 : bx + 1);
    SntMV mv3 = x3 < blk_cols ? mv_field[i - blk_cols + (x3 - bx)] : zero;
    return snt_mv_predict(&mv1, &mv2, &mv3);
}

/* Compute SAD (Sum of Absolute Differences) for motion search. */
uint32_t snt_sad(const uint8_t *a, int a_stride,
                  const uint8_t *b, int b_stride,
                  int w, int h);

#endif /* SNT_MOTION_EST_H */
