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

/* MV prediction: median of spatial neighbors (AMVP-style).
 * left_mv, above_mv, above_right_mv: neighbor MVs (NULL if unavailable).
 * Returns predicted MV. */
SntMV snt_mv_predict(const SntMV *left_mv, const SntMV *above_mv,
                      const SntMV *above_right_mv);

/* Median MV predictor for block (bx, by) from its left, above and
 * above-right neighbours in mv_field, each used only if coded[] says it
 * has already been coded (blocks are coded in macroblock Z order, so the
 * above-right neighbour may not be). */
static inline SntMV snt_mv_predict_at(const SntMV *mv_field, const uint8_t *coded,
                                      int blk_cols, int bx, int by) {
    int i = by * blk_cols + bx;
    const SntMV *l  = (bx > 0 && coded[i - 1]) ? &mv_field[i - 1] : NULL;
    const SntMV *a  = (by > 0 && coded[i - blk_cols]) ? &mv_field[i - blk_cols] : NULL;
    const SntMV *ar = (by > 0 && bx + 1 < blk_cols && coded[i - blk_cols + 1]) ?
                      &mv_field[i - blk_cols + 1] : NULL;
    return snt_mv_predict(l, a, ar);
}

/* Compute SAD (Sum of Absolute Differences) for motion search. */
uint32_t snt_sad(const uint8_t *a, int a_stride,
                  const uint8_t *b, int b_stride,
                  int w, int h);

#endif /* SNT_MOTION_EST_H */
