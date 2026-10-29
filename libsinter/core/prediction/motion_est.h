/*
 * Sinter Codec — Motion Estimation (Encoder-Side)
 * Implements inter prediction search for optimal motion vectors.
 */
#ifndef SNT_MOTION_EST_H
#define SNT_MOTION_EST_H

#include <stdint.h>
#include "../sinter_types.h"

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
 * predictions from snt_mc_luma. cur and ref are frame_w x frame_h planes
 * with stride frame_w. All MVs in quarter-pel units. Returns the refined
 * MV; *sad_out (if non-NULL) receives its SAD. */
SntMV snte_subpel_refine(const uint8_t *cur, const uint8_t *ref,
                          int frame_w, int frame_h, int bx, int by,
                          int blk_size, SntMV mv, SntMV pred,
                          double lambda, uint32_t *sad_out);

/* MV prediction: median of spatial neighbors (AMVP-style).
 * left_mv, above_mv, above_right_mv: neighbor MVs (NULL if unavailable).
 * Returns predicted MV. */
SntMV snt_mv_predict(const SntMV *left_mv, const SntMV *above_mv,
                      const SntMV *above_right_mv);

/* Compute SAD (Sum of Absolute Differences) for motion search. */
uint32_t snt_sad(const uint8_t *a, int a_stride,
                  const uint8_t *b, int b_stride,
                  int w, int h);

#endif /* SNT_MOTION_EST_H */
