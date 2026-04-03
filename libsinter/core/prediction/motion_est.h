/*
 * Sinter Codec — Motion Estimation (Encoder-Side)
 * Implements inter prediction search for optimal motion vectors.
 */
#ifndef SNT_MOTION_EST_H
#define SNT_MOTION_EST_H

#include <stdint.h>
#include "../sinter_types.h"

/* Full-search integer-pel motion estimation.
 * Searches [-search_range, +search_range] for the best matching block.
 * cur: current frame block. ref: reference frame. stride: row stride.
 * bx, by: block position. blk_size: block dimension.
 * search_range: half-width of search window in pixels.
 * Returns: best MV in 1/8-pel units. */
SntMV snte_motion_search(const uint8_t *cur, const uint8_t *ref,
                          int stride, int bx, int by,
                          int blk_size, int search_range,
                          int frame_w, int frame_h);

/* Refine integer-pel MV to 1/8-pel using the existing sub-pixel filters.
 * Returns refined MV. */
SntMV snte_subpel_refine(const uint8_t *cur, const uint8_t *ref,
                          int stride, int bx, int by,
                          int blk_size, SntMV int_mv,
                          int frame_w, int frame_h);

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
