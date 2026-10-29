/*
 * Sinter Codec — Global motion estimation (encoder-only)
 *
 * The format carries one global motion vector per reference and frame
 * (translational global motion compensation, as MPEG-4 Part 2's GMC with
 * one warping point); skipped blocks are predicted with it. How the
 * encoder finds it is not normative.
 */
#ifndef SNT_GLOBAL_MOTION_H
#define SNT_GLOBAL_MOTION_H

#include <stdint.h>
#include "../sinter_types.h"
#include "spline.h"

/* Estimate the dominant translation of cur (w x h, stride w) relative to
 * ref, in quarter-pel units: Lucas-Kanade on a three-level image pyramid
 * with outlier rejection, then the quarter-pel candidate near it (or
 * zero) under which the most 8x8 blocks predict well, by the sum over
 * blocks of min(SSE, skip_sse). */
SntMV snte_global_motion(const uint8_t *cur, const SntSplinePlane *ref, int w, int h,
                         double skip_sse);

#endif /* SNT_GLOBAL_MOTION_H */
