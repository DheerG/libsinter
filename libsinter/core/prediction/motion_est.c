/*
 * Sinter Codec — Motion Estimation
 * Implements Section 1.5: Inter Prediction (encoder-side search)
 *
 * Full-search integer-pel ME with 1/8-pel refinement.
 * Diamond search for speed presets > 0 (future optimization).
 * MV prediction via median of spatial neighbors (AMVP-style).
 */

#include <stdlib.h>
#include <limits.h>
#include "motion_est.h"

/* ── SAD computation ──────────────────────────────────────────────── */
uint32_t snt_sad(const uint8_t *a, int a_stride,
                  const uint8_t *b, int b_stride,
                  int w, int h) {
    uint32_t sad = 0;
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            int diff = (int)a[y * a_stride + x] - (int)b[y * b_stride + x];
            sad += (uint32_t)(diff < 0 ? -diff : diff);
        }
    }
    return sad;
}

/* ── Full-search integer-pel motion estimation ────────────────────── */
SntMV snte_motion_search(const uint8_t *cur, const uint8_t *ref,
                          int stride, int bx, int by,
                          int blk_size, int search_range,
                          int frame_w, int frame_h) {
    SntMV best_mv = {0, 0};
    uint32_t best_sad = UINT32_MAX;

    const uint8_t *cur_blk = &cur[by * stride + bx];

    for (int dy = -search_range; dy <= search_range; dy++) {
        for (int dx = -search_range; dx <= search_range; dx++) {
            int rx = bx + dx;
            int ry = by + dy;

            /* Bounds check */
            if (rx < 0 || ry < 0 ||
                rx + blk_size > frame_w || ry + blk_size > frame_h)
                continue;

            const uint8_t *ref_blk = &ref[ry * stride + rx];
            uint32_t sad = snt_sad(cur_blk, stride, ref_blk, stride,
                                    blk_size, blk_size);

            if (sad < best_sad) {
                best_sad = sad;
                best_mv.x = (int16_t)(dx * 8); /* convert to 1/8-pel */
                best_mv.y = (int16_t)(dy * 8);
            }
        }
    }

    return best_mv;
}

/* ── Sub-pixel refinement (half-pel and quarter-pel) ──────────────── */
SntMV snte_subpel_refine(const uint8_t *cur, const uint8_t *ref,
                          int stride, int bx, int by,
                          int blk_size, SntMV int_mv,
                          int frame_w, int frame_h) {
    /* For the reference encoder, just return the integer-pel MV.
     * Full sub-pixel refinement would interpolate at fractional
     * positions using snt_mc_luma_ref and pick the best SAD.
     * This is a placeholder — the sub-pixel filters exist in obmc.c
     * and would be called here for each candidate fractional position. */
    (void)cur; (void)ref; (void)stride; (void)bx; (void)by;
    (void)blk_size; (void)frame_w; (void)frame_h;
    return int_mv;
}

/* ── MV prediction (median of spatial neighbors) ──────────────────── */
static int16_t median3(int16_t a, int16_t b, int16_t c) {
    if (a > b) { int16_t t = a; a = b; b = t; }
    if (b > c) { int16_t t = b; b = c; c = t; }
    if (a > b) { int16_t t = a; a = b; b = t; }
    return b;
}

SntMV snt_mv_predict(const SntMV *left_mv, const SntMV *above_mv,
                      const SntMV *above_right_mv) {
    SntMV pred = {0, 0};
    int count = 0;
    int16_t mvx[3] = {0, 0, 0};
    int16_t mvy[3] = {0, 0, 0};

    if (left_mv) { mvx[count] = left_mv->x; mvy[count] = left_mv->y; count++; }
    if (above_mv) { mvx[count] = above_mv->x; mvy[count] = above_mv->y; count++; }
    if (above_right_mv) { mvx[count] = above_right_mv->x; mvy[count] = above_right_mv->y; count++; }

    if (count == 0) return pred;
    if (count == 1) { pred.x = mvx[0]; pred.y = mvy[0]; return pred; }
    if (count == 2) { pred.x = (mvx[0] + mvx[1] + 1) >> 1; pred.y = (mvy[0] + mvy[1] + 1) >> 1; return pred; }

    pred.x = median3(mvx[0], mvx[1], mvx[2]);
    pred.y = median3(mvy[0], mvy[1], mvy[2]);
    return pred;
}
