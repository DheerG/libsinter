/*
 * Sinter Codec — Motion Estimation
 * Implements Section 1.5: Inter Prediction (encoder-side search)
 *
 * Full-search integer-pel ME with quarter-pel refinement.
 * Diamond search for speed presets > 0 (future optimization).
 * MV prediction: median of three neighbours (H.263, see motion_est.h).
 */

#include <stdlib.h>
#include <limits.h>
#include "motion_est.h"
#include "mc.h"

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
int snte_mvd_bits(int d) {
    uint32_t u = (d >= 0) ? 2u * (uint32_t)d : (uint32_t)(-2 * d - 1);
    int prefix = 0;
    for (uint32_t c = u + 1; c > 1; c >>= 1) prefix++;
    return 2 * prefix + 1;
}

SntMV snte_motion_search(const uint8_t *cur, const uint8_t *ref,
                          int stride, int bx, int by,
                          int blk_size, int search_range,
                          int frame_w, int frame_h,
                          SntMV pred, double lambda) {
    SntMV best_mv = {0, 0};
    double best_cost = 1e300;

    const uint8_t *cur_blk = &cur[by * stride + bx];

    /* Candidate 0 is the predicted MV (rounded to integer-pel), so a
     * far-away predictor is still considered; the window then covers
     * [-search_range, +search_range]. */
    for (int i = -1; i < (2 * search_range + 1) * (2 * search_range + 1); i++) {
        int dx = (i < 0) ? pred.x / 4 : i % (2 * search_range + 1) - search_range;
        int dy = (i < 0) ? pred.y / 4 : i / (2 * search_range + 1) - search_range;
        int rx = bx + dx;
        int ry = by + dy;

        /* Bounds check */
        if (rx < 0 || ry < 0 ||
            rx + blk_size > frame_w || ry + blk_size > frame_h)
            continue;

        const uint8_t *ref_blk = &ref[ry * stride + rx];
        double cost = snt_sad(cur_blk, stride, ref_blk, stride, blk_size, blk_size)
                    + lambda * (snte_mvd_bits(4 * dx - pred.x) + snte_mvd_bits(4 * dy - pred.y));

        if (cost < best_cost) {
            best_cost = cost;
            best_mv.x = (int16_t)(dx * 4); /* quarter-pel units */
            best_mv.y = (int16_t)(dy * 4);
        }
    }

    return best_mv;
}

/* ── Sub-pixel refinement (half-pel and quarter-pel) ──────────────── */
/* SAD of the in-frame part of the block against its MC prediction */
static uint32_t mc_sad(const uint8_t *cur, const SntSplinePlane *ref, int w, int h,
                       int bx, int by, int blk, SntMV mv) {
    uint8_t pred[SNT_MC_MAX_BLK * SNT_MC_MAX_BLK];
    snt_mc_luma(ref, bx, by, mv, blk, blk, pred);
    uint32_t sad = 0;
    for (int y = 0; y < blk && by + y < h; y++)
        for (int x = 0; x < blk && bx + x < w; x++) {
            int d = (int)cur[(by + y) * w + bx + x] - (int)pred[y * blk + x];
            sad += (uint32_t)(d < 0 ? -d : d);
        }
    return sad;
}

SntMV snte_subpel_refine(const uint8_t *cur, const SntSplinePlane *ref,
                          int frame_w, int frame_h, int bx, int by,
                          int blk_size, SntMV mv, SntMV pred,
                          double lambda, uint32_t *sad_out) {
    uint32_t best_sad = mc_sad(cur, ref, frame_w, frame_h, bx, by, blk_size, mv);
    double best_cost = best_sad + lambda * (snte_mvd_bits(mv.x - pred.x) +
                                            snte_mvd_bits(mv.y - pred.y));
    for (int step = 2; step >= 1; step /= 2) {
        SntMV center = mv;
        for (int dy = -step; dy <= step; dy += step)
            for (int dx = -step; dx <= step; dx += step) {
                if (dx == 0 && dy == 0) continue;
                SntMV c = { (int16_t)(center.x + dx), (int16_t)(center.y + dy) };
                uint32_t sad = mc_sad(cur, ref, frame_w, frame_h, bx, by, blk_size, c);
                double cost = sad + lambda * (snte_mvd_bits(c.x - pred.x) +
                                              snte_mvd_bits(c.y - pred.y));
                if (cost < best_cost) {
                    best_cost = cost;
                    best_sad = sad;
                    mv = c;
                }
            }
    }
    if (sad_out) *sad_out = best_sad;
    return mv;
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
    SntMV pred;
    pred.x = median3(left_mv->x, above_mv->x, above_right_mv->x);
    pred.y = median3(left_mv->y, above_mv->y, above_right_mv->y);
    return pred;
}
