/*
 * Sinter Codec — Sub-pixel motion compensation (normative)
 *
 * Luma: separable 6-tap Lanczos-3 interpolation (windowed sinc) at
 * quarter-pel positions. Taps are Q6, rounded from the kernel, and each
 * phase sums to 64. The horizontal pass keeps its Q6 result, so rounding
 * happens once, after the vertical pass.
 *
 * Chroma: bilinear interpolation at 1/8-pel positions.
 *
 * Both are textbook interpolators, not a codec-specific filter design.
 * The integer-pel case reproduces the reference sample exactly.
 */

#include "mc.h"
#include "../../util/fixed_point.h"

/* Lanczos-3 weights at t = 0, 1/4, 1/2, 3/4 for samples -2 .. +3 */
static const int16_t qpel_taps[4][6] = {
    { 0,  0, 64,  0,  0, 0 },
    { 2, -9, 57, 17, -4, 1 },
    { 2, -9, 39, 39, -9, 2 },
    { 1, -4, 17, 57, -9, 2 },
};

static inline int sample(const uint8_t *ref, int w, int h, int x, int y) {
    x = x < 0 ? 0 : (x >= w ? w - 1 : x);
    y = y < 0 ? 0 : (y >= h ? h - 1 : y);
    return ref[y * w + x];
}

/* Split an MV component into integer and fraction (floor division, so
 * negative MVs need no implementation-defined shifts). */
static inline void split_mv(int v, int denom, int *ip, int *frac) {
    *frac = ((v % denom) + denom) % denom;
    *ip = (v - *frac) / denom;
}

void snt_mc_luma(const uint8_t *ref, int w, int h, int x0, int y0,
                 SntMV mv, int bw, int bh, uint8_t *dst) {
    int ix, fx, iy, fy;
    split_mv(mv.x, 4, &ix, &fx);
    split_mv(mv.y, 4, &iy, &fy);

    /* Horizontal pass over bh + 5 rows (two above, three below), Q6 */
    int16_t tmp[(SNT_MC_MAX_BLK + 5) * SNT_MC_MAX_BLK];
    for (int y = 0; y < bh + 5; y++) {
        int sy = y0 + iy + y - 2;
        for (int x = 0; x < bw; x++) {
            int sx = x0 + ix + x - 2;
            int32_t s = 0;
            for (int t = 0; t < 6; t++)
                s += qpel_taps[fx][t] * sample(ref, w, h, sx + t, sy);
            tmp[y * bw + x] = (int16_t)s;
        }
    }

    /* Vertical pass, Q12 -> pixel */
    for (int y = 0; y < bh; y++)
        for (int x = 0; x < bw; x++) {
            int32_t s = 0;
            for (int t = 0; t < 6; t++)
                s += qpel_taps[fy][t] * tmp[(y + t) * bw + x];
            dst[y * bw + x] = (s < 0) ? 0 : snt_clip_pixel((s + 2048) >> 12);
        }
}

void snt_mc_chroma(const uint8_t *ref, int w, int h, int x0, int y0,
                   SntMV mv, int bw, int bh, uint8_t *dst) {
    int ix, fx, iy, fy;
    split_mv(mv.x, 8, &ix, &fx);
    split_mv(mv.y, 8, &iy, &fy);

    for (int y = 0; y < bh; y++)
        for (int x = 0; x < bw; x++) {
            int sx = x0 + ix + x, sy = y0 + iy + y;
            int a = sample(ref, w, h, sx, sy),     b = sample(ref, w, h, sx + 1, sy);
            int c = sample(ref, w, h, sx, sy + 1), d = sample(ref, w, h, sx + 1, sy + 1);
            dst[y * bw + x] = (uint8_t)(((8 - fx) * (8 - fy) * a + fx * (8 - fy) * b +
                                         (8 - fx) * fy * c + fx * fy * d + 32) >> 6);
        }
}

void snt_pred_luma8(const uint8_t *ref_f, const uint8_t *ref_b, int w, int h,
                    int x0, int y0, int dir, SntMV mvf, SntMV mvb, uint8_t *pred) {
    if (dir == SNT_DIR_BWD) {
        snt_mc_luma(ref_b, w, h, x0, y0, mvb, 8, 8, pred);
        return;
    }
    snt_mc_luma(ref_f, w, h, x0, y0, mvf, 8, 8, pred);
    if (dir == SNT_DIR_BI) {
        uint8_t b[64];
        snt_mc_luma(ref_b, w, h, x0, y0, mvb, 8, 8, b);
        for (int i = 0; i < 64; i++)
            pred[i] = (uint8_t)((pred[i] + b[i] + 1) >> 1);
    }
}

void snt_pred_chroma4(const uint8_t *ref_f, const uint8_t *ref_b, int w, int h,
                      int x0, int y0, int dir, SntMV mvf, SntMV mvb, uint8_t *pred) {
    if (dir == SNT_DIR_BWD) {
        snt_mc_chroma(ref_b, w, h, x0, y0, mvb, 4, 4, pred);
        return;
    }
    snt_mc_chroma(ref_f, w, h, x0, y0, mvf, 4, 4, pred);
    if (dir == SNT_DIR_BI) {
        uint8_t b[16];
        snt_mc_chroma(ref_b, w, h, x0, y0, mvb, 4, 4, b);
        for (int i = 0; i < 16; i++)
            pred[i] = (uint8_t)((pred[i] + b[i] + 1) >> 1);
    }
}
