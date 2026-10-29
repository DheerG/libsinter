/*
 * Sinter Codec — Sub-pixel motion compensation (normative)
 *
 * Luma: cubic B-spline interpolation at quarter-pel positions (spline.c).
 *
 * Chroma: bilinear interpolation at 1/8-pel positions, as in MPEG-4
 * Part 2's sprite and global motion warping (1999), which interpolates
 * bilinearly at up to 1/16 sample. The chroma vector is the luma vector
 * halved, kept at 1/8-sample precision. See PROVENANCE.md.
 */

#include "mc.h"
#include "../../util/fixed_point.h"

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

void snt_mc_luma(const SntSplinePlane *ref, int x0, int y0,
                 SntMV mv, int bw, int bh, uint8_t *dst) {
    snt_spline_interp(ref, 4 * x0 + mv.x, 4 * y0 + mv.y, bw, bh, dst);
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

void snt_pred_luma8(const SntSplinePlane *ref_f, const SntSplinePlane *ref_b,
                    int x0, int y0, int dir, SntMV mvf, SntMV mvb, uint8_t *pred) {
    if (dir == SNT_DIR_BWD) {
        snt_mc_luma(ref_b, x0, y0, mvb, 8, 8, pred);
        return;
    }
    snt_mc_luma(ref_f, x0, y0, mvf, 8, 8, pred);
    if (dir == SNT_DIR_BI) {
        uint8_t b[64];
        snt_mc_luma(ref_b, x0, y0, mvb, 8, 8, b);
        for (int i = 0; i < 64; i++)
            pred[i] = (uint8_t)((pred[i] + b[i] + 1) >> 1);
    }
}

uint32_t snt_pred_sad8(const SntSplinePlane *ref_f, const SntSplinePlane *ref_b, int x0, int y0,
                       int dira, SntMV af, SntMV ab, int dirb, SntMV bf, SntMV bb) {
    uint8_t a[64], b[64];
    snt_pred_luma8(ref_f, ref_b, x0, y0, dira, af, ab, a);
    snt_pred_luma8(ref_f, ref_b, x0, y0, dirb, bf, bb, b);
    uint32_t sad = 0;
    for (int i = 0; i < 64; i++)
        sad += (uint32_t)(a[i] > b[i] ? a[i] - b[i] : b[i] - a[i]);
    return sad;
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
