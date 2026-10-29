/*
 * Sinter Codec — Cubic B-spline interpolation of reference frames (normative)
 * Sources and changes: see spline.h.
 */

#include <stdlib.h>
#include "spline.h"
#include "mc.h"

/* Pole z = sqrt(3) - 2 of the cubic B-spline prefilter, and z / (z^2 - 1)
 * for the anticausal pass's starting value (mirror boundary), both Q16 */
#define Z_Q16   (-17560)
#define AC_Q16  18918
/* Samples summed for the causal pass's starting value (|z|^12 < 2^-22) */
#define INIT_TAPS 12

/* Floor division, defined for negative numerators */
static inline int64_t fdiv(int64_t a, int64_t d) {
    return a >= 0 ? a / d : -((-a + d - 1) / d);
}

static inline int64_t mul_q16(int64_t a, int64_t q16) {
    return fdiv(a * q16 + 32768, 65536);
}

static inline int clampi(int v, int lo, int hi) {
    return v < lo ? lo : (v > hi ? hi : v);
}

/* In-place cubic B-spline prefilter of x[0..n-1] (n >= 2): gain 6,
 * then a causal and an anticausal first-order recursion with pole z */
static void prefilter_1d(int64_t *x, int n) {
    for (int k = 0; k < n; k++) x[k] *= 6;
    int64_t acc = 0, zk = 65536;
    int taps = n < INIT_TAPS ? n : INIT_TAPS;
    for (int k = 0; k < taps; k++) {
        acc += zk * x[k];
        zk = mul_q16(zk, Z_Q16);
    }
    x[0] = fdiv(acc + 32768, 65536);
    for (int k = 1; k < n; k++)
        x[k] += mul_q16(x[k - 1], Z_Q16);
    x[n - 1] = mul_q16(x[n - 1] + mul_q16(x[n - 2], Z_Q16), AC_Q16);
    for (int k = n - 2; k >= 0; k--)
        x[k] = mul_q16(x[k + 1] - x[k], Z_Q16);
}

int32_t *snt_spline_coefs(const uint8_t *pix, int w, int h) {
    int wp = w + 2 * SNT_SPLINE_PAD, hp = h + 2 * SNT_SPLINE_PAD;
    int32_t *coef = (int32_t *)malloc((size_t)wp * hp * sizeof(int32_t));
    int64_t *line = (int64_t *)malloc((size_t)(wp > hp ? wp : hp) * sizeof(int64_t));
    if (!coef || !line) { free(coef); free(line); return NULL; }

    /* Rows of the edge-extended plane, samples in Q8 */
    for (int y = 0; y < hp; y++) {
        const uint8_t *row = pix + clampi(y - SNT_SPLINE_PAD, 0, h - 1) * w;
        for (int x = 0; x < wp; x++)
            line[x] = (int64_t)row[clampi(x - SNT_SPLINE_PAD, 0, w - 1)] << 8;
        prefilter_1d(line, wp);
        for (int x = 0; x < wp; x++)
            coef[y * wp + x] = (int32_t)line[x];
    }
    /* Then columns */
    for (int x = 0; x < wp; x++) {
        for (int y = 0; y < hp; y++)
            line[y] = coef[y * wp + x];
        prefilter_1d(line, hp);
        for (int y = 0; y < hp; y++)
            coef[y * wp + x] = (int32_t)line[y];
    }
    free(line);
    return coef;
}

/* Cubic B-spline weights for samples -1 .. +2 at t = 0, 1/4, 1/2, 3/4,
 * in 1/384ths: ((1 - t)^3, 4 - 6t^2 + 3t^3, 1 + 3t + 3t^2 - 3t^3, t^3) / 6 */
static const int32_t bspline_w[4][4] = {
    { 64, 256,  64,  0 },
    { 27, 235, 121,  1 },
    {  8, 184, 184,  8 },
    {  1, 121, 235, 27 },
};

void snt_spline_interp(const SntSplinePlane *p, int qx, int qy, int bw, int bh,
                       uint8_t *dst) {
    int fx = ((qx % 4) + 4) % 4, fy = ((qy % 4) + 4) % 4;
    int x0 = (qx - fx) / 4, y0 = (qy - fy) / 4;

    if (fx == 0 && fy == 0) {
        /* Whole-sample position: the spline equals the samples */
        for (int y = 0; y < bh; y++)
            for (int x = 0; x < bw; x++)
                dst[y * bw + x] = p->pix[clampi(y0 + y, 0, p->h - 1) * p->w +
                                         clampi(x0 + x, 0, p->w - 1)];
        return;
    }

    int wp = p->w + 2 * SNT_SPLINE_PAD, hp = p->h + 2 * SNT_SPLINE_PAD;
    const int32_t *wx = bspline_w[fx], *wy = bspline_w[fy];

    /* Horizontal pass over bh + 3 coefficient rows (one above, two
     * below), kept in Q8 */
    int32_t tmp[(SNT_MC_MAX_BLK + 3) * SNT_MC_MAX_BLK];
    for (int r = 0; r < bh + 3; r++) {
        const int32_t *row = p->coef +
            (size_t)clampi(y0 - 1 + r + SNT_SPLINE_PAD, 0, hp - 1) * wp;
        for (int x = 0; x < bw; x++) {
            int cx = x0 + x - 1 + SNT_SPLINE_PAD;
            int64_t s = 0;
            for (int t = 0; t < 4; t++)
                s += (int64_t)wx[t] * row[clampi(cx + t, 0, wp - 1)];
            tmp[r * bw + x] = (int32_t)fdiv(s + 192, 384);
        }
    }

    /* Vertical pass, Q8 -> pixel */
    for (int y = 0; y < bh; y++)
        for (int x = 0; x < bw; x++) {
            int64_t s = 0;
            for (int t = 0; t < 4; t++)
                s += (int64_t)wy[t] * tmp[(y + t) * bw + x];
            int64_t v = fdiv(s + 384 * 128, 384 * 256);
            dst[y * bw + x] = (uint8_t)(v < 0 ? 0 : (v > 255 ? 255 : v));
        }
}
