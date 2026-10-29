/*
 * Sinter Codec — Global motion estimation (encoder-only)
 *
 * B. D. Lucas and T. Kanade, "An iterative image registration technique
 * with an application to stereo vision", IJCAI 1981: Gauss-Newton steps
 * on the squared difference between the frame and the shifted reference,
 * coarse to fine over an image pyramid. Pixels whose residual is far
 * above the typical one (moving foreground objects) are left out of each
 * step, so the estimate follows the background.
 */

#include <stdlib.h>
#include <math.h>
#include "global_motion.h"

#define LEVELS  3
#define ITERS   10

typedef struct { float *p; int w, h; } Img;

static Img img_from_u8(const uint8_t *src, int w, int h) {
    Img m = { malloc(sizeof(float) * w * h), w, h };
    for (int i = 0; i < w * h; i++) m.p[i] = src[i];
    return m;
}

static Img img_half(const Img *a) {
    Img m = { malloc(sizeof(float) * (a->w / 2) * (a->h / 2)), a->w / 2, a->h / 2 };
    for (int y = 0; y < m.h; y++)
        for (int x = 0; x < m.w; x++) {
            const float *s = a->p + 2 * y * a->w + 2 * x;
            m.p[y * m.w + x] = 0.25f * (s[0] + s[1] + s[a->w] + s[a->w + 1]);
        }
    return m;
}

/* Bilinear sample with edge clamping */
static float sample(const Img *m, float x, float y) {
    if (x < 0) x = 0;
    if (y < 0) y = 0;
    if (x > m->w - 1.001f) x = m->w - 1.001f;
    if (y > m->h - 1.001f) y = m->h - 1.001f;
    int ix = (int)x, iy = (int)y;
    float fx = x - ix, fy = y - iy;
    const float *s = m->p + iy * m->w + ix;
    return (1 - fy) * ((1 - fx) * s[0] + fx * s[1]) + fy * ((1 - fx) * s[m->w] + fx * s[m->w + 1]);
}

/* Lucas-Kanade translation of cur against ref at one pyramid level,
 * starting from (*dx, *dy) */
static void lk_level(const Img *cur, const Img *ref, const Img *gx, const Img *gy,
                     double *dx, double *dy) {
    int margin = 2 + cur->w / 32;
    for (int it = 0; it < ITERS; it++) {
        /* Typical residual, for rejecting outliers */
        double sum_abs = 0;
        long n = 0;
        for (int y = margin; y < cur->h - margin; y++)
            for (int x = margin; x < cur->w - margin; x++) {
                float e = cur->p[y * cur->w + x] - sample(ref, x + (float)*dx, y + (float)*dy);
                sum_abs += fabsf(e);
                n++;
            }
        if (n == 0) return;
        double limit = 2.5 * (sum_abs / n) + 1.0;

        double a = 0, b = 0, c = 0, ex = 0, ey = 0;
        for (int y = margin; y < cur->h - margin; y++)
            for (int x = margin; x < cur->w - margin; x++) {
                float wx = x + (float)*dx, wy = y + (float)*dy;
                float e = cur->p[y * cur->w + x] - sample(ref, wx, wy);
                if (fabsf(e) > limit) continue;
                float u = sample(gx, wx, wy), v = sample(gy, wx, wy);
                a += u * u; b += u * v; c += v * v;
                ex += u * e; ey += v * e;
            }
        double det = a * c - b * b;
        if (det < 1e-6) return;
        double ux = (c * ex - b * ey) / det, uy = (a * ey - b * ex) / det;
        *dx += ux;
        *dy += uy;
        if (fabs(ux) < 0.01 && fabs(uy) < 0.01) return;
    }
}

static Img gradient(const Img *m, int horiz) {
    Img g = { calloc((size_t)m->w * m->h, sizeof(float)), m->w, m->h };
    for (int y = 1; y < m->h - 1; y++)
        for (int x = 1; x < m->w - 1; x++) {
            const float *s = m->p + y * m->w + x;
            g.p[y * m->w + x] = horiz ? 0.5f * (s[1] - s[-1]) : 0.5f * (s[m->w] - s[-m->w]);
        }
    return g;
}

/* Sum over full 8x8 blocks of min(SSE of the block predicted at mv, cap) */
static double robust_cost(const uint8_t *cur, const SntSplinePlane *ref, int w, int h,
                          SntMV mv, double cap) {
    double total = 0;
    uint8_t pred[64];
    for (int by = 0; by + 8 <= h; by += 8)
        for (int bx = 0; bx + 8 <= w; bx += 8) {
            snt_spline_interp(ref, 4 * bx + mv.x, 4 * by + mv.y, 8, 8, pred);
            double sse = 0;
            for (int y = 0; y < 8; y++)
                for (int x = 0; x < 8; x++) {
                    double d = (double)cur[(by + y) * w + bx + x] - pred[y * 8 + x];
                    sse += d * d;
                }
            total += sse < cap ? sse : cap;
        }
    return total;
}

SntMV snte_global_motion(const uint8_t *cur, const SntSplinePlane *ref, int w, int h,
                         double skip_sse) {
    SntMV zero = {0, 0};
    if (!ref || !ref->coef || w < 64 || h < 64) return zero;

    Img pc[LEVELS], pr[LEVELS];
    pc[0] = img_from_u8(cur, w, h);
    pr[0] = img_from_u8(ref->pix, w, h);
    for (int l = 1; l < LEVELS; l++) {
        pc[l] = img_half(&pc[l - 1]);
        pr[l] = img_half(&pr[l - 1]);
    }
    double dx = 0, dy = 0;
    for (int l = LEVELS - 1; l >= 0; l--) {
        Img gx = gradient(&pr[l], 1), gy = gradient(&pr[l], 0);
        lk_level(&pc[l], &pr[l], &gx, &gy, &dx, &dy);
        free(gx.p); free(gy.p);
        if (l > 0) { dx *= 2; dy *= 2; }
    }
    for (int l = 0; l < LEVELS; l++) { free(pc[l].p); free(pr[l].p); }

    /* Quarter-pel candidates around the estimate, and zero */
    SntMV best = zero;
    double best_cost = robust_cost(cur, ref, w, h, zero, skip_sse);
    if (!isfinite(dx) || !isfinite(dy) || fabs(dx) > 64 || fabs(dy) > 64) return best;
    int qx = (int)lround(4 * dx), qy = (int)lround(4 * dy);
    for (int oy = -1; oy <= 1; oy++)
        for (int ox = -1; ox <= 1; ox++) {
            SntMV c = { (int16_t)(qx + ox), (int16_t)(qy + oy) };
            if (c.x == 0 && c.y == 0) continue;
            double cost = robust_cost(cur, ref, w, h, c, skip_sse);
            if (cost < best_cost) { best_cost = cost; best = c; }
        }
    return best;
}
