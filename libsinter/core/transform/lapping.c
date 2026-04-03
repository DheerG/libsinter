/*
 * Sinter Codec — Lapped Transform Pre/Post Filters
 * Implements Section 1.1: Lapped Transforms (TDLT)
 *
 * Time-Domain Lapped Transforms use invertible pre/post filters at block
 * boundaries to decorrelate energy across blocks, structurally eliminating
 * blocking artifacts without any deblocking filter.
 *
 * Implemented as lifting steps for bit-exact reversibility:
 *   4-point: 2 lifting steps
 *   8-point: 4 lifting steps
 *
 * Pre-filter decorrelates before DCT; post-filter reconstructs after IDCT.
 * post_filter(pre_filter(x)) == x (bit-exact for integer input).
 *
 * Processing order: raster scan, horizontal before vertical.
 * Frame edges: identity (no lapping).
 * Tile boundaries: 4-pixel overlap with raised cosine blend (handled in tile.c).
 *
 * Reference: Malvar, "Signal Processing with Lapped Transforms" (1992).
 * Daala/Xiph implementation (BSD licensed).
 */

#include <stdint.h>
#include <string.h>
#include "lapping.h"
#include "../../util/fixed_point.h"

/* ── 4-point pre-filter (2 lifting steps) ─────────────────────────
 *
 * Operates on 4 samples straddling a block boundary:
 *   [a, b | c, d]  where | is the boundary
 *
 * Lifting step 1 (butterfly):
 *   b' = b + ((c * T0 + (1 << 14)) >> 15)
 *   c' = c - ((b' * T1 + (1 << 14)) >> 15)
 *
 * Lifting step 2 (scale):
 *   a' = a + ((b' * S0 + (1 << 14)) >> 15)
 *   d' = d + ((c' * S1 + (1 << 14)) >> 15)
 *
 * Coefficients derived from Type-IV DST factorization:
 *   T0 = 11585  (≈ 0.3536 * 32768, cos(pi/4)/sqrt(2))
 *   T1 = 11585
 *   S0 = 6270   (≈ 0.1914 * 32768)
 *   S1 = 6270
 */

/* Q15 lifting coefficients for 4-point filter.
 * Derived from Type-IV DST factorization (Malvar, 1992).
 * These are the theoretically-optimal values for maximum
 * cross-boundary decorrelation. Lifting guarantees bit-exact
 * reconstruction regardless of coefficient magnitude. */
/* Q15 lifting coefficients: 60% of Malvar theoretical values.
 * Full-strength Malvar (T0=11585) creates boundary energy that
 * the Q14 integer DCT amplifies through the post-filter chain,
 * limiting PSNR to ~22 dB even at low QP. At 60% strength,
 * we get meaningful deblocking while keeping the post-filter
 * error amplification manageable. */
#define LAP4_T0  6951    /* 0.60 * 11585 ≈ 0.212 in Q15 */
#define LAP4_T1  6951
#define LAP4_S0  3762    /* 0.60 * 6270 ≈ 0.115 in Q15 */
#define LAP4_S1  3762

/* 8-point filter at 60% of Timothy's reviewed values */
#define LAP8_T0  7724    /* 0.60 * 12873 */
#define LAP8_T1  7724
#define LAP8_T2  5320    /* 0.60 * 8867 */
#define LAP8_T3  5320
#define LAP8_S0  3459    /* 0.60 * 5765 */
#define LAP8_S1  3459
#define LAP8_S2  1881    /* 0.60 * 3135 */
#define LAP8_S3  1881

/* ── Q15 lifting step helper ─────────────────────────────────────── */
static inline int16_t lift(int16_t x, int16_t coeff) {
    /* Symmetric rounding (ties away from zero) for both positive and
     * negative products. Avoids asymmetric drift on negative values.
     * Both pre-filter and post-filter call this same function, so
     * bit-exact invertibility is preserved regardless of rounding mode. */
    return (int16_t)snt_rshift_round_signed((int32_t)x * coeff, 15);
}

/* ── 4-point pre-filter (forward, applied before DCT) ────────────── */
void snt_lap_pre_filter_4_ref(int16_t *data, int stride) {
    /* data points to first of 4 samples across the boundary:
     * data[0], data[stride], data[2*stride], data[3*stride]
     * Samples 0,1 are from the left block; 2,3 from the right block. */
    int16_t a = data[0];
    int16_t b = data[stride];
    int16_t c = data[2 * stride];
    int16_t d = data[3 * stride];

    /* Lifting step 1: butterfly */
    b = b + lift(c, LAP4_T0);
    c = c - lift(b, LAP4_T1);

    /* Lifting step 2: scale */
    a = a + lift(b, LAP4_S0);
    d = d + lift(c, LAP4_S1);

    data[0]          = a;
    data[stride]     = b;
    data[2 * stride] = c;
    data[3 * stride] = d;
}

/* ── 4-point post-filter (inverse, applied after IDCT) ───────────── */
void snt_lap_post_filter_4_ref(int16_t *data, int stride) {
    /* Exact inverse of pre-filter: undo lifting steps in reverse order */
    int16_t a = data[0];
    int16_t b = data[stride];
    int16_t c = data[2 * stride];
    int16_t d = data[3 * stride];

    /* Undo lifting step 2 (scale) */
    d = d - lift(c, LAP4_S1);
    a = a - lift(b, LAP4_S0);

    /* Undo lifting step 1 (butterfly) */
    c = c + lift(b, LAP4_T1);
    b = b - lift(c, LAP4_T0);

    data[0]          = a;
    data[stride]     = b;
    data[2 * stride] = c;
    data[3 * stride] = d;
}

/* ── 8-point pre-filter (4 lifting steps) ────────────────────────── */
void snt_lap_pre_filter_8_ref(int16_t *data, int stride) {
    /* 8 samples straddling boundary: [a,b,c,d | e,f,g,h] */
    int16_t s[8];
    for (int i = 0; i < 8; i++)
        s[i] = data[i * stride];

    /* Lifting step 1: inner butterfly */
    s[3] = s[3] + lift(s[4], LAP8_T0);
    s[4] = s[4] - lift(s[3], LAP8_T1);

    /* Lifting step 2: mid butterfly */
    s[2] = s[2] + lift(s[5], LAP8_T2);
    s[5] = s[5] - lift(s[2], LAP8_T3);

    /* Lifting step 3: outer scale */
    s[1] = s[1] + lift(s[3], LAP8_S0);
    s[6] = s[6] + lift(s[4], LAP8_S1);

    /* Lifting step 4: edge scale */
    s[0] = s[0] + lift(s[2], LAP8_S2);
    s[7] = s[7] + lift(s[5], LAP8_S3);

    for (int i = 0; i < 8; i++)
        data[i * stride] = s[i];
}

/* ── 8-point post-filter (inverse of pre-filter) ─────────────────── */
void snt_lap_post_filter_8_ref(int16_t *data, int stride) {
    int16_t s[8];
    for (int i = 0; i < 8; i++)
        s[i] = data[i * stride];

    /* Undo lifting step 4 */
    s[7] = s[7] - lift(s[5], LAP8_S3);
    s[0] = s[0] - lift(s[2], LAP8_S2);

    /* Undo lifting step 3 */
    s[6] = s[6] - lift(s[4], LAP8_S1);
    s[1] = s[1] - lift(s[3], LAP8_S0);

    /* Undo lifting step 2 */
    s[5] = s[5] + lift(s[2], LAP8_T3);
    s[2] = s[2] - lift(s[5], LAP8_T2);

    /* Undo lifting step 1 */
    s[4] = s[4] + lift(s[3], LAP8_T1);
    s[3] = s[3] - lift(s[4], LAP8_T0);

    for (int i = 0; i < 8; i++)
        data[i * stride] = s[i];
}

/* ── 16-point pre-filter (applies two cascaded 8-point filters) ───── */
/* For blocks larger than 8, we apply the 8-point filter at the boundary.
 * The lapping window size matches the smaller of block_size and 8.
 * Blocks >= 16 still use 8-point filters at their boundaries
 * (consistent with Daala's approach where filter size = min(block_size, 8)). */

/* ── Apply pre-filter to a row of block boundaries (horizontal) ──── */
void snt_lap_pre_filter_row_ref(int16_t *row, int width, int blk_size) {
    /* Apply pre-filter at each block boundary along the row.
     * Filter size = min(blk_size, 8). Larger blocks still use 8-point
     * filters — the lapping window doesn't grow beyond 8 samples. */
    int filter_size = blk_size < 8 ? 4 : 8;
    int half = filter_size / 2;

    for (int b = blk_size; b < width; b += blk_size) {
        int16_t *start = &row[b - half];
        if (filter_size == 4) {
            snt_lap_pre_filter_4_ref(start, 1);
        } else {
            snt_lap_pre_filter_8_ref(start, 1);
        }
    }
}

/* ── Apply post-filter to a row of block boundaries (horizontal) ── */
void snt_lap_post_filter_row_ref(int16_t *row, int width, int blk_size) {
    int filter_size = blk_size < 8 ? 4 : 8;
    int half = filter_size / 2;

    for (int b = blk_size; b < width; b += blk_size) {
        int16_t *start = &row[b - half];
        if (filter_size == 4) {
            snt_lap_post_filter_4_ref(start, 1);
        } else {
            snt_lap_post_filter_8_ref(start, 1);
        }
    }
}

/* ── Apply lapping to a 2D block grid (horizontal then vertical) ── */
void snt_lap_pre_filter_2d_ref(int16_t *buf, int width, int height,
                                int stride, int blk_size) {
    /* Step 1: horizontal pre-filter (all rows) */
    for (int y = 0; y < height; y++) {
        snt_lap_pre_filter_row_ref(&buf[y * stride], width, blk_size);
    }

    /* Step 2: vertical pre-filter (all columns at block boundaries) */
    int filter_size = blk_size < 8 ? 4 : 8;
    int half = filter_size / 2;

    for (int b = blk_size; b < height; b += blk_size) {
        for (int x = 0; x < width; x++) {
            int16_t *col_start = &buf[(b - half) * stride + x];
            if (filter_size == 4) {
                snt_lap_pre_filter_4_ref(col_start, stride);
            } else {
                snt_lap_pre_filter_8_ref(col_start, stride);
            }
        }
    }
}

void snt_lap_post_filter_2d_ref(int16_t *buf, int width, int height,
                                 int stride, int blk_size) {
    /* Inverse order: vertical post-filter first, then horizontal.
     * This is the exact inverse of the 2D pre-filter sequence. */
    int filter_size = blk_size < 8 ? 4 : 8;
    int half = filter_size / 2;

    /* Step 1: vertical post-filter */
    for (int b = blk_size; b < height; b += blk_size) {
        for (int x = 0; x < width; x++) {
            int16_t *col_start = &buf[(b - half) * stride + x];
            if (filter_size == 4) {
                snt_lap_post_filter_4_ref(col_start, stride);
            } else {
                snt_lap_post_filter_8_ref(col_start, stride);
            }
        }
    }

    /* Step 2: horizontal post-filter */
    for (int y = 0; y < height; y++) {
        snt_lap_post_filter_row_ref(&buf[y * stride], width, blk_size);
    }
}
