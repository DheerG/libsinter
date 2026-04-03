/*
 * Sinter Codec — Fixed-Point Arithmetic Utilities
 * Implements Section 3: Fixed-Point Precision Spec
 *
 * Global rounding convention: round-half-up
 *   (x + (1 << (shift-1))) >> shift
 * This convention is used everywhere in the codec. No exceptions.
 */
#ifndef SNT_FIXED_POINT_H
#define SNT_FIXED_POINT_H

#include <stdint.h>

/* ── Fixed-point format tags ──────────────────────────────────────────
 * Q15: 16-bit fractional (filter coefficients, probabilities)
 * Q8:  8-bit fractional  (pixel values, blending weights)
 * Q7:  7-bit fractional  (interpolation filter coefficients, sum to 128)
 * Q6:  6-bit fractional  (OBMC blending weights, sum to 64)
 * Q16: 16-bit fractional (affine MV computation)
 * ──────────────────────────────────────────────────────────────────── */
typedef int16_t  q15_t;
typedef int16_t  q8_t;
typedef int16_t  q7_t;
typedef int16_t  q6_t;
typedef int32_t  q16_t;

/* ── Rounding macro (the ONE rounding rule) ───────────────────────── */
/* Round-half-up right shift for non-negative values.
 * For signed transform coefficients, use snt_rshift_round_signed. */
#define SNT_RSHIFT_ROUND(x, shift) \
    (((x) + (1 << ((shift) - 1))) >> (shift))

/* Symmetric rounding for signed values (Gary's note: avoid asymmetry
 * on negative transform coefficients). Rounds toward nearest, ties
 * away from zero. */
static inline int32_t snt_rshift_round_signed(int32_t x, int shift) {
    if (x >= 0) {
        return (x + (1 << (shift - 1))) >> shift;
    } else {
        return -(((-x) + (1 << (shift - 1))) >> shift);
    }
}

/* ── Fixed-point multiply helpers ─────────────────────────────────── */

/* Multiply two Q15 values, return Q15 result */
static inline q15_t snt_q15_mul(q15_t a, q15_t b) {
    int32_t product = (int32_t)a * (int32_t)b;
    return (q15_t)SNT_RSHIFT_ROUND(product, 15);
}

/* Multiply int16 coefficient by Q15 weight, return int16 result */
static inline int16_t snt_coeff_mul_q15(int16_t coeff, q15_t weight) {
    int32_t product = (int32_t)coeff * (int32_t)weight;
    return (int16_t)snt_rshift_round_signed(product, 15);
}

/* ── Saturation / clamping ────────────────────────────────────────── */

static inline int16_t snt_clip_i16(int32_t x) {
    if (x > INT16_MAX) return INT16_MAX;
    if (x < INT16_MIN) return INT16_MIN;
    return (int16_t)x;
}

static inline uint8_t snt_clip_pixel(int32_t x) {
    if (x > 255) return 255;
    if (x < 0) return 0;
    return (uint8_t)x;
}

static inline int32_t snt_clip_i32(int64_t x) {
    if (x > INT32_MAX) return INT32_MAX;
    if (x < INT32_MIN) return INT32_MIN;
    return (int32_t)x;
}

/* ── Bit manipulation ─────────────────────────────────────────────── */

/* Integer log2 (floor) for positive values */
static inline int snt_ilog2(uint32_t x) {
    int r = 0;
    while (x >>= 1) r++;
    return r;
}

/* Integer square root (floor) — used for PVQ gain computation.
 * Decoder-safe: no floating point. */
static inline uint32_t snt_isqrt(uint64_t x) {
    if (x == 0) return 0;
    uint64_t r = 1;
    uint64_t next = (r + x / r) >> 1;
    while (next < r) {
        r = next;
        next = (r + x / r) >> 1;
    }
    return (uint32_t)r;
}

#endif /* SNT_FIXED_POINT_H */
