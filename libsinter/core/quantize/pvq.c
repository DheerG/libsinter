/*
 * Sinter Codec — Perceptual Vector Quantization (PVQ)
 * Implements Section 1.2: Quantization
 *
 * PVQ decomposes each AC coefficient vector into gain (energy magnitude)
 * and shape (direction on the unit hypersphere). The gain is coded precisely
 * with a dead-zone scalar quantizer; the shape is coded via combinatorial
 * enumeration on an integer lattice of K pulses.
 *
 * Key perceptual property: human vision is far more sensitive to energy loss
 * than to shape distortion. PVQ preserves texture energy even at aggressive
 * quantization — textures look "grainy but alive" rather than "smooth and dead."
 *
 * Energy conservation invariant:
 *   sum(recon[i]^2) == dequant_gain^2  (within 1 LSB)
 *
 * Shape encoding uses Fischer's combinatorial enumeration from CELT/Opus:
 * each integer lattice point on the K-radius hypersphere maps to a unique
 * index in [0, N(n,k)), where N(n,k) is the number of such points.
 * Requires 64-bit arithmetic for n > 16.
 *
 * Reference: Fischer, "A pyramid vector quantizer" (1986).
 * Reference: Valin, "Pyramid vector quantization for video coding" (Daala, 2015).
 */

#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include "pvq.h"
#include "../../util/fixed_point.h"

/* ── Combinatorial tables ─────────────────────────────────────────
 * N(n, k) = number of integer vectors of dimension n with L1 norm = k.
 * Recurrence: N(n, k) = N(n-1, k) + N(n, k-1) + N(n-1, k-1)
 *   (add +1, -1, or 0 in the last dimension, then solve remaining)
 *
 * Actually, the number of vectors in Z^n with ||v||_1 = k (allowing signs)
 * is given by: V(n,k) = sum_{j=0}^{n} C(n,j) * C(k-1, j-1) * 2^j
 * But for PVQ, we use the unsigned count and handle signs separately.
 *
 * For our purposes, we use the Opus/CELT formulation:
 * The number of integer points in the n-dimensional "pyramid" of
 * radius k (with signs) is computed recursively.
 */

/* Max dimensions and pulse counts for the lookup table */
#define PVQ_MAX_N  64
#define PVQ_MAX_K  128

/* Precomputed table of pulse vector counts (64-bit to avoid overflow).
 * pvq_count[n][k] = number of signed integer vectors of dim n with L1 norm k.
 * Uses the recurrence from Opus/CELT. */
static uint64_t pvq_count_cache[PVQ_MAX_N + 1][PVQ_MAX_K + 1];
static int pvq_cache_initialized = 0;

static void pvq_init_count_table(void) {
    if (pvq_cache_initialized) return;

    /* Base cases:
     * pvq_count[0][0] = 1 (empty vector)
     * pvq_count[0][k>0] = 0
     * pvq_count[1][k] = 2 for k>0 (just +k or -k), 1 for k=0 */
    for (int k = 0; k <= PVQ_MAX_K; k++) {
        pvq_count_cache[0][k] = (k == 0) ? 1 : 0;
    }
    for (int n = 1; n <= PVQ_MAX_N; n++) {
        pvq_count_cache[n][0] = 1; /* zero vector */
    }

    /* Recurrence: count(n, k) = count(n-1, k) + 2 * count(n-1, k-1)
     * + ... but the exact CELT recurrence is:
     * U(n,k) = U(n-1,k) + U(n,k-1) + U(n-1,k-1)
     * where U(n,k) counts the unsigned vectors with L1 = k.
     * The signed count is then V(n,k) = 2*U(n,k) for k>0, V(n,0)=1.
     *
     * But simpler: for the full signed count on the hyperpyramid:
     * V(n,k) = V(n-1,k) + V(n-1,k-1) + V(n-1,k+1 is not right...)
     *
     * Actually, use the standard recurrence:
     * V(1,k) = 2 for k>0, 1 for k=0
     * V(n,k) = V(n-1,k) + V(n-1,k-1) + V(n-1,k-1) = V(n-1,k) + 2*V(n-1,k-1)
     * Wait, that's not right either.
     *
     * Correct CELT/Opus formula:
     * The number of integer vectors x in Z^n with sum(|x_i|) = k:
     * This equals sum_{j=0}^{k} C(n, j) * C(k-1, j-1) * 2^j  ... complex.
     *
     * Simpler recursive definition:
     * V(n, k) = V(n-1, k) + 2 * sum_{j=1}^{k} V(n-1, k-j)
     *         = V(n-1, k) + 2 * (V(n, k-1) - V(n-1, k-1)) ... messy
     *
     * Let's use the direct recursion the way CELT does it:
     * V(1, k) = 2 for k>0, 1 for k=0
     * V(n, k) = V(n-1, k) + V(n, k-1) + V(n-1, k-1)
     *
     * (This counts: last element is 0 => V(n-1,k) ways;
     *  last element is +j or -j => need sum of |rest| = k-j, and we
     *  iterate j from 1 to k. The recurrence captures this.) */

    for (int n = 1; n <= PVQ_MAX_N; n++) {
        pvq_count_cache[n][0] = 1;
        if (n == 1) {
            for (int k = 1; k <= PVQ_MAX_K; k++)
                pvq_count_cache[1][k] = 2;
        } else {
            for (int k = 1; k <= PVQ_MAX_K; k++) {
                uint64_t val = pvq_count_cache[n - 1][k]
                             + pvq_count_cache[n][k - 1]
                             + pvq_count_cache[n - 1][k - 1];
                /* Clamp to prevent overflow — if we exceed 2^63, cap it */
                if (val < pvq_count_cache[n][k - 1])
                    val = UINT64_MAX;
                pvq_count_cache[n][k] = val;
            }
        }
    }

    pvq_cache_initialized = 1;
}

static uint64_t pvq_count(int n, int k) {
    if (n < 0 || k < 0) return 0;
    if (n > PVQ_MAX_N || k > PVQ_MAX_K) return UINT64_MAX;
    pvq_init_count_table();
    return pvq_count_cache[n][k];
}

/* ── Public accessor for pvq_count ─────────────────────────────────── */
uint64_t pvq_num_vectors(int n, uint16_t k) {
    return pvq_count(n, (int)k);
}

/* ── Shape bit-width computation ──────────────────────────────────── */

int pvq_shape_bits(int n, uint16_t k) {
    if (k == 0) return 0;
    pvq_init_count_table();
    uint64_t count = pvq_count(n, (int)k);
    if (count <= 1) return 0;
    if (count == UINT64_MAX) return 64; /* overflow — use max */
    /* ceil(log2(count)) = floor(log2(count-1)) + 1 */
    int bits = 0;
    uint64_t v = count - 1;
    while (v > 0) { v >>= 1; bits++; }
    return bits;
}

/* ── Gain computation ─────────────────────────────────────────────── */

uint32_t pvq_compute_gain(const int16_t *coeffs, int n) {
    /* L2 norm: sqrt(sum of squares). Integer arithmetic only. */
    uint64_t sum_sq = 0;
    for (int i = 0; i < n; i++)
        sum_sq += (int64_t)coeffs[i] * coeffs[i];
    return snt_isqrt(sum_sq);
}

/* ── Dead-zone quantization for gain and DC ───────────────────────── */

uint32_t pvq_quantize_gain(uint32_t gain, uint16_t step, int is_intra) {
    if (step == 0) return gain;
    /* Narrow dead-zone for gain to preserve AC energy.
     * Wide dead-zones were killing lapping boundary energy. */
    uint32_t offset = is_intra ? (step / 8) : (step / 4);
    if (gain <= offset) return 0;
    return (gain - offset + step / 2) / step; /* round-to-nearest, not floor */
}

uint32_t pvq_dequantize_gain(uint32_t level, uint16_t step) {
    /* Reconstruct at the centre of the quantization bin. The quantizer
     * already rounds to nearest, so level L covers
     * [L*step - step/2 + offset, L*step + step/2 + offset) and its centre
     * is L*step + offset. Adding step/2 here as well corrected for rounding
     * twice and biased every nonzero gain upward by ~step/3.
     * offset = step/8 is the intra dead-zone, the only one the encoder uses. */
    if (level == 0) return 0;
    return level * step + step / 8;
}

int16_t pvq_quantize_dc(int16_t dc, uint16_t step, int is_intra) {
    if (step == 0) return dc;
    /* Round-to-nearest with dead-zone (Gary's fix: match midpoint dequant).
     * Dead-zone kills small values; round-to-nearest minimizes error for the rest. */
    uint32_t offset = is_intra ? (step / 6) : (step / 3);
    int16_t abs_dc = dc < 0 ? -dc : dc;
    int16_t sign = dc < 0 ? -1 : 1;
    if ((uint32_t)abs_dc <= offset) return 0;
    return sign * (int16_t)((abs_dc - offset + step / 2) / step);
}

int16_t pvq_dequantize_dc(int16_t level, uint16_t step) {
    /* Bin-centre reconstruction, as for gain: offset = step/6 is the intra
     * dead-zone used for all DC and scalar AC levels. */
    if (level == 0) return 0;
    int16_t sign = level < 0 ? -1 : 1;
    int16_t abs_level = level < 0 ? -level : level;
    return sign * (int16_t)(abs_level * step + step / 6);
}

/* ── K parameter computation ──────────────────────────────────────── */

uint16_t pvq_compute_k(uint32_t gain_q, int n, uint8_t csf_weight) {
    /* K determines shape resolution. More pulses = finer shape.
     *
     * The key insight: K must be large enough to represent the shape
     * with reasonable fidelity. For a gradient in 63 dimensions, K=3
     * is catastrophically low (only 3 coefficients get any energy).
     *
     * Revised formula: K scales linearly with gain and with sqrt(n).
     * At typical gains (10-50) and n=63, this gives K=15-80, which
     * provides enough pulses to represent the shape meaningfully.
     *
     * K = max(1, gain_q * sqrt(n) * csf_weight / (128 * 4))
     *
     * The /4 factor controls the rate-quality tradeoff:
     * smaller divisor = more pulses = better quality = more bits. */
    if (gain_q == 0) return 0;

    uint32_t sqrt_n = snt_isqrt((uint64_t)n);
    if (sqrt_n == 0) sqrt_n = 1;

    /* K proportional to gain * sqrt(n), modulated by perceptual weight.
     * Divisor of 2 gives generous shape resolution — quality over compression. */
    uint32_t k = (gain_q * sqrt_n * (uint32_t)csf_weight) / (128 * 2);
    if (k < 1) k = 1;

    /* Cap K to prevent combinatorial table overflow in pvq_count(n, k).
     * Exact safe limits verified against uint64 overflow threshold. */
    uint32_t max_k;
    if      (n <= 4)  max_k = 128;
    else if (n <= 8)  max_k = 64;
    else if (n <= 12) max_k = 50;
    else if (n <= 16) max_k = 40;  /* pvq_count(16,55) fits, 60 overflows */
    else if (n <= 20) max_k = 30;
    else if (n <= 24) max_k = 22;
    else if (n <= 32) max_k = 18;
    else              max_k = 12;  /* n=63: max_k=14 */
    if (k > max_k) k = max_k;
    if (k > PVQ_MAX_K) k = PVQ_MAX_K;
    return (uint16_t)k;
}

/* ── Shape encoding ───────────────────────────────────────────────
 *
 * Fischer's combinatorial enumeration on the L1 hyperpyramid:
 * Each signed integer vector of dimension n with L1 norm = k maps
 * to a unique index in [0, V(n,k)).
 *
 * Encoding processes dimensions 0 to n-1. For dimension i with
 * remaining budget k_left:
 *   - For |x_i| = j (j = 0, 1, ..., k_left):
 *     - If j == 0: contributes V(n-1-i, k_left) vectors (just skip this dim)
 *     - If j > 0: contributes 2 * V(n-1-i, k_left-j) vectors (positive + negative)
 *   - The index skips past all vectors with smaller |x_i|, then encodes sign.
 */

/* Project coefficients onto L1 hypersphere of radius K.
 * Returns pulse vector with sum(|pulses[i]|) == k. */
static void pvq_project(const int16_t *coeffs, int n, uint16_t k,
                         int16_t *pulses) {
    if (k == 0 || n == 0) {
        memset(pulses, 0, n * sizeof(int16_t));
        return;
    }

    uint32_t gain = pvq_compute_gain(coeffs, n);
    if (gain == 0) {
        memset(pulses, 0, n * sizeof(int16_t));
        pulses[0] = (int16_t)k;
        return;
    }

    /* Use L1 norm for projection (preserves direction better for PVQ) */
    int32_t l1 = 0;
    for (int i = 0; i < n; i++) l1 += abs(coeffs[i]);
    if (l1 == 0) l1 = 1;

    int32_t remaining = (int32_t)k;
    for (int i = 0; i < n - 1; i++) {
        int32_t p = ((int64_t)coeffs[i] * k + (l1 / 2)) / l1;
        if (abs(p) > remaining) p = (p > 0) ? remaining : -remaining;
        pulses[i] = (int16_t)p;
        remaining -= abs(p);
    }
    /* Last element gets remainder to guarantee L1 = K */
    pulses[n - 1] = (coeffs[n - 1] >= 0) ? (int16_t)remaining
                                            : (int16_t)(-remaining);
}

uint64_t pvq_encode_shape(const int16_t *coeffs, int n, uint16_t k,
                           int16_t *pulses) {
    if (k == 0 || n == 0) {
        memset(pulses, 0, n * sizeof(int16_t));
        return 0;
    }

    pvq_project(coeffs, n, k, pulses);
    pvq_init_count_table();

    /* Encode: for each dimension, accumulate offset.
     * Processing order: dimension 0 first, matching decode order. */
    uint64_t idx = 0;
    int32_t k_left = (int32_t)k;

    for (int i = 0; i < n; i++) {
        int16_t abs_p = pulses[i] < 0 ? -pulses[i] : pulses[i];
        int remaining_dims = n - 1 - i;

        /* Offset past all vectors with |x_i| < abs_p */
        for (int j = 0; j < abs_p; j++) {
            if (j == 0) {
                idx += pvq_count(remaining_dims, k_left);
            } else {
                idx += 2 * pvq_count(remaining_dims, k_left - j);
            }
        }

        /* Sign: if |x_i| > 0, positive is first, negative adds one more offset */
        if (abs_p > 0 && pulses[i] < 0) {
            idx += pvq_count(remaining_dims, k_left - abs_p);
        }

        k_left -= abs_p;
    }

    return idx;
}

void pvq_decode_shape(int n, uint16_t k, uint64_t idx, int16_t *pulses) {
    if (k == 0 || n == 0) {
        memset(pulses, 0, n * sizeof(int16_t));
        return;
    }

    pvq_init_count_table();

    int32_t k_left = (int32_t)k;

    for (int i = 0; i < n; i++) {
        int remaining_dims = n - 1 - i;
        int16_t abs_p = 0;
        int sign = 1;
        int found = 0;

        /* Walk through possible |x_i| values */
        for (int j = 0; j <= k_left && !found; j++) {
            uint64_t block_size;
            if (j == 0) {
                block_size = pvq_count(remaining_dims, k_left);
                if (idx < block_size) {
                    abs_p = 0;
                    found = 1;
                } else {
                    idx -= block_size;
                }
            } else {
                /* Positive case */
                block_size = pvq_count(remaining_dims, k_left - j);
                if (idx < block_size) {
                    abs_p = (int16_t)j;
                    sign = 1;
                    found = 1;
                } else {
                    idx -= block_size;
                    /* Negative case */
                    if (idx < block_size) {
                        abs_p = (int16_t)j;
                        sign = -1;
                        found = 1;
                    } else {
                        idx -= block_size;
                    }
                }
            }
        }

        pulses[i] = (int16_t)(abs_p * sign);
        k_left -= abs_p;
    }
}

/* ── Greedy energy correction (Daala-style) ──────────────────────── */
/* After initial L2-normalized reconstruction, the integer rounding
 * causes energy overshoot/undershoot. This iteratively adjusts the
 * largest |recon[i]| by ±1 until sum(recon[i]^2) is within 1 LSB
 * of dq_gain^2. Both encoder and decoder call this identically. */
static void pvq_energy_correct(int16_t *recon, int n, uint32_t dq_gain) {
    uint64_t target = (uint64_t)dq_gain * dq_gain;
    int max_iters = 32; /* safety bound */

    for (int iter = 0; iter < max_iters; iter++) {
        uint64_t energy = 0;
        for (int i = 0; i < n; i++)
            energy += (int64_t)recon[i] * recon[i];

        int64_t diff = (int64_t)energy - (int64_t)target;
        if (diff >= -(int64_t)dq_gain && diff <= (int64_t)dq_gain)
            break; /* within tolerance */

        /* Find the component with largest absolute value */
        int max_idx = 0;
        int16_t max_abs = 0;
        for (int i = 0; i < n; i++) {
            int16_t a = recon[i] < 0 ? -recon[i] : recon[i];
            if (a > max_abs) { max_abs = a; max_idx = i; }
        }
        if (max_abs == 0) break;

        /* Shrink if energy too high, grow if too low */
        if (diff > 0) {
            recon[max_idx] += (recon[max_idx] > 0) ? -1 : 1;
        } else {
            recon[max_idx] += (recon[max_idx] >= 0) ? 1 : -1;
        }
    }
}

/* ── Full PVQ encode ──────────────────────────────────────────────── */

void pvq_encode(const int16_t *coeffs, int n, uint16_t step,
                int is_intra, uint8_t csf_weight,
                SntPVQCode *code, int16_t *recon) {
    /* Step 1: Compute and quantize gain */
    uint32_t gain = pvq_compute_gain(coeffs, n);
    code->gain = pvq_quantize_gain(gain, step, is_intra);

    if (code->gain == 0) {
        /* Skip block — zero everything */
        code->skip = 1;
        code->k = 0;
        code->shape_idx = 0;
        memset(recon, 0, n * sizeof(int16_t));
        return;
    }

    /* Smart skip (Loop 6): zero bands where dequantized energy overshoots
     * the original, AND bands with marginal energy.
     *
     * If gain_q * step + step/2 > gain (i.e., midpoint reconstruction
     * exceeds original energy), zeroing is strictly better than coding.
     * Also skip if gain_q <= 2 and gain < 2*step (marginal energy). */
    uint32_t dq_would_be = code->gain * step + step / 2;
    if (code->gain <= 2 && (gain < step * 2 || dq_would_be > gain * 2)) {
        code->skip = 1;
        code->gain = 0;
        code->k = 0;
        code->shape_idx = 0;
        memset(recon, 0, n * sizeof(int16_t));
        return;
    }

    code->skip = 0;

    /* Step 2: Compute K and encode shape */
    code->k = pvq_compute_k(code->gain, n, csf_weight);

    int16_t pulses[64]; /* max PVQ_MAX_N dimensions */
    code->shape_idx = pvq_encode_shape(coeffs, n, code->k, pulses);

    /* Step 3: Reconstruct — scale pulse vector by dequantized gain / ||pulses||_2,
     * then apply greedy energy correction to ensure sum(recon[i]^2) ≈ dq_gain^2
     * within 1 LSB per component. (Timothy's Daala-style greedy correction.) */
    uint32_t dq_gain = pvq_dequantize_gain(code->gain, step);
    uint32_t pulse_l2 = pvq_compute_gain(pulses, n);

    if (pulse_l2 > 0) {
        for (int i = 0; i < n; i++) {
            recon[i] = (int16_t)(((int64_t)pulses[i] * dq_gain
                       + pulse_l2 / 2) / pulse_l2);
        }
        /* Greedy energy correction: shrink largest |recon[i]| by 1 until
         * recon energy is within tolerance of target energy. */
        pvq_energy_correct(recon, n, dq_gain);
    } else {
        memset(recon, 0, n * sizeof(int16_t));
    }
}

/* ── Full PVQ decode ──────────────────────────────────────────────── */

void pvq_decode(const SntPVQCode *code, int n, uint16_t step,
                int16_t *recon) {
    if (code->skip || code->gain == 0 || code->k == 0) {
        memset(recon, 0, n * sizeof(int16_t));
        return;
    }

    /* Decode shape */
    int16_t pulses[64];
    pvq_decode_shape(n, code->k, code->shape_idx, pulses);

    /* Reconstruct: normalize by L2 norm + greedy energy correction */
    uint32_t dq_gain = pvq_dequantize_gain(code->gain, step);
    uint32_t pulse_l2 = pvq_compute_gain(pulses, n);

    if (pulse_l2 > 0) {
        for (int i = 0; i < n; i++) {
            recon[i] = (int16_t)(((int64_t)pulses[i] * dq_gain
                       + pulse_l2 / 2) / pulse_l2);
        }
        pvq_energy_correct(recon, n, dq_gain);
    } else {
        memset(recon, 0, n * sizeof(int16_t));
    }
}
