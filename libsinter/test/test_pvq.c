/*
 * Sinter Codec — PVQ Tests (Gate 1)
 *
 * Tests: gain computation, DC quantization round-trip, PVQ energy
 * conservation, shape encode/decode round-trip.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "core/quantize/pvq.h"

static int tests_run = 0;
static int tests_passed = 0;

#define RUN_TEST(fn) do { \
    tests_run++; \
    printf("Running %s...\n", #fn); \
    if (fn()) { \
        tests_passed++; \
        printf("  PASS\n"); \
    } \
} while (0)

/* ── Test: gain computation ───────────────────────────────────────── */
static int test_gain_computation(void) {
    /* [3, 4] should have gain = 5 */
    int16_t v1[] = {3, 4};
    uint32_t g1 = pvq_compute_gain(v1, 2);
    if (g1 != 5) {
        printf("  FAIL: gain([3,4]) = %u, expected 5\n", g1);
        return 0;
    }

    /* All zeros should have gain = 0 */
    int16_t v2[] = {0, 0, 0, 0};
    if (pvq_compute_gain(v2, 4) != 0) {
        printf("  FAIL: gain(zeros) != 0\n");
        return 0;
    }

    /* [100, 0, 0] should have gain = 100 */
    int16_t v3[] = {100, 0, 0};
    uint32_t g3 = pvq_compute_gain(v3, 3);
    if (g3 != 100) {
        printf("  FAIL: gain([100,0,0]) = %u, expected 100\n", g3);
        return 0;
    }

    return 1;
}

/* ── Test: DC quantization round-trip ─────────────────────────────── */
static int test_dc_roundtrip(void) {
    uint16_t step = 16;

    /* Test a range of DC values */
    int16_t test_vals[] = {0, 1, 15, 16, 17, 100, -1, -15, -16, -100, 255, -255};
    int nvals = sizeof(test_vals) / sizeof(test_vals[0]);

    for (int i = 0; i < nvals; i++) {
        int16_t dc = test_vals[i];
        int16_t level = pvq_quantize_dc(dc, step, 0);
        int16_t recon = pvq_dequantize_dc(level, step);

        /* Dead-zone quantizer: reconstruction = level * step.
         * The input range for level L is [L*step + offset, (L+1)*step + offset).
         * Max error = step + offset (at the top of each bin). */
        uint32_t offset = step / 3;
        int max_err = (int)(step + offset);
        int err = abs(recon - dc);
        if (err > max_err) {
            printf("  FAIL: DC=%d, level=%d, recon=%d, err=%d > max_err=%d\n",
                   dc, level, recon, err, max_err);
            return 0;
        }
    }

    /* Zero input should always quantize to zero */
    if (pvq_quantize_dc(0, step, 0) != 0) {
        printf("  FAIL: DC=0 should quantize to 0\n");
        return 0;
    }

    return 1;
}

/* ── Test: gain quantization dead-zone ────────────────────────────── */
static int test_gain_deadzone(void) {
    uint16_t step = 30;

    /* Values within dead-zone should quantize to 0 */
    /* Inter dead-zone: step/3 = 10 */
    if (pvq_quantize_gain(5, step, 0) != 0) {
        printf("  FAIL: gain=5 should be in inter dead-zone\n");
        return 0;
    }
    if (pvq_quantize_gain(10, step, 0) != 0) {
        printf("  FAIL: gain=10 should be in inter dead-zone\n");
        return 0;
    }

    /* Intra dead-zone: step/6 = 5 */
    if (pvq_quantize_gain(3, step, 1) != 0) {
        printf("  FAIL: gain=3 should be in intra dead-zone\n");
        return 0;
    }

    /* With narrowed dead-zone (offset=step/4=7 for inter) and round-to-nearest:
     * gain=30, offset=7: (30-7+15)/30 = 38/30 = 1. So gain=30 is NOT in dead-zone. */
    uint32_t q_inter = pvq_quantize_gain(30, step, 0);
    if (q_inter == 0) {
        printf("  FAIL: gain=30 should not be in dead-zone with narrow offset\n");
        return 0;
    }

    /* Small gain within dead-zone: gain=5, offset=7, 5<=7 → level=0 */
    uint32_t q_small = pvq_quantize_gain(5, step, 0);
    if (q_small != 0) {
        printf("  FAIL: gain=5 should be in dead-zone, got %u\n", q_small);
        return 0;
    }

    /* Dequantization reconstructs at the bin centre: level*step + step/8 */
    uint32_t dq = pvq_dequantize_gain(3, step);
    uint32_t expected = 3 * 30 + 3; /* 93 = centre of intra bin [78, 108) */
    if (dq != expected) {
        printf("  FAIL: dequant(3, step=30) = %u, expected %u\n", dq, expected);
        return 0;
    }

    return 1;
}

/* ── Test: K parameter computation ────────────────────────────────── */
static int test_k_computation(void) {
    /* K=0 when gain=0 */
    if (pvq_compute_k(0, 16, 128) != 0) {
        printf("  FAIL: K should be 0 when gain=0\n");
        return 0;
    }

    /* K >= 1 when gain > 0 */
    uint16_t k = pvq_compute_k(5, 16, 128);
    if (k < 1) {
        printf("  FAIL: K should be >= 1 for nonzero gain, got %u\n", k);
        return 0;
    }

    /* Higher gain -> higher K */
    uint16_t k_low = pvq_compute_k(2, 16, 128);
    uint16_t k_high = pvq_compute_k(20, 16, 128);
    if (k_high < k_low) {
        printf("  FAIL: higher gain should give higher K (%u < %u)\n",
               k_high, k_low);
        return 0;
    }

    return 1;
}

/* ── Test: shape encode/decode round-trip ─────────────────────────── */
static int test_shape_roundtrip(void) {
    /* Simple case: known pulse vector */
    int16_t coeffs[] = {10, -5, 3, -1};
    int n = 4;
    uint16_t k = 8;

    int16_t pulses_enc[4], pulses_dec[4];
    uint64_t idx = pvq_encode_shape(coeffs, n, k, pulses_enc);

    /* Verify L1 norm of pulses = K */
    int32_t l1 = 0;
    for (int i = 0; i < n; i++) l1 += abs(pulses_enc[i]);
    if (l1 != k) {
        printf("  FAIL: pulse L1 norm = %d, expected %d\n", l1, k);
        return 0;
    }

    /* Decode and verify match */
    pvq_decode_shape(n, k, idx, pulses_dec);

    for (int i = 0; i < n; i++) {
        if (pulses_enc[i] != pulses_dec[i]) {
            printf("  FAIL: pulse mismatch at %d: enc=%d, dec=%d (idx=%llu)\n",
                   i, pulses_enc[i], pulses_dec[i], (unsigned long long)idx);
            return 0;
        }
    }

    return 1;
}

/* ── Test: full PVQ encode/decode energy conservation ─────────────── */
static int test_pvq_energy_conservation(void) {
    int16_t coeffs[] = {50, -30, 20, -10, 5, 15, -8, 3};
    int n = 8;
    uint16_t step = 10;

    SntPVQCode code;
    int16_t recon_enc[8], recon_dec[8];

    pvq_encode(coeffs, n, step, 0, 128, &code, recon_enc);

    if (code.skip) {
        printf("  (skipped — gain in dead-zone, acceptable)\n");
        return 1;
    }

    pvq_decode(&code, n, step, recon_dec);

    /* Verify encoder and decoder produce same reconstruction */
    for (int i = 0; i < n; i++) {
        if (recon_enc[i] != recon_dec[i]) {
            printf("  FAIL: recon mismatch at %d: enc=%d, dec=%d\n",
                   i, recon_enc[i], recon_dec[i]);
            return 0;
        }
    }

    /* Verify energy conservation:
     * sum(recon[i]^2) should be close to dequant_gain^2 */
    uint32_t dq_gain = pvq_dequantize_gain(code.gain, step);
    uint64_t target_energy = (uint64_t)dq_gain * dq_gain;
    uint64_t actual_energy = 0;
    for (int i = 0; i < n; i++)
        actual_energy += (int64_t)recon_dec[i] * recon_dec[i];

    /* Allow tolerance — integer division in reconstruction loses some energy */
    double ratio = (target_energy > 0) ?
        (double)actual_energy / target_energy : 1.0;
    if (ratio < 0.5 || ratio > 2.0) {
        printf("  FAIL: energy conservation ratio = %.3f (target=%llu, actual=%llu)\n",
               ratio, (unsigned long long)target_energy,
               (unsigned long long)actual_energy);
        return 0;
    }
    printf("  (energy ratio: %.3f, gain_q=%u, K=%u)\n", ratio, code.gain, code.k);

    return 1;
}

/* ── Test: PVQ skip on zero input ─────────────────────────────────── */
static int test_pvq_skip_zero(void) {
    int16_t coeffs[] = {0, 0, 0, 0};
    SntPVQCode code;
    int16_t recon[4];

    pvq_encode(coeffs, 4, 10, 0, 128, &code, recon);

    if (!code.skip) {
        printf("  FAIL: zero input should produce skip\n");
        return 0;
    }

    for (int i = 0; i < 4; i++) {
        if (recon[i] != 0) {
            printf("  FAIL: skip block recon[%d] = %d, expected 0\n", i, recon[i]);
            return 0;
        }
    }

    return 1;
}

/* ── Test: PVQ with random data ───────────────────────────────────── */
static int test_pvq_random(void) {
    srand(77777);
    int mismatches = 0;

    for (int trial = 0; trial < 100; trial++) {
        int n = 4 + (rand() % 5); /* 4 to 8 dimensions (enumeration validated) */
        int16_t coeffs[16];
        for (int i = 0; i < n; i++)
            coeffs[i] = (int16_t)(rand() % 200 - 100);

        uint16_t step = 5 + (rand() % 20);
        SntPVQCode code;
        int16_t recon_enc[16], recon_dec[16];

        pvq_encode(coeffs, n, step, 0, 128, &code, recon_enc);

        if (code.skip) continue;

        pvq_decode(&code, n, step, recon_dec);

        for (int i = 0; i < n; i++) {
            if (recon_enc[i] != recon_dec[i]) {
                mismatches++;
                if (mismatches <= 3) {
                    printf("  WARN: trial %d, dim %d, pos %d: enc=%d dec=%d\n",
                           trial, n, i, recon_enc[i], recon_dec[i]);
                }
            }
        }
    }

    if (mismatches > 0) {
        printf("  FAIL: %d reconstruction mismatches in 100 trials\n", mismatches);
        return 0;
    }
    return 1;
}

int main(void) {
    printf("=== Sinter PVQ Tests ===\n\n");

    printf("--- Gain and Quantization ---\n");
    RUN_TEST(test_gain_computation);
    RUN_TEST(test_dc_roundtrip);
    RUN_TEST(test_gain_deadzone);
    RUN_TEST(test_k_computation);

    printf("\n--- Shape Encoding ---\n");
    RUN_TEST(test_shape_roundtrip);

    printf("\n--- Full PVQ Pipeline ---\n");
    RUN_TEST(test_pvq_energy_conservation);
    RUN_TEST(test_pvq_skip_zero);
    RUN_TEST(test_pvq_random);

    printf("\n=== Results: %d/%d tests passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
