/*
 * Sinter Codec — Transform Tests (Gate 0 / Gate 1)
 *
 * Gate 0: Lossless round-trip (post_filter(pre_filter(x)) == x)
 * Gate 1: DCT forward/inverse round-trip, 2D transform round-trip
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "core/transform/lapping.h"
#include "core/transform/dct.h"

static int tests_run = 0;
static int tests_passed = 0;

#define ASSERT_EQ(a, b, msg) do { \
    if ((a) != (b)) { \
        printf("  FAIL: %s (got %d, expected %d)\n", msg, (int)(a), (int)(b)); \
        return 0; \
    } \
} while (0)

/* ── Test: 4-point lapping pre/post filter inverse ────────────────── */
static int test_lap4_roundtrip(void) {
    /* Test with random-ish data in [−200, 200] */
    int16_t data[] = {42, -117, 88, -63};
    int16_t orig[4];
    memcpy(orig, data, sizeof(data));

    snt_lap_pre_filter_4_ref(data, 1);
    snt_lap_post_filter_4_ref(data, 1);

    for (int i = 0; i < 4; i++)
        ASSERT_EQ(data[i], orig[i], "lap4 roundtrip mismatch");
    return 1;
}

/* ── Test: 8-point lapping pre/post filter inverse ────────────────── */
static int test_lap8_roundtrip(void) {
    int16_t data[] = {42, -117, 88, -63, 127, -22, 56, -99};
    int16_t orig[8];
    memcpy(orig, data, sizeof(data));

    snt_lap_pre_filter_8_ref(data, 1);
    snt_lap_post_filter_8_ref(data, 1);

    for (int i = 0; i < 8; i++)
        ASSERT_EQ(data[i], orig[i], "lap8 roundtrip mismatch");
    return 1;
}

/* ── Test: lapping with stride (vertical application) ─────────────── */
static int test_lap4_strided(void) {
    int16_t buf[4 * 8]; /* 4 rows, stride 8 */
    memset(buf, 0, sizeof(buf));
    buf[0 * 8] = 42;
    buf[1 * 8] = -117;
    buf[2 * 8] = 88;
    buf[3 * 8] = -63;

    int16_t orig[4] = {42, -117, 88, -63};

    snt_lap_pre_filter_4_ref(&buf[0], 8);
    snt_lap_post_filter_4_ref(&buf[0], 8);

    for (int i = 0; i < 4; i++)
        ASSERT_EQ(buf[i * 8], orig[i], "lap4 strided roundtrip mismatch");
    return 1;
}

/* ── Test: 4-point lapping with random data (1000 trials) ─────────── */
static int test_lap4_random(void) {
    srand(12345);
    for (int trial = 0; trial < 1000; trial++) {
        int16_t data[4], orig[4];
        for (int i = 0; i < 4; i++) {
            data[i] = (int16_t)(rand() % 512 - 256);
            orig[i] = data[i];
        }
        snt_lap_pre_filter_4_ref(data, 1);
        snt_lap_post_filter_4_ref(data, 1);
        for (int i = 0; i < 4; i++) {
            if (data[i] != orig[i]) {
                printf("  FAIL: lap4 random trial %d, index %d: got %d expected %d\n",
                       trial, i, data[i], orig[i]);
                return 0;
            }
        }
    }
    return 1;
}

/* ── Test: 8-point lapping with random data (1000 trials) ─────────── */
static int test_lap8_random(void) {
    srand(54321);
    for (int trial = 0; trial < 1000; trial++) {
        int16_t data[8], orig[8];
        for (int i = 0; i < 8; i++) {
            data[i] = (int16_t)(rand() % 512 - 256);
            orig[i] = data[i];
        }
        snt_lap_pre_filter_8_ref(data, 1);
        snt_lap_post_filter_8_ref(data, 1);
        for (int i = 0; i < 8; i++) {
            if (data[i] != orig[i]) {
                printf("  FAIL: lap8 random trial %d, index %d: got %d expected %d\n",
                       trial, i, data[i], orig[i]);
                return 0;
            }
        }
    }
    return 1;
}

/* ── Test: 2D lapping round-trip ──────────────────────────────────── */
static int test_lap_2d_roundtrip(void) {
    /* 16x16 buffer with block size 8 — has one internal boundary each way */
    const int W = 16, H = 16, S = 16, BS = 8;
    int16_t buf[16 * 16], orig[16 * 16];

    srand(99999);
    for (int i = 0; i < W * H; i++) {
        buf[i] = (int16_t)(rand() % 512 - 256);
        orig[i] = buf[i];
    }

    snt_lap_pre_filter_2d_ref(buf, W, H, S, BS);
    snt_lap_post_filter_2d_ref(buf, W, H, S, BS);

    for (int i = 0; i < W * H; i++) {
        if (buf[i] != orig[i]) {
            int y = i / S, x = i % S;
            printf("  FAIL: 2D lap roundtrip at (%d,%d): got %d expected %d\n",
                   x, y, buf[i], orig[i]);
            return 0;
        }
    }
    return 1;
}

/* ── Test: 4-point DCT forward/inverse round-trip ─────────────────── */
static int test_dct4_roundtrip(void) {
    srand(11111);
    int max_err = 0;
    for (int trial = 0; trial < 1000; trial++) {
        int16_t in[4], coeffs[4], out[4];
        for (int i = 0; i < 4; i++)
            in[i] = (int16_t)(rand() % 512 - 256);

        snt_dct_fwd_4_ref(in, coeffs);
        snt_dct_inv_4_ref(coeffs, out);

        for (int i = 0; i < 4; i++) {
            int err = abs(out[i] - in[i]);
            if (err > max_err) max_err = err;
            /* Allow rounding tolerance of 1 LSB for integer DCT */
            if (err > 1) {
                printf("  FAIL: dct4 trial %d, index %d: got %d expected %d (err %d)\n",
                       trial, i, out[i], in[i], err);
                return 0;
            }
        }
    }
    printf("  (dct4 max roundtrip error: %d)\n", max_err);
    return 1;
}

/* ── Test: 8-point DCT forward/inverse round-trip ─────────────────── */
static int test_dct8_roundtrip(void) {
    srand(22222);
    int max_err = 0;
    for (int trial = 0; trial < 1000; trial++) {
        int16_t in[8], coeffs[8], out[8];
        for (int i = 0; i < 8; i++)
            in[i] = (int16_t)(rand() % 512 - 256);

        snt_dct_fwd_8_ref(in, coeffs);
        snt_dct_inv_8_ref(coeffs, out);

        for (int i = 0; i < 8; i++) {
            int err = abs(out[i] - in[i]);
            if (err > max_err) max_err = err;
            if (err > 2) {
                printf("  FAIL: dct8 trial %d, index %d: got %d expected %d (err %d)\n",
                       trial, i, out[i], in[i], err);
                return 0;
            }
        }
    }
    printf("  (dct8 max roundtrip error: %d)\n", max_err);
    return 1;
}

/* ── Test: 2D DCT round-trip (4x4) ───────────────────────────────── */
static int test_dct_2d_4x4_roundtrip(void) {
    srand(33333);
    int max_err = 0;
    for (int trial = 0; trial < 100; trial++) {
        int16_t block[4 * 4], coeffs[16], out[4 * 4];
        for (int i = 0; i < 16; i++)
            block[i] = (int16_t)(rand() % 512 - 256);

        snt_dct_fwd_2d_ref(block, 4, 4, coeffs);
        snt_dct_inv_2d_ref(coeffs, 4, out, 4);

        for (int i = 0; i < 16; i++) {
            int err = abs(out[i] - block[i]);
            if (err > max_err) max_err = err;
            if (err > 2) {
                printf("  FAIL: 2d dct4x4 trial %d, index %d: got %d expected %d\n",
                       trial, i, out[i], block[i]);
                return 0;
            }
        }
    }
    printf("  (2d dct4x4 max roundtrip error: %d)\n", max_err);
    return 1;
}

/* ── Test: 2D DCT round-trip (8x8) ───────────────────────────────── */
static int test_dct_2d_8x8_roundtrip(void) {
    srand(44444);
    int max_err = 0;
    for (int trial = 0; trial < 100; trial++) {
        int16_t block[8 * 8], coeffs[64], out[8 * 8];
        for (int i = 0; i < 64; i++)
            block[i] = (int16_t)(rand() % 512 - 256);

        snt_dct_fwd_2d_ref(block, 8, 8, coeffs);
        snt_dct_inv_2d_ref(coeffs, 8, out, 8);

        for (int i = 0; i < 64; i++) {
            int err = abs(out[i] - block[i]);
            if (err > max_err) max_err = err;
            if (err > 3) {
                printf("  FAIL: 2d dct8x8 trial %d, index %d: got %d expected %d\n",
                       trial, i, out[i], block[i]);
                return 0;
            }
        }
    }
    printf("  (2d dct8x8 max roundtrip error: %d)\n", max_err);
    return 1;
}

/* ── Test: flat gray (128) through DCT produces only DC ──────────── */
static int test_dct_flat_gray(void) {
    int16_t block[4 * 4], coeffs[16];
    for (int i = 0; i < 16; i++) block[i] = 128;

    snt_dct_fwd_2d_ref(block, 4, 4, coeffs);

    /* DC should be nonzero, all AC should be zero */
    if (coeffs[0] == 0) {
        printf("  FAIL: flat gray DC is zero\n");
        return 0;
    }
    for (int i = 1; i < 16; i++) {
        if (coeffs[i] != 0) {
            printf("  FAIL: flat gray AC[%d] = %d (expected 0)\n", i, coeffs[i]);
            return 0;
        }
    }
    return 1;
}

/* ── Runner ───────────────────────────────────────────────────────── */
#define RUN_TEST(fn) do { \
    tests_run++; \
    printf("Running %s...\n", #fn); \
    if (fn()) { \
        tests_passed++; \
        printf("  PASS\n"); \
    } \
} while (0)

int main(void) {
    printf("=== Sinter Transform Tests ===\n\n");

    /* Gate 0: Lapped transform perfect reconstruction */
    printf("--- Gate 0: Lapped Transform Round-Trip ---\n");
    RUN_TEST(test_lap4_roundtrip);
    RUN_TEST(test_lap8_roundtrip);
    RUN_TEST(test_lap4_strided);
    RUN_TEST(test_lap4_random);
    RUN_TEST(test_lap8_random);
    RUN_TEST(test_lap_2d_roundtrip);

    /* Gate 1: DCT round-trip */
    printf("\n--- Gate 1: DCT Round-Trip ---\n");
    RUN_TEST(test_dct4_roundtrip);
    RUN_TEST(test_dct8_roundtrip);
    RUN_TEST(test_dct_2d_4x4_roundtrip);
    RUN_TEST(test_dct_2d_8x8_roundtrip);

    /* Gate 3: Reference signals */
    printf("\n--- Gate 3: Reference Signal Tests ---\n");
    RUN_TEST(test_dct_flat_gray);

    printf("\n=== Results: %d/%d tests passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
