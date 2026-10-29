/*
 * Sinter Codec — Range coder tests
 *
 * Tests: model initialization and halving, encode/decode round trips
 * (uniform, skewed, many interleaved models, long random runs that force
 * carries), trailing-zero trimming, determinism.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "core/entropy/rangecoder.h"

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

/* Deterministic pseudo-random numbers (LCG) */
static uint32_t rng_state = 1;
static uint32_t rng(void) {
    rng_state = rng_state * 1103515245u + 12345u;
    return rng_state >> 8;
}

/* Encode n symbols with one model per context ctx[i] (or model 0 when
 * ctx is NULL), decode them, and compare. Returns the byte count, or -1
 * on a mismatch. */
static int roundtrip(const int *syms, const int *ctx, int n, int nsym, int nctx,
                     int inc, uint32_t limit) {
    SntModel *em = malloc(sizeof(SntModel) * nctx), *dm = malloc(sizeof(SntModel) * nctx);
    for (int c = 0; c < nctx; c++) {
        snt_model_init(&em[c], nsym, 1, inc, limit);
        snt_model_init(&dm[c], nsym, 1, inc, limit);
    }
    uint32_t cap = (uint32_t)n * 2 + 64;
    uint8_t *buf = malloc(cap);
    SntRangeEnc e;
    snt_rc_enc_init(&e, buf, cap);
    for (int i = 0; i < n; i++)
        snt_rc_encode_sym(&e, &em[ctx ? ctx[i] : 0], syms[i]);
    uint32_t len = snt_rc_enc_finish(&e);

    SntRangeDec d;
    snt_rc_dec_init(&d, buf, len);
    int ok = 1;
    for (int i = 0; i < n && ok; i++) {
        int s = snt_rc_decode_sym(&d, &dm[ctx ? ctx[i] : 0]);
        if (s != syms[i]) {
            printf("  FAIL: symbol %d decoded as %d, expected %d\n", i, s, syms[i]);
            ok = 0;
        }
    }
    free(em); free(dm); free(buf);
    return ok ? (int)len : -1;
}

/* ── Test: model initialization ───────────────────────────────────── */
static int test_model_init(void) {
    SntModel m;
    snt_model_init(&m, 4, 1, SNT_RC_INC_COEFF, SNT_RC_LIMIT_COEFF);
    if (m.total != 4) { printf("  FAIL: total=%u\n", m.total); return 0; }
    for (int i = 0; i < 4; i++)
        if (m.freq[i] != 1) { printf("  FAIL: freq[%d]=%d\n", i, m.freq[i]); return 0; }
    if (snt_model_cum(&m, 3) != 3) { printf("  FAIL: cum(3)\n"); return 0; }
    return 1;
}

/* ── Test: halving keeps counts positive and the total within limit ── */
static int test_model_halving(void) {
    SntModel m;
    snt_model_init(&m, 256, 1, 24, 8192);
    for (int i = 0; i < 5000; i++) {
        snt_model_update(&m, i % 3 == 0 ? 7 : 0);
        if (m.total > 8192) { printf("  FAIL: total %u over limit\n", m.total); return 0; }
    }
    uint32_t sum = 0;
    for (int i = 0; i < 256; i++) {
        if (m.freq[i] == 0) { printf("  FAIL: freq[%d] = 0\n", i); return 0; }
        sum += m.freq[i];
    }
    if (sum != m.total) { printf("  FAIL: sum %u != total %u\n", sum, m.total); return 0; }
    if (m.freq[0] <= m.freq[7] || m.freq[7] <= m.freq[1]) {
        printf("  FAIL: counts do not follow the statistics\n"); return 0;
    }
    return 1;
}

/* ── Test: uniform round trip ─────────────────────────────────────── */
static int test_rc_uniform(void) {
    int syms[4000];
    for (int i = 0; i < 4000; i++) syms[i] = (int)(rng() % 16);
    int len = roundtrip(syms, NULL, 4000, 16, 1, 24, 8192);
    if (len < 0) return 0;
    /* 4 bits per symbol: about 2000 bytes */
    if (len < 1900 || len > 2100) { printf("  FAIL: %d bytes for 2000 expected\n", len); return 0; }
    return 1;
}

/* ── Test: skewed source compresses well ──────────────────────────── */
static int test_rc_skewed(void) {
    int syms[20000];
    for (int i = 0; i < 20000; i++) {
        uint32_t r = rng() % 1000;
        syms[i] = r < 900 ? 0 : (r < 980 ? 1 : (int)(2 + r % 254));
    }
    int len = roundtrip(syms, NULL, 20000, 256, 1, 24, 8192);
    if (len < 0) return 0;
    /* Entropy about 0.75 bits/symbol: well under 1 bit */
    if (len > 20000 / 8) { printf("  FAIL: %d bytes, not compressed\n", len); return 0; }
    printf("  (%d bytes for 20000 symbols)\n", len);
    return 1;
}

/* ── Test: many interleaved models (as the codec's contexts) ──────── */
static int test_rc_contexts(void) {
    enum { N = 30000 };
    static int syms[N], ctx[N];
    for (int i = 0; i < N; i++) {
        ctx[i] = (int)(rng() % 5);
        syms[i] = ctx[i] == 0 ? (int)(rng() % 2) : (int)((rng() % 3) * ctx[i]);
    }
    return roundtrip(syms, ctx, N, 256, 5, 24, 2048) >= 0;
}

/* ── Test: long runs of near-certain symbols force carries ────────── */
static int test_rc_carries(void) {
    enum { N = 200000 };
    static int syms[N];
    for (int t = 0; t < 8; t++) {
        for (int i = 0; i < N; i++)
            syms[i] = (rng() % 4096) == 0 ? (int)(1 + rng() % 255) : 255 * (t & 1);
        if (roundtrip(syms, NULL, N, 256, 1, 32, 65000) < 0) return 0;
    }
    return 1;
}

/* ── Test: a stream that ends in likely symbols drops zero bytes ──── */
static int test_rc_trailing(void) {
    int syms[3] = {0, 0, 0};
    int len = roundtrip(syms, NULL, 3, 2, 1, 24, 2048);
    if (len < 0) return 0;
    if (len > 1) { printf("  FAIL: %d bytes for 3 likely symbols\n", len); return 0; }
    int one[1] = {1};
    return roundtrip(one, NULL, 1, 256, 1, 24, 2048) >= 0 &&
           roundtrip(one, NULL, 0, 256, 1, 24, 2048) == 0;
}

/* ── Test: identical input gives identical bytes ──────────────────── */
static int test_rc_determinism(void) {
    uint8_t a[4096], b[4096];
    uint32_t la = 0, lb = 0;
    for (int pass = 0; pass < 2; pass++) {
        SntModel m;
        snt_model_init(&m, 256, 1, 24, 8192);
        SntRangeEnc e;
        snt_rc_enc_init(&e, pass ? b : a, 4096);
        rng_state = 99;
        for (int i = 0; i < 3000; i++)
            snt_rc_encode_sym(&e, &m, (int)(rng() % 7));
        *(pass ? &lb : &la) = snt_rc_enc_finish(&e);
    }
    if (la != lb || memcmp(a, b, la)) { printf("  FAIL: outputs differ\n"); return 0; }
    return 1;
}

int main(void) {
    printf("=== Sinter Range Coder Tests ===\n\n");

    printf("--- Models ---\n");
    RUN_TEST(test_model_init);
    RUN_TEST(test_model_halving);

    printf("\n--- Round Trips ---\n");
    RUN_TEST(test_rc_uniform);
    RUN_TEST(test_rc_skewed);
    RUN_TEST(test_rc_contexts);

    printf("\n--- Stress Tests ---\n");
    RUN_TEST(test_rc_carries);
    RUN_TEST(test_rc_trailing);
    RUN_TEST(test_rc_determinism);

    printf("\n=== Results: %d/%d tests passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
