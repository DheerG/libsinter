/*
 * Sinter Codec — rANS Entropy Coder Tests (Gate 1)
 *
 * Tests: CDF operations, encode/decode round-trip, dual interleaved mode,
 * adaptation, adversarial distributions.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "core/entropy/rans.h"

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

/* ── Test: CDF initialization ─────────────────────────────────────── */
static int test_cdf_init(void) {
    RansCDF cdf;
    rans_cdf_init_uniform(&cdf, 4, RANS_ADAPT_RATE_COEFF);

    /* Should be uniform: each symbol gets 32768/4 = 8192 */
    if (cdf.cdf[0] != 0) { printf("  FAIL: cdf[0]=%d\n", cdf.cdf[0]); return 0; }
    if (cdf.cdf[4] != RANS_PROB_TOTAL) {
        printf("  FAIL: cdf[4]=%d\n", cdf.cdf[4]); return 0;
    }

    /* All frequencies should be equal */
    for (int i = 0; i < 4; i++) {
        uint16_t freq = cdf.cdf[i + 1] - cdf.cdf[i];
        if (freq != 8192) {
            printf("  FAIL: freq[%d]=%d, expected 8192\n", i, freq);
            return 0;
        }
    }
    return 1;
}

/* ── Test: CDF lookup ─────────────────────────────────────────────── */
static int test_cdf_lookup(void) {
    RansCDF cdf;
    rans_cdf_init_uniform(&cdf, 4, RANS_ADAPT_RATE_COEFF);

    /* Cumulative 0 -> symbol 0, 8192 -> symbol 1, etc. */
    if (rans_cdf_lookup(&cdf, 0) != 0) { printf("  FAIL: lookup(0)\n"); return 0; }
    if (rans_cdf_lookup(&cdf, 8191) != 0) { printf("  FAIL: lookup(8191)\n"); return 0; }
    if (rans_cdf_lookup(&cdf, 8192) != 1) { printf("  FAIL: lookup(8192)\n"); return 0; }
    if (rans_cdf_lookup(&cdf, 32767) != 3) { printf("  FAIL: lookup(32767)\n"); return 0; }
    return 1;
}

/* ── Test: single-state rANS round-trip ───────────────────────────── */
static int test_rans_single_roundtrip(void) {
    RansCDF cdf;
    rans_cdf_init_uniform(&cdf, 4, RANS_ADAPT_RATE_COEFF);

    /* Encode 8 symbols */
    uint16_t symbols[] = {0, 1, 2, 3, 3, 2, 1, 0};
    int nsym = 8;

    uint16_t rev_buf[256];
    RansEncoder enc;
    rans_enc_init(&enc, rev_buf, 256, 0); /* single state */

    /* CRITICAL: encode in REVERSE order (LIFO) */
    for (int i = nsym - 1; i >= 0; i--) {
        uint16_t cdf_lo, freq;
        rans_cdf_symbol_range(&cdf, symbols[i], &cdf_lo, &freq);
        rans_enc_put(&enc, cdf_lo, freq, (uint32_t)i);
    }
    rans_enc_flush(&enc);

    /* Get encoded data */
    uint32_t enc_len;
    const uint16_t *enc_data = rans_enc_data(&enc, &enc_len);

    /* Convert to byte stream for decoder */
    uint8_t byte_buf[512];
    for (uint32_t i = 0; i < enc_len; i++) {
        byte_buf[2 * i]     = (uint8_t)(enc_data[i] & 0xFF);
        byte_buf[2 * i + 1] = (uint8_t)(enc_data[i] >> 8);
    }

    /* Decode in FORWARD order */
    RansDecoder dec;
    rans_dec_init(&dec, byte_buf, enc_len * 2, 0);

    /* Reset CDF for decode (must match encoder's initial state) */
    RansCDF cdf_dec;
    rans_cdf_init_uniform(&cdf_dec, 4, RANS_ADAPT_RATE_COEFF);

    for (int i = 0; i < nsym; i++) {
        uint16_t sym = rans_dec_symbol(&dec, &cdf_dec);
        if (sym != symbols[i]) {
            printf("  FAIL: symbol %d: decoded %d, expected %d\n", i, sym, symbols[i]);
            return 0;
        }
    }
    return 1;
}

/* ── Test: dual interleaved rANS round-trip ───────────────────────── */
static int test_rans_dual_roundtrip(void) {
    RansCDF cdf;
    rans_cdf_init_uniform(&cdf, 8, RANS_ADAPT_RATE_COEFF);

    uint16_t symbols[] = {0, 7, 3, 5, 2, 6, 1, 4, 0, 7, 3, 5};
    int nsym = 12;

    uint16_t rev_buf[256];
    RansEncoder enc;
    rans_enc_init(&enc, rev_buf, 256, 1); /* dual interleaved */

    for (int i = nsym - 1; i >= 0; i--) {
        uint16_t cdf_lo, freq;
        rans_cdf_symbol_range(&cdf, symbols[i], &cdf_lo, &freq);
        rans_enc_put(&enc, cdf_lo, freq, (uint32_t)i);
    }
    rans_enc_flush(&enc);

    uint32_t enc_len;
    const uint16_t *enc_data = rans_enc_data(&enc, &enc_len);

    uint8_t byte_buf[512];
    for (uint32_t i = 0; i < enc_len; i++) {
        byte_buf[2 * i]     = (uint8_t)(enc_data[i] & 0xFF);
        byte_buf[2 * i + 1] = (uint8_t)(enc_data[i] >> 8);
    }

    RansDecoder dec;
    rans_dec_init(&dec, byte_buf, enc_len * 2, 1);

    RansCDF cdf_dec;
    rans_cdf_init_uniform(&cdf_dec, 8, RANS_ADAPT_RATE_COEFF);

    for (int i = 0; i < nsym; i++) {
        uint16_t sym = rans_dec_symbol(&dec, &cdf_dec);
        if (sym != symbols[i]) {
            printf("  FAIL: symbol %d: decoded %d, expected %d\n", i, sym, symbols[i]);
            return 0;
        }
    }
    return 1;
}

/* ── Test: rANS with adaptive CDF ─────────────────────────────────── */
static int test_rans_adaptive(void) {
    uint16_t symbols[] = {0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 3, 0, 0, 0};
    int nsym = 16;

    /* Encode with adaptation */
    RansCDF cdf_enc;
    rans_cdf_init_uniform(&cdf_enc, 4, RANS_ADAPT_RATE_COEFF);

    uint16_t rev_buf[512];
    RansEncoder enc;
    rans_enc_init(&enc, rev_buf, 512, 0);

    /* Record CDF states and encode in reverse */
    uint16_t cdf_lo_arr[16], freq_arr[16];
    RansCDF cdf_states[16]; /* save CDF before each symbol */

    /* First pass: compute CDF state for each symbol (forward order) */
    {
        RansCDF cdf_tmp;
        rans_cdf_init_uniform(&cdf_tmp, 4, RANS_ADAPT_RATE_COEFF);
        for (int i = 0; i < nsym; i++) {
            memcpy(&cdf_states[i], &cdf_tmp, sizeof(RansCDF));
            rans_cdf_symbol_range(&cdf_tmp, symbols[i], &cdf_lo_arr[i], &freq_arr[i]);
            rans_cdf_adapt(&cdf_tmp, symbols[i]);
        }
    }

    /* Encode in reverse using saved CDF states */
    for (int i = nsym - 1; i >= 0; i--) {
        rans_enc_put(&enc, cdf_lo_arr[i], freq_arr[i], (uint32_t)i);
    }
    rans_enc_flush(&enc);

    uint32_t enc_len;
    const uint16_t *enc_data = rans_enc_data(&enc, &enc_len);

    uint8_t byte_buf[1024];
    for (uint32_t i = 0; i < enc_len; i++) {
        byte_buf[2 * i]     = (uint8_t)(enc_data[i] & 0xFF);
        byte_buf[2 * i + 1] = (uint8_t)(enc_data[i] >> 8);
    }

    /* Decode with adaptation (forward order) */
    RansDecoder dec;
    rans_dec_init(&dec, byte_buf, enc_len * 2, 0);

    RansCDF cdf_dec;
    rans_cdf_init_uniform(&cdf_dec, 4, RANS_ADAPT_RATE_COEFF);

    for (int i = 0; i < nsym; i++) {
        uint16_t sym = rans_dec_symbol(&dec, &cdf_dec);
        rans_cdf_adapt(&cdf_dec, sym);
        if (sym != symbols[i]) {
            printf("  FAIL: symbol %d: decoded %d, expected %d\n", i, sym, symbols[i]);
            return 0;
        }
    }
    return 1;
}

/* ── Test: adversarial distribution (one symbol nearly certain) ──── */
static int test_rans_skewed(void) {
    RansCDF cdf;
    rans_cdf_init_uniform(&cdf, 4, RANS_ADAPT_RATE_COEFF);

    /* Heavily adapt toward symbol 0 */
    for (int i = 0; i < 200; i++)
        rans_cdf_adapt(&cdf, 0);

    /* Verify symbol 0 has high probability but all symbols still encodable */
    uint16_t cdf_lo, freq;
    rans_cdf_symbol_range(&cdf, 0, &cdf_lo, &freq);
    if (freq < RANS_PROB_TOTAL / 2) {
        printf("  FAIL: symbol 0 should dominate after adaptation, freq=%d\n", freq);
        return 0;
    }

    /* All symbols should still have freq >= 1 */
    for (int i = 0; i < 4; i++) {
        rans_cdf_symbol_range(&cdf, (uint16_t)i, &cdf_lo, &freq);
        if (freq == 0) {
            printf("  FAIL: symbol %d has zero frequency\n", i);
            return 0;
        }
    }

    /* Encode and decode a mixed sequence with skewed CDF */
    uint16_t symbols[] = {0, 0, 0, 0, 1, 0, 0, 2, 0, 3};
    int nsym = 10;

    uint16_t rev_buf[256];
    RansEncoder enc;
    rans_enc_init(&enc, rev_buf, 256, 0);

    for (int i = nsym - 1; i >= 0; i--) {
        rans_cdf_symbol_range(&cdf, symbols[i], &cdf_lo, &freq);
        rans_enc_put(&enc, cdf_lo, freq, (uint32_t)i);
    }
    rans_enc_flush(&enc);

    uint32_t enc_len;
    const uint16_t *enc_data = rans_enc_data(&enc, &enc_len);

    uint8_t byte_buf[512];
    for (uint32_t i = 0; i < enc_len; i++) {
        byte_buf[2 * i]     = (uint8_t)(enc_data[i] & 0xFF);
        byte_buf[2 * i + 1] = (uint8_t)(enc_data[i] >> 8);
    }

    RansDecoder dec;
    rans_dec_init(&dec, byte_buf, enc_len * 2, 0);

    for (int i = 0; i < nsym; i++) {
        uint16_t cum = rans_dec_get(&dec);
        uint16_t sym = rans_cdf_lookup(&cdf, cum);
        rans_cdf_symbol_range(&cdf, sym, &cdf_lo, &freq);
        rans_dec_advance(&dec, cdf_lo, freq);
        if (sym != symbols[i]) {
            printf("  FAIL: symbol %d: decoded %d, expected %d\n", i, sym, symbols[i]);
            return 0;
        }
    }
    return 1;
}

/* ── Test: random symbols round-trip (stress test) ────────────────── */
static int test_rans_random(void) {
    srand(42424);
    int num_symbols = 16;
    int nsym = 200;

    RansCDF cdf;
    rans_cdf_init_uniform(&cdf, (uint16_t)num_symbols, RANS_ADAPT_RATE_COEFF);

    uint16_t symbols[200];
    for (int i = 0; i < nsym; i++)
        symbols[i] = (uint16_t)(rand() % num_symbols);

    uint16_t rev_buf[2048];
    RansEncoder enc;
    rans_enc_init(&enc, rev_buf, 2048, 1); /* dual */

    for (int i = nsym - 1; i >= 0; i--) {
        uint16_t cdf_lo, freq;
        rans_cdf_symbol_range(&cdf, symbols[i], &cdf_lo, &freq);
        rans_enc_put(&enc, cdf_lo, freq, (uint32_t)i);
    }
    rans_enc_flush(&enc);

    uint32_t enc_len;
    const uint16_t *enc_data = rans_enc_data(&enc, &enc_len);

    uint8_t byte_buf[4096];
    for (uint32_t i = 0; i < enc_len; i++) {
        byte_buf[2 * i]     = (uint8_t)(enc_data[i] & 0xFF);
        byte_buf[2 * i + 1] = (uint8_t)(enc_data[i] >> 8);
    }

    RansDecoder dec;
    rans_dec_init(&dec, byte_buf, enc_len * 2, 1);

    /* Must use same CDF as encoder */
    RansCDF cdf_dec;
    rans_cdf_init_uniform(&cdf_dec, (uint16_t)num_symbols, RANS_ADAPT_RATE_COEFF);

    for (int i = 0; i < nsym; i++) {
        uint16_t sym = rans_dec_symbol(&dec, &cdf_dec);
        if (sym != symbols[i]) {
            printf("  FAIL: random symbol %d: decoded %d, expected %d\n",
                   i, sym, symbols[i]);
            return 0;
        }
    }

    printf("  (encoded %d symbols in %u 16-bit words = %.2f bits/symbol)\n",
           nsym, enc_len, (float)(enc_len * 16) / nsym);
    return 1;
}

/* ── Test: bitstream determinism ──────────────────────────────────── */
static int test_rans_determinism(void) {
    srand(99999);
    uint16_t symbols[50];
    for (int i = 0; i < 50; i++)
        symbols[i] = (uint16_t)(rand() % 8);

    RansCDF cdf;
    rans_cdf_init_uniform(&cdf, 8, RANS_ADAPT_RATE_MODE);

    /* Encode twice */
    uint16_t rev1[256], rev2[256];
    uint32_t len1, len2;

    for (int pass = 0; pass < 2; pass++) {
        uint16_t *rev = (pass == 0) ? rev1 : rev2;
        RansEncoder enc;
        rans_enc_init(&enc, rev, 256, 0);

        for (int i = 49; i >= 0; i--) {
            uint16_t cdf_lo, freq;
            rans_cdf_symbol_range(&cdf, symbols[i], &cdf_lo, &freq);
            rans_enc_put(&enc, cdf_lo, freq, (uint32_t)i);
        }
        rans_enc_flush(&enc);

        uint32_t len;
        rans_enc_data(&enc, &len);
        if (pass == 0) len1 = len; else len2 = len;
    }

    if (len1 != len2) {
        printf("  FAIL: different lengths: %u vs %u\n", len1, len2);
        return 0;
    }

    /* Compare the actual data regions */
    uint32_t off1 = 256 - len1;
    uint32_t off2 = 256 - len2;
    if (memcmp(&rev1[off1], &rev2[off2], len1 * sizeof(uint16_t)) != 0) {
        printf("  FAIL: bitstreams differ\n");
        return 0;
    }

    return 1;
}

int main(void) {
    printf("=== Sinter rANS Entropy Tests ===\n\n");

    printf("--- CDF Operations ---\n");
    RUN_TEST(test_cdf_init);
    RUN_TEST(test_cdf_lookup);

    printf("\n--- Single-State rANS ---\n");
    RUN_TEST(test_rans_single_roundtrip);

    printf("\n--- Dual Interleaved rANS ---\n");
    RUN_TEST(test_rans_dual_roundtrip);

    printf("\n--- Adaptive Coding ---\n");
    RUN_TEST(test_rans_adaptive);
    RUN_TEST(test_rans_skewed);

    printf("\n--- Stress Tests ---\n");
    RUN_TEST(test_rans_random);
    RUN_TEST(test_rans_determinism);

    printf("\n=== Results: %d/%d tests passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
