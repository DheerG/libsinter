/*
 * Sinter Codec — rANS Entropy Coder
 * Implements Section 1.3: Entropy Coding
 *
 * Range Asymmetric Numeral Systems (rANS) with:
 * - 32-bit state, 16-bit renormalization
 * - 15-bit CDF probability resolution (totals sum to 32768)
 * - Normative dual interleaved rANS (Giesen's approach)
 * - PARA: rate=5 for coefficients, rate=4 for modes
 *
 * CRITICAL: ANS is LIFO (stack-based).
 *   Encoder processes symbols in REVERSE order.
 *   Decoder processes symbols in FORWARD order.
 * Getting this wrong produces plausible-looking but incorrect output.
 *
 * Reference: Duda, "Asymmetric numeral systems" (2009). Public domain.
 * Reference: Giesen, "Interleaved entropy coders" (2014).
 */

#include <stdint.h>
#include <string.h>
#include "rans.h"

/* ── CDF operations ───────────────────────────────────────────────── */

void rans_cdf_init_uniform(RansCDF *cdf, uint16_t num_symbols, uint8_t adapt_rate) {
    cdf->num_symbols = num_symbols;
    cdf->adapt_rate = adapt_rate;

    /* Uniform: each symbol gets RANS_PROB_TOTAL / num_symbols probability.
     * Ensure no symbol has zero probability (minimum 1). */
    uint16_t base = RANS_PROB_TOTAL / num_symbols;
    uint16_t remainder = RANS_PROB_TOTAL - base * num_symbols;

    cdf->cdf[0] = 0;
    for (uint16_t i = 0; i < num_symbols; i++) {
        uint16_t freq = base + (i < remainder ? 1 : 0);
        cdf->cdf[i + 1] = cdf->cdf[i] + freq;
    }
}

void rans_cdf_adapt(RansCDF *cdf, uint16_t symbol) {
    /* PARA exponential moving average update.
     * Move the CDF toward a distribution with all probability on 'symbol'.
     * Rate controls adaptation speed (lower = faster).
     *
     * For each CDF entry: cdf[i] += ((target[i] - cdf[i]) >> rate)
     * where target[i] = i <= symbol ? 0 : RANS_PROB_TOTAL
     * (i.e., target is a step function at 'symbol'). */
    uint8_t rate = cdf->adapt_rate;
    uint16_t nsym = cdf->num_symbols;

    /* cdf[0] is always 0, cdf[nsym] is always RANS_PROB_TOTAL. */
    for (uint16_t i = 1; i < nsym; i++) {
        uint16_t target = (i <= symbol) ? 0 : RANS_PROB_TOTAL;
        /* Rounding offset for faster CDF convergence. The spec is v0.1 —
         * this is the normative adaptation formula. Both encoder and decoder
         * use identical code, so interop is guaranteed. (Timothy: faster
         * convergence measurably improves coding efficiency.) */
        int32_t delta = ((int32_t)target - cdf->cdf[i] + (1 << (rate - 1))) >> rate;
        cdf->cdf[i] = (uint16_t)((int32_t)cdf->cdf[i] + delta);
    }

    /* Ensure monotonicity and bounds.
     * Forward pass: enforce cdf[i] >= cdf[i-1] + 1 (min probability of 1).
     * Backward pass: enforce cdf[i] <= cdf[i+1] - 1 (all entries < PROB_TOTAL).
     * Both passes are needed because the forward pass can push entries above
     * PROB_TOTAL when many symbols cluster at the same value. */
    cdf->cdf[nsym] = RANS_PROB_TOTAL;
    for (uint16_t i = 1; i < nsym; i++) {
        if (cdf->cdf[i] <= cdf->cdf[i - 1])
            cdf->cdf[i] = cdf->cdf[i - 1] + 1;
    }
    for (int i = (int)nsym - 1; i >= 1; i--) {
        if (cdf->cdf[i] >= cdf->cdf[i + 1])
            cdf->cdf[i] = cdf->cdf[i + 1] - 1;
    }
}

uint16_t rans_cdf_lookup(const RansCDF *cdf, uint16_t cum_val) {
    /* Linear search (reference impl; binary search for production) */
    for (uint16_t i = 0; i < cdf->num_symbols; i++) {
        if (cum_val < cdf->cdf[i + 1])
            return i;
    }
    return cdf->num_symbols - 1;
}

/* ── Encoder ──────────────────────────────────────────────────────── */

void rans_enc_init(RansEncoder *enc, uint16_t *rev_buf, uint32_t rev_cap, int dual) {
    enc->state0 = RANS_RENORM_LO; /* initial state in [RENORM_LO, 2*RENORM_LO) */
    enc->state1 = RANS_RENORM_LO;
    enc->dual = dual;
    enc->rev_buf = rev_buf;
    enc->rev_cap = rev_cap;
    enc->rev_pos = rev_cap; /* write from end */
}

/* Renormalize: output 16-bit words until state is in valid range.
 * The valid range for encoding symbol with freq f is:
 *   state in [RENORM_LO * f, RENORM_LO * f * (1 << 16))
 * We output 16-bit words to bring state into this range. */
static void rans_enc_renorm(RansEncoder *enc, RansState *state, uint16_t freq) {
    /* Max state before encoding: freq * (RENORM_LO << 16) / PROB_TOTAL
     * If state >= this, output a 16-bit word. */
    uint32_t max_state = ((uint32_t)freq << RANS_RENORM_BITS) *
                         (RANS_RENORM_LO >> RANS_PROB_BITS);
    /* Simplified: max_state = freq << (RENORM_BITS + RENORM_BITS - PROB_BITS)
     * = freq << (16 + 16 - 15) = freq << 17 */
    /* Use 64-bit to prevent overflow when freq is large */
    uint64_t max_state64 = (uint64_t)freq << (RANS_STATE_BITS - RANS_PROB_BITS);
    if (max_state64 > UINT32_MAX) max_state64 = UINT32_MAX; /* clamp: no renorm needed */
    max_state = (uint32_t)max_state64;

    while (*state >= max_state) {
        if (enc->rev_pos > 0) {
            enc->rev_pos--;
            enc->rev_buf[enc->rev_pos] = (uint16_t)(*state & 0xFFFF);
        }
        *state >>= 16;
    }
}

void rans_enc_put(RansEncoder *enc, uint16_t cdf_lo, uint16_t freq, uint32_t fwd_idx) {
    if (freq == 0) return; /* safety: shouldn't happen with proper CDFs */

    /* In dual mode: even forward indices use state0, odd use state1.
     * This matches the decoder which uses sym_count % 2. */
    RansState *state;
    if (enc->dual) {
        state = (fwd_idx % 2 == 0) ? &enc->state0 : &enc->state1;
    } else {
        state = &enc->state0;
    }

    /* Renormalize before encoding */
    rans_enc_renorm(enc, state, freq);

    /* Encode: state' = (state / freq) * PROB_TOTAL + (state % freq) + cdf_lo */
    *state = ((*state / freq) << RANS_PROB_BITS) + (*state % freq) + cdf_lo;
}

void rans_enc_flush(RansEncoder *enc) {
    /* Write final state(s) to buffer as 32-bit values (two 16-bit words each) */
    if (enc->dual) {
        /* State1 first (will be read last by decoder) */
        if (enc->rev_pos >= 2) {
            enc->rev_pos -= 2;
            enc->rev_buf[enc->rev_pos]     = (uint16_t)(enc->state1 >> 16);
            enc->rev_buf[enc->rev_pos + 1] = (uint16_t)(enc->state1 & 0xFFFF);
        }
    }
    /* State0 */
    if (enc->rev_pos >= 2) {
        enc->rev_pos -= 2;
        enc->rev_buf[enc->rev_pos]     = (uint16_t)(enc->state0 >> 16);
        enc->rev_buf[enc->rev_pos + 1] = (uint16_t)(enc->state0 & 0xFFFF);
    }
}

const uint16_t *rans_enc_data(const RansEncoder *enc, uint32_t *len_words) {
    *len_words = enc->rev_cap - enc->rev_pos;
    return &enc->rev_buf[enc->rev_pos];
}

/* ── Decoder ──────────────────────────────────────────────────────── */

static uint16_t read16(RansDecoder *dec) {
    if (dec->pos + 1 < dec->buf_len) {
        uint16_t val = (uint16_t)dec->buf[dec->pos]
                     | ((uint16_t)dec->buf[dec->pos + 1] << 8);
        dec->pos += 2;
        return val;
    }
    return 0;
}

void rans_dec_init(RansDecoder *dec, const uint8_t *buf, uint32_t buf_len, int dual) {
    dec->buf = buf;
    dec->buf_len = buf_len;
    dec->pos = 0;
    dec->dual = dual;
    dec->sym_count = 0;

    /* Read initial state0 (32-bit = two 16-bit words, high then low) */
    uint16_t hi = read16(dec);
    uint16_t lo = read16(dec);
    dec->state0 = ((uint32_t)hi << 16) | lo;

    if (dual) {
        hi = read16(dec);
        lo = read16(dec);
        dec->state1 = ((uint32_t)hi << 16) | lo;
    } else {
        dec->state1 = 0;
    }
}

uint16_t rans_dec_get(const RansDecoder *dec) {
    /* Extract cumulative frequency from current state.
     * In dual mode: even sym_count uses state0, odd uses state1. */
    const RansState *state;
    if (dec->dual) {
        state = (dec->sym_count % 2 == 0) ? &dec->state0 : &dec->state1;
    } else {
        state = &dec->state0;
    }
    return (uint16_t)(*state & (RANS_PROB_TOTAL - 1));
}

void rans_dec_advance(RansDecoder *dec, uint16_t cdf_lo, uint16_t freq) {
    RansState *state;
    if (dec->dual) {
        state = (dec->sym_count % 2 == 0) ? &dec->state0 : &dec->state1;
        dec->sym_count++;
    } else {
        state = &dec->state0;
    }

    /* Decode: state' = freq * (state >> PROB_BITS) + (state & PROB_MASK) - cdf_lo */
    *state = (uint32_t)freq * (*state >> RANS_PROB_BITS)
           + (*state & (RANS_PROB_TOTAL - 1)) - cdf_lo;

    /* Renormalize: read 16-bit words while state is too small */
    while (*state < RANS_RENORM_LO) {
        *state = (*state << 16) | read16(dec);
    }
}

uint16_t rans_dec_symbol(RansDecoder *dec, RansCDF *cdf) {
    uint16_t cum = rans_dec_get(dec);
    uint16_t sym = rans_cdf_lookup(cdf, cum);
    uint16_t cdf_lo, freq;
    rans_cdf_symbol_range(cdf, sym, &cdf_lo, &freq);
    rans_dec_advance(dec, cdf_lo, freq);
    return sym;
}
