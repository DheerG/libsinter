/*
 * Sinter Codec — Daala-Style Per-Dimension PVQ Shape Coding via rANS
 *
 * The combinatorial enumeration in pvq.c computes pvq_count(n, k) at
 * each dimension during encode/decode. These counts ARE the exact
 * probability distribution for each dimension's pulse value. This
 * module feeds those counts directly into rANS as per-symbol
 * probabilities, achieving theoretically optimal shape coding.
 *
 * Estimated improvement: 10-20% bitrate reduction on AC data compared
 * to the flat combinatorial index approach.
 */

#include <string.h>
#include "pvq_shape_rans.h"
#include "pvq.h"

/* Get pvq_count through the public accessor */
extern uint64_t pvq_num_vectors(int n, uint16_t k);

/* Scale a combinatorial count to a 15-bit rANS probability.
 * total: total count for this dimension. sub: sub-count for this symbol.
 * Returns probability in [1, RANS_PROB_TOTAL-1], ensuring no zero-probability symbols. */
static uint16_t scale_prob(uint64_t sub, uint64_t total) {
    if (total == 0) return RANS_PROB_TOTAL / 2;
    uint64_t p = (sub * RANS_PROB_TOTAL) / total;
    if (p < 1) p = 1;
    if (p >= RANS_PROB_TOTAL) p = RANS_PROB_TOTAL - 1;
    return (uint16_t)p;
}

void pvq_shape_encode_rans(const int16_t *pulses, int n, uint16_t k,
                            RansEncoder *enc, uint32_t fwd_idx) {
    if (k == 0 || n == 0) return;

    int32_t k_left = (int32_t)k;

    /* Process dimensions in forward order for CDF synchronization.
     * rANS encoding is done in reverse by the caller;
     * we record ranges and let the caller handle ordering. */
    for (int i = 0; i < n; i++) {
        int16_t abs_p = pulses[i] < 0 ? -pulses[i] : pulses[i];
        int remaining = n - 1 - i;

        /* Compute total count for this dimension */
        uint64_t total = pvq_num_vectors(remaining + 1, (uint16_t)k_left);
        if (total == 0) total = 1;

        /* Compute CDF position for this pulse value */
        uint64_t cdf_lo_64 = 0;
        for (int j = 0; j < abs_p; j++) {
            if (j == 0) {
                cdf_lo_64 += pvq_num_vectors(remaining, (uint16_t)k_left);
            } else {
                cdf_lo_64 += 2 * pvq_num_vectors(remaining, (uint16_t)(k_left - j));
            }
        }
        /* Add sign offset for negative */
        if (abs_p > 0 && pulses[i] < 0) {
            cdf_lo_64 += pvq_num_vectors(remaining, (uint16_t)(k_left - abs_p));
        }

        /* Compute frequency for this symbol */
        uint64_t freq_64;
        if (abs_p == 0) {
            freq_64 = pvq_num_vectors(remaining, (uint16_t)k_left);
        } else {
            freq_64 = pvq_num_vectors(remaining, (uint16_t)(k_left - abs_p));
        }

        /* Scale to 15-bit rANS probabilities */
        uint16_t cdf_lo = scale_prob(cdf_lo_64, total);
        uint16_t freq = scale_prob(freq_64, total);
        if (freq < 1) freq = 1;
        /* Ensure cdf_lo + freq <= RANS_PROB_TOTAL */
        if (cdf_lo + freq > RANS_PROB_TOTAL)
            freq = RANS_PROB_TOTAL - cdf_lo;

        rans_enc_put(enc, cdf_lo, freq, fwd_idx + (uint32_t)i);

        k_left -= abs_p;
    }
}

void pvq_shape_decode_rans(int16_t *pulses, int n, uint16_t k,
                            RansDecoder *dec) {
    if (k == 0 || n == 0) {
        memset(pulses, 0, n * sizeof(int16_t));
        return;
    }

    int32_t k_left = (int32_t)k;

    for (int i = 0; i < n; i++) {
        int remaining = n - 1 - i;
        uint64_t total = pvq_num_vectors(remaining + 1, (uint16_t)k_left);
        if (total == 0) { pulses[i] = 0; continue; }

        /* Get cumulative value from rANS state */
        uint16_t cum = rans_dec_get(dec);

        /* Map cumulative back to pulse value by walking the CDF */
        uint64_t cum_scaled = ((uint64_t)cum * total) / RANS_PROB_TOTAL;
        int16_t abs_p = 0;
        int sign = 1;

        uint64_t running = 0;
        int found = 0;
        for (int j = 0; j <= k_left && !found; j++) {
            uint64_t count;
            if (j == 0) {
                count = pvq_num_vectors(remaining, (uint16_t)k_left);
                if (cum_scaled < running + count) {
                    abs_p = 0; found = 1;
                } else {
                    running += count;
                }
            } else {
                count = pvq_num_vectors(remaining, (uint16_t)(k_left - j));
                /* Positive */
                if (cum_scaled < running + count) {
                    abs_p = (int16_t)j; sign = 1; found = 1;
                } else {
                    running += count;
                    /* Negative */
                    if (cum_scaled < running + count) {
                        abs_p = (int16_t)j; sign = -1; found = 1;
                    } else {
                        running += count;
                    }
                }
            }
        }

        pulses[i] = (int16_t)(abs_p * sign);

        /* Advance rANS state */
        uint16_t cdf_lo = scale_prob(running, total);
        uint64_t freq_64 = (abs_p == 0) ?
            pvq_num_vectors(remaining, (uint16_t)k_left) :
            pvq_num_vectors(remaining, (uint16_t)(k_left - abs_p));
        uint16_t freq = scale_prob(freq_64, total);
        if (freq < 1) freq = 1;
        if (cdf_lo + freq > RANS_PROB_TOTAL)
            freq = RANS_PROB_TOTAL - cdf_lo;

        rans_dec_advance(dec, cdf_lo, freq);

        k_left -= abs_p;
    }
}
