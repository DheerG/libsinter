/*
 * Sinter Codec — Frequency-Domain Intra Prediction
 * Implements Section 1.4: Intra Prediction
 *
 * Predicts DCT coefficients from neighboring blocks' DCT coefficients
 * directly in the frequency domain. This avoids the pixel-domain angular
 * prediction modes used in H.264/H.265 and their associated patent clusters.
 *
 * Two-pass scheme:
 *   - DC and low-AC bands: predicted from pre-lapped neighbor coefficients
 *   - High-AC bands: zero prediction (R-D gated — skip when predictor
 *     residual would exceed original)
 *
 * Read from reconstructed buffer, not original input (causal context).
 */

#include <string.h>
#include "intra.h"

/* Number of low-frequency coefficients to predict (diagonal scan).
 * For 4x4: predict DC + 3 lowest AC = 4 coefficients.
 * For 8x8: predict DC + 5 lowest AC = 6 coefficients. */
#define LOW_FREQ_COUNT_4 4
#define LOW_FREQ_COUNT_8 6

void snt_intra_predict_freq_ref(const int16_t *left_coeffs,
                                 const int16_t *above_coeffs,
                                 int size, int16_t *pred) {
    int n = size * size;
    memset(pred, 0, n * sizeof(int16_t));

    int low_count = (size == 4) ? LOW_FREQ_COUNT_4 : LOW_FREQ_COUNT_8;

    /* DC prediction: average of left and above DC (or whichever is available) */
    if (left_coeffs && above_coeffs) {
        pred[0] = (int16_t)(((int32_t)left_coeffs[0] + above_coeffs[0] + 1) >> 1);
    } else if (left_coeffs) {
        pred[0] = left_coeffs[0];
    } else if (above_coeffs) {
        pred[0] = above_coeffs[0];
    }
    /* else: pred[0] = 0 (no neighbors available) */

    /* Low-AC prediction: simple average of corresponding coefficients
     * from left and above neighbors. Only predict the lowest frequencies
     * (diagonal band indices < low_count). */
    for (int i = 1; i < low_count && i < n; i++) {
        int row = i / size;
        int col = i % size;
        /* Only predict if coefficient is in a low-frequency band
         * (row + col < threshold) */
        if (row + col < 3) {
            int32_t p = 0;
            int count = 0;
            if (left_coeffs) { p += left_coeffs[i]; count++; }
            if (above_coeffs) { p += above_coeffs[i]; count++; }
            if (count > 0)
                pred[i] = (int16_t)((p + count / 2) / count);
        }
    }
    /* High-AC: zero prediction (the R-D gating in the encoder will
     * decide whether to use this prediction or not). */
}

void snt_intra_reconstruct_freq_ref(const int16_t *residual,
                                     const int16_t *pred,
                                     int size, int16_t *coeffs) {
    int n = size * size;
    for (int i = 0; i < n; i++)
        coeffs[i] = residual[i] + pred[i];
}
