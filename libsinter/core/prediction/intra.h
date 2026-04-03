/*
 * Sinter Codec — Frequency-Domain Intra Prediction
 * Implements Section 1.4: Intra Prediction
 */
#ifndef SNT_INTRA_H
#define SNT_INTRA_H

#include <stdint.h>

/* Predict DC and low-AC coefficients from neighboring blocks' coefficients.
 * coeffs: current block's transform coefficients (to be predicted).
 * left_coeffs: left neighbor's reconstructed coefficients (NULL if unavailable).
 * above_coeffs: above neighbor's reconstructed coefficients (NULL if unavailable).
 * size: block size (4 or 8).
 * pred: output prediction (size*size coefficients).
 * residual: output = coeffs - pred. */
void snt_intra_predict_freq_ref(const int16_t *left_coeffs,
                                 const int16_t *above_coeffs,
                                 int size, int16_t *pred);

/* Reconstruct: add prediction back to residual. */
void snt_intra_reconstruct_freq_ref(const int16_t *residual,
                                     const int16_t *pred,
                                     int size, int16_t *coeffs);

#endif /* SNT_INTRA_H */
