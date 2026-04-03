/*
 * Sinter Codec — Perceptual Transform Distortion (PTD)
 * Implements perceptual optimization: CSF-weighted transform-domain MSE
 *
 * PTD replaces MSE in the RDO inner loop with a perceptually-weighted
 * distortion metric. It exploits the HVS's reduced sensitivity to high
 * spatial frequencies via the Contrast Sensitivity Function (CSF).
 *
 * PTD = sum_i (csf_weight[i] * (orig[i] - recon[i])^2)
 *
 * Near-zero marginal cost over MSE. Structurally different from SSIM
 * (no patent risk — Wang et al. patents cover SSIM specifically).
 *
 * Invariant: for any input, CSF-weighted distortion <= unweighted MSE
 * (because all CSF weights are <= 1 after normalization).
 */

#include <stdint.h>
#include "ptd.h"
#include "csf_weights.h"

float psy_ptd_distortion(const int16_t *orig, const int16_t *recon,
                          int size) {
    const float *csf = psy_csf_table(size);
    if (!csf) return 0.0f;

    int n = size * size;
    float dist = 0.0f;

    for (int i = 0; i < n; i++) {
        float diff = (float)(orig[i] - recon[i]);
        dist += csf[i] * diff * diff;
    }
    return dist;
}

float psy_mse_distortion(const int16_t *orig, const int16_t *recon, int n) {
    float dist = 0.0f;
    for (int i = 0; i < n; i++) {
        float diff = (float)(orig[i] - recon[i]);
        dist += diff * diff;
    }
    return dist / n;
}
