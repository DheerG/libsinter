/*
 * Sinter Codec — Activity Masking / Per-CTU QP Offset
 * Implements perceptual optimization: spatial activity masking
 *
 * The HVS is less sensitive to distortion in textured/busy regions
 * and more sensitive in smooth/flat regions. We compute a per-CTU
 * QP offset based on local AC variance:
 *
 *   delta_QP = alpha * log2(activity / mean_activity), clamp [-6, +6]
 *
 * where alpha=0.5 (configurable via PerceptualConfig).
 * Activity = variance of AC coefficients (exclude DC).
 * Must be computed BEFORE the RDO loop, not during.
 */

#include <math.h>
#include "activity_mask.h"

float psy_compute_activity(const int16_t *coeffs, int n) {
    /* Activity = variance of AC coefficients (skip DC at index 0).
     * var = E[x^2] - E[x]^2, computed over AC coefficients. */
    if (n <= 1) return 0.0f;

    int64_t sum = 0, sum_sq = 0;
    for (int i = 1; i < n; i++) {
        sum += coeffs[i];
        sum_sq += (int64_t)coeffs[i] * coeffs[i];
    }
    int ac_count = n - 1;
    float mean = (float)sum / ac_count;
    float var = (float)sum_sq / ac_count - mean * mean;
    return var > 0.0f ? var : 0.0f;
}

int8_t psy_activity_qp_offset(float activity, float mean_activity, float alpha) {
    if (mean_activity < 1.0f) return 0;
    if (activity < 1.0f) activity = 1.0f;

    float delta = alpha * log2f(activity / mean_activity);

    /* Clamp to [-6, +6] */
    if (delta > 6.0f) delta = 6.0f;
    if (delta < -6.0f) delta = -6.0f;

    return (int8_t)(delta > 0 ? delta + 0.5f : delta - 0.5f);
}
