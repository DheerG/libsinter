/*
 * Sinter Codec — Contrast Sensitivity Function Weight Table
 * Implements perceptual optimization: CSF-weighted distortion
 *
 * The CSF describes the human visual system's sensitivity to different
 * spatial frequencies. Peak sensitivity is at 4-8 cycles/degree,
 * rolling off at both low and high frequencies (bandpass shape).
 *
 * Weights are precomputed as a LUT indexed by (block_size, coefficient_position).
 * float32 minimum — float16 causes visible banding in smooth gradients.
 * Normalized so DC-only error equals unweighted distortion.
 */

#include <math.h>
#include <string.h>
#include "csf_weights.h"

/* CSF model parameters (updated for modern displays, ~96 DPI, 60cm viewing) */
#define CSF_PEAK_FREQ  6.0f    /* cycles/degree at peak sensitivity */
#define CSF_BANDWIDTH  2.5f    /* octave bandwidth of the bandpass */

/* Precomputed weight tables */
static float csf_table_4x4[16];
static float csf_table_8x8[64];
static int csf_initialized = 0;

/* Simple CSF model: bandpass approximation.
 * w(f) = exp(-0.5 * ((log2(f/f_peak))^2) / bw^2)
 * where f is spatial frequency in cycles/degree. */
static float csf_weight_for_freq(float freq_cpd) {
    if (freq_cpd < 0.1f) return 1.0f; /* DC: full weight */
    float log_ratio = logf(freq_cpd / CSF_PEAK_FREQ) / logf(2.0f);
    float w = expf(-0.5f * (log_ratio * log_ratio) / (CSF_BANDWIDTH * CSF_BANDWIDTH));
    /* Clamp to [0.01, 1.0] — never zero (every coefficient matters a little) */
    if (w < 0.01f) w = 0.01f;
    if (w > 1.0f) w = 1.0f;
    return w;
}

static void csf_init_table(float *table, int size) {
    /* Map each coefficient position to a spatial frequency.
     * For an NxN block at ~96 DPI viewed at 60cm:
     * Nyquist frequency ≈ 30 cycles/degree.
     * Coefficient (u,v) corresponds to frequency sqrt(u^2+v^2) * Nyquist/N. */
    float nyquist_cpd = 30.0f;

    for (int v = 0; v < size; v++) {
        for (int u = 0; u < size; u++) {
            float fu = (float)u * nyquist_cpd / (float)size;
            float fv = (float)v * nyquist_cpd / (float)size;
            float freq = sqrtf(fu * fu + fv * fv);
            table[v * size + u] = csf_weight_for_freq(freq);
        }
    }

    /* Normalize: DC weight = 1.0 (CSF-weighted DC error = unweighted DC error) */
    float dc_weight = table[0];
    if (dc_weight > 0.0f) {
        for (int i = 0; i < size * size; i++)
            table[i] /= dc_weight;
    }
}

void psy_csf_init(void) {
    if (csf_initialized) return;
    csf_init_table(csf_table_4x4, 4);
    csf_init_table(csf_table_8x8, 8);
    csf_initialized = 1;
}

float psy_csf_weight(int size, int coeff_idx) {
    psy_csf_init();
    if (size == 4 && coeff_idx < 16) return csf_table_4x4[coeff_idx];
    if (size == 8 && coeff_idx < 64) return csf_table_8x8[coeff_idx];
    return 1.0f; /* fallback */
}

const float *psy_csf_table(int size) {
    psy_csf_init();
    if (size == 4) return csf_table_4x4;
    if (size == 8) return csf_table_8x8;
    return NULL;
}
