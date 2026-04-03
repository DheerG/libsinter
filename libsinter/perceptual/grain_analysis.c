/*
 * Sinter Codec — Film Grain Analysis (Encoder-Side)
 * Implements grain parameter estimation from source video.
 *
 * Estimates AR model on the original source BEFORE encoding.
 * Placeholder for full implementation — the reference encoder
 * uses default grain parameters for now.
 */

#include <stdint.h>
#include "../core/sinter_types.h"

void psy_estimate_grain(const uint8_t *frame_y, int width, int height,
                         int stride, SntGrainParams *out) {
    /* Placeholder: estimate grain strength from high-frequency energy.
     * A full implementation would:
     * 1. Apply mild Wiener denoise
     * 2. Compute HF residual
     * 3. Fit AR(3)/AR(4) model to residual
     * For now, set reasonable defaults. */
    (void)frame_y; (void)width; (void)height; (void)stride;

    out->enabled = 0; /* disabled by default until full estimation is implemented */
    out->ar_coeff_y[0] = 128; /* Q8: 0.5 */
    out->ar_coeff_y[1] = 64;  /* Q8: 0.25 */
    out->ar_coeff_y[2] = 32;  /* Q8: 0.125 */
    out->ar_coeff_u[0] = 96;
    out->ar_coeff_u[1] = 48;
    out->ar_coeff_u[2] = 24;
    out->ar_coeff_u[3] = 64;  /* luma correlation tap */
    out->ar_coeff_v[0] = 96;
    out->ar_coeff_v[1] = 48;
    out->ar_coeff_v[2] = 24;
    out->ar_coeff_v[3] = 64;
    out->strength_y = 32;
    out->strength_u = 16;
    out->strength_v = 16;
    out->seed = 0;
}
