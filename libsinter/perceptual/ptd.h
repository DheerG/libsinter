#ifndef PSY_PTD_H
#define PSY_PTD_H

#include <stdint.h>

/* Compute CSF-weighted transform-domain distortion (PTD).
 * orig: original transform coefficients. recon: reconstructed.
 * size: block size (4 or 8). */
float psy_ptd_distortion(const int16_t *orig, const int16_t *recon, int size);

/* Compute plain MSE for comparison. */
float psy_mse_distortion(const int16_t *orig, const int16_t *recon, int n);

#endif
