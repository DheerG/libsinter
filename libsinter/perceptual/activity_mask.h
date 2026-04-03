#ifndef PSY_ACTIVITY_MASK_H
#define PSY_ACTIVITY_MASK_H

#include <stdint.h>

/* Compute local activity (AC variance) for a block's coefficients.
 * coeffs: transform coefficients. n: number of coefficients. */
float psy_compute_activity(const int16_t *coeffs, int n);

/* Compute per-CTU QP offset from activity masking.
 * activity: this block's activity. mean_activity: frame average.
 * alpha: masking strength (default 0.5). */
int8_t psy_activity_qp_offset(float activity, float mean_activity, float alpha);

#endif
