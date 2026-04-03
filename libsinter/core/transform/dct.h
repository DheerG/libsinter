/*
 * Sinter Codec — Integer DCT Forward/Inverse
 * See dct.c for implementation details.
 */
#ifndef SNT_DCT_H
#define SNT_DCT_H

#include <stdint.h>

/* Forward and inverse 4-point DCT (1D).
 * in/out: arrays of 4 int16_t values.
 * Operates in-place if in == out. */
void snt_dct_fwd_4_ref(const int16_t *in, int16_t *out);
void snt_dct_inv_4_ref(const int16_t *in, int16_t *out);

/* Forward and inverse 8-point DCT (1D). */
void snt_dct_fwd_8_ref(const int16_t *in, int16_t *out);
void snt_dct_inv_8_ref(const int16_t *in, int16_t *out);

/* Generic N-point forward/inverse 1D DCT (supports 4, 8, 16, 32). */
void snt_dct_fwd_n_ref(const int16_t *in, int16_t *out, int size);
void snt_dct_inv_n_ref(const int16_t *in, int16_t *out, int size);

/* Forward and inverse 2D DCT on a block.
 * block: pointer to top-left of block in a contiguous buffer.
 * stride: row stride of the buffer.
 * size: block width/height (must be 4, 8, 16, or 32).
 * coeffs: output coefficient buffer (size*size int16_t, contiguous). */
void snt_dct_fwd_2d_ref(const int16_t *block, int stride, int size,
                         int16_t *coeffs);
void snt_dct_inv_2d_ref(const int16_t *coeffs, int size,
                         int16_t *block, int stride);

#endif /* SNT_DCT_H */
