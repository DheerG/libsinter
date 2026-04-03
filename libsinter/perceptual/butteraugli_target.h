#ifndef PSY_BUTTERAUGLI_TARGET_H
#define PSY_BUTTERAUGLI_TARGET_H

#include <stdint.h>

/* Compute perceptual distance between original and reconstructed frame.
 * Returns score: 0 = identical, ~1.0 = just noticeable, higher = worse. */
float psy_perceptual_distance(const uint8_t *orig, const uint8_t *recon,
                               int width, int height, int stride);

/* Callback for encoding at a specific QP. */
typedef void (*PsyEncodeFunc)(void *ctx, const uint8_t *input,
                               int width, int height, int stride,
                               uint8_t qp, uint8_t *recon,
                               uint8_t *bitstream, uint32_t bs_cap,
                               uint32_t *bs_len);

/* QP bisection for constant-quality encoding.
 * Returns optimal QP that achieves target perceptual quality. */
uint8_t psy_butteraugli_bisect(const uint8_t *orig, int width, int height,
                                int stride, float target,
                                PsyEncodeFunc encode_func, void *enc_ctx,
                                uint8_t *recon_buf, uint8_t *bs_buf,
                                uint32_t bs_cap);

#endif
