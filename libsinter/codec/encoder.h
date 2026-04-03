#ifndef SNTE_ENCODER_H
#define SNTE_ENCODER_H

#include <stdint.h>
#include "../core/sinter_types.h"

/* Convert QP to quantization step size. */
uint16_t snte_qp_to_step(uint8_t qp);

/* Encode a single luma frame.
 * ctx: encoder context (must have header and psy_cfg set).
 * input_y: input luma plane. width/height/stride: frame dimensions.
 * bitstream: output buffer. bs_capacity: buffer size.
 * bs_len: output, number of bytes written.
 * Returns 0 on success. */
int snte_encode_frame(SntEncoderCtx *ctx, const uint8_t *input_y,
                       int width, int height, int stride,
                       uint8_t *bitstream, uint32_t bs_capacity,
                       uint32_t *bs_len);

#endif
