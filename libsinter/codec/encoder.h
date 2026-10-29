#ifndef SNTE_ENCODER_H
#define SNTE_ENCODER_H

#include <stdint.h>
#include "../core/sinter_types.h"

/* Encode a single YUV420 frame (luma + chroma).
 * ctx: encoder context (must have header and psy_cfg set).
 * input_y: input luma plane. input_u/input_v: chroma planes (may be NULL for luma-only).
 * width/height: frame dimensions. stride_y: luma row stride. stride_uv: chroma row stride.
 * bitstream: output buffer. bs_capacity: buffer size.
 * bs_len: output, number of bytes written.
 * Returns 0 on success. */
int snte_encode_frame(SntEncoderCtx *ctx, const uint8_t *input_y,
                       const uint8_t *input_u, const uint8_t *input_v,
                       int width, int height, int stride_y, int stride_uv,
                       uint8_t *bitstream, uint32_t bs_capacity,
                       uint32_t *bs_len);

#ifdef SNT_TRAIN_STATS
#include <stdio.h>
/* Training builds: write the symbol counts collected so far */
void snte_train_dump(FILE *f);
/* Training builds: ideal bits per model set and model, and other bits */
void snte_bits_dump(FILE *f);
#endif

#endif
