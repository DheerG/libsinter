#ifndef SNTD_DECODER_H
#define SNTD_DECODER_H

#include <stdint.h>
#include "../core/sinter_types.h"

/* Decode a single frame from bitstream.
 * ctx: decoder context. bitstream/bs_len: input data.
 * output_y: decoded luma output. out_width/out_height: frame dims.
 * stride: output row stride.
 * Returns 0 on success. */
int sntd_decode_frame(SntDecoderCtx *ctx, const uint8_t *bitstream,
                       uint32_t bs_len, uint8_t *output_y,
                       int *out_width, int *out_height, int stride);

#endif
