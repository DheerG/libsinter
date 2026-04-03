#ifndef SNTD_DECODER_H
#define SNTD_DECODER_H

#include <stdint.h>
#include "../core/sinter_types.h"

/* Decode a single YUV420 frame from bitstream.
 * ctx: decoder context. bitstream/bs_len: input data.
 * output_y: decoded luma. output_u/output_v: decoded chroma (may be NULL).
 * out_width/out_height: frame dims. stride_y/stride_uv: row strides.
 * Returns 0 on success. */
int sntd_decode_frame(SntDecoderCtx *ctx, const uint8_t *bitstream,
                       uint32_t bs_len, uint8_t *output_y,
                       uint8_t *output_u, uint8_t *output_v,
                       int *out_width, int *out_height,
                       int stride_y, int stride_uv);

#endif
