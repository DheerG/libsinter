#ifndef SNT_GRAIN_SYNTH_H
#define SNT_GRAIN_SYNTH_H

#include <stdint.h>
#include "../sinter_types.h"

/* Generate grain LUTs from parameters. */
void snt_grain_generate_lut(const SntGrainParams *params,
                             int16_t *lut_y, int16_t *lut_u, int16_t *lut_v);

/* Apply grain synthesis to a decoded frame (last post-processing step). */
void snt_grain_apply(const SntGrainParams *params,
                      uint8_t *frame_y, uint8_t *frame_u, uint8_t *frame_v,
                      int width, int height, int stride_y, int stride_uv);

#endif
