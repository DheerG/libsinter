/*
 * Sinter Codec — In-loop deblocking filter (normative)
 */
#ifndef SNT_DEBLOCK_H
#define SNT_DEBLOCK_H

#include <stdint.h>
#include "../sinter_types.h"

/* Filter the internal 8x8 block edges of a w x h luma plane in place:
 * all vertical edges first, then all horizontal edges.
 * step: the frame's quantizer step (thresholds scale with it).
 * blk_skip / mv_field: per-block skip flags and MVs (blk_cols blocks per
 * row); NULL for key frames. An edge between two skipped blocks with the
 * same MV is not filtered. */
void snt_deblock_luma(uint8_t *plane, int w, int h, int stride, uint16_t step,
                      const int *blk_skip, const SntMV *mv_field, int blk_cols);

/* The same filter on the 4x4 block edges of a chroma plane, with the
 * chroma quantizer step. Chroma block (i, j) takes the skip flag and MV of
 * luma block (i, j), which covers the same area. */
void snt_deblock_chroma(uint8_t *plane, int w, int h, int stride, uint16_t step,
                        const int *blk_skip, const SntMV *mv_field, int blk_cols);

#endif /* SNT_DEBLOCK_H */
