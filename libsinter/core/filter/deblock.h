/*
 * Sinter Codec — In-loop deblocking filter (normative)
 * ITU-T H.263 Annex J; see deblock.c.
 */
#ifndef SNT_DEBLOCK_H
#define SNT_DEBLOCK_H

#include <stdint.h>

/* Filter the internal 8x8 block edges of a w x h luma plane in place:
 * all horizontal edges first, then all vertical edges.
 * step: the frame's quantizer step.
 * blk_skip: per-block skip flags (blk_cols blocks per row); NULL for key
 * frames. An edge between two skipped blocks is not filtered. */
void snt_deblock_luma(uint8_t *plane, int w, int h, int stride, uint16_t step,
                      const int *blk_skip, int blk_cols);

/* The same filter on the 4x4 block edges of a chroma plane, with the
 * chroma quantizer step. Chroma block (i, j) takes the skip flag of luma
 * block (i, j), which covers the same area. */
void snt_deblock_chroma(uint8_t *plane, int w, int h, int stride, uint16_t step,
                        const int *blk_skip, int blk_cols);

#endif /* SNT_DEBLOCK_H */
