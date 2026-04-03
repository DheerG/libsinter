/*
 * Sinter Codec — Lapped Transform Pre/Post Filters
 * See lapping.c for implementation details.
 */
#ifndef SNT_LAPPING_H
#define SNT_LAPPING_H

#include <stdint.h>

/* 4-point pre/post filter pair (2 lifting steps).
 * data: pointer to first of 4 samples
 * stride: distance between samples (1 for horizontal, row_stride for vertical) */
void snt_lap_pre_filter_4_ref(int16_t *data, int stride);
void snt_lap_post_filter_4_ref(int16_t *data, int stride);

/* 8-point pre/post filter pair (4 lifting steps). */
void snt_lap_pre_filter_8_ref(int16_t *data, int stride);
void snt_lap_post_filter_8_ref(int16_t *data, int stride);

/* Apply pre/post filter along a row at every block boundary.
 * row: pointer to start of row. width: row length. blk_size: block dimension.
 * Frame edges (0 and width) are skipped (identity filter). */
void snt_lap_pre_filter_row_ref(int16_t *row, int width, int blk_size);
void snt_lap_post_filter_row_ref(int16_t *row, int width, int blk_size);

/* Apply 2D lapping to an entire tile/frame buffer.
 * Pre: horizontal first, then vertical.
 * Post: vertical first, then horizontal (exact inverse). */
void snt_lap_pre_filter_2d_ref(int16_t *buf, int width, int height,
                                int stride, int blk_size);
void snt_lap_post_filter_2d_ref(int16_t *buf, int width, int height,
                                 int stride, int blk_size);

#endif /* SNT_LAPPING_H */
