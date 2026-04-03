/*
 * Sinter Codec — Tile Management
 * Implements multi-tile parallelism with 4-pixel overlap blending.
 *
 * Each tile has independent entropy state (CDF reset at tile boundary).
 * The 4-pixel overlap zone at tile boundaries is decoded redundantly
 * by both adjacent tiles and blended with raised cosine weights.
 */
#ifndef SNT_TILE_H
#define SNT_TILE_H

#include <stdint.h>
#include "../sinter_types.h"

#define SNT_TILE_OVERLAP 4

/* Compute tile grid from frame dimensions and tile counts. */
void snt_tile_grid(int frame_w, int frame_h,
                    int num_cols, int num_rows,
                    SntTileInfo *tiles);

/* Blend overlap zones between horizontally adjacent tiles. */
void snt_tile_blend_horizontal(int16_t *recon, int stride,
                                const SntTileInfo *left, const SntTileInfo *right);

/* Blend overlap zones between vertically adjacent tiles. */
void snt_tile_blend_vertical(int16_t *recon, int stride,
                              const SntTileInfo *above, const SntTileInfo *below);

/* Get total number of tiles. */
static inline int snt_tile_count(int num_cols, int num_rows) {
    return num_cols * num_rows;
}

#endif /* SNT_TILE_H */
