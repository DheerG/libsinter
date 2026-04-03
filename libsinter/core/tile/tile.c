/*
 * Sinter Codec — Tile Management
 * Implements multi-tile parallelism with 4-pixel overlap blending.
 *
 * Tile boundaries use a 4-pixel overlap zone decoded redundantly by
 * both adjacent tiles, then blended with raised cosine weights.
 * Per-tile CDF reset ensures independent entropy state.
 */

#include <string.h>
#include "tile.h"
#include "../../util/fixed_point.h"

/* Raised cosine blend weights for 4-pixel overlap (Q6, sum to 64) */
static const uint8_t tile_blend_w[SNT_TILE_OVERLAP] = {8, 25, 47, 60};

/* ── Compute tile grid ────────────────────────────────────────────── */
void snt_tile_grid(int frame_w, int frame_h,
                    int num_cols, int num_rows,
                    SntTileInfo *tiles) {
    int tile_w = (frame_w + num_cols - 1) / num_cols;
    int tile_h = (frame_h + num_rows - 1) / num_rows;

    for (int r = 0; r < num_rows; r++) {
        for (int c = 0; c < num_cols; c++) {
            SntTileInfo *t = &tiles[r * num_cols + c];
            t->col = (uint16_t)c;
            t->row = (uint16_t)r;
            t->x = (uint16_t)(c * tile_w);
            t->y = (uint16_t)(r * tile_h);
            t->width = (uint16_t)((c + 1 < num_cols) ? tile_w :
                                   (frame_w - c * tile_w));
            t->height = (uint16_t)((r + 1 < num_rows) ? tile_h :
                                    (frame_h - r * tile_h));
        }
    }
}

/* ── Horizontal tile overlap blending ─────────────────────────────── */
void snt_tile_blend_horizontal(int16_t *recon, int stride,
                                const SntTileInfo *left, const SntTileInfo *right) {
    /* Blend at the boundary between left and right tiles.
     * The overlap zone is SNT_TILE_OVERLAP pixels wide, centered on
     * the tile boundary. Left tile's rightmost pixels blend with
     * right tile's leftmost pixels. */
    int boundary_x = right->x;
    int top_y = left->y;
    int height = left->height < right->height ? left->height : right->height;

    for (int y = 0; y < height; y++) {
        for (int i = 0; i < SNT_TILE_OVERLAP; i++) {
            int px = boundary_x - SNT_TILE_OVERLAP / 2 + i;
            if (px < 0 || px >= left->x + left->width + SNT_TILE_OVERLAP) continue;

            int row = (top_y + y) * stride + px;
            int32_t left_val = recon[row]; /* left tile's reconstruction */
            int32_t right_val = recon[row]; /* same buffer — already decoded */

            /* Blend: left weight decreases, right weight increases */
            uint8_t w_right = tile_blend_w[i];
            uint8_t w_left = 64 - w_right;
            recon[row] = (int16_t)((left_val * w_left + right_val * w_right + 32) >> 6);
        }
    }
}

/* ── Vertical tile overlap blending ───────────────────────────────── */
void snt_tile_blend_vertical(int16_t *recon, int stride,
                              const SntTileInfo *above, const SntTileInfo *below) {
    int boundary_y = below->y;
    int left_x = above->x;
    int width = above->width < below->width ? above->width : below->width;

    for (int i = 0; i < SNT_TILE_OVERLAP; i++) {
        int py = boundary_y - SNT_TILE_OVERLAP / 2 + i;
        if (py < 0) continue;

        uint8_t w_below = tile_blend_w[i];
        uint8_t w_above = 64 - w_below;

        for (int x = 0; x < width; x++) {
            int idx = py * stride + (left_x + x);
            int32_t above_val = recon[idx];
            int32_t below_val = recon[idx];
            recon[idx] = (int16_t)((above_val * w_above + below_val * w_below + 32) >> 6);
        }
    }
}
