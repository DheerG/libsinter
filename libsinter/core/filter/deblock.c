/*
 * Sinter Codec — In-loop deblocking filter (normative)
 *
 * Smooths block edges (8x8 luma, 4x4 chroma) where the step across the edge is small enough
 * to be a quantization artifact rather than image content. It runs on the
 * reconstructed planes before the frame is output and stored as the
 * reference, so encoder and decoder predict from the same filtered frame.
 *
 * For samples p1 p0 | q0 q1 across an edge, filter when
 *   |p0 - q0| < alpha, |p1 - p0| < beta and |q1 - q0| < beta,
 * with delta = (4 (q0 - p0) + (p1 - q1)) / 8 clipped to [-tc, tc]
 * (the classic 4-tap edge filter, H.263 Annex J lineage):
 *   p0 += delta, q0 -= delta.
 * alpha, beta and tc scale with the quantizer step.
 */

#include "deblock.h"

/* Thresholds in 1/16ths of the quantizer step, tuned by BD-rate on Xiph
 * CIF clips; stronger settings start smoothing real texture. */
#define DB_ALPHA 24
#define DB_BETA   6
#define DB_TC     4

static inline int iabs(int v) { return v < 0 ? -v : v; }

/* Round-half-away division by 8 (no implementation-defined shifts) */
static inline int div8(int v) { return v >= 0 ? (v + 4) / 8 : -((-v + 4) / 8); }

static inline uint8_t clip8(int v) { return (uint8_t)(v < 0 ? 0 : (v > 255 ? 255 : v)); }

/* Filter one edge position; s points at q0, d is the step across the edge */
static inline void filter_edge(uint8_t *s, int d, int alpha, int beta, int tc) {
    int p1 = s[-2 * d], p0 = s[-d], q0 = s[0], q1 = s[d];
    if (iabs(p0 - q0) >= alpha || iabs(p1 - p0) >= beta || iabs(q1 - q0) >= beta)
        return;
    int delta = div8(4 * (q0 - p0) + (p1 - q1));
    if (delta > tc) delta = tc;
    if (delta < -tc) delta = -tc;
    s[-d] = clip8(p0 + delta);
    s[0] = clip8(q0 - delta);
}

/* Is the edge between blocks a and b filtered? */
static inline int edge_active(const int *blk_skip, const SntMV *mv, int a, int b) {
    if (!blk_skip) return 1;
    return !(blk_skip[a] && blk_skip[b] && mv[a].x == mv[b].x && mv[a].y == mv[b].y);
}

/* Filter a plane whose blocks are bs x bs */
static void deblock_plane(uint8_t *plane, int w, int h, int stride, uint16_t step,
                          int bs, const int *blk_skip, const SntMV *mv_field, int blk_cols) {
    int alpha = step * DB_ALPHA / 16, beta = step * DB_BETA / 16, tc = step * DB_TC / 16;
    if (alpha < 1 || beta < 1 || tc < 1)
        return;  /* fine quantizer: nothing to smooth */

    /* Vertical edges (filter across x); q1 must lie inside the plane */
    for (int x = bs; x + 1 < w; x += bs)
        for (int y = 0; y < h; y++) {
            int blk = (y / bs) * blk_cols + x / bs;
            if (edge_active(blk_skip, mv_field, blk - 1, blk))
                filter_edge(&plane[y * stride + x], 1, alpha, beta, tc);
        }

    /* Horizontal edges (filter across y) */
    for (int y = bs; y + 1 < h; y += bs)
        for (int x = 0; x < w; x++) {
            int blk = (y / bs) * blk_cols + x / bs;
            if (edge_active(blk_skip, mv_field, blk - blk_cols, blk))
                filter_edge(&plane[y * stride + x], stride, alpha, beta, tc);
        }
}

void snt_deblock_luma(uint8_t *plane, int w, int h, int stride, uint16_t step,
                      const int *blk_skip, const SntMV *mv_field, int blk_cols) {
    deblock_plane(plane, w, h, stride, step, 8, blk_skip, mv_field, blk_cols);
}

void snt_deblock_chroma(uint8_t *plane, int w, int h, int stride, uint16_t step,
                        const int *blk_skip, const SntMV *mv_field, int blk_cols) {
    deblock_plane(plane, w, h, stride, step, 4, blk_skip, mv_field, blk_cols);
}
