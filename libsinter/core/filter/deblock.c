/*
 * Sinter Codec — In-loop deblocking filter (normative)
 *
 * ITU-T H.263 Annex J, "Deblocking filter mode" (1998). For the four
 * samples A B | C D across a block edge (A, B in the block above or to
 * the left):
 *   d  = (A - 4B + 4C - D) / 8
 *   d1 = UpDownRamp(d, STRENGTH)
 *   d2 = clipd1((A - D) / 4, d1 / 2)
 *   B' = clip(B + d1), C' = clip(C - d1), A' = A - d2, D' = D + d2
 * where UpDownRamp(x, s) = sign(x) * max(0, |x| - max(0, 2 (|x| - s)))
 * passes small steps, fades out between s and 2s, and leaves larger
 * (real) edges alone, and STRENGTH comes from the quantizer (Table J.2).
 * "/" truncates toward zero, as in H.263. An edge is filtered when at
 * least one of its two blocks is coded. Horizontal edges are filtered
 * before vertical ones.
 *
 * Changes from the source (see PROVENANCE.md): the QUANT that indexes
 * Table J.2 is a fifth of Sinter's quantizer step, rounded up. H.263's
 * quantizer step is 2 QUANT, but its dead zone is much wider than
 * Sinter's, so a given step leaves smaller edge errors here; QUANT =
 * step / 2 smoothed real detail and cost 6% (by BD-rate, a fifth and
 * a sixth were best, within 0.1%, of 1/6 .. 3/2). Table J.2 (QUANT 1..31) is continued past 31 with
 * its own pattern, STRENGTH = (QUANT + 7) / 3. Chroma is filtered on its
 * 4x4 transform-block edges.
 */

#include "deblock.h"

/* H.263 Table J.2: STRENGTH for QUANT = 1 .. 31 */
static const uint8_t strength_j2[32] = {
    0,  1,  1,  2,  2,  3,  3,  4,  4,  4,  5,  5,  6,  6,  7,  7,
    7,  8,  8,  8,  9,  9,  9, 10, 10, 10, 11, 11, 11, 12, 12, 12
};

static int strength_for_step(uint16_t step) {
    int quant = (step + 4) / 5;
    if (quant < 1) quant = 1;
    return quant <= 31 ? strength_j2[quant] : (quant + 7) / 3;
}

static inline int iabs(int v) { return v < 0 ? -v : v; }
static inline int imax(int a, int b) { return a > b ? a : b; }
static inline uint8_t clip8(int v) { return (uint8_t)(v < 0 ? 0 : (v > 255 ? 255 : v)); }

static inline int up_down_ramp(int x, int s) {
    int m = imax(0, iabs(x) - imax(0, 2 * (iabs(x) - s)));
    return x < 0 ? -m : m;
}

/* Filter one edge position; p points at C, step is the distance across
 * the edge */
static inline void filter_edge(uint8_t *p, int step, int strength) {
    int a = p[-2 * step], b = p[-step], c = p[0], d = p[step];
    int dd = (a - 4 * b + 4 * c - d) / 8;
    int d1 = up_down_ramp(dd, strength);
    if (d1 == 0) return;
    int lim = iabs(d1 / 2);
    int d2 = (a - d) / 4;
    d2 = d2 < -lim ? -lim : (d2 > lim ? lim : d2);
    p[-2 * step] = clip8(a - d2);
    p[-step] = clip8(b + d1);
    p[0] = clip8(c - d1);
    p[step] = clip8(d + d2);
}

/* Is the edge between blocks a and b filtered? Not if both are skipped. */
static inline int edge_active(const int *blk_skip, int a, int b) {
    return !blk_skip || !(blk_skip[a] && blk_skip[b]);
}

/* Filter a plane whose blocks are bs x bs */
static void deblock_plane(uint8_t *plane, int w, int h, int stride, uint16_t step,
                          int bs, const int *blk_skip, int blk_cols) {
    int s = strength_for_step(step);

    /* Horizontal edges (filter across y); D must lie inside the plane */
    for (int y = bs; y + 1 < h; y += bs)
        for (int x = 0; x < w; x++) {
            int blk = (y / bs) * blk_cols + x / bs;
            if (edge_active(blk_skip, blk - blk_cols, blk))
                filter_edge(&plane[y * stride + x], stride, s);
        }

    /* Vertical edges (filter across x) */
    for (int x = bs; x + 1 < w; x += bs)
        for (int y = 0; y < h; y++) {
            int blk = (y / bs) * blk_cols + x / bs;
            if (edge_active(blk_skip, blk - 1, blk))
                filter_edge(&plane[y * stride + x], 1, s);
        }
}

void snt_deblock_luma(uint8_t *plane, int w, int h, int stride, uint16_t step,
                      const int *blk_skip, int blk_cols) {
    deblock_plane(plane, w, h, stride, step, 8, blk_skip, blk_cols);
}

void snt_deblock_chroma(uint8_t *plane, int w, int h, int stride, uint16_t step,
                        const int *blk_skip, int blk_cols) {
    deblock_plane(plane, w, h, stride, step, 4, blk_skip, blk_cols);
}
