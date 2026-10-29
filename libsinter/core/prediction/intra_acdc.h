/*
 * Sinter Codec — Intra AC/DC coefficient prediction (normative)
 *
 * ITU-T H.263 Annex I, "Advanced INTRA coding mode" (1998). An intra
 * block's DCT coefficients are predicted from its neighbours' in one of
 * three modes, chosen by the encoder and sent with the block:
 *   DC   - the DC from the mean of the left and above blocks' DC
 *   VERT - the DC and the first coefficient row from the block above
 *   HORZ - the DC and the first coefficient column from the block left
 * The same family of prediction is in MPEG-4 Part 2 (1999).
 *
 * Changes from the source (see PROVENANCE.md): the mode is sent per 8x8
 * luma block rather than per macroblock; the zigzag scan is kept for all
 * modes (Annex I switches to alternate scans); a neighbour coded in
 * inter mode contributes the coefficients of its reconstruction instead
 * of counting as unavailable; in DC mode with one neighbour outside the
 * frame, the other is used alone.
 */
#ifndef SNT_INTRA_ACDC_H
#define SNT_INTRA_ACDC_H

#include <stdint.h>
#include <string.h>

enum {
    SNT_IPRED_DC   = 0,
    SNT_IPRED_VERT = 1,
    SNT_IPRED_HORZ = 2,
    SNT_NUM_IPRED  = 3
};

/* Per-block state kept for the neighbours' prediction: the DC, the first
 * row's AC (raster indices 1..7) and the first column's AC (8, 16 .. 56)
 * of the block's reconstruction, on the orthonormal 8x8 DCT's scale. */
#define SNT_ACDC_STATE 15

/* A neighbour outside the frame: mid-grey DC, no AC */
#define SNT_ACDC_GREY_DC 1024

static inline void snt_acdc_store(const int16_t *coeffs, int16_t *st) {
    st[0] = coeffs[0];
    for (int i = 1; i < 8; i++) {
        st[i] = coeffs[i];
        st[7 + i] = coeffs[8 * i];
    }
}

/* Coefficient prediction (64, raster) for an intra block in `mode`, from
 * the left and above neighbours' states (NULL when outside the frame). */
static inline void snt_acdc_predict(int mode, const int16_t *left, const int16_t *above,
                                    int16_t *pred) {
    memset(pred, 0, 64 * sizeof(int16_t));
    if (mode == SNT_IPRED_VERT) {
        pred[0] = above ? above[0] : SNT_ACDC_GREY_DC;
        if (above)
            for (int i = 1; i < 8; i++) pred[i] = above[i];
    } else if (mode == SNT_IPRED_HORZ) {
        pred[0] = left ? left[0] : SNT_ACDC_GREY_DC;
        if (left)
            for (int i = 1; i < 8; i++) pred[8 * i] = left[7 + i];
    } else if (left && above) {
        int s = left[0] + above[0] + 1;  /* halve, rounding down */
        pred[0] = (int16_t)(s >= 0 ? s / 2 : -((-s + 1) / 2));
    } else {
        pred[0] = left ? left[0] : (above ? above[0] : SNT_ACDC_GREY_DC);
    }
}

/* MPEG-4 Part 2's choice of prediction direction (ISO/IEC 14496-2, 1999,
 * intra DC and AC prediction): with the DC of the left (A), above-left (B)
 * and above (C) blocks, predict from above if |A - B| < |B - C|, else from
 * the left. Returns 1 (above, as SNT_IPRED_VERT), 2 (left, as
 * SNT_IPRED_HORZ), or 0 if a neighbour is outside the frame. Sinter sends
 * the mode; this only selects the mode's model (cdf_set.h). */
static inline int snt_acdc_direction(const int16_t *left, const int16_t *aboveleft,
                                     const int16_t *above) {
    if (!left || !aboveleft || !above) return 0;
    int ab = left[0] - aboveleft[0], bc = aboveleft[0] - above[0];
    return (ab < 0 ? -ab : ab) < (bc < 0 ? -bc : bc) ? 1 : 2;
}

#endif /* SNT_INTRA_ACDC_H */
