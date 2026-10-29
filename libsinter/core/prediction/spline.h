/*
 * Sinter Codec — Cubic B-spline interpolation of reference frames (normative)
 *
 * M. Unser, A. Aldroubi and M. Eden, "Fast B-spline transforms for
 * continuous image representation and interpolation", IEEE Trans. PAMI
 * 13(3), 1991; and "B-spline signal processing, Part II: efficient design
 * and applications", IEEE Trans. Signal Processing 41(2), 1993.
 *
 * A reference frame is turned once into cubic B-spline coefficients by a
 * recursive prefilter (one causal and one anticausal first-order pass per
 * row, then per column, pole sqrt(3) - 2). The B-spline through those
 * coefficients passes exactly through the frame's samples, and any
 * fractional position is then a 4 x 4 weighted sum of coefficients.
 *
 * Changes from the sources (see PROVENANCE.md): fixed-point arithmetic
 * (pole and boundary constants in Q16, coefficients in Q8) so encoder
 * and decoder agree exactly; the frame is first extended by repeating
 * its edge samples SNT_SPLINE_PAD pixels on each side, and the mirror
 * boundary of the sources applies at the edge of that extension; a
 * motion vector with no fractional part copies the reference samples,
 * which is what the exact B-spline gives there.
 */
#ifndef SNT_SPLINE_H
#define SNT_SPLINE_H

#include <stdint.h>

#define SNT_SPLINE_PAD 32  /* edge extension on each side, pixels */

/* A reference plane and its B-spline coefficients */
typedef struct {
    const uint8_t *pix;   /* w x h samples, stride w */
    const int32_t *coef;  /* Q8 coefficients of the extended plane,
                             (w + 2 PAD) x (h + 2 PAD) */
    int w, h;
} SntSplinePlane;

/* Coefficients of the w x h plane pix (stride w), in a new buffer the
 * caller frees with free(). NULL if out of memory. */
int32_t *snt_spline_coefs(const uint8_t *pix, int w, int h);

/* Interpolate a bw x bh block at quarter-pel position (qx, qy) of its
 * top-left sample (quarter-pel units). dst: stride bw. */
void snt_spline_interp(const SntSplinePlane *p, int qx, int qy, int bw, int bh,
                       uint8_t *dst);

#endif /* SNT_SPLINE_H */
