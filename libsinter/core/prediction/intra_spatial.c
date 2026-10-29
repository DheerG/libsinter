/*
 * Sinter Codec — Spatial intra prediction (normative)
 *
 * Pixel-domain predictors for 8x8 intra blocks from the reconstructed
 * neighbours (before the deblocking filter, which runs on whole frames).
 * Vertical and horizontal extrapolation and their mean are textbook
 * predictors (DPCM-style edge extension).
 */

#include "intra_spatial.h"
#include "../../util/fixed_point.h"

void snt_intra_spatial_pred8(const int16_t *recon, int stride, int x0, int y0,
                             int mode, uint8_t *pred) {
    uint8_t above[8], left[8];
    for (int i = 0; i < 8; i++) {
        above[i] = (y0 > 0) ? snt_clip_pixel(recon[(y0 - 1) * stride + x0 + i]) : 128;
        left[i]  = (x0 > 0) ? snt_clip_pixel(recon[(y0 + i) * stride + x0 - 1]) : 128;
    }
    for (int y = 0; y < 8; y++)
        for (int x = 0; x < 8; x++) {
            int v;
            switch (mode) {
            case SNT_IPRED_V: v = above[x]; break;
            case SNT_IPRED_H: v = left[y]; break;
            default:          v = (above[x] + left[y] + 1) >> 1; break;
            }
            pred[y * 8 + x] = (uint8_t)v;
        }
}
