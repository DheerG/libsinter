/*
 * Sinter Codec — Spatial intra prediction (normative)
 */
#ifndef SNT_INTRA_SPATIAL_H
#define SNT_INTRA_SPATIAL_H

#include <stdint.h>

/* Intra modes. Mode 0 is the transform-domain DC/low-AC prediction done
 * by the codec itself; the others predict pixels from the reconstructed
 * row above and column left of the block. */
enum {
    SNT_IPRED_TD  = 0,  /* transform-domain DC/AC prediction */
    SNT_IPRED_V   = 1,  /* copy the row above downwards */
    SNT_IPRED_H   = 2,  /* copy the left column across */
    SNT_IPRED_AVG = 3,  /* mean of the V and H predictions */
    SNT_NUM_IPRED = 4
};

/* Build the 8x8 prediction for spatial mode `mode` (not SNT_IPRED_TD) of
 * the block at (x0, y0). recon: reconstructed plane so far (stride),
 * clipped to 0..255 on read; neighbours outside the frame read as 128.
 * pred: 64 samples, stride 8. */
void snt_intra_spatial_pred8(const int16_t *recon, int stride, int x0, int y0,
                             int mode, uint8_t *pred);

#endif /* SNT_INTRA_SPATIAL_H */
