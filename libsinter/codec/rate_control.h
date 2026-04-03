/*
 * Sinter Codec — Rate Control
 * CRF (Constant Rate Factor) and CBR modes.
 */
#ifndef SNTE_RATE_CONTROL_H
#define SNTE_RATE_CONTROL_H

#include <stdint.h>

typedef struct {
    float   crf;           /* CRF target (0=lossless, 51=worst) */
    float   bitrate;       /* target bitrate in kbps (0=CRF mode) */
    float   fps;           /* frame rate for bitrate calculation */
    /* Internal state */
    float   qp_smooth;     /* temporally smoothed QP */
    float   complexity_avg; /* running average frame complexity */
    int     frame_count;
} SntRateCtrl;

/* Initialize rate control. crf=0 means use fixed QP. bitrate=0 means CRF. */
void snte_rc_init(SntRateCtrl *rc, float crf, float bitrate, float fps);

/* Compute QP for the next frame given its estimated complexity.
 * complexity: proxy for frame difficulty (e.g., mean activity).
 * Returns QP in [0..63]. */
uint8_t snte_rc_get_qp(SntRateCtrl *rc, float complexity);

/* Update rate control after encoding a frame.
 * actual_bits: number of bits used. */
void snte_rc_update(SntRateCtrl *rc, uint32_t actual_bits, float complexity);

#endif
