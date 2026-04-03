/*
 * Sinter Codec — Rate Control
 * Implements CRF and CBR rate control for the encoder.
 *
 * CRF mode: QP = CRF + log2(complexity / avg_complexity) * scale
 * CBR mode: QP adjusted per frame to hit target bitrate using
 *           a leaky bucket model.
 *
 * Temporal QP smoothing prevents flicker from frame-to-frame
 * QP oscillation (Jens-Rainer's requirement).
 */

#include <math.h>
#include "rate_control.h"
#include "../core/quantize/qp_table.h"

#define RC_SMOOTH_FACTOR 0.8f  /* temporal QP smoothing (0.8 = 80% old + 20% new) */

void snte_rc_init(SntRateCtrl *rc, float crf, float bitrate, float fps) {
    rc->crf = crf;
    rc->bitrate = bitrate;
    rc->fps = (fps > 0) ? fps : 30.0f;
    rc->qp_smooth = crf;
    rc->complexity_avg = 1.0f;
    rc->frame_count = 0;
}

uint8_t snte_rc_get_qp(SntRateCtrl *rc, float complexity) {
    float target_qp;

    if (rc->bitrate > 0) {
        /* CBR mode: adjust QP to hit target bits per frame.
         * bits_per_frame = bitrate * 1000 / fps.
         * If complexity > avg, need more bits -> lower QP.
         * QP = base_qp + alpha * log2(complexity / avg) */
        float bits_per_frame = rc->bitrate * 1000.0f / rc->fps;
        (void)bits_per_frame; /* TODO: actual buffer model */
        float ratio = (rc->complexity_avg > 0.1f) ?
                       complexity / rc->complexity_avg : 1.0f;
        target_qp = rc->crf + 3.0f * log2f(ratio + 0.001f);
    } else {
        /* CRF mode: constant quality factor with complexity adaptation.
         * Higher complexity -> slightly higher QP to maintain visual quality. */
        float ratio = (rc->complexity_avg > 0.1f) ?
                       complexity / rc->complexity_avg : 1.0f;
        target_qp = rc->crf + 2.0f * log2f(ratio + 0.001f);
    }

    /* Temporal smoothing: 80% previous + 20% current (prevents flicker) */
    if (rc->frame_count > 0) {
        rc->qp_smooth = RC_SMOOTH_FACTOR * rc->qp_smooth +
                        (1.0f - RC_SMOOTH_FACTOR) * target_qp;
    } else {
        rc->qp_smooth = target_qp;
    }

    /* Clamp to valid range */
    int qp = (int)(rc->qp_smooth + 0.5f);
    if (qp < 0) qp = 0;
    if (qp > 63) qp = 63;

    return (uint8_t)qp;
}

void snte_rc_update(SntRateCtrl *rc, uint32_t actual_bits, float complexity) {
    (void)actual_bits; /* used in full CBR model */

    /* Update running complexity average (exponential moving average) */
    if (rc->frame_count == 0) {
        rc->complexity_avg = complexity;
    } else {
        rc->complexity_avg = 0.9f * rc->complexity_avg + 0.1f * complexity;
    }
    rc->frame_count++;
}
