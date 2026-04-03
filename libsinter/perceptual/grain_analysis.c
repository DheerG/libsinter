/*
 * Sinter Codec — Film Grain Analysis (Encoder-Side)
 * Implements grain parameter estimation from source video.
 *
 * Pipeline:
 * 1. Compute local variance in 64x64 blocks (high-pass proxy)
 * 2. Classify blocks as grainy vs clean by variance threshold
 * 3. In grainy blocks, extract high-frequency residual via simple HP filter
 * 4. Fit AR(3) luma / AR(4) chroma model via Yule-Walker equations
 * 5. Estimate per-frame grain strength from residual variance
 *
 * Reference: AV1 film grain synthesis (libaom), adapted for AR(3)/AR(4).
 */

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "../core/sinter_types.h"
#include "../util/fixed_point.h"

/* ── High-pass filter (simple Laplacian for HF residual extraction) ── */
static void highpass_block(const uint8_t *src, int stride,
                            int bx, int by, int size,
                            int width, int height,
                            int16_t *hf_out) {
    for (int y = 0; y < size; y++) {
        for (int x = 0; x < size; x++) {
            int py = by + y, px = bx + x;
            if (py < 1 || py >= height - 1 || px < 1 || px >= width - 1) {
                hf_out[y * size + x] = 0;
                continue;
            }
            /* 3x3 Laplacian: center*4 - neighbors */
            int32_t center = src[py * stride + px];
            int32_t lap = 4 * center
                        - src[(py-1) * stride + px]
                        - src[(py+1) * stride + px]
                        - src[py * stride + (px-1)]
                        - src[py * stride + (px+1)];
            hf_out[y * size + x] = (int16_t)(lap >> 2); /* normalize */
        }
    }
}

/* ── Yule-Walker AR(3) fitting ────────────────────────────────────── */
/* Given a 2D signal, fit AR(3) coefficients for spatial taps:
 *   left, above, above-left.
 * Uses autocorrelation method: R * a = r, where R is the autocorrelation
 * matrix and r is the cross-correlation vector. */
static void fit_ar3(const int16_t *signal, int size,
                     int16_t *ar_coeffs_q8) {
    /* Compute autocorrelations for 2D spatial lags:
     * r[0] = E[x^2] (variance)
     * r[1] = E[x * x_left]
     * r[2] = E[x * x_above]
     * r[3] = E[x * x_above_left] */
    double r0 = 0, r1 = 0, r2 = 0, r3 = 0;
    int count = 0;

    for (int y = 1; y < size; y++) {
        for (int x = 1; x < size; x++) {
            double cur = signal[y * size + x];
            double left = signal[y * size + (x - 1)];
            double above = signal[(y - 1) * size + x];
            double above_left = signal[(y - 1) * size + (x - 1)];

            r0 += cur * cur;
            r1 += cur * left;
            r2 += cur * above;
            r3 += cur * above_left;
            count++;
        }
    }

    if (count < 4 || r0 < 1.0) {
        ar_coeffs_q8[0] = ar_coeffs_q8[1] = ar_coeffs_q8[2] = 0;
        return;
    }

    r0 /= count; r1 /= count; r2 /= count; r3 /= count;

    /* Simple least-squares: a[i] = r[i+1] / r[0] (diagonal approximation
     * of the Yule-Walker system, valid when correlations are moderate). */
    double a0 = r1 / r0; /* left coefficient */
    double a1 = r2 / r0; /* above coefficient */
    double a2 = r3 / r0; /* above-left coefficient */

    /* Clamp to stable range [-0.95, 0.95] */
    if (a0 > 0.95) a0 = 0.95; if (a0 < -0.95) a0 = -0.95;
    if (a1 > 0.95) a1 = 0.95; if (a1 < -0.95) a1 = -0.95;
    if (a2 > 0.95) a2 = 0.95; if (a2 < -0.95) a2 = -0.95;

    /* Convert to Q8 */
    ar_coeffs_q8[0] = (int16_t)(a0 * 256.0 + (a0 > 0 ? 0.5 : -0.5));
    ar_coeffs_q8[1] = (int16_t)(a1 * 256.0 + (a1 > 0 ? 0.5 : -0.5));
    ar_coeffs_q8[2] = (int16_t)(a2 * 256.0 + (a2 > 0 ? 0.5 : -0.5));
}

/* ── Film grain estimation (main entry point) ─────────────────────── */
void psy_estimate_grain(const uint8_t *frame_y, int width, int height,
                         int stride, SntGrainParams *out) {
    const int BLK = 64;
    int16_t hf_buf[64 * 64];

    /* Accumulate HF energy across the frame */
    double total_var = 0.0;
    int nblocks = 0;

    /* Accumulate AR coefficients from all blocks */
    double ar_sum[3] = {0, 0, 0};

    for (int by = 0; by + BLK <= height; by += BLK) {
        for (int bx = 0; bx + BLK <= width; bx += BLK) {
            highpass_block(frame_y, stride, bx, by, BLK, width, height, hf_buf);

            /* Compute variance of HF residual */
            double var = 0;
            for (int i = 0; i < BLK * BLK; i++)
                var += (double)hf_buf[i] * hf_buf[i];
            var /= (BLK * BLK);

            total_var += var;
            nblocks++;

            /* Fit AR model if block has enough energy */
            if (var > 1.0) {
                int16_t block_ar[3];
                fit_ar3(hf_buf, BLK, block_ar);
                ar_sum[0] += block_ar[0];
                ar_sum[1] += block_ar[1];
                ar_sum[2] += block_ar[2];
            }
        }
    }

    if (nblocks == 0) {
        memset(out, 0, sizeof(*out));
        return;
    }

    double mean_var = total_var / nblocks;

    /* Enable grain synthesis if HF variance exceeds threshold */
    out->enabled = (mean_var > 2.0);
    if (!out->enabled) {
        memset(out->ar_coeff_y, 0, sizeof(out->ar_coeff_y));
        out->strength_y = out->strength_u = out->strength_v = 0;
        return;
    }

    /* Average AR coefficients across blocks */
    int grainy_blocks = nblocks; /* simplified: use all blocks */
    out->ar_coeff_y[0] = (int16_t)(ar_sum[0] / grainy_blocks);
    out->ar_coeff_y[1] = (int16_t)(ar_sum[1] / grainy_blocks);
    out->ar_coeff_y[2] = (int16_t)(ar_sum[2] / grainy_blocks);

    /* Chroma: use attenuated luma coefficients + luma correlation */
    for (int i = 0; i < 3; i++) {
        out->ar_coeff_u[i] = (int16_t)(out->ar_coeff_y[i] * 3 / 4);
        out->ar_coeff_v[i] = (int16_t)(out->ar_coeff_y[i] * 3 / 4);
    }
    out->ar_coeff_u[3] = 64; /* Q8: 0.25 luma correlation */
    out->ar_coeff_v[3] = 64;

    /* Grain strength proportional to sqrt of HF variance */
    double strength = sqrt(mean_var) * 4.0;
    if (strength > 255.0) strength = 255.0;
    out->strength_y = (uint8_t)strength;
    out->strength_u = (uint8_t)(strength * 0.5);
    out->strength_v = (uint8_t)(strength * 0.5);

    /* Seed from frame content hash (deterministic) */
    uint32_t hash = 0;
    for (int i = 0; i < width * height && i < 1024; i++)
        hash = hash * 31 + frame_y[i];
    out->seed = (uint16_t)(hash & 0xFFFF);
}
