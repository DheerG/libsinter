/*
 * Sinter Codec — Butteraugli Constant-Quality Encoding Mode
 * Implements perceptual optimization: outer-loop QP bisection
 *
 * Encodes a frame, measures perceptual distance using a simplified
 * Butteraugli-like metric, adjusts QP, and re-encodes until the
 * target quality is met. Maximum 3 iterations (diminishing returns).
 *
 * The perceptual distance metric is a simplified multi-scale SSIM
 * variant computed in a perceptually uniform space. Full Butteraugli
 * would use the Jyrki Alakuijala formulation; this is a practical
 * approximation for the reference encoder.
 */

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "butteraugli_target.h"

/* Simplified perceptual distance: weighted combination of
 * mean absolute error and edge-preserving structural difference.
 * Returns a score where 0 = identical, 1.0 ≈ just noticeable,
 * higher = more visible distortion. */
float psy_perceptual_distance(const uint8_t *orig, const uint8_t *recon,
                               int width, int height, int stride) {
    if (width < 8 || height < 8) return 0.0f;

    /* Compute block-wise differences at 8x8 scale */
    float total_dist = 0.0f;
    int nblocks = 0;

    for (int by = 0; by + 8 <= height; by += 8) {
        for (int bx = 0; bx + 8 <= width; bx += 8) {
            float mae = 0.0f;  /* mean absolute error */
            float orig_var = 0.0f, recon_var = 0.0f;
            float orig_mean = 0.0f, recon_mean = 0.0f;

            for (int y = 0; y < 8; y++) {
                for (int x = 0; x < 8; x++) {
                    float o = (float)orig[(by + y) * stride + bx + x];
                    float r = (float)recon[(by + y) * stride + bx + x];
                    mae += fabsf(o - r);
                    orig_mean += o;
                    recon_mean += r;
                }
            }
            mae /= 64.0f;
            orig_mean /= 64.0f;
            recon_mean /= 64.0f;

            for (int y = 0; y < 8; y++) {
                for (int x = 0; x < 8; x++) {
                    float o = (float)orig[(by + y) * stride + bx + x] - orig_mean;
                    float r = (float)recon[(by + y) * stride + bx + x] - recon_mean;
                    orig_var += o * o;
                    recon_var += r * r;
                }
            }
            orig_var /= 64.0f;
            recon_var /= 64.0f;

            /* Structural difference: penalize variance loss (blur) more
             * than variance gain (noise). This is the key perceptual insight. */
            float var_ratio = (recon_var + 1.0f) / (orig_var + 1.0f);
            float struct_penalty = (var_ratio < 1.0f) ? (1.0f - var_ratio) * 2.0f : 0.0f;

            /* Combined distance: MAE + structural penalty */
            float block_dist = mae / 255.0f + struct_penalty;
            total_dist += block_dist;
            nblocks++;
        }
    }

    return (nblocks > 0) ? total_dist / nblocks : 0.0f;
}

/* QP bisection for constant-quality encoding.
 * Returns the QP that achieves the target Butteraugli score.
 * encode_func: callback that encodes at a given QP and returns the
 *              reconstructed frame and bitstream length. */
uint8_t psy_butteraugli_bisect(const uint8_t *orig, int width, int height,
                                int stride, float target,
                                PsyEncodeFunc encode_func, void *enc_ctx,
                                uint8_t *recon_buf, uint8_t *bs_buf,
                                uint32_t bs_cap) {
    uint8_t lo_qp = 0, hi_qp = 63;
    uint8_t best_qp = 32;
    int max_iters = 3; /* diminishing returns beyond 3 */

    for (int iter = 0; iter < max_iters; iter++) {
        uint8_t mid_qp = (lo_qp + hi_qp) / 2;
        uint32_t bs_len;

        /* Encode at mid_qp */
        encode_func(enc_ctx, orig, width, height, stride,
                    mid_qp, recon_buf, bs_buf, bs_cap, &bs_len);

        /* Measure perceptual distance */
        float dist = psy_perceptual_distance(orig, recon_buf, width, height, stride);

        if (dist > target) {
            /* Quality too low — decrease QP */
            hi_qp = mid_qp;
        } else {
            /* Quality sufficient — can increase QP for smaller bitstream */
            lo_qp = mid_qp;
        }
        best_qp = (lo_qp + hi_qp) / 2;
    }

    return best_qp;
}
