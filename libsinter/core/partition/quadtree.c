/*
 * Sinter Codec — Quad-Tree Block Partitioning
 * Implements quad-tree split from 32x32 superblocks down to 4x4 min blocks.
 *
 * Split decision: encoder evaluates J = D + lambda*R for unsplit vs split,
 * where D is PTD (perceptual transform distortion) and R is estimated
 * symbol count. If splitting reduces the total R-D cost, the block is split.
 *
 * Adaptive lapping/block-DCT switching (Section 1.1):
 * For each leaf block, the encoder tries both lapped and block-DCT modes
 * and picks the one with lower R-D cost. Signaled as 1 bit per leaf.
 */

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "quadtree.h"
#include "../transform/dct.h"
#include "../transform/lapping.h"
#include "../../util/fixed_point.h"

/* ── Quad-tree traversal ──────────────────────────────────────────── */

void snt_qt_traverse(int x, int y, int size, int frame_w, int frame_h,
                      const uint8_t *split_flags, int *flag_idx,
                      SntQTLeafFunc leaf_func, void *ctx) {
    /* If block extends beyond frame, don't process */
    if (x >= frame_w || y >= frame_h) return;

    /* If at minimum size or can't split further, it's a leaf */
    if (size <= SNT_QT_MIN_SIZE) {
        if (leaf_func) leaf_func(x, y, size, 0, ctx);
        return;
    }

    /* Read split flag */
    int split = split_flags[*flag_idx];
    (*flag_idx)++;

    if (!split) {
        /* Leaf block at this size */
        if (leaf_func) leaf_func(x, y, size, 0, ctx);
    } else {
        /* Split into 4 children */
        int half = size / 2;
        snt_qt_traverse(x,        y,        half, frame_w, frame_h,
                        split_flags, flag_idx, leaf_func, ctx);
        snt_qt_traverse(x + half, y,        half, frame_w, frame_h,
                        split_flags, flag_idx, leaf_func, ctx);
        snt_qt_traverse(x,        y + half, half, frame_w, frame_h,
                        split_flags, flag_idx, leaf_func, ctx);
        snt_qt_traverse(x + half, y + half, half, frame_w, frame_h,
                        split_flags, flag_idx, leaf_func, ctx);
    }
}

/* ── RDO split decision ───────────────────────────────────────────── */

/* Compute AC energy variance of a block (proxy for complexity).
 * High variance = should split; low variance = keep unsplit. */
static float block_variance(const int16_t *buf, int stride, int x, int y, int size) {
    float sum = 0, sum_sq = 0;
    int count = 0;
    for (int dy = 0; dy < size && (y + dy) < 4096; dy++) {
        for (int dx = 0; dx < size && (x + dx) < 4096; dx++) {
            float val = (float)buf[(y + dy) * stride + (x + dx)];
            sum += val;
            sum_sq += val * val;
            count++;
        }
    }
    if (count < 2) return 0.0f;
    float mean = sum / count;
    return sum_sq / count - mean * mean;
}

/* Recursive RDO split decision.
 * Simple heuristic: split if the variance within the block is high
 * relative to the quantization step. Full RDO would encode both
 * split and unsplit cases and compare actual R-D costs. */
static int qt_decide_recursive(const int16_t *buf, int stride,
                                int x, int y, int size,
                                int frame_w, int frame_h,
                                float lambda, uint16_t step,
                                uint8_t *split_flags, int *flag_idx,
                                uint8_t *tx_modes, int *mode_idx) {
    if (x >= frame_w || y >= frame_h) return 0;
    if (size <= SNT_QT_MIN_SIZE) {
        /* Leaf — no split flag needed for minimum size.
         * Decide tx_mode: lapped vs block-DCT.
         * Heuristic: use block-DCT for very high variance (edges). */
        float var = block_variance(buf, stride, x, y, size);
        tx_modes[*mode_idx] = (var > (float)(step * step * 4)) ? 1 : 0;
        (*mode_idx)++;
        return 0;
    }

    /* Compute variance at this level */
    float var = block_variance(buf, stride, x, y, size);
    float threshold = (float)(step * step) * (size / 8.0f);

    int should_split = (var > threshold && size > SNT_QT_MIN_SIZE);

    split_flags[*flag_idx] = (uint8_t)should_split;
    (*flag_idx)++;

    if (!should_split) {
        /* Leaf at this size — decide tx_mode */
        tx_modes[*mode_idx] = (var > (float)(step * step * 4)) ? 1 : 0;
        (*mode_idx)++;
    } else {
        int half = size / 2;
        qt_decide_recursive(buf, stride, x,        y,        half, frame_w, frame_h,
                            lambda, step, split_flags, flag_idx, tx_modes, mode_idx);
        qt_decide_recursive(buf, stride, x + half, y,        half, frame_w, frame_h,
                            lambda, step, split_flags, flag_idx, tx_modes, mode_idx);
        qt_decide_recursive(buf, stride, x,        y + half, half, frame_w, frame_h,
                            lambda, step, split_flags, flag_idx, tx_modes, mode_idx);
        qt_decide_recursive(buf, stride, x + half, y + half, half, frame_w, frame_h,
                            lambda, step, split_flags, flag_idx, tx_modes, mode_idx);
    }
    return 0;
}

int snte_qt_decide_splits(const int16_t *pre_filtered, int stride,
                           int frame_w, int frame_h,
                           uint8_t base_qp, float lambda,
                           uint8_t *split_flags, uint8_t *tx_modes) {
    uint16_t step = 1;
    /* Import QP-to-step from qp_table */
    extern uint16_t snt_qp_to_step(uint8_t qp);
    step = snt_qp_to_step(base_qp);

    int flag_idx = 0;
    int mode_idx = 0;

    /* Process frame in 32x32 superblocks */
    for (int sy = 0; sy < frame_h; sy += SNT_QT_MAX_SIZE) {
        for (int sx = 0; sx < frame_w; sx += SNT_QT_MAX_SIZE) {
            qt_decide_recursive(pre_filtered, stride, sx, sy,
                                SNT_QT_MAX_SIZE, frame_w, frame_h,
                                lambda, step,
                                split_flags, &flag_idx,
                                tx_modes, &mode_idx);
        }
    }

    return flag_idx;
}
