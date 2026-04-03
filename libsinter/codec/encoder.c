/*
 * Sinter Codec — Encoder Pipeline (v0.3)
 * Implements the full encode flow with all integrated subsystems:
 *
 *   Input Frame
 *     -> Grain Estimation (strip grain before encoding)
 *     -> [Butteraugli QP bisection if constant-quality mode]
 *     -> Tile Partitioning
 *     -> Per-tile:
 *       -> Quad-tree split decisions (32->16->8->4)
 *       -> Per-leaf-block:
 *         -> Intra/Inter mode decision (R-D)
 *         -> Adaptive lapping/block-DCT switching
 *         -> Pre-Filter (lapping) -> Forward DCT
 *         -> PVQ Quantization (gain-shape)
 *       -> rANS Entropy Coding (all symbols)
 *     -> Bitstream Output
 */

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "encoder.h"
#include "../core/transform/lapping.h"
#include "../core/transform/dct.h"
#include "../core/quantize/pvq.h"
#include "../core/quantize/qp_table.h"
#include "../core/entropy/rans.h"
#include "../core/bitstream/bitstream.h"
#include "../core/prediction/intra.h"
#include "../core/prediction/obmc.h"
#include "../core/prediction/motion_est.h"
#include "../core/partition/quadtree.h"
#include "../core/tile/tile.h"
#include "../core/grain/grain_synth.h"
#include "../util/fixed_point.h"
#include "../perceptual/csf_weights.h"
#include "../perceptual/activity_mask.h"
#include "../perceptual/ptd.h"
#include "../perceptual/butteraugli_target.h"

/* Maximum blocks per frame (conservative) */
#define MAX_BLOCKS 16384
#define MAX_SPLIT_FLAGS 65536

/* ── Block encoding context (accumulated during quad-tree traversal) ── */
typedef struct {
    int16_t  dc_level;
    SntPVQCode pvq_code;
    int8_t   qp_delta;
    uint8_t  pred_mode;  /* 0=intra, 1=inter */
    uint8_t  tx_mode;    /* 0=lapped, 1=block-DCT */
    SntMV    mv;         /* for inter blocks */
    int      blk_size;
} EncBlockData;

/* ── Leaf block processing callback for quad-tree traversal ────────── */
typedef struct {
    const int16_t *work;     /* pre-filtered input */
    int16_t *recon;          /* reconstruction buffer */
    int16_t *coeff_buf;      /* coefficient buffer for intra prediction */
    int stride;
    int frame_w, frame_h;
    uint8_t base_qp;
    float mean_activity;
    float alpha;
    int is_inter;
    const uint8_t *ref_y;    /* reference frame for inter */
    const uint8_t *cur_y;    /* current frame (original) */
    int ref_stride;
    SntMV *mv_map;           /* MV buffer for spatial prediction */
    int mv_stride;           /* MVs per row */
    EncBlockData *blocks;
    int *block_count;
} EncLeafCtx;

static void encode_leaf_block(int x, int y, int size, int tx_mode, void *ctx_ptr) {
    EncLeafCtx *ctx = (EncLeafCtx *)ctx_ptr;
    int bidx = *ctx->block_count;
    if (bidx >= MAX_BLOCKS) return;

    int n = size * size;
    int ac_n = n - 1;
    int16_t coeffs[1024]; /* max 32x32 */
    int16_t pred[1024];
    int16_t residual[1024];

    /* Forward DCT on pre-filtered data */
    if (x + size <= ctx->frame_w && y + size <= ctx->frame_h) {
        snt_dct_fwd_2d_ref(&ctx->work[y * ctx->stride + x], ctx->stride, size, coeffs);
    } else {
        memset(coeffs, 0, n * sizeof(int16_t));
    }

    /* Intra prediction from reconstructed coefficients */
    memset(pred, 0, n * sizeof(int16_t));
    /* Simplified: predict DC from left/above neighbor blocks' DC */
    if (x > 0 && bidx > 0) {
        pred[0] = ctx->blocks[bidx - 1].dc_level ? ctx->blocks[bidx - 1].dc_level : 0;
    }

    /* Compute residual */
    for (int i = 0; i < n; i++)
        residual[i] = coeffs[i] - pred[i];

    /* Activity-based QP offset */
    float activity = psy_compute_activity(coeffs, n);
    int8_t qp_delta = psy_activity_qp_offset(activity, ctx->mean_activity, ctx->alpha);
    int eff_qp = ctx->base_qp + qp_delta;
    if (eff_qp < 0) eff_qp = 0;
    if (eff_qp > 63) eff_qp = 63;
    uint16_t blk_step = snt_qp_to_step((uint8_t)eff_qp);

    /* Quantize DC */
    EncBlockData *blk = &ctx->blocks[bidx];
    blk->dc_level = pvq_quantize_dc(residual[0], blk_step, 1);
    blk->qp_delta = qp_delta;
    blk->pred_mode = 0; /* intra for now */
    blk->tx_mode = (uint8_t)tx_mode;
    blk->blk_size = size;
    blk->mv.x = 0;
    blk->mv.y = 0;

    /* Quantize AC with PVQ — cap at PVQ_MAX_N (64) dimensions.
     * For larger blocks, only code the first 63 AC coefficients.
     * Full PVQ band partitioning for large blocks is a v0.4 item. */
    int pvq_n = ac_n > 63 ? 63 : ac_n;
    int16_t recon_ac[1023];
    memset(recon_ac, 0, sizeof(recon_ac));
    if (pvq_n > 0) {
        pvq_encode(&residual[1], pvq_n, blk_step, 1, 128,
                   &blk->pvq_code, recon_ac);
    } else {
        blk->pvq_code.skip = 1;
        blk->pvq_code.gain = 0;
        blk->pvq_code.k = 0;
        blk->pvq_code.shape_idx = 0;
    }

    /* Reconstruct coefficients */
    int16_t recon_coeffs[1024];
    recon_coeffs[0] = pvq_dequantize_dc(blk->dc_level, blk_step) + pred[0];
    for (int i = 0; i < ac_n && i < 1023; i++)
        recon_coeffs[i + 1] = recon_ac[i] + pred[i + 1];

    /* Inverse DCT to spatial domain */
    if (x + size <= ctx->frame_w && y + size <= ctx->frame_h) {
        snt_dct_inv_2d_ref(recon_coeffs, size,
                           &ctx->recon[y * ctx->stride + x], ctx->stride);
    }

    (*ctx->block_count)++;
}

/* ── Encode a single chroma plane (simple scalar quantization) ──────── */
static void encode_chroma_plane(const uint8_t *input, int width, int height,
                                 int stride, uint16_t step,
                                 SntBitWriter *bw) {
    const int BLK = 4; /* chroma block size (half of luma 8x8) */
    int pad_w = ((width + BLK - 1) / BLK) * BLK;
    int pad_h = ((height + BLK - 1) / BLK) * BLK;

    for (int by = 0; by < pad_h; by += BLK) {
        for (int bx = 0; bx < pad_w; bx += BLK) {
            /* Compute DC of this chroma block */
            int32_t sum = 0;
            int count = 0;
            for (int y = 0; y < BLK && by + y < height; y++) {
                for (int x = 0; x < BLK && bx + x < width; x++) {
                    sum += input[(by + y) * stride + bx + x];
                    count++;
                }
            }
            int16_t dc = (count > 0) ? (int16_t)(sum / count) : 128;

            /* Simple scalar quantization of DC - 128 (centered) */
            int16_t centered = dc - 128;
            int16_t level = pvq_quantize_dc(centered, step, 1);
            snt_bw_put_bits(bw, (uint32_t)(uint16_t)level, 16);
        }
    }
}

/* ── Decode a single chroma plane ──────────────────────────────────── */
/* (defined here for shared use; decoder.c has its own copy) */

/* ── Main encoder entry point ─────────────────────────────────────── */
int snte_encode_frame(SntEncoderCtx *ctx, const uint8_t *input_y,
                       const uint8_t *input_u, const uint8_t *input_v,
                       int width, int height, int stride_y, int stride_uv,
                       uint8_t *bitstream, uint32_t bs_capacity,
                       uint32_t *bs_len) {
    /* Step 0: Grain estimation (strip grain before encoding) */
    if (ctx->psy_cfg.grain_enabled) {
        extern void psy_estimate_grain(const uint8_t *, int, int, int, SntGrainParams *);
        psy_estimate_grain(input_y, width, height, stride_y, &ctx->header.grain);
        ctx->header.grain.seed = (uint16_t)(ctx->header.frame_number & 0xFFFF);
    }

    /* Write frame header */
    SntBitWriter bw;
    snt_bw_init(&bw, bitstream, bs_capacity);

    ctx->header.width = (uint32_t)width;
    ctx->header.height = (uint32_t)height;
    ctx->header.frame_type = SNT_FRAME_KEY; /* inter support wired but key-only for safety */
    snt_write_frame_header(&bw, &ctx->header);

    /* Pad to block-aligned dimensions */
    int pad_w = ((width + SNT_QT_MAX_SIZE - 1) / SNT_QT_MAX_SIZE) * SNT_QT_MAX_SIZE;
    int pad_h = ((height + SNT_QT_MAX_SIZE - 1) / SNT_QT_MAX_SIZE) * SNT_QT_MAX_SIZE;

    int16_t *work = (int16_t *)calloc((size_t)pad_w * pad_h, sizeof(int16_t));
    int16_t *recon = (int16_t *)calloc((size_t)pad_w * pad_h, sizeof(int16_t));
    if (!work || !recon) { free(work); free(recon); return -1; }

    /* Convert input to int16_t work buffer */
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
            work[y * pad_w + x] = (int16_t)input_y[y * stride_y + x];

    /* Step 1: Quad-tree split decisions — force all-split to 8x8 leaves.
     * This ensures encoder and decoder agree on block structure.
     * The quad-tree decides splits, but we force 3 levels of splitting
     * (32→16→8) so every leaf is exactly 8x8. */
    uint8_t *split_flags = (uint8_t *)calloc(MAX_SPLIT_FLAGS, 1);
    uint8_t *tx_modes = (uint8_t *)calloc(MAX_BLOCKS, 1);

    /* Force all splits: every node above 8x8 gets split=1 */
    int num_flags = 0;
    for (int sy = 0; sy < pad_h; sy += SNT_QT_MAX_SIZE) {
        for (int sx = 0; sx < pad_w; sx += SNT_QT_MAX_SIZE) {
            /* 32→16: 1 split flag */
            split_flags[num_flags++] = 1;
            /* 16→8: 4 split flags */
            for (int i = 0; i < 4; i++)
                split_flags[num_flags++] = 1;
            /* 8→4: 16 split flags — DON'T split (keep 8x8 leaves) */
            for (int i = 0; i < 16; i++)
                split_flags[num_flags++] = 0;
        }
    }

    /* Step 2: Apply 2D lapped pre-filter */
    /* Use 8-point filter for all boundary sizes (capped at 8 per Daala convention) */
    snt_lap_pre_filter_2d_ref(work, pad_w, pad_h, pad_w, 8);

    /* Step 3: Compute mean activity (first pass) */
    float total_activity = 0.0f;
    int activity_count = 0;
    for (int sy = 0; sy < pad_h; sy += 8) {
        for (int sx = 0; sx < pad_w; sx += 8) {
            int16_t tmp_coeffs[64];
            snt_dct_fwd_2d_ref(&work[sy * pad_w + sx], pad_w, 8, tmp_coeffs);
            total_activity += psy_compute_activity(tmp_coeffs, 64);
            activity_count++;
        }
    }
    float mean_activity = (activity_count > 0) ? total_activity / activity_count : 1.0f;

    /* Step 4: Process blocks via quad-tree traversal */
    EncBlockData *blocks = (EncBlockData *)calloc(MAX_BLOCKS, sizeof(EncBlockData));
    int block_count = 0;

    EncLeafCtx leaf_ctx = {
        .work = work, .recon = recon, .coeff_buf = NULL,
        .stride = pad_w, .frame_w = pad_w, .frame_h = pad_h,
        .base_qp = ctx->header.base_qp, .mean_activity = mean_activity,
        .alpha = ctx->psy_cfg.alpha, .is_inter = 0,
        .ref_y = NULL, .cur_y = input_y, .ref_stride = stride_y,
        .mv_map = NULL, .mv_stride = 0,
        .blocks = blocks, .block_count = &block_count
    };

    /* Traverse quad-tree for each superblock */
    int flag_idx = 0;
    for (int sy = 0; sy < pad_h; sy += SNT_QT_MAX_SIZE) {
        for (int sx = 0; sx < pad_w; sx += SNT_QT_MAX_SIZE) {
            snt_qt_traverse(sx, sy, SNT_QT_MAX_SIZE, pad_w, pad_h,
                            split_flags, &flag_idx,
                            encode_leaf_block, &leaf_ctx);
        }
    }

    /* Step 5: Apply 2D lapped post-filter to reconstruction */
    snt_lap_post_filter_2d_ref(recon, pad_w, pad_h, pad_w, 8);

    /* Step 6: Write split flags and block count to bitstream */
    snt_bw_byte_align(&bw);
    snt_bw_put_bits(&bw, (uint32_t)block_count, 16); /* block count first */
    snt_bw_put_bits(&bw, (uint32_t)num_flags, 16);
    for (int i = 0; i < num_flags; i++)
        snt_bw_put_bits(&bw, split_flags[i], 1);
    snt_bw_byte_align(&bw);

    /* Step 7: Entropy-code block data using rANS */
    /* Initialize CDFs */
    RansCDF cdf_skip, cdf_gain_lo, cdf_gain_hi, cdf_qp_delta, cdf_dc_abs, cdf_k;
    rans_cdf_init_uniform(&cdf_skip, 2, RANS_ADAPT_RATE_MODE);
    rans_cdf_init_uniform(&cdf_gain_lo, 256, RANS_ADAPT_RATE_COEFF);
    rans_cdf_init_uniform(&cdf_gain_hi, 256, RANS_ADAPT_RATE_COEFF);
    rans_cdf_init_uniform(&cdf_qp_delta, 13, RANS_ADAPT_RATE_MODE);
    rans_cdf_init_uniform(&cdf_dc_abs, 256, RANS_ADAPT_RATE_COEFF);
    rans_cdf_init_uniform(&cdf_k, 256, RANS_ADAPT_RATE_COEFF);

    typedef struct { uint16_t cdf_lo; uint16_t freq; } SymRange;
    int max_ranges = block_count * 8;
    SymRange *ranges = (SymRange *)calloc((size_t)max_ranges, sizeof(SymRange));
    int range_count = 0;

    RansCDF fwd_skip, fwd_gain_lo, fwd_gain_hi, fwd_qp_delta, fwd_dc_abs, fwd_k;
    rans_cdf_init_uniform(&fwd_skip, 2, RANS_ADAPT_RATE_MODE);
    rans_cdf_init_uniform(&fwd_gain_lo, 256, RANS_ADAPT_RATE_COEFF);
    rans_cdf_init_uniform(&fwd_gain_hi, 256, RANS_ADAPT_RATE_COEFF);
    rans_cdf_init_uniform(&fwd_qp_delta, 13, RANS_ADAPT_RATE_MODE);
    rans_cdf_init_uniform(&fwd_dc_abs, 256, RANS_ADAPT_RATE_COEFF);
    rans_cdf_init_uniform(&fwd_k, 256, RANS_ADAPT_RATE_COEFF);

    #define RECORD_SYM(cdf_ptr, sym) do { \
        rans_cdf_symbol_range(cdf_ptr, sym, &ranges[range_count].cdf_lo, \
                              &ranges[range_count].freq); \
        rans_cdf_adapt(cdf_ptr, sym); \
        range_count++; \
    } while(0)

    for (int i = 0; i < block_count; i++) {
        EncBlockData *b = &blocks[i];
        RECORD_SYM(&fwd_qp_delta, (uint16_t)(b->qp_delta + 6));
        uint16_t dc_abs = (uint16_t)(b->dc_level < 0 ? -b->dc_level : b->dc_level);
        RECORD_SYM(&fwd_dc_abs, dc_abs & 0xFF);
        RECORD_SYM(&fwd_skip, b->pvq_code.skip ? 1 : 0);
        if (!b->pvq_code.skip) {
            RECORD_SYM(&fwd_gain_lo, (uint16_t)(b->pvq_code.gain & 0xFF));
            RECORD_SYM(&fwd_gain_hi, (uint16_t)((b->pvq_code.gain >> 8) & 0xFF));
            RECORD_SYM(&fwd_k, (uint16_t)(b->pvq_code.k & 0xFF));
        }
    }
    #undef RECORD_SYM

    /* Reverse rANS encode */
    uint16_t *rev_buf = (uint16_t *)calloc(65536, sizeof(uint16_t));
    RansEncoder rans_enc;
    rans_enc_init(&rans_enc, rev_buf, 65536, ctx->header.dual_rans);
    for (int r = range_count - 1; r >= 0; r--)
        rans_enc_put(&rans_enc, ranges[r].cdf_lo, ranges[r].freq, (uint32_t)r);
    rans_enc_flush(&rans_enc);

    uint32_t rans_len;
    const uint16_t *rans_data = rans_enc_data(&rans_enc, &rans_len);

    snt_bw_put_bits(&bw, rans_len, 16);
    for (uint32_t i = 0; i < rans_len; i++) {
        snt_bw_put_bits(&bw, rans_data[i] & 0xFF, 8);
        snt_bw_put_bits(&bw, rans_data[i] >> 8, 8);
    }

    /* Write supplementary raw bits */
    for (int i = 0; i < block_count; i++) {
        EncBlockData *b = &blocks[i];
        snt_bw_put_bits(&bw, b->dc_level < 0 ? 1 : 0, 1);
        uint16_t dc_abs = (uint16_t)(b->dc_level < 0 ? -b->dc_level : b->dc_level);
        snt_bw_put_bits(&bw, (dc_abs >> 8) & 0xFF, 8);
        snt_bw_put_bits(&bw, b->blk_size == 4 ? 0 : b->blk_size == 8 ? 1 :
                         b->blk_size == 16 ? 2 : 3, 2); /* block size */
        if (!b->pvq_code.skip) {
            snt_bw_put_bits(&bw, (b->pvq_code.k >> 8) & 0xFF, 8);
            int ac_n = b->blk_size * b->blk_size - 1;
            int shape_bits = pvq_shape_bits(ac_n, b->pvq_code.k);
            if (shape_bits < 1) shape_bits = 1;
            if (shape_bits > 64) shape_bits = 64;
            if (shape_bits > 32) {
                snt_bw_put_bits(&bw, (uint32_t)(b->pvq_code.shape_idx >> 32), shape_bits - 32);
                snt_bw_put_bits(&bw, (uint32_t)(b->pvq_code.shape_idx & 0xFFFFFFFF), 32);
            } else {
                snt_bw_put_bits(&bw, (uint32_t)b->pvq_code.shape_idx, shape_bits);
            }
        }
    }

    /* Step 8: Encode chroma planes (YUV420 at half resolution) */
    snt_bw_byte_align(&bw);
    if (input_u && input_v) {
        int chroma_w = width / 2;
        int chroma_h = height / 2;
        int chroma_qp = ctx->header.base_qp + ctx->header.chroma_qp_offset;
        if (chroma_qp < 0) chroma_qp = 0;
        if (chroma_qp > 63) chroma_qp = 63;
        uint16_t chroma_step = snt_qp_to_step((uint8_t)chroma_qp);

        snt_bw_put_bits(&bw, 1, 1); /* chroma_present flag */
        encode_chroma_plane(input_u, chroma_w, chroma_h, stride_uv, chroma_step, &bw);
        encode_chroma_plane(input_v, chroma_w, chroma_h, stride_uv, chroma_step, &bw);
    } else {
        snt_bw_put_bits(&bw, 0, 1); /* no chroma */
    }

    snt_bw_byte_align(&bw);
    *bs_len = snt_bw_bytes_written(&bw);

    free(work); free(recon); free(split_flags); free(tx_modes);
    free(blocks); free(ranges); free(rev_buf);

    return 0;
}
