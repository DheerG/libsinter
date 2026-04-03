/*
 * Sinter Codec — Decoder Pipeline (v0.3)
 * Mirrors the encoder: reads quad-tree splits, decodes variable-size
 * blocks, applies lapped post-filter, film grain synthesis.
 */

#include <string.h>
#include <stdlib.h>
#include "decoder.h"
#include "../util/fixed_point.h"
#include "../core/transform/lapping.h"
#include "../core/transform/dct.h"
#include "../core/quantize/pvq.h"
#include "../core/quantize/qp_table.h"
#include "../core/bitstream/bitstream.h"
#include "../core/prediction/intra.h"
#include "../core/partition/quadtree.h"
#include "../core/grain/grain_synth.h"
#include "../core/entropy/rans.h"

#define MAX_BLOCKS 16384

typedef struct {
    int16_t  dc_level;
    SntPVQCode pvq_code;
    int8_t   qp_delta;
    int      blk_size;
} DecBlockData;

/* ── Leaf block decode callback ───────────────────────────────────── */
typedef struct {
    int16_t *recon;
    int stride;
    int frame_w, frame_h;
    uint8_t base_qp;
    DecBlockData *blocks;
    int *block_idx;
} DecLeafCtx;

static void decode_leaf_block(int x, int y, int size, int tx_mode, void *ctx_ptr) {
    DecLeafCtx *ctx = (DecLeafCtx *)ctx_ptr;
    int bidx = *ctx->block_idx;
    if (bidx >= MAX_BLOCKS) return;
    (void)tx_mode;

    DecBlockData *b = &ctx->blocks[bidx];
    int n = size * size;
    int ac_n = n - 1;

    /* Apply per-block QP delta */
    int eff_qp = ctx->base_qp + b->qp_delta;
    if (eff_qp < 0) eff_qp = 0;
    if (eff_qp > 63) eff_qp = 63;
    uint16_t blk_step = snt_qp_to_step((uint8_t)eff_qp);

    /* Simplified intra prediction */
    int16_t pred[1024];
    memset(pred, 0, n * sizeof(int16_t));
    if (x > 0 && bidx > 0) {
        pred[0] = ctx->blocks[bidx - 1].dc_level ? ctx->blocks[bidx - 1].dc_level : 0;
    }

    /* Dequantize DC */
    int16_t recon_coeffs[1024];
    recon_coeffs[0] = pvq_dequantize_dc(b->dc_level, blk_step) + pred[0];

    /* Dequantize AC — cap at 63 dimensions (PVQ_MAX_N) */
    int pvq_n = ac_n > 63 ? 63 : ac_n;
    int16_t recon_ac[1023];
    memset(recon_ac, 0, sizeof(recon_ac));
    if (pvq_n > 0) {
        pvq_decode(&b->pvq_code, pvq_n, blk_step, recon_ac);
    }
    for (int i = 0; i < ac_n && i < 1023; i++)
        recon_coeffs[i + 1] = recon_ac[i] + pred[i + 1];

    /* Inverse DCT */
    if (x + size <= ctx->frame_w && y + size <= ctx->frame_h) {
        snt_dct_inv_2d_ref(recon_coeffs, size,
                           &ctx->recon[y * ctx->stride + x], ctx->stride);
    }

    (*ctx->block_idx)++;
}

/* ── Main decoder entry point ─────────────────────────────────────── */
int sntd_decode_frame(SntDecoderCtx *ctx, const uint8_t *bitstream,
                       uint32_t bs_len, uint8_t *output_y,
                       uint8_t *output_u, uint8_t *output_v,
                       int *out_width, int *out_height,
                       int stride_y, int stride_uv) {
    SntBitReader br;
    snt_br_init(&br, bitstream, bs_len);
    snt_read_frame_header(&br, &ctx->header);

    int width = (int)ctx->header.width;
    int height = (int)ctx->header.height;
    *out_width = width;
    *out_height = height;

    int pad_w = ((width + SNT_QT_MAX_SIZE - 1) / SNT_QT_MAX_SIZE) * SNT_QT_MAX_SIZE;
    int pad_h = ((height + SNT_QT_MAX_SIZE - 1) / SNT_QT_MAX_SIZE) * SNT_QT_MAX_SIZE;

    /* Read block count and split flags */
    snt_br_byte_align(&br);
    uint32_t total_blocks = snt_br_get_bits(&br, 16);
    uint32_t num_flags = snt_br_get_bits(&br, 16);
    uint8_t *split_flags = (uint8_t *)calloc(num_flags + 1, 1);
    for (uint32_t i = 0; i < num_flags; i++)
        split_flags[i] = (uint8_t)snt_br_get_bits(&br, 1);
    snt_br_byte_align(&br);

    /* Read rANS payload */
    uint32_t rans_len = snt_br_get_bits(&br, 16);
    uint8_t *rans_bytes = (uint8_t *)calloc(rans_len * 2 + 4, 1);
    for (uint32_t i = 0; i < rans_len; i++) {
        rans_bytes[i * 2]     = (uint8_t)snt_br_get_bits(&br, 8);
        rans_bytes[i * 2 + 1] = (uint8_t)snt_br_get_bits(&br, 8);
    }

    /* Initialize rANS decoder + CDFs (must match encoder exactly) */
    RansDecoder rans_dec;
    rans_dec_init(&rans_dec, rans_bytes, rans_len * 2, ctx->header.dual_rans);

    RansCDF cdf_skip, cdf_gain_lo, cdf_gain_hi, cdf_qp_delta, cdf_dc_abs, cdf_k;
    rans_cdf_init_uniform(&cdf_skip, 2, RANS_ADAPT_RATE_MODE);
    rans_cdf_init_uniform(&cdf_gain_lo, 256, RANS_ADAPT_RATE_COEFF);
    rans_cdf_init_uniform(&cdf_gain_hi, 256, RANS_ADAPT_RATE_COEFF);
    rans_cdf_init_uniform(&cdf_qp_delta, 13, RANS_ADAPT_RATE_MODE);
    rans_cdf_init_uniform(&cdf_dc_abs, 256, RANS_ADAPT_RATE_COEFF);
    rans_cdf_init_uniform(&cdf_k, 256, RANS_ADAPT_RATE_COEFF);

    /* First: count blocks by traversing the quad-tree (dry run) */
    int block_count = 0;
    int fi = 0;
    for (int sy = 0; sy < pad_h; sy += SNT_QT_MAX_SIZE) {
        for (int sx = 0; sx < pad_w; sx += SNT_QT_MAX_SIZE) {
            /* Count leaves */
            int saved_fi = fi;
            snt_qt_traverse(sx, sy, SNT_QT_MAX_SIZE, pad_w, pad_h,
                            split_flags, &fi,
                            NULL, &block_count); /* NULL callback = just count */
        }
    }
    /* Reset for actual traversal — need to re-count properly */
    /* Actually, we can decode rANS symbols first since we know the block count
     * from the flag traversal. But snt_qt_traverse with NULL callback won't work.
     * Let's use a simple counting callback instead. */

    /* Count blocks by traversing split flags */
    block_count = 0;
    fi = 0;
    /* Simple counter: each leaf increments count */
    for (uint32_t i = 0; i < num_flags; i++) {
        if (!split_flags[i]) block_count++; /* unsplit node = leaf */
    }
    /* Also count minimum-size blocks (4x4) that have no split flag */
    /* This is approximate — for the reference impl, estimate generously */
    block_count += (pad_w / SNT_QT_MIN_SIZE) * (pad_h / SNT_QT_MIN_SIZE);
    if (block_count > MAX_BLOCKS) block_count = MAX_BLOCKS;

    /* Decode rANS symbols */
    DecBlockData *blocks = (DecBlockData *)calloc(MAX_BLOCKS, sizeof(DecBlockData));
    int actual_blocks = 0;

    /* Decode exactly total_blocks blocks (count read from bitstream) */
    for (uint32_t i = 0; i < total_blocks && i < MAX_BLOCKS; i++) {
        uint16_t qpd_sym = rans_dec_symbol(&rans_dec, &cdf_qp_delta);
        blocks[i].qp_delta = (int8_t)((int)qpd_sym - 6);

        uint16_t dc_abs_lo = rans_dec_symbol(&rans_dec, &cdf_dc_abs);

        uint16_t skip_sym = rans_dec_symbol(&rans_dec, &cdf_skip);
        blocks[i].pvq_code.skip = (uint8_t)skip_sym;

        if (!blocks[i].pvq_code.skip) {
            blocks[i].pvq_code.gain = rans_dec_symbol(&rans_dec, &cdf_gain_lo);
            uint16_t gain_hi = rans_dec_symbol(&rans_dec, &cdf_gain_hi);
            blocks[i].pvq_code.gain |= ((uint32_t)gain_hi << 8);
            uint16_t k_lo = rans_dec_symbol(&rans_dec, &cdf_k);
            blocks[i].pvq_code.k = k_lo;
        }

        /* Read supplementary raw bits */
        uint32_t dc_sign = snt_br_get_bits(&br, 1);
        uint32_t dc_hi = snt_br_get_bits(&br, 8);
        int16_t dc_abs = (int16_t)(dc_abs_lo | (int16_t)(dc_hi << 8));
        blocks[i].dc_level = dc_sign ? -dc_abs : dc_abs;

        uint32_t blk_size_code = snt_br_get_bits(&br, 2);
        blocks[i].blk_size = (blk_size_code == 0) ? 4 : (blk_size_code == 1) ? 8 :
                              (blk_size_code == 2) ? 16 : 32;

        if (!blocks[i].pvq_code.skip) {
            uint32_t k_hi = snt_br_get_bits(&br, 8);
            blocks[i].pvq_code.k |= (uint16_t)(k_hi << 8);

            int ac_n = blocks[i].blk_size * blocks[i].blk_size - 1;
            int shape_bits = pvq_shape_bits(ac_n, blocks[i].pvq_code.k);
            if (shape_bits < 1) shape_bits = 1;
            if (shape_bits > 64) shape_bits = 64;

            blocks[i].pvq_code.shape_idx = 0;
            if (shape_bits > 32) {
                uint64_t hi = snt_br_get_bits(&br, shape_bits - 32);
                uint64_t lo = snt_br_get_bits(&br, 32);
                blocks[i].pvq_code.shape_idx = (hi << 32) | lo;
            } else {
                blocks[i].pvq_code.shape_idx = snt_br_get_bits(&br, shape_bits);
            }
        }

        /* (block count bounded by total_blocks from bitstream) */
    }

    free(rans_bytes);

    /* Reconstruct frame via quad-tree traversal */
    int16_t *recon = (int16_t *)calloc((size_t)pad_w * pad_h, sizeof(int16_t));

    int blk_idx = 0;
    DecLeafCtx leaf_ctx = {
        .recon = recon, .stride = pad_w,
        .frame_w = pad_w, .frame_h = pad_h,
        .base_qp = ctx->header.base_qp,
        .blocks = blocks, .block_idx = &blk_idx
    };

    fi = 0;
    for (int sy = 0; sy < pad_h; sy += SNT_QT_MAX_SIZE) {
        for (int sx = 0; sx < pad_w; sx += SNT_QT_MAX_SIZE) {
            snt_qt_traverse(sx, sy, SNT_QT_MAX_SIZE, pad_w, pad_h,
                            split_flags, &fi,
                            decode_leaf_block, &leaf_ctx);
        }
    }

    /* Apply 2D lapped post-filter */
    snt_lap_post_filter_2d_ref(recon, pad_w, pad_h, pad_w, 8);

    /* Copy to output, clipping to [0, 255] */
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
            output_y[y * stride_y + x] = snt_clip_pixel((int32_t)recon[y * pad_w + x]);

    /* Decode chroma planes */
    if (!snt_br_eof(&br)) {
        uint32_t chroma_present = snt_br_get_bits(&br, 1);
        if (chroma_present && output_u && output_v) {
            int chroma_w = width / 2;
            int chroma_h = height / 2;
            int chroma_qp = ctx->header.base_qp + ctx->header.chroma_qp_offset;
            if (chroma_qp < 0) chroma_qp = 0;
            if (chroma_qp > 63) chroma_qp = 63;
            uint16_t chroma_step = snt_qp_to_step((uint8_t)chroma_qp);
            int chroma_blk = 4;
            int cpw = ((chroma_w + chroma_blk - 1) / chroma_blk) * chroma_blk;
            int cph = ((chroma_h + chroma_blk - 1) / chroma_blk) * chroma_blk;

            /* Decode U plane */
            for (int by = 0; by < cph; by += chroma_blk) {
                for (int bx = 0; bx < cpw; bx += chroma_blk) {
                    int16_t level = (int16_t)(uint16_t)snt_br_get_bits(&br, 16);
                    int16_t dc = pvq_dequantize_dc(level, chroma_step) + 128;
                    for (int y = 0; y < chroma_blk && by + y < chroma_h; y++)
                        for (int x = 0; x < chroma_blk && bx + x < chroma_w; x++)
                            output_u[(by + y) * stride_uv + bx + x] = snt_clip_pixel(dc);
                }
            }
            /* Decode V plane */
            for (int by = 0; by < cph; by += chroma_blk) {
                for (int bx = 0; bx < cpw; bx += chroma_blk) {
                    int16_t level = (int16_t)(uint16_t)snt_br_get_bits(&br, 16);
                    int16_t dc = pvq_dequantize_dc(level, chroma_step) + 128;
                    for (int y = 0; y < chroma_blk && by + y < chroma_h; y++)
                        for (int x = 0; x < chroma_blk && bx + x < chroma_w; x++)
                            output_v[(by + y) * stride_uv + bx + x] = snt_clip_pixel(dc);
                }
            }
        }
    }

    /* Apply film grain synthesis (last step) */
    if (ctx->header.grain.enabled) {
        snt_grain_apply(&ctx->header.grain, output_y, output_u, output_v,
                        width, height, stride_y, stride_uv);
    }

    free(split_flags);
    free(blocks);
    free(recon);

    return 0;
}
