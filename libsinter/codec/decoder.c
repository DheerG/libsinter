/*
 * Sinter Codec — Simple Decoder (matches encoder_simple.c)
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
#include "../core/grain/grain_synth.h"

#define BLK 8

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

    int pad_w = ((width + BLK - 1) / BLK) * BLK;
    int pad_h = ((height + BLK - 1) / BLK) * BLK;
    int blk_cols = pad_w / BLK;
    int blk_rows = pad_h / BLK;

    uint16_t step = snt_qp_to_step(ctx->header.base_qp);

    snt_br_byte_align(&br);
    uint32_t total_blocks = snt_br_get_bits(&br, 16);
    (void)total_blocks; /* we know it's blk_cols * blk_rows */

    int16_t *recon = (int16_t *)calloc((size_t)pad_w * pad_h, sizeof(int16_t));

    /* Decode each block */
    for (int by = 0; by < blk_rows; by++) {
        for (int bx = 0; bx < blk_cols; bx++) {
            int16_t recon_coeffs[64];

            /* DC */
            int16_t dc_level = (int16_t)(uint16_t)snt_br_get_bits(&br, 16);
            recon_coeffs[0] = pvq_dequantize_dc(dc_level, step);

            /* AC via band-partitioned PVQ (4 bands of ~15 dims) */
            #define NUM_BANDS 4
            static const int band_bounds[NUM_BANDS + 1] = {0, 15, 31, 47, 63};

            for (int b = 0; b < NUM_BANDS; b++) {
                int start = band_bounds[b];
                int n = band_bounds[b + 1] - start;

                uint8_t skip = (uint8_t)snt_br_get_bits(&br, 1);
                if (!skip) {
                    SntPVQCode code;
                    code.skip = 0;
                    code.gain = snt_br_get_bits(&br, 16);
                    code.k = (uint16_t)snt_br_get_bits(&br, 16);
                    int shape_bits = (int)snt_br_get_bits(&br, 8);
                    if (shape_bits > 64) shape_bits = 64;
                    if (shape_bits > 32) {
                        uint64_t hi = snt_br_get_bits(&br, shape_bits - 32);
                        uint64_t lo = snt_br_get_bits(&br, 32);
                        code.shape_idx = (hi << 32) | lo;
                    } else {
                        code.shape_idx = snt_br_get_bits(&br, shape_bits);
                    }

                    int16_t band_recon[16];
                    pvq_decode(&code, n, step, band_recon);
                    for (int i = 0; i < n; i++)
                        recon_coeffs[1 + start + i] = band_recon[i];
                } else {
                    for (int i = 0; i < n; i++)
                        recon_coeffs[1 + start + i] = 0;
                }
            }

            snt_dct_inv_2d_ref(recon_coeffs, BLK,
                               &recon[by * BLK * pad_w + bx * BLK], pad_w);
        }
    }

    /* Post-filter */
    snt_lap_post_filter_2d_ref(recon, pad_w, pad_h, pad_w, BLK);

    /* Copy to output */
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
            output_y[y * stride_y + x] = snt_clip_pixel((int32_t)recon[y * pad_w + x]);

    /* Decode chroma */
    snt_br_byte_align(&br);
    if (!snt_br_eof(&br)) {
        uint32_t chroma_present = snt_br_get_bits(&br, 1);
        if (chroma_present && output_u && output_v) {
            int cw = width / 2, ch = height / 2;
            int cqp = ctx->header.base_qp + ctx->header.chroma_qp_offset;
            if (cqp < 0) cqp = 0; if (cqp > 63) cqp = 63;
            uint16_t cstep = snt_qp_to_step((uint8_t)cqp);

            for (int y = 0; y < ch; y += 4)
                for (int x = 0; x < cw; x += 4) {
                    int16_t level = (int16_t)(uint16_t)snt_br_get_bits(&br, 16);
                    int16_t dc = pvq_dequantize_dc(level, cstep) + 128;
                    for (int dy = 0; dy < 4 && y+dy < ch; dy++)
                        for (int dx = 0; dx < 4 && x+dx < cw; dx++)
                            output_u[(y+dy)*stride_uv+x+dx] = snt_clip_pixel(dc);
                }
            for (int y = 0; y < ch; y += 4)
                for (int x = 0; x < cw; x += 4) {
                    int16_t level = (int16_t)(uint16_t)snt_br_get_bits(&br, 16);
                    int16_t dc = pvq_dequantize_dc(level, cstep) + 128;
                    for (int dy = 0; dy < 4 && y+dy < ch; dy++)
                        for (int dx = 0; dx < 4 && x+dx < cw; dx++)
                            output_v[(y+dy)*stride_uv+x+dx] = snt_clip_pixel(dc);
                }
        }
    }

    /* Grain synthesis */
    if (ctx->header.grain.enabled && output_u && output_v) {
        snt_grain_apply(&ctx->header.grain, output_y, output_u, output_v,
                        width, height, stride_y, stride_uv);
    }

    free(recon);
    return 0;
}
