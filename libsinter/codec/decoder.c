/*
 * Sinter Codec — Decoder Pipeline
 * Implements the top-level decode flow:
 *
 *   Bitstream -> rANS Entropy Decoding -> PVQ Dequantization
 *   -> Inverse DCT -> Post-Filter (lapping) -> Film Grain -> Output
 */

#include <string.h>
#include <stdlib.h>
#include "decoder.h"
#include "../util/fixed_point.h"
#include "../core/transform/lapping.h"
#include "../core/transform/dct.h"
#include "../core/quantize/pvq.h"
#include "../core/bitstream/bitstream.h"
#include "../core/prediction/intra.h"
#include "../core/grain/grain_synth.h"

/* QP to step (same table as encoder — must match exactly) */
extern uint16_t snte_qp_to_step(uint8_t qp);

int sntd_decode_frame(SntDecoderCtx *ctx, const uint8_t *bitstream,
                       uint32_t bs_len, uint8_t *output_y,
                       int *out_width, int *out_height, int stride) {
    const int BLK = 8;

    /* Read frame header */
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
    int total_blocks = blk_cols * blk_rows;

    uint16_t step = snte_qp_to_step(ctx->header.base_qp);

    /* Read block data */
    snt_br_byte_align(&br);

    int16_t *dc_levels = (int16_t *)calloc(total_blocks, sizeof(int16_t));
    SntPVQCode *pvq_codes = (SntPVQCode *)calloc(total_blocks, sizeof(SntPVQCode));
    int16_t *recon = (int16_t *)calloc((size_t)pad_w * pad_h, sizeof(int16_t));
    int16_t *coeff_buf = (int16_t *)calloc((size_t)total_blocks * BLK * BLK, sizeof(int16_t));
    if (!dc_levels || !pvq_codes || !recon || !coeff_buf) {
        free(dc_levels); free(pvq_codes); free(recon); free(coeff_buf);
        return -1;
    }

    /* Read DC levels */
    for (int i = 0; i < total_blocks; i++) {
        dc_levels[i] = (int16_t)(uint16_t)snt_br_get_bits(&br, 16);
    }

    /* Read PVQ codes */
    for (int i = 0; i < total_blocks; i++) {
        pvq_codes[i].skip = (uint8_t)snt_br_get_bits(&br, 1);
        if (!pvq_codes[i].skip) {
            pvq_codes[i].gain = snt_br_get_bits(&br, 16);
            pvq_codes[i].k = (uint16_t)snt_br_get_bits(&br, 16);
            pvq_codes[i].shape_idx = snt_br_get_bits(&br, 32);
        }
    }

    /* Decode blocks */
    for (int by = 0; by < blk_rows; by++) {
        for (int bx = 0; bx < blk_cols; bx++) {
            int bidx = by * blk_cols + bx;

            /* Intra prediction from reconstructed coefficients (frequency domain) */
            int16_t pred[64];
            int left_bidx = by * blk_cols + (bx - 1);
            int above_bidx = (by - 1) * blk_cols + bx;
            const int16_t *left_c = (bx > 0) ? &coeff_buf[left_bidx * BLK * BLK] : NULL;
            const int16_t *above_c = (by > 0) ? &coeff_buf[above_bidx * BLK * BLK] : NULL;
            snt_intra_predict_freq_ref(left_c, above_c, BLK, pred);

            /* Dequantize DC */
            int16_t *recon_coeffs = &coeff_buf[bidx * BLK * BLK];
            recon_coeffs[0] = pvq_dequantize_dc(dc_levels[bidx], step) + pred[0];

            /* Dequantize AC with PVQ */
            int16_t recon_ac[63];
            pvq_decode(&pvq_codes[bidx], 63, step, recon_ac);
            for (int i = 0; i < 63; i++)
                recon_coeffs[i + 1] = recon_ac[i] + pred[i + 1];

            /* Inverse DCT */
            snt_dct_inv_2d_ref(recon_coeffs, BLK,
                               &recon[by * BLK * pad_w + bx * BLK], pad_w);
        }
    }

    /* Apply 2D lapped post-filter */
    snt_lap_post_filter_2d_ref(recon, pad_w, pad_h, pad_w, BLK);

    /* Copy to output, clipping to [0, 255] */
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            output_y[y * stride + x] = snt_clip_pixel((int32_t)recon[y * pad_w + x]);
        }
    }

    /* Apply film grain synthesis (last step) */
    if (ctx->header.grain.enabled) {
        snt_grain_apply(&ctx->header.grain, output_y, NULL, NULL,
                        width, height, stride, 0);
    }

    free(dc_levels);
    free(pvq_codes);
    free(recon);
    free(coeff_buf);

    return 0;
}
