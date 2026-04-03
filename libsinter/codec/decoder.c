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
#include "../core/quantize/qp_table.h"
#include "../core/entropy/rans.h"

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

    uint16_t step = snt_qp_to_step(ctx->header.base_qp);

    /* Read block data: rANS-coded symbols + raw supplementary data */
    snt_br_byte_align(&br);

    int16_t *dc_levels = (int16_t *)calloc(total_blocks, sizeof(int16_t));
    SntPVQCode *pvq_codes = (SntPVQCode *)calloc(total_blocks, sizeof(SntPVQCode));
    int8_t *qp_deltas = (int8_t *)calloc(total_blocks, sizeof(int8_t));
    int16_t *recon = (int16_t *)calloc((size_t)pad_w * pad_h, sizeof(int16_t));
    int16_t *coeff_buf = (int16_t *)calloc((size_t)total_blocks * BLK * BLK, sizeof(int16_t));
    if (!dc_levels || !pvq_codes || !recon || !coeff_buf || !qp_deltas) {
        free(dc_levels); free(pvq_codes); free(recon);
        free(coeff_buf); free(qp_deltas);
        return -1;
    }

    /* Read rANS payload */
    uint32_t rans_len = snt_br_get_bits(&br, 16);
    uint8_t *rans_bytes = (uint8_t *)calloc(rans_len * 2, 1);
    for (uint32_t i = 0; i < rans_len; i++) {
        uint8_t lo_byte = (uint8_t)snt_br_get_bits(&br, 8);
        uint8_t hi_byte = (uint8_t)snt_br_get_bits(&br, 8);
        rans_bytes[i * 2]     = lo_byte;
        rans_bytes[i * 2 + 1] = hi_byte;
    }

    /* Initialize rANS decoder */
    RansDecoder rans_dec;
    rans_dec_init(&rans_dec, rans_bytes, rans_len * 2, ctx->header.dual_rans);

    /* Initialize CDFs (must match encoder's initial state exactly) */
    RansCDF cdf_skip, cdf_gain_lo, cdf_gain_hi, cdf_qp_delta, cdf_dc_abs, cdf_k;
    rans_cdf_init_uniform(&cdf_skip, 2, RANS_ADAPT_RATE_MODE);
    rans_cdf_init_uniform(&cdf_gain_lo, 256, RANS_ADAPT_RATE_COEFF);
    rans_cdf_init_uniform(&cdf_gain_hi, 256, RANS_ADAPT_RATE_COEFF);
    rans_cdf_init_uniform(&cdf_qp_delta, 13, RANS_ADAPT_RATE_MODE);
    rans_cdf_init_uniform(&cdf_dc_abs, 256, RANS_ADAPT_RATE_COEFF);
    rans_cdf_init_uniform(&cdf_k, 256, RANS_ADAPT_RATE_COEFF);

    /* Decode rANS symbols (forward order, must match encoder's symbol order:
     * qp_delta, dc_abs_lo, skip, [gain_lo, gain_hi, k_lo if !skip]) */
    for (int i = 0; i < total_blocks; i++) {
        /* QP delta */
        uint16_t qpd_sym = rans_dec_symbol(&rans_dec, &cdf_qp_delta);
        qp_deltas[i] = (int8_t)((int)qpd_sym - 6);

        /* DC abs low byte */
        uint16_t dc_abs_lo = rans_dec_symbol(&rans_dec, &cdf_dc_abs);
        dc_levels[i] = (int16_t)dc_abs_lo; /* sign + high byte from raw section */

        /* Skip flag */
        uint16_t skip_sym = rans_dec_symbol(&rans_dec, &cdf_skip);
        pvq_codes[i].skip = (uint8_t)skip_sym;

        /* PVQ gain (full) and K if not skip */
        if (!pvq_codes[i].skip) {
            pvq_codes[i].gain = rans_dec_symbol(&rans_dec, &cdf_gain_lo);
            uint16_t gain_hi = rans_dec_symbol(&rans_dec, &cdf_gain_hi);
            pvq_codes[i].gain |= ((uint32_t)gain_hi << 8);
            uint16_t k_lo = rans_dec_symbol(&rans_dec, &cdf_k);
            pvq_codes[i].k = k_lo; /* high byte from raw section */
        }
    }

    free(rans_bytes);

    /* Read supplementary raw bits (overflow bytes, signs, shape) */
    for (int i = 0; i < total_blocks; i++) {
        /* DC sign + high byte */
        uint32_t dc_sign = snt_br_get_bits(&br, 1);
        uint32_t dc_hi = snt_br_get_bits(&br, 8);
        int16_t dc_abs = (int16_t)(dc_levels[i] | (int16_t)(dc_hi << 8));
        dc_levels[i] = dc_sign ? -dc_abs : dc_abs;

        if (!pvq_codes[i].skip) {
            /* K high byte */
            uint32_t k_hi = snt_br_get_bits(&br, 8);
            pvq_codes[i].k |= (uint16_t)(k_hi << 8);
            /* Shape index: variable-width, same computation as encoder */
            int shape_bits = pvq_shape_bits(63, pvq_codes[i].k);
            if (shape_bits < 1) shape_bits = 1;
            if (shape_bits > 64) shape_bits = 64;

            pvq_codes[i].shape_idx = 0;
            if (shape_bits > 32) {
                uint64_t hi = snt_br_get_bits(&br, shape_bits - 32);
                uint64_t lo = snt_br_get_bits(&br, 32);
                pvq_codes[i].shape_idx = (hi << 32) | lo;
            } else {
                pvq_codes[i].shape_idx = snt_br_get_bits(&br, shape_bits);
            }
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

            /* Apply per-block QP delta for dequantization */
            int eff_qp = ctx->header.base_qp + qp_deltas[bidx];
            if (eff_qp < 0) eff_qp = 0;
            if (eff_qp > 63) eff_qp = 63;
            uint16_t blk_step = snt_qp_to_step((uint8_t)eff_qp);

            /* Dequantize DC */
            int16_t *recon_coeffs = &coeff_buf[bidx * BLK * BLK];
            recon_coeffs[0] = pvq_dequantize_dc(dc_levels[bidx], blk_step) + pred[0];

            /* Dequantize AC with PVQ */
            int16_t recon_ac[63];
            pvq_decode(&pvq_codes[bidx], 63, blk_step, recon_ac);
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
    free(qp_deltas);
    free(recon);
    free(coeff_buf);

    return 0;
}
