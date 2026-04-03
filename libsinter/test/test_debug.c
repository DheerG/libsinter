/* Quick debug test to trace the encode/decode pipeline values */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "core/bitstream/bitstream.h"
#include "core/transform/lapping.h"
#include "core/transform/dct.h"
#include "core/quantize/pvq.h"
#include "util/fixed_point.h"

int main(void) {
    /* Test 1: bitstream round-trip */
    printf("--- Bitstream round-trip ---\n");
    uint8_t buf[100];
    SntBitWriter bw;
    snt_bw_init(&bw, buf, 100);
    snt_bw_put_bits(&bw, 0xABCD, 16);
    snt_bw_put_bits(&bw, 0x1234, 16);
    snt_bw_put_bits(&bw, (uint32_t)(uint16_t)(int16_t)-42, 16);

    SntBitReader br;
    snt_br_init(&br, buf, 100);
    uint32_t v1 = snt_br_get_bits(&br, 16);
    uint32_t v2 = snt_br_get_bits(&br, 16);
    int16_t v3 = (int16_t)(uint16_t)snt_br_get_bits(&br, 16);
    printf("  Written: 0xABCD, 0x1234, -42\n");
    printf("  Read:    0x%X, 0x%X, %d\n", v1, v2, v3);
    printf("  %s\n\n", (v1==0xABCD && v2==0x1234 && v3==-42) ? "PASS" : "FAIL");

    /* Test 2: single 8x8 block encode/decode (no lapping) */
    printf("--- Single block encode/decode (flat 128, no lapping) ---\n");
    int16_t block[64], coeffs[64], recon_coeffs[64], recon_block[64];
    for (int i = 0; i < 64; i++) block[i] = 128;

    snt_dct_fwd_2d_ref(block, 8, 8, coeffs);
    printf("  DC coeff: %d, AC[1]: %d, AC[2]: %d\n", coeffs[0], coeffs[1], coeffs[2]);

    uint16_t step = 22; /* QP=20 */
    int16_t dc_level = pvq_quantize_dc(coeffs[0], step, 1);
    int16_t dc_recon = pvq_dequantize_dc(dc_level, step);
    printf("  DC: orig=%d, level=%d, recon=%d, err=%d\n",
           coeffs[0], dc_level, dc_recon, coeffs[0] - dc_recon);

    /* PVQ on AC */
    SntPVQCode pvq_code;
    int16_t recon_ac[63];
    pvq_encode(&coeffs[1], 63, step, 1, 128, &pvq_code, recon_ac);
    printf("  PVQ: skip=%d, gain=%u, K=%u\n", pvq_code.skip, pvq_code.gain, pvq_code.k);

    /* Reconstruct */
    recon_coeffs[0] = dc_recon;
    for (int i = 0; i < 63; i++) recon_coeffs[i+1] = recon_ac[i];

    snt_dct_inv_2d_ref(recon_coeffs, 8, recon_block, 8);
    printf("  Recon pixels: [%d, %d, %d, %d, ...]\n",
           recon_block[0], recon_block[1], recon_block[2], recon_block[3]);

    double mse = 0;
    for (int i = 0; i < 64; i++) {
        double d = block[i] - recon_block[i];
        mse += d*d;
    }
    mse /= 64;
    printf("  Block MSE: %.2f (PSNR: %.2f dB)\n\n", mse, 10*log10(255*255/mse));

    /* Test 3: 16x16 with lapping */
    printf("--- 16x16 flat 128 with lapping ---\n");
    int16_t frame[16*16], recon_frame[16*16];
    for (int i = 0; i < 256; i++) frame[i] = 128;

    /* Pre-filter */
    int16_t work[16*16];
    memcpy(work, frame, sizeof(frame));
    snt_lap_pre_filter_2d_ref(work, 16, 16, 16, 8);
    printf("  After pre-filter: [%d, %d, ..., %d, %d] (row 0)\n",
           work[0], work[1], work[14], work[15]);
    printf("  After pre-filter: [%d, %d, ..., %d, %d] (row 7)\n",
           work[7*16], work[7*16+1], work[7*16+14], work[7*16+15]);

    /* DCT + quantize + IDCT for each 8x8 block */
    for (int by = 0; by < 2; by++) {
        for (int bx = 0; bx < 2; bx++) {
            int16_t c[64], rc[64];
            snt_dct_fwd_2d_ref(&work[by*8*16 + bx*8], 16, 8, c);

            /* Quantize */
            rc[0] = pvq_dequantize_dc(pvq_quantize_dc(c[0], step, 1), step);
            SntPVQCode code;
            int16_t rac[63];
            pvq_encode(&c[1], 63, step, 1, 128, &code, rac);
            for (int i = 0; i < 63; i++) rc[i+1] = rac[i];

            if (by == 0 && bx == 0)
                printf("  Block(0,0): DC=%d->%d, skip=%d, gain=%u\n",
                       c[0], rc[0], code.skip, code.gain);

            snt_dct_inv_2d_ref(rc, 8, &recon_frame[by*8*16 + bx*8], 16);
        }
    }

    /* Post-filter */
    snt_lap_post_filter_2d_ref(recon_frame, 16, 16, 16, 8);

    /* Check quality */
    mse = 0;
    for (int i = 0; i < 256; i++) {
        int32_t val = snt_clip_pixel((int32_t)recon_frame[i]);
        double d = 128.0 - val;
        mse += d*d;
    }
    mse /= 256;
    printf("  Post-filter recon: [%d, %d, ..., %d, %d]\n",
           (int)snt_clip_pixel(recon_frame[0]),
           (int)snt_clip_pixel(recon_frame[1]),
           (int)snt_clip_pixel(recon_frame[14]),
           (int)snt_clip_pixel(recon_frame[15]));
    printf("  Frame MSE: %.2f (PSNR: %.2f dB)\n", mse,
           mse > 0.001 ? 10*log10(255*255/mse) : 99.99);

    return 0;
}
