/*
 * Sinter Codec — Rate-Distortion Optimization (Encoder-Only)
 *
 * Placeholder for full RDO. The reference encoder uses a simplified
 * approach — PTD-based distortion with fixed lambda. Production
 * encoders would search over block partitions, prediction modes,
 * and quantizer settings.
 */

#include <stdint.h>

/* Lambda for RDO: controls quality-bitrate trade-off.
 * Higher lambda = more compression, lower quality. */
float snte_rdo_lambda(uint8_t qp) {
    /* Exponential lambda model (common in video codecs) */
    static const float lambda_table[64] = {
        0.01f, 0.02f, 0.03f, 0.04f, 0.06f, 0.08f, 0.10f, 0.13f,
        0.16f, 0.20f, 0.25f, 0.32f, 0.40f, 0.50f, 0.63f, 0.80f,
        1.0f,  1.3f,  1.6f,  2.0f,  2.5f,  3.2f,  4.0f,  5.0f,
        6.3f,  8.0f, 10.0f, 12.5f, 16.0f, 20.0f, 25.0f, 32.0f,
        40.0f, 50.0f, 63.0f, 80.0f, 100.f, 125.f, 160.f, 200.f,
        250.f, 320.f, 400.f, 500.f, 630.f, 800.f, 1000.f,1250.f,
        1600.f,2000.f,2500.f,3200.f,4000.f,5000.f,6300.f,8000.f,
        10000.f,12500.f,16000.f,20000.f,25000.f,32000.f,40000.f,50000.f
    };
    if (qp > 63) qp = 63;
    return lambda_table[qp];
}
