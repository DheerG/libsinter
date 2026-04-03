/*
 * Sinter Codec — Film Grain Synthesis (Decoder-Side)
 * Implements Section 1.6: Film Grain Synthesis
 *
 * AR(3) luma / AR(4) chroma model.
 * Grain amplitude scales with local luma (photon noise model).
 * Seed derived from frame_number for deterministic seeking.
 *
 * The decoder synthesizes grain and adds it to the reconstructed frame
 * as the LAST post-processing step.
 */

#include <string.h>
#include "grain_synth.h"
#include "../../util/fixed_point.h"

/* Grain LUT dimensions (tiled across the frame with random offsets) */
#define GRAIN_LUT_SIZE 64

/* Simple LCG PRNG — deterministic, reproducible, adequate for grain */
static uint32_t grain_prng(uint32_t *state) {
    *state = *state * 1103515245u + 12345u;
    return (*state >> 16) & 0x7FFF;
}

/* Generate a Gaussian-ish sample from uniform PRNG (CLT approximation) */
static int16_t grain_gaussian(uint32_t *state) {
    int32_t sum = 0;
    for (int i = 0; i < 4; i++)
        sum += (int32_t)(grain_prng(state) & 0xFF) - 128;
    return (int16_t)(sum >> 1); /* range approximately [-256, 256] */
}

void snt_grain_generate_lut(const SntGrainParams *params,
                             int16_t *lut_y, int16_t *lut_u, int16_t *lut_v) {
    /* Seed from frame number for deterministic seeking */
    uint32_t state = (uint32_t)params->seed * 7919u + 1;

    /* Generate white noise base */
    int16_t white_y[GRAIN_LUT_SIZE * GRAIN_LUT_SIZE];
    int16_t white_u[GRAIN_LUT_SIZE * GRAIN_LUT_SIZE];
    int16_t white_v[GRAIN_LUT_SIZE * GRAIN_LUT_SIZE];

    for (int i = 0; i < GRAIN_LUT_SIZE * GRAIN_LUT_SIZE; i++) {
        white_y[i] = grain_gaussian(&state);
        white_u[i] = grain_gaussian(&state);
        white_v[i] = grain_gaussian(&state);
    }

    /* Apply 2D spatial AR(3) filter to luma.
     * Taps reference spatially adjacent samples: left, above, above-left.
     * This produces proper 2D spatial correlation instead of 1D raster streaking.
     * (Jens-Rainer's review: 1D raster AR creates horizontal streaking artifacts.) */
    #define G(y,x) ((y) * GRAIN_LUT_SIZE + (x))
    memcpy(lut_y, white_y, GRAIN_LUT_SIZE * GRAIN_LUT_SIZE * sizeof(int16_t));
    for (int gy = 1; gy < GRAIN_LUT_SIZE; gy++) {
        for (int gx = 1; gx < GRAIN_LUT_SIZE; gx++) {
            int32_t val = white_y[G(gy,gx)];
            val += ((int32_t)lut_y[G(gy, gx-1)]   * params->ar_coeff_y[0]) >> 8; /* left */
            val += ((int32_t)lut_y[G(gy-1, gx)]   * params->ar_coeff_y[1]) >> 8; /* above */
            val += ((int32_t)lut_y[G(gy-1, gx-1)] * params->ar_coeff_y[2]) >> 8; /* above-left */
            lut_y[G(gy,gx)] = snt_clip_i16(val);
        }
    }

    /* Apply 2D spatial AR(4) filter to chroma.
     * 4th tap = co-located luma grain (chroma-correlated). */
    memcpy(lut_u, white_u, GRAIN_LUT_SIZE * GRAIN_LUT_SIZE * sizeof(int16_t));
    memcpy(lut_v, white_v, GRAIN_LUT_SIZE * GRAIN_LUT_SIZE * sizeof(int16_t));
    for (int gy = 1; gy < GRAIN_LUT_SIZE; gy++) {
        for (int gx = 1; gx < GRAIN_LUT_SIZE; gx++) {
            int32_t val_u = white_u[G(gy,gx)];
            val_u += ((int32_t)lut_u[G(gy, gx-1)]   * params->ar_coeff_u[0]) >> 8;
            val_u += ((int32_t)lut_u[G(gy-1, gx)]   * params->ar_coeff_u[1]) >> 8;
            val_u += ((int32_t)lut_u[G(gy-1, gx-1)] * params->ar_coeff_u[2]) >> 8;
            val_u += ((int32_t)lut_y[G(gy,gx)]      * params->ar_coeff_u[3]) >> 8;
            lut_u[G(gy,gx)] = snt_clip_i16(val_u);

            int32_t val_v = white_v[G(gy,gx)];
            val_v += ((int32_t)lut_v[G(gy, gx-1)]   * params->ar_coeff_v[0]) >> 8;
            val_v += ((int32_t)lut_v[G(gy-1, gx)]   * params->ar_coeff_v[1]) >> 8;
            val_v += ((int32_t)lut_v[G(gy-1, gx-1)] * params->ar_coeff_v[2]) >> 8;
            val_v += ((int32_t)lut_y[G(gy,gx)]      * params->ar_coeff_v[3]) >> 8;
            lut_v[G(gy,gx)] = snt_clip_i16(val_v);
        }
    }
    #undef G

    /* Scale by strength */
    for (int i = 0; i < GRAIN_LUT_SIZE * GRAIN_LUT_SIZE; i++) {
        lut_y[i] = (int16_t)(((int32_t)lut_y[i] * params->strength_y) >> 8);
        lut_u[i] = (int16_t)(((int32_t)lut_u[i] * params->strength_u) >> 8);
        lut_v[i] = (int16_t)(((int32_t)lut_v[i] * params->strength_v) >> 8);
    }
}

void snt_grain_apply(const SntGrainParams *params,
                      uint8_t *frame_y, uint8_t *frame_u, uint8_t *frame_v,
                      int width, int height, int stride_y, int stride_uv) {
    if (!params->enabled) return;

    /* Generate grain LUTs */
    int16_t lut_y[GRAIN_LUT_SIZE * GRAIN_LUT_SIZE];
    int16_t lut_u[GRAIN_LUT_SIZE * GRAIN_LUT_SIZE];
    int16_t lut_v[GRAIN_LUT_SIZE * GRAIN_LUT_SIZE];
    snt_grain_generate_lut(params, lut_y, lut_u, lut_v);

    /* Tile grain LUT across the frame with per-block random offsets.
     * Grain amplitude scales with local luma (photon noise model). */
    uint32_t offset_state = (uint32_t)params->seed * 3571u + 7;

    for (int by = 0; by < height; by += GRAIN_LUT_SIZE) {
        for (int bx = 0; bx < width; bx += GRAIN_LUT_SIZE) {
            /* Random offset into the LUT for this tile */
            int ox = grain_prng(&offset_state) % GRAIN_LUT_SIZE;
            int oy = grain_prng(&offset_state) % GRAIN_LUT_SIZE;

            for (int y = 0; y < GRAIN_LUT_SIZE && by + y < height; y++) {
                for (int x = 0; x < GRAIN_LUT_SIZE && bx + x < width; x++) {
                    int gy = (y + oy) % GRAIN_LUT_SIZE;
                    int gx = (x + ox) % GRAIN_LUT_SIZE;
                    int gi = gy * GRAIN_LUT_SIZE + gx;

                    int py = by + y;
                    int px = bx + x;

                    /* Luma: scale grain by local luma value (photon noise model).
                     * Brighter regions get more grain. */
                    uint8_t luma = frame_y[py * stride_y + px];
                    int32_t grain_scaled = ((int32_t)lut_y[gi] * luma) >> 8;
                    int32_t val = (int32_t)luma + grain_scaled;
                    frame_y[py * stride_y + px] = snt_clip_pixel(val);
                }
            }

            /* Chroma (half resolution) */
            int chroma_h = height / 2;
            int chroma_w = width / 2;
            int cbx = bx / 2;
            int cby = by / 2;
            int clut = GRAIN_LUT_SIZE / 2;

            for (int y = 0; y < clut && cby + y < chroma_h; y++) {
                for (int x = 0; x < clut && cbx + x < chroma_w; x++) {
                    int gi = ((y * 2 + oy) % GRAIN_LUT_SIZE) * GRAIN_LUT_SIZE +
                             ((x * 2 + ox) % GRAIN_LUT_SIZE);
                    int py = cby + y;
                    int px = cbx + x;

                    /* Scale chroma grain by co-located luma (photon noise model).
                     * Jim's review: spec says all planes scale by luma. */
                    uint8_t coloc_luma = frame_y[(py * 2) * stride_y + (px * 2)];
                    int32_t val_u = (int32_t)frame_u[py * stride_uv + px]
                                  + (((int32_t)lut_u[gi] * coloc_luma) >> 8);
                    int32_t val_v = (int32_t)frame_v[py * stride_uv + px]
                                  + (((int32_t)lut_v[gi] * coloc_luma) >> 8);
                    frame_u[py * stride_uv + px] = snt_clip_pixel(val_u);
                    frame_v[py * stride_uv + px] = snt_clip_pixel(val_v);
                }
            }
        }
    }
}
