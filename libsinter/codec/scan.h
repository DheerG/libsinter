/*
 * Sinter Codec — Coefficient scan order (normative)
 *
 * PVQ/scalar bands are consecutive runs of AC coefficients in zigzag
 * order, so each band holds similar spatial frequencies and energy falls
 * off band by band. (Raster order put column-0 coefficients, the low
 * vertical frequencies, into every other band.)
 */
#ifndef SNT_SCAN_H
#define SNT_SCAN_H

#include <stdint.h>

/* snt_zigzag8[k] = raster index of the k-th coefficient in zigzag order */
static const uint8_t snt_zigzag8[64] = {
     0,  1,  8, 16,  9,  2,  3, 10, 17, 24, 32, 25, 18, 11,  4,  5,
    12, 19, 26, 33, 40, 48, 41, 34, 27, 20, 13,  6,  7, 14, 21, 28,
    35, 42, 49, 56, 57, 50, 43, 36, 29, 22, 15, 23, 30, 37, 44, 51,
    58, 59, 52, 45, 38, 31, 39, 46, 53, 60, 61, 54, 47, 55, 62, 63,
};

#endif /* SNT_SCAN_H */
