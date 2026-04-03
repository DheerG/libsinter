/*
 * Sinter Codec — QP to Quantization Step Mapping (Normative)
 */
#include "qp_table.h"

static const uint16_t qp_to_step[64] = {
     1,  1,  2,  2,  3,  3,  4,  4,  5,  6,  7,  8,  9, 10, 11, 12,
    14, 16, 18, 20, 22, 24, 28, 32, 36, 40, 44, 48, 56, 64, 72, 80,
    88, 96,112,128,144,160,176,192,224,256,288,320,352,384,448,512,
    576,640,704,768,896,1024,1152,1280,1408,1536,1792,2048,2304,2560,2816,3072
};

uint16_t snt_qp_to_step(uint8_t qp) {
    if (qp > 63) qp = 63;
    return qp_to_step[qp];
}
