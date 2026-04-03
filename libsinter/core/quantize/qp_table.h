/*
 * Sinter Codec — QP to Quantization Step Mapping (Normative)
 *
 * This table is used by BOTH encoder and decoder.
 * Must be identical in both — it's part of the bitstream spec.
 */
#ifndef SNT_QP_TABLE_H
#define SNT_QP_TABLE_H

#include <stdint.h>

/* Convert QP [0..63] to quantization step size. Normative. */
uint16_t snt_qp_to_step(uint8_t qp);

#endif
