/*
 * Sinter Codec — Daala-Style Per-Dimension PVQ Shape Coding
 *
 * Instead of encoding a flat combinatorial index, encode each dimension's
 * pulse value as an rANS symbol using the exact combinatorial counts
 * as the probability distribution. This achieves theoretically optimal
 * compression for PVQ shape vectors.
 *
 * At each dimension i with remaining budget k_left:
 *   P(|x_i| = j) = pvq_count(remaining_dims, k_left - j) / pvq_count(n-i, k_left)
 *   for j = 0: unsigned count
 *   for j > 0: count * 2 (positive + negative)
 */
#ifndef SNT_PVQ_SHAPE_RANS_H
#define SNT_PVQ_SHAPE_RANS_H

#include <stdint.h>
#include "../entropy/rans.h"

/* Encode PVQ shape vector per-dimension through rANS.
 * pulses: the pulse vector to encode (from pvq_project).
 * n: vector dimension. k: total pulses.
 * enc: rANS encoder to write symbols to.
 * fwd_idx: starting forward index for dual rANS state selection. */
void pvq_shape_encode_rans(const int16_t *pulses, int n, uint16_t k,
                            RansEncoder *enc, uint32_t fwd_idx);

/* Decode PVQ shape vector per-dimension from rANS.
 * pulses: output pulse vector.
 * n: vector dimension. k: total pulses.
 * dec: rANS decoder to read symbols from. */
void pvq_shape_decode_rans(int16_t *pulses, int n, uint16_t k,
                            RansDecoder *dec);

#endif
