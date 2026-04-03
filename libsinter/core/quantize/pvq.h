/*
 * Sinter Codec — Perceptual Vector Quantization (PVQ)
 * See pvq.c for implementation details.
 */
#ifndef SNT_PVQ_H
#define SNT_PVQ_H

#include <stdint.h>

/* PVQ encoded representation of a coefficient vector */
typedef struct {
    uint32_t gain;      /* quantized gain (L2 norm of AC coefficients) */
    uint16_t k;         /* number of pulses (shape resolution) */
    uint64_t shape_idx; /* combinatorial index encoding the shape vector */
    uint8_t  skip;      /* 1 if gain quantized to zero (skip shape) */
} SntPVQCode;

/* ── Gain-shape decomposition ─────────────────────────────────────── */

/* Compute the L2 norm (gain) of a coefficient vector.
 * Integer-only: returns floor(sqrt(sum of squares)).
 * coeffs: AC coefficient vector (DC excluded).
 * n: vector length. */
uint32_t pvq_compute_gain(const int16_t *coeffs, int n);

/* Quantize PVQ gain with dead-zone.
 * gain: unquantized gain from pvq_compute_gain.
 * step: quantization step size (from QP).
 * is_intra: 1 for intra blocks (narrower dead-zone).
 * Returns: quantized gain level. */
uint32_t pvq_quantize_gain(uint32_t gain, uint16_t step, int is_intra);

/* Dequantize PVQ gain. Exact: result = level * step. */
uint32_t pvq_dequantize_gain(uint32_t level, uint16_t step);

/* ── Shape encoding (Fischer's combinatorial enumeration) ─────────── */

/* Compute K (number of pulses) from quantized gain and CSF weight.
 * Higher gain and higher perceptual importance -> more pulses -> finer shape.
 * gain_q: quantized gain level.
 * n: vector dimension.
 * csf_weight: perceptual importance [0..255], 128 = neutral. */
uint16_t pvq_compute_k(uint32_t gain_q, int n, uint8_t csf_weight);

/* Encode shape: project AC coefficients onto hypersphere of radius K,
 * return the integer pulse vector and its combinatorial index.
 * coeffs: AC coefficients (will be normalized internally).
 * n: vector length.
 * k: number of pulses.
 * pulses: output pulse vector (n int16_t values summing to ±K).
 * Returns: combinatorial index of the pulse vector. */
uint64_t pvq_encode_shape(const int16_t *coeffs, int n, uint16_t k,
                           int16_t *pulses);

/* Decode shape: reconstruct pulse vector from combinatorial index.
 * n: vector length. k: number of pulses. idx: combinatorial index.
 * pulses: output pulse vector. */
void pvq_decode_shape(int n, uint16_t k, uint64_t idx, int16_t *pulses);

/* ── Full PVQ encode/decode ───────────────────────────────────────── */

/* Encode a block's AC coefficients using PVQ.
 * coeffs: AC coefficients (DC excluded). n: vector length.
 * step: quantization step. is_intra: intra block flag.
 * csf_weight: perceptual weight for K computation.
 * code: output PVQ code.
 * recon: output reconstructed coefficients (gain * normalized_shape). */
void pvq_encode(const int16_t *coeffs, int n, uint16_t step,
                int is_intra, uint8_t csf_weight,
                SntPVQCode *code, int16_t *recon);

/* Decode a block's AC coefficients from PVQ code.
 * code: PVQ code from bitstream. n: vector length. step: quant step.
 * recon: output reconstructed coefficients. */
void pvq_decode(const SntPVQCode *code, int n, uint16_t step,
                int16_t *recon);

/* ── Public accessor for pvq_count ─────────────────────────────────── */
uint64_t pvq_num_vectors(int n, uint16_t k);

/* ── Shape bit-width computation ──────────────────────────────────── */

/* Compute exact number of bits needed for a shape index with
 * dimension n and pulse count k: ceil(log2(pvq_count(n, k))).
 * Both encoder and decoder compute this identically from (n, k). */
int pvq_shape_bits(int n, uint16_t k);

/* ── DC scalar quantization (separate from PVQ) ──────────────────── */

/* Quantize DC coefficient with dead-zone. */
int16_t pvq_quantize_dc(int16_t dc, uint16_t step, int is_intra);

/* Dequantize DC. Exact: result = level * step. */
int16_t pvq_dequantize_dc(int16_t level, uint16_t step);

#endif /* SNT_PVQ_H */
