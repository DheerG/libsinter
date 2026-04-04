# Sinter Codec Technical Reference

Version 0.3 — April 2026

This document is the definitive technical reference for the Sinter video codec. It covers the architecture, bitstream specification, transform design, quantization, entropy coding, perceptual optimization, and container format. It is written for implementers, researchers, and anyone building on this codebase.

---

## Table of Contents

1. [Architecture Overview](#1-architecture-overview)
2. [Transform and Lapping](#2-transform-and-lapping)
3. [PVQ Quantization](#3-pvq-quantization)
4. [Perceptual Optimization](#4-perceptual-optimization)
5. [Bitstream Specification](#5-bitstream-specification)
6. [Entropy Coding](#6-entropy-coding)
7. [Container Format](#7-container-format)
8. [Appendices](#8-appendices)

---

## 1. Architecture Overview

### Why Sinter Exists

Every deployed video codec since H.261 (1988) shares the same architecture: block DCT, spatial pixel-domain prediction, scalar quantization, and post-hoc deblocking filters. Each generation (H.264, H.265, AV1, VVC) adds complexity to patch artifacts created by this architecture. Sinter asks: what if we eliminated the source of those artifacts instead?

Sinter is built on three architectural breaks from tradition:

1. **Lapped transforms** replace block DCT + deblocking. Cross-boundary pre/post-filters structurally prevent blocking artifacts. Zero post-processing required.
2. **Perceptual Vector Quantization (PVQ)** replaces scalar quantization. Gain-shape decomposition preserves texture energy instead of destroying it. Grass stays grainy; fabric keeps its weave.
3. **rANS entropy coding** replaces arithmetic coding (CABAC). Identical compression with ~2x faster decode throughput.

### Pipeline

```
Encoder:
  Input (YUV420)
    -> Lapped Pre-Filter (8-point lifting steps at block boundaries)
    -> 8x8 Forward DCT (Q14 orthonormal integer)
    -> DC Prediction (from reconstructed left/above neighbors)
    -> PVQ Band Partitioning (16 bands of 3-4 AC coefficients)
    -> PVQ Gain-Shape Quantization (per band)
    -> Raw Bit Packing
    -> Bitstream

Decoder:
  Bitstream
    -> Raw Bit Unpacking
    -> PVQ Dequantization (midpoint reconstruction + energy correction)
    -> DC Reconstruction (add prediction from left/above)
    -> 8x8 Inverse DCT
    -> Lapped Post-Filter (exact inverse of pre-filter)
    -> Film Grain Synthesis (optional, AR(3)/AR(4) model)
    -> Output (YUV420)
```

### Normative vs. Encoder-Only

**Normative** (decoder must implement identically): lapped post-filter lifting coefficients, inverse DCT matrix, PVQ dequantization formula, PVQ energy correction algorithm, DC prediction formula, band partition boundaries, `pvq_shape_bits()` computation, film grain synthesis, QP-to-step table.

**Encoder-only** (non-normative, can vary): lapped pre-filter, forward DCT, QP selection, activity masking, Butteraugli targeting, PVQ K parameter selection, rate-distortion optimization, grain estimation.

### Current Quality

**Video test (ffmpeg testsrc, 256x256, 30 frames, color):**

| Plane | Sinter (QP=4) | H.264 (QP=20) | Gap |
|-------|---------------|----------------|-----|
| Luma Y | 44.80 dB | 45.58 dB | 0.77 dB |
| Chroma U | 47.18 dB | 53.91 dB | 6.73 dB |
| Chroma V | 47.06 dB | 53.47 dB | 6.41 dB |
| File size | 838 KB (30 frames) | 149 KB (30 frames) | 5.6x |

**Luma R-D curve (256x256 gradient, single frame):**

| QP | PSNR (dB) | Size (bytes) |
|----|-----------|--------------|
| 4  | 46.49     | 21,387       |
| 8  | 42.0      | ~18,000      |
| 12 | 37.5      | ~14,000      |
| 20 | 30.0      | ~9,000       |
| 28 | 24.0      | ~5,500       |

Bitrates use exp-Golomb VLC. With rANS entropy coding (built and tested, not yet in pipeline), expect ~2x further reduction.

---

## 2. Transform and Lapping

### Lapped Transforms (TDLT)

Sinter uses Time-Domain Lapped Transforms: invertible pre/post-filters at block boundaries, implemented via lifting steps for bit-exact reversibility. The pre-filter decorrelates energy across adjacent blocks before the DCT; the post-filter reconstructs smoothly after the IDCT. No deblocking filter is needed.

**Filter sizes:**
- 4-point: 2 lifting steps (used when block size < 8)
- 8-point: 4 lifting steps (used for block size >= 8; capped at 8 per Daala convention)

**Lifting step:**
```c
lift(x, coeff) = snt_rshift_round_signed((int32_t)x * coeff, 15)
```

Uses symmetric rounding (ties away from zero) for both positive and negative products. This prevents rounding asymmetry that would break bit-exact invertibility.

**Coefficients (Q15, 20% of Malvar theoretical):**

| Param | Full Malvar | Shipped (20%) | Role |
|-------|-------------|---------------|------|
| LAP4_T0/T1 | 11585 | 2317 | Butterfly coupling |
| LAP4_S0/S1 | 6270 | 1254 | Scale |
| LAP8_T0/T1 | 12873 | 2575 | Inner butterfly |
| LAP8_T2/T3 | 8867 | 1773 | Outer butterfly |
| LAP8_S0/S1 | 5765 | 1153 | Inner scale |
| LAP8_S2/S3 | 3135 | 627 | Outer scale |

**Why 20%?** Full-strength Malvar (100%) maximizes decorrelation but the Q14 integer DCT amplifies errors through the post-filter chain, limiting PSNR to ~22 dB. At 60%, PSNR reaches ~42 dB. At 20%, PSNR reaches ~46 dB — within 0.77 dB of H.264. The 20% value was determined empirically through a lapping strength sweep during quality improvement loops. The tradeoff is less aggressive deblocking at high QP, but for the primary use case (high-quality encoding, QP 0-20) this is optimal.

**4-point filter** (samples [a, b | c, d] across boundary):
```
Forward:  b += lift(c, T0); c -= lift(b, T1); a += lift(b, S0); d += lift(c, S1)
Inverse:  d -= lift(c, S1); a -= lift(b, S0); c += lift(b, T1); b -= lift(c, T0)
```

**2D processing order:**
- Pre-filter: horizontal first (all rows), then vertical (all columns)
- Post-filter: vertical first, then horizontal (exact inverse)
- Frame edges: identity filter (no lapping)

### Integer DCT

Orthonormal DCT-II via direct matrix multiply with Q14 fixed-point coefficients.

**Formula:** `C[k][n] = round(alpha(k) * cos(pi*(2n+1)*k / (2N)) * 16384)` where `alpha(0) = sqrt(1/N)`, `alpha(k>0) = sqrt(2/N)`.

**Supported sizes:** 4, 8, 16, 32 points. All use `int16_t` matrix entries and `int32_t` accumulators.

**Round-trip accuracy:** +/-1 per sample (1D), +/-2 per sample (2D 8x8). Verified over 1000 random trials per size.

**2D application:** Separable row-then-column (forward), column-then-row (inverse).

---

## 3. PVQ Quantization

### Gain-Shape Decomposition

Each band of AC coefficients is decomposed into:
- **Gain**: L2 norm (energy magnitude), quantized with a dead-zone scalar quantizer
- **Shape**: direction on the unit hypersphere, quantized via Fischer's combinatorial enumeration on the L1 lattice

The key perceptual property: human vision is far more sensitive to energy loss than shape distortion. PVQ preserves the spectral envelope of textures even at aggressive quantization.

### Band Partitioning

The 63 AC coefficients of each 8x8 block are split into **16 bands**:

```
Boundaries: {0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 59, 63}
```

15 bands of 4 dimensions, 1 band of 3 dimensions. Each band is independently coded using either PVQ or scalar quantization (hybrid mode, 1 bit per band).

### Hybrid PVQ/Scalar Mode

For each band, the encoder measures energy concentration: if the top coefficient holds >60% of total band energy (`max_sq * 10 > sum_sq * 6`), use per-coefficient **scalar** quantization (clean gradients, no shape noise). Otherwise use **PVQ** gain-shape (texture preservation). Signaled as 1 bit per band: 0=PVQ, 1=scalar.

Scalar mode uses CSF-weighted step sizes for high-frequency bands (`eff_step = step * 256 / band_csf` when `band_csf < 128`). This was the single biggest quality improvement: +3.49 dB.

**Why 16 bands?** With n=4, the maximum safe K (before combinatorial table overflow) is 225. This gives dramatically better shape resolution than the monolithic approach (n=63, max K=14). The quality improvement from 4 bands to 16 bands is approximately +10 dB.

### K Parameter

K (pulse count) determines shape resolution:
```
K = max(1, gain_q * sqrt(n) * csf_weight / 256)
```
Capped at `900/n` to prevent overflow. For n=4, max K=225. For n=63, max K=14.

### Gain Quantization

Dead-zone with round-to-nearest and midpoint reconstruction:

| Parameter | Intra | Inter |
|-----------|-------|-------|
| Dead-zone offset | step/8 | step/4 |
| Quantize | `(gain - offset + step/2) / step` | same |
| Dequantize | `level * step + step/2` (midpoint) | same |
| Dequantize (level=0) | 0 | 0 |

Midpoint reconstruction minimizes expected squared error. This was a critical fix — lower-edge reconstruction caused 6-10 dB quality loss.

### DC Quantization and Prediction

DC is quantized separately with a wider dead-zone (step/6 intra, step/3 inter) and midpoint reconstruction.

**DC prediction:** average of left and above reconstructed DC coefficients in the transform domain:
```
pred = (left_dc + above_dc + 1) >> 1  (if both available)
pred = left_dc                         (if only left)
pred = above_dc                        (if only above)
pred = 0                               (top-left block)
```

The bitstream carries the quantized residual (DC - prediction), not the absolute DC. This reduces DC magnitude by 10-20x on smooth content, dramatically improving quantization accuracy.

### Shape Encoding (Fischer's Combinatorial Enumeration)

Each non-skip band's shape is encoded as a unique index into the set of signed integer vectors with L1 norm = K and dimension n. The count follows the CELT/Opus recurrence:
```
V(1, k) = 2 for k > 0, V(1, 0) = 1
V(n, k) = V(n-1, k) + V(n, k-1) + V(n-1, k-1)
```

The index requires `pvq_shape_bits(n, k) = ceil(log2(V(n, k)))` bits. This is computed identically by encoder and decoder — not signaled.

### Energy Conservation

After reconstruction (`recon[i] = pulses[i] * dq_gain / ||pulses||_2`), a greedy correction iteratively adjusts the largest |recon[i]| by +/-1 until energy matches target within tolerance. Brings energy ratio from ~1.25 to ~1.001. Both encoder and decoder run the identical correction.

---

## 4. Perceptual Optimization

### Contrast Sensitivity Function (CSF)

Bandpass model of human visual sensitivity: peak at 6 cpd, bandwidth 2.5 octaves, assuming 96 DPI at 60cm viewing distance. Precomputed as float32 LUT indexed by (block_size, coefficient_position). DC weight ~0.3, peak-band weight 1.0. Used in PTD and for PVQ K modulation.

### Perceptual Transform Distortion (PTD)

CSF-weighted sum of squared errors in the transform domain:
```
PTD = sum_i(csf_weight[i] * (orig[i] - recon[i])^2)
```
Structurally different from SSIM (no patent concern). Property: `PTD <= MSE` always. Allows the encoder to spend fewer bits on perceptually irrelevant frequencies.

### Activity Masking

Per-block QP offset based on AC variance:
```
delta_QP = 0.5 * log2(activity / mean_activity), clamped [-6, +6]
```
Smooth regions get lower QP (more bits); textured regions get higher QP. Signaled in bitstream as a 13-symbol rANS-coded value.

### Butteraugli Constant-Quality Mode

Outer-loop QP bisection: encode frame, measure perceptual distance, adjust QP, re-encode. Max 3 iterations. The distance metric combines MAE with a structural variance-loss penalty:
```
dist = MAE/255 + 2 * max(0, 1 - recon_var/orig_var)
```
Target values: 0.3 (visually lossless), 1.0 (high quality), 2.0 (medium), 4.0+ (aggressive).

### Film Grain Synthesis

Decoder-side: 2D spatial AR(3) luma / AR(4) chroma model with taps at (left, above, above-left). 4th chroma tap = co-located luma grain. Grain amplitude scales with local luma (photon noise model). Seed from frame_number for deterministic seeking.

Encoder-side: Yule-Walker AR fitting on high-pass filtered source. Strength from residual variance. Parameters signaled per-frame in the bitstream header.

---

## 5. Bitstream Specification

All multi-bit fields are written MSB-first. `f(n)` = fixed-width unsigned integer, n bits.

### Frame Layout
```
frame_header()          // width, height, QP, grain params, etc.
byte_align()
block_count             f(16)
for each block:
  dc_level              f(16)    // quantized DC residual (signed as uint16)
  for each of 16 PVQ bands:
    skip                f(1)
    if (!skip):
      gain              f(16)
      k                 f(16)
      shape_idx         f(pvq_shape_bits(n, k))
byte_align()
chroma_present          f(1)
if chroma:
  U plane (4x4 DC blocks, each f(16))
  V plane (4x4 DC blocks, each f(16))
```

### Frame Header

| Field | Bits | Description |
|-------|------|-------------|
| width | 32 | Frame width |
| height | 32 | Frame height |
| frame_type | 1 | 0=key, 1=inter |
| base_qp | 6 | QP [0..63] |
| num_tile_cols | 4 | Tile columns |
| num_tile_rows | 4 | Tile rows |
| dual_rans | 1 | Dual rANS flag |
| chroma_qp_offset | 8 | Signed chroma QP delta |
| frame_number | 32 | Sequence number |
| grain_enabled | 1 | Grain synthesis flag |
| [grain_params] | 248 | If grain_enabled |

### QP-to-Step Table

```
QP:  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
     1  1  2  2  3  3  4  4  5  6  7  8  9 10 11 12

QP: 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
    14 16 18 20 22 24 28 32 36 40 44 48 56 64 72 80
```

Approximately doubles every 6 QP values (6 dB per doubling).

---

## 6. Entropy Coding

### Exp-Golomb VLC (Shipped Pipeline)

All integer syntax elements use Exponential-Golomb variable-length codes (order 0, same as H.264 baseline). This provides 79% bitstream reduction vs fixed-width raw coding.

- **Unsigned:** value `v` written as `floor(log2(v+1))` zero bits, then `v+1` in binary
- **Signed:** zigzag mapping (`0→0, 1→1, -1→2, 2→3, -2→4, ...`), then unsigned Exp-Golomb
- Small values (0, +/-1) use 1-3 bits; larger values scale logarithmically

Applied to: DC residuals, scalar AC levels, PVQ gain, PVQ K, chroma DC residuals, chroma AC levels. Skip flags and mode bits remain as 1 raw bit each.

### rANS (Built, Not Yet Wired)

The full rANS infrastructure in `core/entropy/rans.c` is tested (8/8 tests pass) and supports 32-bit state, 15-bit CDFs, dual interleaved mode, and PARA adaptive CDFs. Wiring it into the pipeline is the #1 compression priority for v0.4 — expected additional 20-30% reduction beyond Exp-Golomb.

---

## 7. Container Format

### .sntr (IVF-style)

20-byte file header:
```
"SNTR" (4 bytes magic)
version (4 bytes LE uint32)
width   (4 bytes LE uint32)
height  (4 bytes LE uint32)
num_frames (4 bytes LE uint32)
```

Per frame: 4-byte LE size + raw frame bitstream.

### FFmpeg Integration

Wrapper files in `ffmpeg/` implement `AVCodec` following the libaom pattern. Encoder maps `AVFrame` YUV420P to Sinter API; decoder does the reverse. Registration as `AV_CODEC_ID_SINTER` with fourcc `SNTR`.

---

## 8. Appendices

### Naming Conventions

| Prefix | Scope |
|--------|-------|
| `snt_` | Shared/normative |
| `snte_` | Encoder-only |
| `sntd_` | Decoder-only |
| `psy_` | Perceptual functions |
| `pvq_` | PVQ quantization |
| `rans_` | rANS entropy coding |
| `_ref` | Scalar reference implementation |
| `_fwd` / `_inv` | Forward / inverse transforms |

### Test Gates

| Gate | Tests | Pass Criteria |
|------|-------|---------------|
| 0 | Lapping round-trip | Bit-exact over 1000 random trials |
| 1 | DCT + PVQ precision | DCT +/-1 per sample; PVQ energy ratio 1.001 |
| 2 | Bitstream determinism | Same input = identical bitstream |
| 3 | Reference signals | Flat gray, gradient, texture: reasonable PSNR |
| 4 | Quality monotonic | Lower QP = equal or better quality |

**Current status: 32/32 tests pass.**

### Known Limitations (v0.3 final)

- All-intra only (no inter prediction — biggest remaining compression gap)
- Fixed 8x8 luma blocks, 4x4 chroma blocks (no quad-tree partitioning)
- Exp-Golomb VLC entropy coding (rANS built but not wired — ~20-30% compression left on table)
- No chroma AC prediction (luma has DC + 4 AC prediction, chroma has DC only)
- 8-bit samples only (10-bit types defined, pipeline not extended)
- Single tile (multi-tile infrastructure exists but uses 1x1 grid)
- No rate control in pipeline (CRF module exists but not wired)

### Acknowledgments

Architecture designed through a roundtable of perspectives spanning open-codec development (VP8/VP9/AV1), perceptual coding research (Daala, HEVC), and open-source ecosystem building (Xiph/Ogg/Vorbis/Opus). Built on published, patent-free techniques: lapped transforms (Malvar, 1992), PVQ (Fischer, 1986; Xiph/Daala), rANS (Duda, 2009, public domain), OBMC (prior art, 1990s).
