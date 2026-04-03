# Sinter — An Experimental Video Codec

Sinter is a bytecode experiment exploring a novel video codec architecture built on
lapped transforms, perceptual vector quantization, and asymmetric numeral systems.
This is a research prototype — not a production codec.

## Why Sinter Exists

Traditional hybrid video codecs (H.264, H.265, AV1) share the same fundamental
architecture: block-based DCT, spatial intra prediction, scalar quantization, and
post-hoc deblocking filters to hide the artifacts that architecture creates. Sinter
asks: what if we eliminated the source of blocking artifacts instead of patching them?

The answer is **lapped transforms** — overlapping pre/post filters at block boundaries
that structurally prevent blocking. Combined with **perceptual vector quantization** that
preserves texture energy rather than destroying it, and **rANS entropy coding** that
decodes faster than arithmetic coding, Sinter explores a fundamentally different
compression pipeline.

## Architecture Overview

```
Encoder Pipeline:
  Input Frame
    -> Tile Partitioning (independent parallel units)
    -> Block Partitioning (quad-tree, 128x128 down to 4x4)
    -> Prediction (frequency-domain intra / OBMC inter + affine + DMVR)
    -> Residual Computation
    -> Pre-Filter (lapped transform boundary overlap)
    -> Forward DCT
    -> PVQ Quantization (gain-shape decomposition)
    -> rANS Entropy Coding (dual interleaved, LIFO)
    -> Bitstream Output

Decoder Pipeline:
  Bitstream Input
    -> rANS Entropy Decoding (dual interleaved, FIFO)
    -> PVQ Dequantization (gain-shape reconstruction)
    -> Inverse DCT
    -> Post-Filter (lapped transform reconstruction)
    -> Motion Compensation / OBMC Blending
    -> Film Grain Synthesis
    -> Output Frame
```

## Core Design Decisions

### 1. Lapped Transforms (No Deblocking Filter)

Every block-based codec produces blocking artifacts at quantized block boundaries.
H.264 added a deblocking filter. H.265 added SAO on top. AV1 added CDEF and loop
restoration. Each generation adds another post-processing layer to hide the same
fundamental problem.

Sinter uses **Time-Domain Lapped Transforms (TDLT)**: invertible pre/post filters
that overlap across block boundaries, applied via lifting steps for bit-exact
reversibility. The pre-filter decorrelates energy across boundaries before the DCT;
the post-filter reconstructs smoothly after inverse DCT. Blocking artifacts are
structurally impossible.

- 4-point filter: 2 lifting steps
- 8-point filter: 4 lifting steps
- Processing order: raster scan, horizontal before vertical
- Frame edges: identity filter (no lapping)
- Tile boundaries: 4-pixel overlap zones with raised cosine blending
- Switchable: encoder can choose block-DCT (identity filter) per superblock when
  lapping hurts (sharp edges, text, synthetic content)

### 2. Perceptual Vector Quantization (PVQ)

Conventional codecs use scalar quantization: each DCT coefficient is independently
divided by a step size. This destroys texture — at low bitrates, grass becomes plastic,
fabric becomes smooth, and film grain vanishes into mush.

PVQ (from Xiph/CELT/Opus) decomposes each coefficient vector into **gain** (energy
magnitude) and **shape** (energy distribution on the unit hypersphere):

- **Gain** is quantized precisely with a dead-zone scalar quantizer
- **Shape** is quantized coarsely via combinatorial enumeration on an integer lattice

The key insight: human perception is far more sensitive to energy loss than to shape
distortion. PVQ preserves the spectral envelope of textures even at aggressive
quantization. The result looks "grainy but textured" rather than "smooth and blurry."

- K parameter (pulses per band) controls the gain-shape trade-off
- Energy conservation invariant: `sum(decoded[i]^2) == gain_decoded^2` within 1 LSB
- Fischer's algorithm for exact combinatorial enumeration (requires 64-bit math)
- DC coefficients use standard scalar quantization with dead-zone

### 3. rANS Entropy Coding

Arithmetic coding (CABAC in H.264/H.265) processes one bit at a time, creating a
serial bottleneck. Sinter uses **range Asymmetric Numeral Systems (rANS)**, which
achieves identical compression ratios with ~2x faster decoding.

- 32-bit state, 16-bit renormalization
- 15-bit CDF probability resolution (totals sum to 32768)
- Normative dual interleaved rANS states (Giesen's approach) with frame-header flag
- Multi-symbol context-adaptive probability models
- PARA (Probability Adaptation Rate Adjustment): rate=5 for coefficients, rate=4 for modes
- CDF reset per-tile (tiles have independent entropy state)
- **LIFO ordering**: encode in reverse, decode forward

### 4. Frequency-Domain Intra Prediction

Traditional codecs predict pixels from neighboring pixels (spatial domain), creating
serial dependency chains that limit parallelism. Sinter predicts DCT coefficients
from neighboring blocks' DCT coefficients directly:

- DC and low-AC bands: predicted from pre-lapped neighbor coefficients
- High-AC bands: lightweight learned linear predictors (R-D gated — skip when
  predictor residual exceeds original)
- Eliminates the prediction-transform-quantize-inverse-filter pipeline

### 5. OBMC Inter Prediction

Overlapped Block Motion Compensation blends predictions from neighboring blocks'
motion vectors, smoothing prediction residuals at block boundaries. This works
synergistically with lapped transforms — the overlap regions align.

- Quad-tree partition: 128x128 down to 4x4
- Blending weights sum to exactly 64 (6-bit fixed-point) at every pixel
- Checkerboard processing: even superblocks first, odd second
- Affine motion at superblock level with per-4x4 derived MVs
- 1/8-pel sub-pixel precision: 8-tap luma, 4-tap chroma
- DMVR (Decoder-side Motion Vector Refinement): optical-flow-based, zero-bit cost

### 6. Perceptual Optimization (Encoder-Side)

The encoder uses perceptual models to allocate bits where they matter:

- **PTD (Perceptual Transform Distortion)**: CSF-weighted transform-domain distortion
  in the RDO inner loop. Exploits reduced sensitivity to high spatial frequencies
- **Activity masking**: per-CTU QP offset based on local AC variance.
  `delta_QP = 0.5 * log2(activity / mean_activity)`, clamped to [-6, +6]
- **Butteraugli outer loop**: constant-perceptual-quality mode via QP bisection
- **Chroma QP offsets**: +4 to +6 at picture/slice level (HVS is ~4x less sensitive
  to chroma error)
- **VMAF-proxy lambda modulation**: adjusts RDO lambda for perceptual gains

### 7. Film Grain Synthesis

Rather than wasting bits encoding noise, Sinter strips film grain before encoding
and resynthesizes it at the decoder:

- AR(3) luma / AR(4) chroma model (4th chroma tap = co-located luma grain)
- Grain amplitude scales with local luma (photon noise model)
- Parameters signaled per-frame in bitstream header
- Deterministic seed derived from frame number (seeking works correctly)

## Fixed-Point Precision

All decoder-path arithmetic is integer-only for cross-platform bit-exactness:

- Filter coefficients: Q15 format (16-bit fractional)
- Accumulators: 32-bit for 8-bit input, 64-bit for 10-bit
- Rounding: `(x + (1 << (shift-1))) >> shift` (round-half-up), defined once globally
- Data types: `int16_t` for 8-bit coefficients, `int32_t` for 10-bit+

## Parallelism

- **Tile grid**: independent entropy state per tile, full parallel encode/decode
- **Tile overlap**: 4-pixel redundant decode zones with raised cosine blending
- **WPP**: wavefront parallel processing within tiles
- **Checkerboard OBMC**: two independent groups processed in parallel

## Project Structure

```
libsinter/
  core/
    transform/      -- fwd/inv DCT, lapping (lifting steps)
    prediction/     -- freq-domain intra, inter/OBMC, affine
    entropy/        -- rANS core, probability models, CDFs
    quantize/       -- PVQ gain-shape, DC scalar quant, QP mapping
    grain/          -- film grain synthesis (decoder-side)
    bitstream/      -- writer, reader, syntax elements
  codec/
    encoder.c       -- encoder pipeline
    decoder.c       -- decoder pipeline
    rdo.c           -- rate-distortion optimization
  perceptual/
    csf_weights.c   -- CSF weight table
    activity_mask.c -- per-CTU activity and QP offset
    ptd.c           -- Perceptual Transform Distortion
    grain_analysis.c-- film grain estimation (encoder-side)
  util/
    fixed_point.h   -- fixed-point math, rounding macros
    memory.h        -- aligned allocation
    platform.h      -- SIMD detection
  test/
    test_transform.c
    test_entropy.c
    test_pvq.c
    test_obmc.c
    test_roundtrip.c
    test_perceptual.c
```

## Naming Conventions

- `snt_` shared/normative, `snte_` encoder-only, `sntd_` decoder-only
- `psy_` perceptual functions
- `obmc_`, `rans_`, `pvq_` module prefixes
- `_fwd`/`_inv` for transforms, `_encode`/`_decode` for codec pairs
- `_ref` suffix for scalar reference implementations
- Fixed-point precision prefixes: `q15_`, `q8_`

## Testing Strategy

Tests are gated — each gate must pass before proceeding:

| Gate | What | Pass Criteria |
|------|------|---------------|
| 0 | Lossless round-trip | QP=0 output bit-identical to input |
| 1 | Component inverses | Transform, rANS, PVQ energy, OBMC weights |
| 2 | Encoder-decoder match | Internal recon buffer == decoded output |
| 3 | Reference signals | Flat gray, impulse, gradients, checkerboard |
| 4 | Video sequences | Motion, scene cuts, standard test media |
| 5 | Perceptual sanity | CSF, masking, grain, Butteraugli targeting |

## Performance Targets

- ~40-50% bitrate reduction over H.264 at equivalent perceptual quality
- ~15-25% bitrate reduction over H.265 at equivalent perceptual quality
- Speed preset 6: real-time 1080p encoding on modern hardware

## Thought Process

The design emerged from a roundtable of perspectives spanning open-codec development
(VP8/VP9/AV1), perceptual coding research (Daala, HEVC), and open-source ecosystem
building (Xiph/Ogg/Vorbis/Opus). The key convergence:

1. **Lapped transforms** are the single highest-leverage architectural change. They
   eliminate an entire class of artifacts and an entire category of post-processing
   complexity. Daala proved they work; Sinter ships them.

2. **PVQ** is the perceptual quantization strategy that conventional codecs lack.
   Scalar quantization optimizes for PSNR; PVQ optimizes for texture preservation.
   Proven in Opus audio, validated in Daala video research.

3. **rANS** matches arithmetic coding compression with better decode throughput.
   Multi-symbol coding reduces operations per block vs. binary CABAC.

4. **Frequency-domain intra prediction** breaks the serial dependency chain that
   makes conventional intra coding slow, while being architecturally distinct.

5. **Perceptual-first encoder design** puts bits where the eye looks and takes them
   from where it doesn't notice. Not PSNR-optimized with perceptual bolted on —
   perceptual from the ground up.

The codec composes well-understood building blocks in a novel architecture. Each
component has decades of published research behind it. The novelty is the combination
and the perceptual-first philosophy.
