# Sinter Video Codec

Sinter is an experimental video codec built on lapped transforms, perceptual vector quantization, and asymmetric numeral systems. It eliminates blocking artifacts structurally (no deblocking filter), preserves texture energy through gain-shape quantization, and decodes faster than arithmetic-coded codecs.

**Quality:** 37-40 dB PSNR at QP=4 on test content (competitive with H.264 baseline).

**Architecture:** Fundamentally different from the H.261 lineage (H.264/H.265/AV1/VVC). No loop filters, no scalar quantization, no pixel-domain angular prediction.

## Quick Start

```bash
make              # build libsinter.a
make test         # run 32 tests (transform, PVQ, entropy, round-trip)

# Encode a raw YUV file
./build/sinterenc -i input.yuv -w 1920 -h 1080 -q 20 -n 1 -o output.sntr

# Decode back to YUV
./build/sinterdec -i output.sntr -o decoded.yuv
```

## How It Works

```
Encoder:  Pixels -> Lapped Pre-Filter -> 8x8 DCT -> DC Prediction
          -> 16-Band PVQ (gain-shape) -> Bitstream

Decoder:  Bitstream -> PVQ Decode -> DC Reconstruct -> Inverse DCT
          -> Lapped Post-Filter -> Film Grain Synthesis -> Pixels
```

Three architectural breaks from traditional codecs:

1. **Lapped transforms** overlap adjacent blocks via lifting-step pre/post filters, structurally preventing blocking artifacts. No deblocking filter needed.
2. **Perceptual Vector Quantization** decomposes each frequency band into energy (gain) and direction (shape). Preserves texture where scalar quantization destroys it.
3. **rANS entropy coding** matches arithmetic coding compression at ~2x decode speed. Adaptive CDFs with PARA. (Implemented and tested; pipeline integration in progress.)

## Documentation

- **[Technical Reference](docs/TECHNICAL_REFERENCE.md)** — Complete specification: bitstream format, transform math, PVQ design, perceptual pipeline, entropy coding, container format
- **[Public API](libsinter/sinter.h)** — 8-function C API with opaque handles

## Project Structure

```
libsinter/
  sinter.h / sinter.c     Public C API
  core/transform/          Lapped filters (lifting) + Q14 integer DCT (4-32pt)
  core/quantize/           PVQ gain-shape + DC scalar quant + QP table
  core/entropy/            Dual interleaved rANS + adaptive CDFs
  core/bitstream/          Frame header + block data syntax
  core/prediction/         Freq-domain intra + OBMC inter + motion estimation
  core/partition/          Quad-tree block partitioning with RDO
  core/grain/              Film grain AR(3)/AR(4) synthesis
  core/tile/               Multi-tile with overlap blending
  codec/                   Encoder + decoder pipelines + rate control
  perceptual/              CSF weights + PTD metric + activity masking + Butteraugli
  test/                    4 test suites, 32 tests
tools/
  sinterenc.c              CLI encoder (YUV -> .sntr)
  sinterdec.c              CLI decoder (.sntr -> YUV)
ffmpeg/
  libsinter_enc.c          FFmpeg encoder wrapper (reference)
  libsinter_dec.c          FFmpeg decoder wrapper (reference)
wasm/
  build.sh                 Emscripten WASM decoder build
  sinter_player.html       Browser playback demo
```

## Quality Results

256x256 diagonal gradient, all-intra, single frame:

| QP | PSNR | Max Error | Bitrate | Compression |
|----|------|-----------|---------|-------------|
| 4  | 39.91 dB | 10 | 9.47 bpp | 1.3x |
| 8  | 36.06 dB | 16 | 8.04 bpp | 1.5x |
| 12 | 31.18 dB | 23 | 6.44 bpp | 1.9x |
| 16 | 27.35 dB | 34 | 5.68 bpp | 2.1x |
| 20 | 24.07 dB | 48 | 4.67 bpp | 2.6x |
| 28 | 19.17 dB | 87 | 2.75 bpp | 4.4x |

Bitrates are high because the shipped pipeline uses raw bit packing (no entropy coding). With rANS, expect ~2x bitrate reduction at the same quality.

## Key Design Decisions

- **60% Malvar lapping coefficients:** Full-strength Malvar (100%) maximizes decorrelation but the Q14 integer DCT amplifies errors through the post-filter chain. 60% balances deblocking quality with reconstruction precision.
- **16-band PVQ (4 dims each):** Smaller bands give dramatically higher K values (max K=225 at n=4 vs K=14 at n=63), improving shape resolution by 10+ dB.
- **Midpoint dequantization:** Reconstructs at the center of each quantization bin, not the lower edge. A textbook fix that every production codec uses.
- **DC prediction:** Average of left and above reconstructed DC coefficients. Reduces DC residual 10-20x on smooth content.

## Current Limitations

- All-intra only (inter prediction module exists but is not wired in)
- Fixed 8x8 block grid (quad-tree exists but not integrated)
- Raw bit packing (rANS validated but not in pipeline — prior integration had bitstream ordering bugs)
- Chroma is DC-only at 4x4 blocks
- 8-bit only (10-bit types defined, pipeline not extended)
- No rate control in pipeline (CRF module exists)

## Ecosystem

- **FFmpeg:** Wrapper files in `ffmpeg/` following libaom pattern. Requires bitstream freeze for upstream submission.
- **WASM:** `wasm/build.sh` compiles decoder to WebAssembly via Emscripten. `sinter_player.html` demonstrates browser playback.
- **Container:** `.sntr` IVF-style format with magic `SNTR`, frame-length-prefixed.
- **pkg-config:** `libsinter.pc.in` template for build system integration.

## Building from Source

Requirements: C11 compiler (GCC or Clang), make, no external dependencies.

```bash
make                    # static library + archive
make test               # all 32 tests
make clean              # remove build artifacts
```

CLI tools are built separately:
```bash
gcc -O2 -I libsinter libsinter/*.c libsinter/**/*.c tools/sinterenc.c -lm -o sinterenc
gcc -O2 -I libsinter libsinter/*.c libsinter/**/*.c tools/sinterdec.c -lm -o sinterdec
```
