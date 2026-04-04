# Sinter Video Codec

A novel video codec built on lapped transforms, perceptual vector quantization, and asymmetric numeral systems. Eliminates blocking artifacts structurally — no deblocking filter needed.

## Results vs H.264

Tested on ffmpeg testsrc (256x256 color bars, 30 frames at 30fps):

| Plane | Sinter (QP=4) | H.264 (QP=20) | Gap |
|-------|---------------|----------------|-----|
| Luma Y | 44.80 dB | 45.58 dB | 0.77 dB |
| Chroma U | 47.18 dB | 53.91 dB | 6.73 dB |
| Chroma V | 47.06 dB | 53.47 dB | 6.41 dB |
| File size | 838 KB | 149 KB | 5.6x |

Luma is within 0.77 dB of H.264. Video is visually comparable on color content.

**v0.4 update:** rANS entropy coding is now fully integrated into the encoder/decoder pipeline. The adaptive rANS coder with 6 CDFs and PARA adaptation replaces the exp-Golomb VLC, providing ~3x compression improvement over the previous bitstream format. Quality is stable across all frame sizes (40 dB at QP=4 from 8x8 to 256x256).

**v0.5 update:** Inter-frame prediction (P-frames) with integer-pel motion search, median MV prediction, and skip mode. Inter frames are 2-6x smaller than keyframes depending on content motion. 10 adaptive CDFs total (6 intra + 4 inter). Lapping bypassed for inter frames (plain block DCT on residual). All patent-safe (H.263-era techniques only).

## Quick Start

```bash
# Build library and tests
make && make test

# Build CLI tools
gcc -O2 -I libsinter libsinter/sinter.c libsinter/core/transform/*.c \
  libsinter/core/quantize/pvq.c libsinter/core/quantize/qp_table.c \
  libsinter/core/quantize/pvq_shape_rans.c libsinter/core/entropy/rans.c \
  libsinter/core/bitstream/bitstream.c libsinter/core/prediction/*.c \
  libsinter/core/partition/quadtree.c libsinter/core/grain/grain_synth.c \
  libsinter/core/tile/tile.c libsinter/perceptual/*.c libsinter/codec/*.c \
  tools/sinterenc.c -lm -o build/sinterenc

gcc -O2 -I libsinter libsinter/sinter.c libsinter/core/transform/*.c \
  libsinter/core/quantize/pvq.c libsinter/core/quantize/qp_table.c \
  libsinter/core/quantize/pvq_shape_rans.c libsinter/core/entropy/rans.c \
  libsinter/core/bitstream/bitstream.c libsinter/core/prediction/*.c \
  libsinter/core/partition/quadtree.c libsinter/core/grain/grain_synth.c \
  libsinter/core/tile/tile.c libsinter/perceptual/*.c libsinter/codec/*.c \
  tools/sinterdec.c -lm -o build/sinterdec
```

## Example: Encode, Decode, and Compare vs H.264

```bash
# Generate a 1-second test video
ffmpeg -f lavfi -i "testsrc=duration=1:size=256x256:rate=30" \
  -pix_fmt yuv420p -f rawvideo test_input.yuv

# Encode with Sinter (QP=4 is highest quality)
./build/sinterenc -i test_input.yuv -w 256 -h 256 -q 4 -n 30 -o test.sntr

# Decode
./build/sinterdec -i test.sntr -o test_decoded.yuv

# Convert to playable MP4
ffmpeg -f rawvideo -pix_fmt yuv420p -s 256x256 -r 30 \
  -i test_decoded.yuv -c:v libx264 -crf 0 sinter_result.mp4

# Encode same source with H.264 for comparison
ffmpeg -f rawvideo -pix_fmt yuv420p -s 256x256 -r 30 \
  -i test_input.yuv -c:v libx264 -qp 20 -g 1 h264_result.mp4

# Watch both side by side
open sinter_result.mp4
open h264_result.mp4
```

## Example: Encode Your Own Video

```bash
# Convert any video to raw YUV
ffmpeg -i my_video.mp4 -f rawvideo -pix_fmt yuv420p -s 256x256 -t 2 input.yuv

# Count frames (file_size / (256*256*1.5))
FRAMES=$(python3 -c "import os; print(os.path.getsize('input.yuv') // (256*256*3//2))")

# Encode (QP: 4=best quality, 16=balanced, 32=small file)
./build/sinterenc -i input.yuv -w 256 -h 256 -q 8 -n $FRAMES -o encoded.sntr

# Decode and play
./build/sinterdec -i encoded.sntr -o decoded.yuv
ffmpeg -f rawvideo -pix_fmt yuv420p -s 256x256 -r 30 \
  -i decoded.yuv -c:v libx264 -crf 0 output.mp4
open output.mp4
```

## Example: Measure PSNR

```bash
# After encoding and decoding, measure quality per plane
python3 -c "
import math
W, H = 256, 256
fy, fuv = W*H, W*H//4

with open('test_input.yuv', 'rb') as f:
    oy, ou, ov = list(f.read(fy)), list(f.read(fuv)), list(f.read(fuv))
with open('test_decoded.yuv', 'rb') as f:
    dy, du, dv = list(f.read(fy)), list(f.read(fuv)), list(f.read(fuv))

for name, o, d in [('Y', oy, dy), ('U', ou, du), ('V', ov, dv)]:
    mse = sum((a-b)**2 for a,b in zip(o,d)) / len(o)
    psnr = 10*math.log10(255**2/mse) if mse > 0 else 99.99
    print(f'{name}: {psnr:.2f} dB')
"
```

## Architecture

```
Luma Pipeline (8x8 blocks):
  Input → Lapped Pre-Filter (20% Malvar lifting) → 8x8 DCT (Q14)
  → DC Prediction (left+above) → AC Prediction (4 lowest coefficients)
  → Per-band mode decision: PVQ (texture) or Scalar (smooth)
  → 16 bands × 4 coefficients, energy-concentration threshold
  → rANS Entropy Coding (6 adaptive CDFs, dual interleaved) → Bitstream

Chroma Pipeline (4x4 blocks, YUV420 half-resolution):
  Input → 4x4 DCT (Q14) → DC Prediction → Scalar Quantization (all 16 coeffs)
  → rANS Entropy Coding (shared CDFs with luma) → Bitstream

Decoder reverses each step identically.
```

### What Makes Sinter Different

| Feature | H.264/H.265/AV1 | Sinter |
|---------|-----------------|--------|
| Blocking | DCT + deblocking filter | Lapped transforms (no filter needed) |
| Quantization | Per-coefficient scalar | Hybrid: PVQ (texture) or scalar (smooth) |
| Entropy coding | CABAC/CAVLC | rANS with PARA-adaptive CDFs |
| Prediction | Pixel-domain angular modes | Frequency-domain coefficient prediction |
| Patents | Multiple patent pools | All techniques public domain |

### Key Innovations

1. **Hybrid PVQ/Scalar per-band** — the encoder measures energy concentration per frequency band. If one coefficient dominates (>60% of band energy), use scalar quantization (clean gradients). If energy is spread, use PVQ gain-shape (preserve texture). 1 bit per band. This was the single biggest quality improvement: +3.49 dB.

2. **16-band PVQ with 4 dimensions each** — splitting 63 AC coefficients into 16 small bands gives max K=225 per band (vs K=14 monolithic). The 16x improvement in shape resolution was worth +10 dB.

3. **Perceptual skip threshold** — PVQ bands with marginal energy are zeroed rather than coded, preventing PVQ from injecting texture noise into smooth content.

4. **CSF-weighted PVQ** — each band receives a Contrast Sensitivity Function weight that modulates the K parameter. High-frequency bands get fewer pulses (perceptually invisible detail), low-frequency bands get more (perceptually critical detail).

## Honest Assessment

**Strengths:**
- Near-H.264 luma quality (0.77 dB gap)
- No blocking artifacts by design
- Patent-free architecture
- Hybrid PVQ/scalar preserves texture while keeping gradients clean

**Weaknesses:**
- File sizes still larger than H.264 (inter + rANS narrows gap significantly)
- Chroma quality 6-7 dB behind H.264, chroma inter not yet implemented
- No sub-pel motion refinement (integer-pel only)
- Fixed block sizes (no quad-tree in pipeline)
- Reference encoder is slow (~100x vs x264)

**Path forward:**
- Sub-pel motion refinement (half-pel bilinear, ~1 dB gain)
- Chroma inter prediction (20-30% overall size reduction)
- Quad-tree partitioning (better content adaptation)
- SIMD optimization for real-time encoding

## Documentation

- **[Technical Reference](docs/TECHNICAL_REFERENCE.md)** — Full specification
- **[Scoreboard](SCOREBOARD.md)** — Quality improvement tracking
- **[Public API](libsinter/sinter.h)** — Stable C API (8 functions)

## Requirements

- C11 compiler (GCC or Clang)
- make
- ffmpeg (for testing and format conversion)
- No external library dependencies
