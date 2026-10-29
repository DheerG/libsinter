# Sinter Video Codec

An experimental, vibecoded video codec built from scratch using [Claude Code](https://claude.ai/code) agent teams. It started as an exploration of lapped transforms, perceptual vector quantization (PVQ), and rANS entropy coding as a patent-free alternative to the H.264/H.265 lineage.

**This was a learning experiment**, not a production codec. The goal was to test one-shot agent team workflows on a domain I had zero prior experience in — and to see how far a simulated expert team could push a novel architecture. Write-up: [One-Shot Wonder](https://dheer.co/one-shot-wonder/) | [Claude Agent Teams](https://dheer.co/claude-agent-teams/)

A second phase then fixed bugs and improved compression, measuring every change against x264 and keeping only the changes that helped.

## Results

On real camera footage, Sinter needs **about 2x the bitrate of H.264 (x264)** for the same luma quality. The original experiment needed 10.5–17.6x.

| Clip | Original | After bug fixes | Now |
|------|----------|-----------------|-----|
| Foreman | 14.5x | 8.3x | **2.0x** |
| Mobile & Calendar | 10.5x | 8.1x | **2.1x** |
| Akiyo | 17.6x | 7.7x | **1.9x** |

Each figure is the bitrate Sinter needs to match x264's luma PSNR (BD-rate), on the Xiph CIF test sequences (352x288, first 150 frames). x264 runs at `-preset medium -tune psnr` with the same 30-frame keyframe interval. Counting chroma (YUV PSNR) gives 2.1x on all three clips; SSIM gives 2.0–2.5x.

The original write-up reported 18.6x from a two-point comparison on a synthetic 256x256 clip. [SCOREBOARD.md](SCOREBOARD.md) and [docs/TECHNICAL_REFERENCE.md](docs/TECHNICAL_REFERENCE.md) record that original design and have not been updated for the changes below.

## What's in It

~6,400 lines of C in the library:

- **8x8 integer DCT**, with AC coefficients grouped into 16 bands in zigzag order
- **Hybrid PVQ/scalar quantization**, chosen per band by rate-distortion cost: PVQ gain-shape where it pays, scalar levels elsewhere
- **Intra prediction**: transform-domain DC/low-AC prediction, plus vertical, horizontal and average spatial predictors
- **Inter prediction** (P-frames): quarter-pel motion compensation (4-tap Catmull-Rom luma, bilinear chroma), rate-aware motion search, median MV prediction, and skip blocks that follow the predicted motion
- **In-loop deblocking filter** on luma and chroma block edges
- **Dual-interleaved rANS** with adaptive contexts per band, block type and plane, count-based adaptation, and cross-frame carry
- **CLI tools** with Y4M I/O (`sinterenc` / `sinterdec`)
- 32 passing tests; encoder and decoder reconstructions match bit-exactly on the test clips

The lapped transform from the original design is still in the tree (`libsinter/core/transform/lapping.c`) but not used: its lifting network was not DC-preserving and added a grid artifact at every block edge.

## Quick Start

```bash
make && make test

# Build CLI
gcc -O2 -I libsinter tools/sinterenc.c -Lbuild -lsinter -lm -o build/sinterenc
gcc -O2 -I libsinter tools/sinterdec.c -Lbuild -lsinter -lm -o build/sinterdec

# Encode/decode
ffmpeg -f lavfi -i "testsrc=duration=1:size=64x64:rate=10" -pix_fmt yuv420p -f yuv4mpegpipe input.y4m -y
./build/sinterenc -i input.y4m -o encoded.sntr -q 12
./build/sinterdec -i encoded.sntr -o decoded.y4m
```

The bitstream has changed since the original experiment, so `.sntr` files from earlier versions do not decode.

## Where the Gap Came From

The original "why the gap exists" estimates were the expert team's guesses. Measured instead:

- With x264 limited to the original Sinter's toolset (no B-frames, one reference frame, full-pel motion, CAVLC, no deblocking), x264 needed 2.2–3.5x its normal bitrate. The original Sinter needed another 2.1–3.5x on top of that. So about half of the original gap was missing tools, and half was bugs and implementation.
- The largest single wins were: choosing PVQ or scalar per band by rate-distortion cost (41%), quarter-pel motion (30%), and skip blocks that follow the predicted motion (10%). Each is measured against the step before it.

What x264 still has that Sinter doesn't: B-frames, variable block sizes, multiple reference frames, mode decisions by trial encoding, and much richer entropy contexts (CABAC).

## Known Issues

- The public C API (`sinter.h` / `sinter.c`) is not built by the Makefile, and `snt_decoder_decode` passes an uninitialized stride. The CLI tools use the internal API and are unaffected.
- The FFmpeg wrapper (`ffmpeg/`) and the WASM build (`wasm/`) depend on that API and are untested.
- The decoder is not hardened against malformed input: a truncated frame can hang it, and a frame larger than the container header says can overflow the output buffer.
- Inter frames assume the luma stride equals the frame width, which is true for the CLI tools.

## Credits

Built with [Claude Code](https://claude.ai/code). The expert team was assembled by asking Claude to name five experts for building an open-source codec competing with H.264/H.265:

- **Jim Bankoski** — Google architect behind VP8 and VP9; open-codec veteran
- **Timothy Terriberry** — Mozilla/Xiph engineer; drove Daala and AV1 research
- **Monty Montgomery** — Xiph.org founder; built Ogg, Vorbis, and Theora ecosystems
- **Gary Sullivan** — Microsoft; co-chaired the H.264 and HEVC standardization efforts
- **Jens-Rainer Ohm** — HEVC co-chair; deep expertise in perceptual coding theory

Plus a principal engineer (facilitator) and patent reviewer. All personas simulated by Claude.

The second phase (bug fixes and compression work) was also done with Claude Code, without personas: each change was measured by BD-rate against x264 on the Xiph clips and kept only if it helped. The git history records each change and its measured effect.

## Requirements

- C11 compiler (GCC or Clang)
- make
- ffmpeg (for Y4M generation)
- No external dependencies
