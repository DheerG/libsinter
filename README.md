# Sinter Video Codec

An experimental, vibecoded video codec built from scratch using [Claude Code](https://claude.ai/code) agent teams. It started as an exploration of lapped transforms, perceptual vector quantization (PVQ), and rANS entropy coding as a patent-free alternative to the H.264/H.265 lineage.

**This was a learning experiment**, not a production codec. The goal was to test one-shot agent team workflows on a domain I had zero prior experience in — and to see how far a simulated expert team could push a novel architecture. Write-up: [One-Shot Wonder](https://dheer.co/one-shot-wonder/) | [Claude Agent Teams](https://dheer.co/claude-agent-teams/)

Later phases then fixed bugs and improved compression, measuring every change against x264 and keeping only the changes that helped.

## Results

On real camera footage, Sinter needs **about 1.4x the bitrate of H.264 (x264)** for the same luma quality. The original experiment needed 10.5–17.6x.

| Clip | Original | After bug fixes | Before B-frames etc. | Now |
|------|----------|-----------------|----------------------|-----|
| Foreman | 14.5x | 8.3x | 2.0x | **1.44x** |
| Mobile & Calendar | 10.5x | 8.1x | 2.1x | **1.38x** |
| Akiyo | 17.6x | 7.7x | 1.9x | **1.45x** |
| Bus (held out) | – | – | 1.76x | **1.43x** |
| News (held out) | – | – | 1.77x | **1.37x** |

Each figure is the bitrate Sinter needs to match x264's luma PSNR (BD-rate), on the Xiph CIF test sequences (352x288, first 150 frames). x264 runs at `-preset medium -tune psnr` (which includes its adaptive B-frames) with the same 30-frame keyframe interval. Bus and News were never used for tuning; they were measured once, at the end. Counting chroma (YUV PSNR) gives 1.35–1.57x; SSIM gives 1.41–1.47x.

The original write-up reported 18.6x from a two-point comparison on a synthetic 256x256 clip. [SCOREBOARD.md](SCOREBOARD.md) and [docs/TECHNICAL_REFERENCE.md](docs/TECHNICAL_REFERENCE.md) record that original design and have not been updated for the changes below.

## What's in It

~6,800 lines of C in the library:

- **16x16 macroblocks**: in P- and B-frames, each is skipped, coded with one motion vector, or split into four 8x8 blocks with their own modes
- **8x8 integer DCT**, with AC coefficients grouped into 16 bands in zigzag order
- **Hybrid PVQ/scalar quantization**, chosen per band by rate-distortion cost: PVQ gain-shape where it pays, scalar levels elsewhere
- **Intra prediction**: transform-domain DC/low-AC prediction, plus vertical, horizontal and average spatial predictors
- **Inter prediction**: quarter-pel motion compensation (6-tap Lanczos-3 luma, bilinear chroma), rate-aware motion search, median MV prediction, and skip blocks that follow the predicted motion
- **B-frames**: two between anchor frames by default, each block predicted forward, backward or from the average of both, in closed GOPs
- **Mode decisions by trial encoding**: each candidate (skip, inter, intra, 16x16 or 8x8) is coded without emitting it, and the lowest distortion + lambda x bits wins
- **In-loop deblocking filter** on luma and chroma block edges
- **Dual-interleaved rANS** with adaptive contexts per band, block type and plane, count-based adaptation, and cross-frame carry
- **CLI tools** with Y4M I/O (`sinterenc` / `sinterdec`); `sinterdec` reorders B-frames to display order
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

`sinterenc` takes `-k` for the keyframe interval (default 30) and `-b` for the number of B-frames between anchors (default 2). The bitstream has changed since the original experiment: `.sntr` files are now container version 3, and files from earlier versions do not decode.

## Where the Gap Came From

The original "why the gap exists" estimates were the expert team's guesses. Measured instead:

- With x264 limited to the original Sinter's toolset (no B-frames, one reference frame, full-pel motion, CAVLC, no deblocking), x264 needed 2.2–3.5x its normal bitrate. The original Sinter needed another 2.1–3.5x on top of that. So about half of the original gap was missing tools, and half was bugs and implementation.
- The largest single wins were: choosing PVQ or scalar per band by rate-distortion cost (41%), quarter-pel motion (30%), B-frames (14%), skip blocks that follow the predicted motion (10%), mode decisions by trial encoding (8%), 16x16 macroblocks (5%) and the sharper 6-tap interpolation filter (5%). Each is measured against the step before it.

What x264 still has that Sinter doesn't: multiple reference frames, adaptive B-frame placement, more partition shapes (16x8, 8x16, 4x4) and transform sizes, weighted prediction, and much richer entropy contexts (CABAC).

## Known Issues

- The public C API (`sinter.h` / `sinter.c`) is not built by the Makefile, and `snt_decoder_decode` passes an uninitialized stride. It also has no B-frame support. The CLI tools use the internal API and are unaffected.
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

The later phases (bug fixes and compression work) were also done with Claude Code, without personas: each change was measured by BD-rate against x264 on the Xiph clips and kept only if it helped. The last phase ran autonomously toward a fixed target (at most 1.5x x264 on each tuning clip). The git history records each change, its measured effect, and the changes tried and dropped.

## Requirements

- C11 compiler (GCC or Clang)
- make
- ffmpeg (for Y4M generation)
- No external dependencies
