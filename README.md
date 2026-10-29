# Sinter Video Codec

An experimental, vibecoded video codec built from scratch using [Claude Code](https://claude.ai/code) agent teams. It started as an exploration of lapped transforms, perceptual vector quantization (PVQ), and rANS entropy coding as a patent-free alternative to the H.264/H.265 lineage.

**This was a learning experiment**, not a production codec. The goal was to test one-shot agent team workflows on a domain I had zero prior experience in — and to see how far a simulated expert team could push a novel architecture. Write-up: [One-Shot Wonder](https://dheer.co/one-shot-wonder/) | [Claude Agent Teams](https://dheer.co/claude-agent-teams/)

In October 2026 the codec was reworked with Claude Code: bugs fixed, compression tools added, and the format rebuilt so that every tool the decoder runs comes from a publication before 2004, or is new and was searched for in patents first (see [Where the tools come from](#where-the-tools-come-from)). Every change was measured against x264.

## Results

On real camera footage, Sinter now needs **about 1.2–1.3x the bitrate of H.264 (x264)** for the same luma quality. The original experiment needed 10.5–17.6x.

| Clip | Original experiment | October 2026 |
|------|---------------------|--------------|
| Foreman | 14.5x | **1.24x** |
| Mobile & Calendar | 10.5x | **1.28x** |
| Akiyo | 17.6x | **1.21x** |
| Bus (held out) | – | **1.32x** |
| News (held out) | – | **1.17x** |

Each figure is the bitrate Sinter needs to match x264's luma PSNR (BD-rate), on the Xiph CIF test sequences (352x288, first 150 frames). x264 runs at `-preset medium -tune psnr` (which includes its adaptive B-frames) with the same 30-frame keyframe interval. Bus and News were never used for tuning. Counting chroma (YUV PSNR) gives 1.26–1.37x; SSIM gives 1.21–1.33x.

The original write-up reported 18.6x from a two-point comparison on a synthetic 256x256 clip. [SCOREBOARD.md](SCOREBOARD.md) and [docs/TECHNICAL_REFERENCE.md](docs/TECHNICAL_REFERENCE.md) record that original design and have not been updated.

## Examples

Left to right: the source, Sinter, and x264 at the same file size. Frame 44 (a B-frame), cropped and zoomed 3x. Here x264 uses its default tuning, 2-pass, at Sinter's bitrate. Each clip also has a side-by-side video of the full 5 seconds, in the same order.

**Foreman**: Sinter 265 kbps, 35.0 dB; x264 264 kbps, 35.4 dB.

![Foreman, frame 44: source, Sinter, x264](docs/images/foreman.png)

[Watch the side-by-side video](docs/videos/foreman.mp4)

**Mobile & Calendar**: Sinter 893 kbps, 31.0 dB; x264 888 kbps, 32.0 dB.

![Mobile & Calendar, frame 44: source, Sinter, x264](docs/images/mobile.png)

[Watch the side-by-side video](docs/videos/mobile.mp4)

**Bus** (held out): Sinter 821 kbps, 32.3 dB; x264 822 kbps, 33.6 dB.

![Bus, frame 44: source, Sinter, x264](docs/images/bus.png)

[Watch the side-by-side video](docs/videos/bus.mp4)

## What's in It

- **16x16 macroblocks**: in P- and B-frames, each is skipped, coded with one motion vector, or split into four 8x8 blocks with their own modes
- **Skipped blocks** repeat the motion of the block to their left (MPEG-1's rule) or follow a **global motion vector** sent per frame (MPEG-4 Part 2's GMC), which the encoder estimates by Lucas-Kanade
- **8x8 integer DCT**, with AC coefficients in zigzag order grouped into 16 bands, each with a coded flag
- **Intra prediction** of the DC and first coefficient row or column from the neighbours (H.263 Annex I)
- **Inter prediction**: quarter-pel motion compensation through a cubic B-spline of each reference frame (luma) and 1/8-pel bilinear (chroma), rate-aware motion search, H.263's median vector prediction
- **B-frames**: two between anchor frames by default, each block predicted forward, backward or from the average of both, in closed GOPs
- **Mode decisions by trial encoding**: each candidate (skip, inter, intra, 16x16 or 8x8) is coded without emitting it, and the lowest distortion + lambda x bits wins
- **In-loop deblocking filter** (H.263 Annex J) on luma and chroma block edges
- **Range coder** with adaptive frequency counts, starting from tables trained on other Xiph clips; levels are coded as binary decisions as in JPEG's arithmetic mode
- **Two new context tools**: the residual's models are chosen from the DCT of its prediction, and the mode flags' models from how well two decoder-side predictions agree ([docs/NOVEL.md](docs/NOVEL.md))
- **CLI tools** with Y4M I/O (`sinterenc` / `sinterdec`); `sinterdec` reorders B-frames to display order
- 32 passing tests; encoder and decoder reconstructions match bit-exactly on every test clip
- About 200 ms to encode and 3–6 ms to decode a CIF frame on a laptop

The lapped transform, PVQ, rANS and the film grain synthesis of the original design are still in the tree, but none of them is part of the format: see [PROVENANCE.md](PROVENANCE.md).

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

`sinterenc` takes `-k` for the keyframe interval (default 30) and `-b` for the number of B-frames between anchors (default 2). The bitstream has changed since the original experiment: `.sntr` files are now container version 4, and files from earlier versions do not decode.

## Where the tools come from

[PROVENANCE.md](PROVENANCE.md) lists every tool the decoder runs, with its source, date and the changes made to it. All but two come from publications between 1974 and 2000: the DCT (1974), Martin's range coder (1979), Witten, Neal & Cleary's adaptive models (1987), Unser's B-splines (1991), JPEG (1992), MPEG-1 (1993), H.263 and its annexes (1996–1998) and MPEG-4 Part 2 (1999–2000). Nothing comes from H.264 or later codecs (VP8 onward, Daala, AV1).

The two new tools are described in [docs/NOVEL.md](docs/NOVEL.md), with the patent searches run before keeping them. Two other new ideas were left out because the searches found them: online logistic mixing of several probability models (Google, EP3677027 / US10887595), which was still built behind an off-by-default compile switch (`-DSNT_EXPERIMENT_MIXING`) and measured at 0.5%, and per-clip probability tables sent as changes from defaults (VP8's frame-header updates), also built behind a switch (`-DSNT_EXPERIMENT_CLIPTABLES`) and measured at 0.1%.

These records are a map of where each idea comes from, not a patent clearance.

The entropy models' starting tables are trained on ten Xiph CIF clips that are never used for evaluation; [tools/train/README.md](tools/train/README.md) shows how to regenerate them.

## Where the Gap Came From

The original "why the gap exists" estimates were the expert team's guesses. Measured instead:

- With x264 limited to the original Sinter's toolset (no B-frames, one reference frame, full-pel motion, CAVLC, no deblocking), x264 needed 2.2–3.5x its normal bitrate. The original Sinter needed another 2.1–3.5x on top of that. So about half of the original gap was missing tools, and half was bugs and implementation.
- The largest single gains in October 2026 came from quarter-pel motion, B-frames, skipped blocks that follow the motion, and mode decisions by trial encoding. The git history gives every change with its measured effect.

What x264 still has that Sinter doesn't: multiple reference frames, adaptive B-frame placement, more partition shapes (16x8, 8x16, 4x4) and transform sizes, weighted prediction, and look-ahead rate control (MB-tree).

## Known Issues

- The public C API (`sinter.h` / `sinter.c`) is not built by the Makefile, and `snt_decoder_decode` passes an uninitialized stride. It also has no B-frame support. The CLI tools use the internal API and are unaffected.
- The FFmpeg wrapper (`ffmpeg/`) and the WASM build (`wasm/`) depend on that API and are untested.
- The decoder is not hardened against malformed input: a damaged frame can decode to garbage, and a frame larger than the container header says can overflow the output buffer.
- Inter frames assume the luma stride equals the frame width, which is true for the CLI tools.
- Changing the entropy model set requires regenerating `libsinter/codec/model_tables.h` (a static assertion catches a mismatch in the number of models).

## Credits

Built with [Claude Code](https://claude.ai/code). The expert team was assembled by asking Claude to name five experts for building an open-source codec competing with H.264/H.265:

- **Jim Bankoski** — Google architect behind VP8 and VP9; open-codec veteran
- **Timothy Terriberry** — Mozilla/Xiph engineer; drove Daala and AV1 research
- **Monty Montgomery** — Xiph.org founder; built Ogg, Vorbis, and Theora ecosystems
- **Gary Sullivan** — Microsoft; co-chaired the H.264 and HEVC standardization efforts
- **Jens-Rainer Ohm** — HEVC co-chair; deep expertise in perceptual coding theory

Plus a principal engineer (facilitator) and patent reviewer. All personas simulated by Claude.

The October 2026 improvements were also done with Claude Code, without personas, much of it autonomously toward fixed targets: each change was measured by BD-rate against x264 on the Xiph clips and kept only if it helped, except replacements made for their sources. The git history records each change, its measured effect, and the changes tried and dropped.

## Requirements

- C11 compiler (GCC or Clang)
- make
- ffmpeg (for Y4M generation)
- No external dependencies (Python 3 only to retrain the entropy tables)
