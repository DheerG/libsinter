# Sinter Video Codec

An experimental, vibecoded video codec built from scratch using [Claude Code](https://claude.ai/code) agent teams. Explores lapped transforms, perceptual vector quantization (PVQ), and rANS entropy coding as a patent-free alternative to the H.264/H.265 lineage.

**This was a learning experiment**, not a production codec. The goal was to test one-shot agent team workflows on a domain I had zero prior experience in — and to see how far a simulated expert team could push a novel architecture. Write-up: [One-Shot Wonder](https://dheer.co/one-shot-wonder/) | [Claude Agent Teams](https://dheer.co/claude-agent-teams/)

## Results

At comparable luma quality (~49 dB): **18.6x larger than H.264.** The architecture produces competitive perceptual quality but cannot match H.264's compression efficiency without adopting the same tools (B-frames, sub-pel MC, CABAC-level entropy coding).

| QP | Sinter PSNR | Sinter Size | H.264 PSNR | H.264 Size |
|----|-------------|-------------|-------------|------------|
| 4  | 49.10 dB    | 136 KB      | 74.59 dB    | 33 KB      |
| 20 | 33.70 dB    | 37 KB       | 55.12 dB    | 11 KB      |
| 28 | 28.41 dB    | 13 KB       | 48.37 dB    | 7 KB       |

(256x256 testsrc, 30 frames. Full BD-rate data in [SCOREBOARD.md](SCOREBOARD.md).)

## What We Built

~5,000 lines of C across 12 improvement loops:

- **Lapped transforms** (TDLT, 20% Malvar lifting) — eliminates blocking artifacts structurally
- **Hybrid PVQ/scalar quantization** — preserves texture where H.264 smooths it away
- **Dual-interleaved rANS** with 10 PARA-adaptive CDFs and cross-frame carry
- **Inter prediction** (P-frames, integer-pel full search, median MV prediction, skip mode)
- **Chroma inter** (luma MV at half-resolution for YUV 4:2:0)
- **CLI tools** with Y4M I/O (`sinterenc` / `sinterdec`)
- 32 passing tests, entirely patent-free

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

## Why the Gap Exists

The expert team's consensus: the 18.6x gap is the product of missing standard tools, each multiplicative:

| Missing Feature | Cost | Notes |
|----------------|------|-------|
| No sub-pel MC | 1.5-2x | Integer-pel can't track sub-pixel motion |
| No B-frames | 1.5-2x | Patent risk on bidirectional prediction |
| PVQ overhead | 1.5-2x | Gain-shape is more verbose than run-level |
| Fewer entropy contexts | 1.2-1.5x | 10 CDFs vs CABAC's hundreds |

Realistic patent-safe ceiling: **4-6x H.264** with half-pel + better contexts. Matching H.264 would require rebuilding most of its toolset, defeating the purpose.

## Credits

Built with [Claude Code](https://claude.ai/code). The expert team was assembled by asking Claude to name five experts for building an open-source codec competing with H.264/H.265:

- **Jim Bankoski** — Google architect behind VP8 and VP9; open-codec veteran
- **Timothy Terriberry** — Mozilla/Xiph engineer; drove Daala and AV1 research
- **Monty Montgomery** — Xiph.org founder; built Ogg, Vorbis, and Theora ecosystems
- **Gary Sullivan** — Microsoft; co-chaired the H.264 and HEVC standardization efforts
- **Jens-Rainer Ohm** — HEVC co-chair; deep expertise in perceptual coding theory

Plus a principal engineer (facilitator) and patent reviewer. All personas simulated by Claude.

## Requirements

- C11 compiler (GCC or Clang)
- make
- ffmpeg (for Y4M generation)
- No external dependencies
