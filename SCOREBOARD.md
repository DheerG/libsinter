# Sinter Quality Improvement Scoreboard

Target: Match or beat H.264 QP=20 on 256x256 test scene
- H.264: **46.71 dB PSNR, 4,382 bytes, max error 43**

## Surgical Improvement Loops

| Loop | Change | PSNR (dB) | Delta | Gap | Size | Max Err | Patent |
|------|--------|-----------|-------|-----|------|---------|--------|
| 0 | Baseline (pre-loop) | 41.26 | — | 5.45 | 60,409 | 16 | — |
| 1 | Wiener filter (reverted) | 41.13 | -0.13 | 5.58 | 60,409 | 16 | CLEAR |
| **2** | **Hybrid PVQ/Scalar per band** | **44.75** | **+3.49** | **1.96** | 120,218 | **13** | **CLEAR** |
| 3 | (next) | | | | | | |

## Rules
- One surgical change per loop
- PE runs roundtable before each change
- Patent reviewer clears each technique
- 20 loops max, then rethink if gap not closed

## Progress
- Session start: 7.63 dB (broken)
- Pre-loop: 41.26 dB
- Current: **44.75 dB** (gap: 1.96 dB to H.264)
- Sinter max error: **13** (H.264: 43) — Sinter is MORE uniform
