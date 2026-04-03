# Sinter Quality Improvement Scoreboard

Target: Match or beat H.264 QP=20 on 256x256 test scene
- H.264: **46.71 dB PSNR, 4,382 bytes, max error 43**

## Loops

| Loop | Change | PSNR (dB) | Size (bytes) | Max Err | Gap to H.264 |
|------|--------|-----------|--------------|---------|---------------|
| 0 | Baseline (broken pipeline) | 7.63 | — | 255 | 39.1 dB |
| — | PVQ band partitioning (4 bands) | 29.95 | 46,794 | — | 16.8 dB |
| — | 9-band PVQ | 37.16 | 67,909 | — | 9.6 dB |
| — | 16-band PVQ | 39.91 | 67,909 | 10 | 6.8 dB |
| — | CSF weights + DC round-to-nearest | 41.18 | 62,146 | 17 | 5.5 dB |
| — | PVQ skip threshold + AC prediction | 41.26 | 60,409 | 16 | 5.4 dB |
| 1 | (pending) | | | | |

## Rules
- One surgical change per loop
- PE runs roundtable before each change
- Patent reviewer clears each technique
- If 20 loops don't close the gap, rethink architecture
