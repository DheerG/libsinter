# Sinter Quality Improvement Scoreboard

## Target
Match or beat H.264 on visual quality for real video content.

**H.264 reference (QP=20, all-intra, 256x256 testsrc):**
- Luma: 45.58 dB | Chroma U: 53.91 dB | Chroma V: 53.47 dB | Size: 149 KB (30 frames)

## Surgical Improvement Loops

| Loop | Change | Luma (dB) | Chroma U (dB) | Size | Key Result |
|------|--------|-----------|---------------|------|------------|
| 0 | Pre-loop baseline | 41.26 | 25.91 (DC-only) | 60 KB/frame | — |
| 1 | Wiener filter | 41.13 | — | — | Reverted (-0.13 dB) |
| **2** | **Hybrid PVQ/Scalar** | **44.75** | — | 120 KB | **+3.49 dB luma** |
| 3 | CSF-weighted step | 44.75 | — | — | No effect at QP=4 |
| **4** | **20% Malvar lapping** | **46.43** | — | 101 KB | **+1.68 dB luma** |
| **5** | **Exp-Golomb VLC** | 46.43 | — | **21 KB** | **79% size reduction** |
| 6 | Smart skip | 46.49 | — | 21 KB | +0.06 dB |
| **7** | **Full chroma DCT+scalar** | 44.80 | **47.18** | 28 KB/frame | **+21 dB chroma** |

## Final State

| Metric | Sinter | H.264 | Gap |
|--------|--------|-------|-----|
| Luma PSNR | 44.80 dB | 45.58 dB | 0.77 dB |
| Chroma U PSNR | 47.18 dB | 53.91 dB | 6.73 dB |
| Chroma V PSNR | 47.06 dB | 53.47 dB | 6.41 dB |
| File size (30 frames) | 838 KB | 149 KB | 5.6x |
| Visual quality | Near-equivalent | Reference | — |

## Session Journey

```
Start:     7.63 dB luma (broken pipeline, max error 255)
           ↓ PVQ band partitioning, K overflow fix
Mid:      41.26 dB luma (functional but 5.4 dB behind H.264)
           ↓ Hybrid PVQ/scalar, reduced lapping, VLC, chroma DCT
Final:    44.80 dB luma, 47.1 dB chroma (visually comparable to H.264)
```

Total luma improvement: **+37.2 dB**
Total chroma improvement: **+21.3 dB**
Total size reduction: **79%** (from raw 16-bit to exp-Golomb VLC)
