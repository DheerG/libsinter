# Compression experiments

This file is the work plan for the `more-compression` branch, the rules for each experiment, and a log of every experiment tried, kept or dropped. Read the log before trying an idea: if it is there, either skip it or change the shape that failed and say why.

## Rules

- **Sources.** A format tool (anything the decoder runs) follows [PROVENANCE.md](../PROVENANCE.md): it comes from a publication before 2004, never from H.264/AVC, its drafts or a later codec (VP8 onward, Daala, AV1), or it is new and a patent search is recorded in [NOVEL.md](NOVEL.md) before it is kept. Encoder-only methods can be patented too, so each one also needs a pre-2004 source in its commit, or a recorded search. If a search finds the same idea in a patent or application, drop or redesign it.
- **Measure** each change by BD-rate against the commit before it (luma, YUV and SSIM), as below. Keep it if the mean luma BD-rate falls and YUV does not rise. Below 0.2%, keep it only if it adds little code.
- Encoder and decoder reconstructions match on every test point, and `make test` passes.
- For a format change, look at stills of P- and B-frames at QP 28 for new artifacts.
- Bump the container version (to 5) at the first format change on this branch.
- Retrain the starting tables ([tools/train](../tools/train/README.md)) when the model set changes. Never train on foreman, mobile, akiyo, bus or news.
- Speed: encoder at most 650 ms per CIF frame, decoder at most 18 ms.
- Commit each kept change on its own, with its saving and source. Log each dropped change below, with its shape, result and reason; don't commit its code. Don't push, and don't change `main`.
- The x264 numbers use `-qp`, which is constant QP: x264's embedded settings read `rc=cqp mbtree=0 aq=0` (and `trellis=1 ref=3`). Any tool that changes QP inside a frame or by look-ahead must also be measured against x264 in CRF mode, which turns MB-tree on.

## Measurement

The harness is in `tools/eval/` (Step 0, done). It reproduces the published ratios against x264 exactly.

```bash
tools/eval/fetch_clips.sh                 # Xiph CIF clips, first 150 frames, into clips/
tools/eval/run.sh base HEAD               # tuning points of a commit (or "." for the working tree)
tools/eval/run.sh new                     # the working tree
python3 tools/eval/bd.py base new         # BD-rates: Y, YUV, SSIM per clip and mean, mismatches, speed
BITS=1 tools/eval/run.sh bits             # bit breakdown ...
python3 tools/eval/bits.py build/eval/bits.bits   # ... by frame type and model group
tools/eval/x264.sh x264                   # x264 points (crf as a second argument: MB-tree on)
CLIPS="foreman mobile akiyo bus news" QPS="4 8 12 16 20 24 28 32" FRAMES=150 tools/eval/run.sh full
python3 tools/eval/bd.py x264 full --ratio
```

- **Tuning:** foreman, mobile and akiyo, first 60 frames, QP 8, 12, 16, 20, 24 and 28 (run.sh's defaults). BD-rate is the classic Bjøntegaard cubic fit over the overlapping range. PSNR is the mean of per-frame PSNR; YUV weights Y:U:V as 6:1:1; SSIM is luma SSIM computed as ffmpeg does, in dB for the BD-rate.
- **Reports against x264:** all five clips, 150 frames, Sinter QP 4 to 32 in steps of 4, and x264 `-preset medium -tune psnr -g 30 -qp 16` to `44`. Bus and news are held out: measure them only at the end, never to tune.
- **Luma against YUV.** A parameter that moves bits between luma and chroma (lambda, chroma QP) shows a luma gain with a YUV loss, or the reverse. That is a trade, not a gain: keep such a change only if both improve.

### Where the bits went at the start (main at 38f3a7d, QP 28, 60 frames)

| Clip | Key frames | Luma AC | Band flags | Luma DC | Modes | MV deltas | Chroma |
|------|-----------|---------|------------|---------|-------|-----------|--------|
| Akiyo | 72% | 32% | 11% | 13% | 15% | 9% | 16% |
| News | 55% | 39% | 11% | 11% | 14% | 12% | 12% |
| Mobile | 45% | 45% | 10% | 5% | 12% | 15% | 13% |
| Foreman | 31% | 25% | 8% | 14% | 21% | 24% | 7% |
| Bus | 20% | 42% | 13% | 9% | 14% | 18% | 4% |

"Key frames" is their share of the file; the other columns are shares of all bits. At QP 16, chroma is 18% of news and 16% of akiyo.

## Queue

Try these roughly in order; the order follows the expected saving and the cost of the work. Gains in brackets are estimates, not measurements. Record each result in the log and update the status here.

| # | Experiment | Source | Status |
|---|------------|--------|--------|
| 1 | **RD-optimized quantization** (encoder). Per coefficient, try the rounded level, one level toward zero, and zero; keep the lowest SSE + lambda x bits from the models' current state. Then a trellis over each band. Luma first, chroma after item 3. Don't copy HEVC's RDOQ (last-position search, coefficient groups). [3-6%] | K. Ramchandran and M. Vetterli, "Rate-distortion optimal fast thresholding with complete JPEG/MPEG decoder compatibility", IEEE TIP 3(5), 1994; J. Wen, M. Luttrell and J. Villasenor, "Trellis-based R-D optimal quantization in H.263+", IEEE TIP 9(8), 2000 | luma kept (-1.5%); chroma with item 3 |
| 2 | **Real model costs in the motion search** (encoder). `motion_est.c` costs vectors with exp-Golomb lengths. Also count chroma SSE and bits in the skip and mode decisions, and refine bi-predicted vectors jointly. [1-2%, more on YUV] | G. J. Sullivan and T. Wiegand, "Rate-distortion optimization for video compression", IEEE Signal Processing Magazine, 1998; S.-W. Wu and A. Gersho, "Joint estimation of forward and backward motion vectors for interpolative prediction of video", IEEE TIP 3(5), 1994 | MV candidates and RD refinement kept (-6.0%); model costs in the search dropped; chroma in decisions next |
| 3 | **Chroma coded like luma.** A coded flag per chroma block, levels as JPEG's binary decisions, and RD zeroing in the encoder. Today every non-skipped 4x4 chroma block codes 16 symbols of 256-symbol alphabets, with no RD decision. [1-3%, more on YUV] | ITU-T H.261 coded block pattern (covers chroma), 1988; ITU-T T.81 (JPEG) arithmetic coding, all components, 1992 | to do |
| 4 | **Luma DC as JPEG's DC decisions**, with JPEG's context: the class of the previous block's DC difference (zero, small or large, by sign). [1-2%] | ITU-T T.81 Annex F.1.4.4.1, 1992 | to do |
| 5 | **MV deltas as binary decisions** (JPEG's zero, sign, category and magnitude bits). Contexts: for y, whether x was zero; block size; direction. Not CABAC's context from the neighbours' MVD sizes (H.264). [1-2%] | ITU-T T.81, 1992; Rissanen and Langdon, 1981 | to do |
| 6 | **Mode contexts from the left and above macroblocks.** The skip/mode map is a binary image, coded with a template of neighbours as JBIG does. Separately, try sending the intra mode once per macroblock, as Annex I does. A neighbour context for intra modes is close to H.264's most probable mode: search before trying it. [1-2%] | ITU-T T.82 (JBIG), 1993; MPEG-1 skipped-macroblock runs, 1993; H.263 Annex I, 1998 | to do |
| 7 | **Multiple reference frames** for P-frames: a reference index per macroblock or 8x8 block, coded with an adaptive model. Check memory: each reference keeps a spline plane of 4 bytes per sample. [2-5%] | ITU-T H.263 Annex U, 2000; T. Wiegand, X. Zhang and B. Girod, "Long-term memory motion-compensated prediction", IEEE TCSVT 9(1), 1999 | to do |
| 8 | **OBMC.** Write it from Annex F: `core/prediction/obmc.c` is from the original design (its own 8-tap filters) and unused. Retune the deblocking QUANT afterwards. [1-4%] | ITU-T H.263 Annex F, 1996; M. T. Orchard and G. J. Sullivan, "Overlapped block motion compensation: an estimation-theoretic approach", IEEE TIP 3(5), 1994 | to do |
| 9 | **Direct mode for B-frames:** the co-located anchor vector scaled by temporal distance, with an optional delta. [<1.5%] | ITU-T H.263 Annex G (PB-frames), 1996; ISO/IEC 14496-2 direct mode, 1999 | to do |
| 10 | **Affine global motion** (2 or 3 warping points), and let coded blocks choose the global motion. [0-2%, foreman and bus] | ISO/IEC 14496-2 version 2 GMC and mcsel, 2000 | to do |
| 11 | **Quintic B-spline** luma interpolation. [0-1%, mobile and bus] | M. Unser et al., 1991; M. Unser, 1993 (the current sources) | to do |
| 12 | **32x32 blocks** for skip and one-vector inter. Overlaps with item 6: try that first. [0-2%, akiyo and news] | G. J. Sullivan and R. L. Baker, "Efficient quadtree coding of images and video", IEEE TIP 3(3), 1994 | to do |
| 13 | **QP per macroblock** with look-ahead allocation: lower QP where many later frames predict from a block. Measure against x264 CQP and CRF. MB-tree (x264, 2009) is later and has later patents nearby: follow the 1994 paper and search first. Also activity masking for SSIM as an option, off by default (it lowers PSNR). [maybe 5-10% on akiyo and news] | ITU-T H.263 DQUANT, 1996; K. Ramchandran, A. Ortega and M. Vetterli, "Bit allocation for dependent quantization with applications to multiresolution and MPEG video coders", IEEE TIP 3(5), 1994; MPEG-2 Test Model 5, 1993 | to do |
| 14 | **Time-domain lapped transform** (pre- and post-filters around the DCT). Large, uncertain on inter frames; last. | T. D. Tran, J. Liang and C. Tu, "Lapped transform via time-domain pre- and post-filtering", IEEE Trans. Signal Processing 51(6), 2003 | to do |

After the queue, pick each next experiment from a fresh bit breakdown.

**Done when** every item in the queue is kept or dropped, and three experiments in a row from the bit breakdown each save under 0.5%. Then:
- measure all five clips against x264;
- update the README, PROVENANCE.md, NOVEL.md (for new tools) and this file;
- make a comparison page like the last one;
- stop.

**Stay away from**, as too close to excluded sources or patents:
- HEVC: merge and MV candidate lists, SAO, ALF.
- H.264: most probable intra mode, CABAC's neighbour contexts, 16x8 and smaller partitions.
- B-frames used as references (hierarchical B).
- VVC: adaptive MV resolution.
- RD decisions based on SSIM.
- Probability mixing (Google's patents; see NOVEL.md).

## Log

### Before this branch (from the git history)

| Experiment | Result | Shape tried, and what would differ next time |
|------------|--------|-----------------------------------------------|
| Coded-band count per block | +2% | A count of coded bands before the band flags. |
| Diagonal intra modes | 0.1%, dropped | Pixel-domain modes, now excluded anyway (H.264-like). |
| MV-delta contexts | 0.1%, dropped | With 256-symbol MV models, which relearn slowly. Binary decisions (queue item 5) are a different shape. |
| Chroma QP offset changes | Traded luma for chroma | Inter +1, key frames +4 kept. |
| AC level models by the previous level's magnitude (0, 1, 2+) | +1.0% | With 256-symbol level models. The binary levels later made such contexts cheap. |
| Motion search range +/-32 | Neutral | Full search at integer pel, then half and quarter. |
| Wider dead zone for inter residuals | +1.5% | A fixed encoder threshold. RD quantization (queue item 1) decides each coefficient from its cost. |
| Lambda 0.7x, 1.4x, 2x of 0.14 step^2 | All worse | For trial-encoding mode decisions. |
| 1 or 3 B-frames between anchors | Within 0.5% of 2 | Fixed GOP; no adaptive placement. |
| B-frame QP offset 0, 2, 5, 6 | Worse than 4 | |
| Deblocking QUANT 3/2, 1, 3/4, 1/3, 1/4, 2/7 of the step | Worse than 1/5 | H.263 Annex J. No deblocking: +3.7%. |
| B-frames left unfiltered | -0.4% YUV, dropped | Strong 8x8 blocking in coarse B-frames. |
| Filtering edges between skipped blocks whose motion differs | Neutral | Annex J's rule kept. |
| Global-motion skip threshold 10 and 100 lambda | Worse than 30 | |
| Trained-table weight 4 to 256 observations | 16 and up within 0.1% | 64 for coefficients, 16 for mode flags. |
| Prediction-class thresholds 2x lower or higher | No better | NOVEL.md, tool 2. |
| Agreement thresholds 2x coarser | Worse | NOVEL.md, tool 4. |
| Logistic mixing of models | -0.5%, left out | Google's patents. Behind `-DSNT_EXPERIMENT_MIXING`; NOVEL.md, tool 1. |
| Per-clip starting tables | -0.1%, left out | VP8's design. Behind `-DSNT_EXPERIMENT_CLIPTABLES`; NOVEL.md, tool 3. |
| PVQ | Removed, 0.3% better without | Used Daala's energy correction. |
| Lapping of the original design | Removed | Not DC-preserving (+/-11 sawtooth at block edges). Queue item 14 is a different design. |

### On this branch

| Experiment | Result | Shape tried, and what would differ next time |
|------------|--------|-----------------------------------------------|
| 1. RD quantization, greedy: each level from {rounded, one toward zero, zero} by SSE + lambda x bits | -1.25% Y, -1.63% YUV, -2.00% SSIM | Superseded by the next row. |
| 1. Plus coded/zero for all bands at once by a two-state Viterbi pass (band flag context = previous band coded), the DC's zero-ness as the first stage | -1.52% Y, -1.96% YUV, -1.76% SSIM; kept | Encoder only. Within a band the only dependency is the skipped Z of the last level, so a fuller trellis has little left to find. |
| 1. Lambda for the level choice 0.7x and 1.4x of the mode lambda | 0.7x: Y -0.6%, YUV +1.2% against 1.0x; 1.4x: Y +3.0% | A luma/chroma trade either way; the level lambda stays equal to the mode lambda. |
| 1. Global lambda 0.12 and 0.165 (from 0.14) with RD quantization | 0.12: Y -0.28%, YUV +0.92%; 0.165: Y +0.90%, YUV -0.60% | A trade again. Revisit after chroma has RD decisions (item 3). |
| 1. Reconstruction offset step/4 and step/10 (from step/6), format change | step/4: +0.92% Y; step/10: -0.24% Y, -0.10% YUV | Too small for a format change. The encoder's rounding offset matched to it (step/8, /10, /14): +0.9% to +3.0% Y, so the step/6 rounding acts as a useful dead zone for the first candidate. |
| 2a. MV-delta costs in the motion search from the MV models (exp-Golomb lengths before) | +0.06% Y | Neutral. Search lambda 0.7x and 1.4x with it: +0.22%, +0.45%. The SAD search only shortlists; the trial encoding decides (2b). |
| 2b. MVs decided by trial encoding among candidates: searched, predicted (no delta), skip motion, left, above and above-right blocks' MVs, the 16x16-searched MV for 8x8 blocks, the 8x8 blocks' MVs for INTER16 | -3.74% Y, -3.54% YUV | The predictor alone: -1.99%; with skip motion, left and above: -2.51%. |
| 2b. Plus bi-directional from the best forward and best backward candidate | -0.17% more | |
| 2b. Plus RD refinement: the 8 quarter-sample neighbours of each MV of the best inter candidate, the other MV held | -2.16% more; all of 2b -5.97% Y, -5.48% YUV, -5.81% SSIM; kept | Encoder at most 453 ms per CIF frame (mobile QP 8). A second refinement pass: -0.20% more for about 30% more encoder time; not kept, worth retrying if time is left. |
