# Compression experiments

This file is the work plan for the `more-compression` branch, the rules for each experiment, and a log of every experiment tried, kept or dropped. Read the log before trying an idea: if it is there, either skip it or change the shape that failed and say why.

**Outcome.** Every item of the queue is kept or dropped, and the last three experiments from the bit breakdown each saved under 0.5%, so the plan is done. Against `main` (38f3a7d), on the tuning clips at 60 frames: luma BD-rate -23.0%, YUV -22.6%, SSIM -27.0% (the full-length results against x264 are in the README). Kept: items 1-9 and 11, the frame-level part of 13, and five context changes from the bit breakdown. Dropped: 10 and 12 after measuring their upper bounds, 13b for patent risk, 14 after measuring it.

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
| 1 | **RD-optimized quantization** (encoder). Per coefficient, try the rounded level, one level toward zero, and zero; keep the lowest SSE + lambda x bits from the models' current state. Then a trellis over each band. Luma first, chroma after item 3. Don't copy HEVC's RDOQ (last-position search, coefficient groups). [3-6%] | K. Ramchandran and M. Vetterli, "Rate-distortion optimal fast thresholding with complete JPEG/MPEG decoder compatibility", IEEE TIP 3(5), 1994; J. Wen, M. Luttrell and J. Villasenor, "Trellis-based R-D optimal quantization in H.263+", IEEE TIP 9(8), 2000 | luma kept (-1.5%); chroma done in item 3 |
| 2 | **Real model costs in the motion search** (encoder). `motion_est.c` costs vectors with exp-Golomb lengths. Also count chroma SSE and bits in the skip and mode decisions, and refine bi-predicted vectors jointly. [1-2%, more on YUV] | G. J. Sullivan and T. Wiegand, "Rate-distortion optimization for video compression", IEEE Signal Processing Magazine, 1998; S.-W. Wu and A. Gersho, "Joint estimation of forward and backward motion vectors for interpolative prediction of video", IEEE TIP 3(5), 1994 | MV candidates and RD refinement kept (-6.0%); chroma in decisions kept (-0.3% Y, -0.6% YUV); model costs in the search dropped |
| 3 | **Chroma coded like luma.** A coded flag per chroma block, levels as JPEG's binary decisions, and RD zeroing in the encoder. Today every non-skipped 4x4 chroma block codes 16 symbols of 256-symbol alphabets, with no RD decision. [1-3%, more on YUV] | ITU-T H.261 coded block pattern (covers chroma), 1988; ITU-T T.81 (JPEG) arithmetic coding, all components, 1992 | kept (-3.4% Y, -2.6% YUV) |
| 4 | **Luma DC as JPEG's DC decisions**, with JPEG's context: the class of the previous block's DC difference (zero, small or large, by sign). [1-2%] | ITU-T T.81 Annex F.1.4.4.1, 1992 | kept (-0.6%) |
| 5 | **MV deltas as binary decisions** (JPEG's zero, sign, category and magnitude bits). Contexts: for y, whether x was zero; block size; direction. Not CABAC's context from the neighbours' MVD sizes (H.264). [1-2%] | ITU-T T.81, 1992; Rissanen and Langdon, 1981 | kept (-0.6%), with the agreement class (new, NOVEL.md 5) |
| 6 | **Mode contexts from the left and above macroblocks.** The skip/mode map is a binary image, coded with a template of neighbours as JBIG does. Separately, try sending the intra mode once per macroblock, as Annex I does. A neighbour context for intra modes is close to H.264's most probable mode: search before trying it. [1-2%] | ITU-T T.82 (JBIG), 1993; MPEG-1 skipped-macroblock runs, 1993; H.263 Annex I, 1998 | skip contexts kept (-0.3%); per-MB intra modes and a B-direction context dropped |
| 7 | **Multiple reference frames** for P-frames: a reference index per macroblock or 8x8 block, coded with an adaptive model. Check memory: each reference keeps a spline plane of 4 bytes per sample. [2-5%] | ITU-T H.263 Annex U, 2000; T. Wiegand, X. Zhang and B. Girod, "Long-term memory motion-compensated prediction", IEEE TCSVT 9(1), 1999 | three references for P-frames kept (-4.1%) |
| 8 | **OBMC.** Write it from Annex F: `core/prediction/obmc.c` is from the original design (its own 8-tap filters) and unused. Retune the deblocking QUANT afterwards. [1-4%] | ITU-T H.263 Annex F, 1996; M. T. Orchard and G. J. Sullivan, "Overlapped block motion compensation: an estimation-theoretic approach", IEEE TIP 3(5), 1994 | kept (-2.0%), above and left neighbours only |
| 9 | **Direct mode for B-frames:** the co-located anchor vector scaled by temporal distance, with an optional delta. [<1.5%] | ITU-T H.263 Annex G (PB-frames), 1996; ISO/IEC 14496-2 direct mode, 1999 | kept (-1.5%), no delta vector |
| 10 | **Affine global motion** (2 or 3 warping points), and let coded blocks choose the global motion. [0-2%, foreman and bus] | ISO/IEC 14496-2 version 2 GMC and mcsel, 2000 | dropped after measuring (global motion predicts under 1% of skipped blocks) |
| 11 | **Quintic B-spline** luma interpolation. [0-1%, mobile and bus] | M. Unser et al., 1991; M. Unser, 1993 (the current sources) | kept (-1.5%, mobile -3.9%) |
| 12 | **32x32 blocks** for skip and one-vector inter. Overlaps with item 6: try that first. [0-2%, akiyo and news] | G. J. Sullivan and R. L. Baker, "Efficient quadtree coding of images and video", IEEE TIP 3(3), 1994 | dropped after measuring (skip runs already cost 0.15 bits per 2x2 group) |
| 13 | **QP per macroblock** with look-ahead allocation: lower QP where many later frames predict from a block. Measure against x264 CQP and CRF. MB-tree (x264, 2009) is later and has later patents nearby: follow the 1994 paper and search first. Also activity masking for SSIM as an option, off by default (it lowers PSNR). [maybe 5-10% on akiyo and news] | ITU-T H.263 DQUANT, 1996; K. Ramchandran, A. Ortega and M. Vetterli, "Bit allocation for dependent quantization with applications to multiresolution and MPEG video coders", IEEE TIP 3(5), 1994; MPEG-2 Test Model 5, 1993 | frame-level part kept (-1.3%); per-macroblock look-ahead QP built, dropped (patent risk) |
| 14 | **Time-domain lapped transform** (pre- and post-filters around the DCT). Large, uncertain on inter frames; last. | T. D. Tran, J. Liang and C. Tu, "Lapped transform via time-domain pre- and post-filtering", IEEE Trans. Signal Processing 51(6), 2003 | key-frame luma version built, dropped (+2.1% at best) |

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
| 2c. Chroma in the skip/inter/intra decisions: each candidate's 4x4 chroma blocks costed as the chroma pass codes them (intra DC prediction estimated from the source) | Weight 1.0: -0.30% Y, -0.72% YUV, +0.32% SSIM (mobile +1.46%); 0.5: -0.61%, -0.50%, +0.12%; 2.0: +0.73%, -0.41%, +1.31% | Deferred until after item 3: every non-skipped block pays 16 chroma symbols today, so chroma costs push towards skip and lose luma detail. Retry once chroma has a coded flag. |
| 3. Chroma coded like luma: a coded flag per 4x4 block (context: block type, luma block coded, plane and U's flag for V), 16 levels in zigzag order as binary decisions, RD levels and coded/zero in the encoder; container version 5 | -3.38% Y, -2.55% YUV, -3.45% SSIM with retrained tables; kept | Chroma lambda 0.25x the luma lambda. Untrained, against untrained: 0 (plain rounding, entropy coding only) -2.28% Y, -2.31% YUV; 0.25x -3.53%, -2.61%; 0.5x -4.92%, -2.48%; 1.0x -6.53%, -0.47%. Above 0.25x the luma gain is bought with chroma (YUV gets worse). Stills of mobile at QP 28: no new artifacts; faint red in small digits fades a little more. |
| 2c, after item 3. Chroma in the decisions, costed with the chroma pass's own RD choice (chroma_rd: levels, coded flag, V's context from U) | Weight on chroma SSE 1.0: -0.30% Y, -0.57% YUV, -0.07% SSIM; kept | 0.5: -0.55%, -0.37%, -0.32%; 0.75: -0.33%, -0.41%, -0.11%; 2.0: +1.12%, +0.05%; 4.0: +7.2%, +4.9%. A heavy chroma weight makes the modes follow chroma at luma's expense. Encoder at most 501 ms per frame. |
| Global lambda 0.12 and 0.165, again after chroma RD (item 3) | 0.12: -0.11% Y, +0.92% YUV; 0.165: +0.77% Y, -0.52% YUV | Still a luma/chroma trade; 0.14 stays. |
| 4. Luma DC as JPEG's DC decisions (Z, S, G by sign, X, M), context = JPEG's class of the left block's DC level | -0.46% Y, -0.47% YUV with retrained tables | Superseded by the next row. |
| 4. Plus whether the above block's DC level is nonzero (ten contexts) | -0.62% Y, -0.63% YUV, -0.77% SSIM with retrained tables; kept | Untrained: the above context alone -0.19%. The DC still takes 12-15% of akiyo's and foreman's bits at QP 28, most of it the zero decision. |
| 5. MV deltas as binary decisions; Z and G by component, block size and (y) whether x is zero | -0.37% Y, -0.37% YUV with retrained tables | Superseded by the next row. |
| 5. Plus the prediction agreement class (tool 4's) for Z and G: new, searched (NOVEL.md 5) | -0.64% Y, -0.59% YUV, -0.68% SSIM with retrained tables; kept | The agreement class alone: -0.35% (untrained on both sides). |
| Escape section removed (nothing used it after items 3-5) | -0.04% | Simpler format. |
| 6. Mode flags (macroblock mode, inter, skip) also by how many of the left and above blocks are skipped (JBIG-style template) | Untrained: -0.10%; retrained tables: -0.30% Y, -0.32% YUV, -0.26% SSIM; kept | Three times the contexts, so it needs trained starting tables to pay. |
| 6. Key frames: one intra mode per macroblock, as H.263 Annex I sends it (encoder tries each mode over the four blocks) | +0.03% Y, +0.02% YUV, +0.37% SSIM | Mode bits saved are lost in residual bits. A neighbour context for intra modes was not tried: it is close to H.264's most probable mode. |
| 6. B-frame direction by the left coded block's direction (4 contexts) | 0.00% (untrained) | The agreement and skip contexts already carry it. |
| 7. P-frames choose between the newest anchor and the one before it (same GOP), flag per inter block or macroblock; skip takes the left block's reference | Untrained: -2.49% Y; retrained tables: -2.41% Y, -2.26% YUV, -2.90% SSIM; kept | No new storage (the older anchor is B-frames' forward reference). For the second P-frame of a GOP the second reference is the key frame, coded 4 QP finer. Encoder at most 582 ms per frame. |
| 7. The GOP's key frame instead of the older anchor as P-frames' second reference (long-term) | -1.05% Y against two references (tables not retrained): foreman -3.56%, mobile +0.85% | The two help different clips, so both are kept: next row. |
| 7. Three references for P-frames: newest anchor, the one before, and from the third anchor on the GOP's key frame; flags "not newest", then "key frame" | Untrained: -1.88% Y; retrained tables: -1.76% Y, -1.59% YUV, -2.14% SSIM against two; kept | The key frame is coded 4 QP finer, and static background predicts well from it. Encoder at most 628 ms per frame (budget 650): a fourth reference needs a faster search first. |
| Encoder speed: SAD with early exit (same result), log2 of counts from a table | -0.04% Y; slowest frame 628 ms -> 484 ms | Needed to stay inside 650 ms with OBMC. |
| 8. OBMC with Annex F's weights, above and left neighbours only (the decoder needs the prediction before the residual), in P- and B-frames, all references | -2.04% Y, -2.02% YUV, -2.08% SSIM with retrained tables; kept | Measured first against a stale baseline as -3.8% (that included the third reference); corrected. P-frames only: +0.40% against P and B. Deblocking QUANT step/4 +0.27%, step/7 +0.01% Y but +0.57% SSIM: step/5 stays. Stills at QP 28 show fewer block edges. Encoder at most 532 ms, decoder 7.9 ms per frame. |
| 9. B-frame direct mode: a fourth B-direction symbol, MVs from the newest anchor's co-located block scaled by TRB/TRD (MPEG-4), for 8x8 blocks and 16x16 macroblocks, no delta | Untrained: -1.51%; retrained tables: -1.51% Y, -1.50% YUV, -1.57% SSIM; kept | Co-located blocks on the second or third reference count as zero motion. Encoder at most 550 ms, decoder 8.0 ms per frame. |
| 9. B-frame skipped blocks fall back to direct motion instead of the global motion (when the left block isn't inter) | -0.04% Y, 0.00% SSIM | Neutral; the global-motion fallback stays. |
| 11. Quintic B-spline luma interpolation (gain 120, two poles, 6 x 6 weights exact in 1/122880ths) instead of cubic | Untrained: -1.39%; retrained tables: -1.49% Y, -1.23% YUV, -1.02% SSIM; kept | Mobile -3.92%, akiyo -0.50%, foreman -0.05% (SSIM +0.56% on foreman). Encoder at most 623 ms, decoder 10.5 ms per frame: little encoder time left in the budget. |
| Encoder speed: plain 8x8 predictions cached per block and motion (OBMC's neighbour predictions repeat for every candidate) | Output identical; slowest frame 623 ms -> 558 ms | |
| 10. Affine global motion; coded blocks using the global motion (MPEG-4 mcsel) | Dropped after measuring: at QP 16-28, skipped blocks predicted with the global motion are 0.3-1% of skipped blocks (foreman 263 of 50909, mobile 448 of 52330, bus 830 of 43259); the rest repeat the left block's motion | Affine warping would improve almost nothing; coded blocks already have the skip motion among their MV candidates. |
| 12. 32x32 blocks for skips (quadtree top level) | Dropped after measuring the upper bound: in 2x2 macroblock groups that are all SKIP16, the mode symbols cost 0.15 bits per group (akiyo QP 28: 759 bits, 1.3% of the file; QP 16: 0.7%; foreman QP 28: 1.2%, but only 28% of groups are all-skip) | A 32x32 flag must be paid in every group, so the JBIG-style skip contexts (item 6) leave almost nothing to gain. |
| 13a. Frame-level dependent quantization: anchors with few frames left after them in the closed GOP coded coarser (+1 QP with 8 or fewer left, +2 with 2, +3 for the last frame); encoder only, needs the caller's keyint | -1.25% Y, -1.44% YUV, -2.00% SSIM with retrained tables; kept | Finer early anchors (-1 QP for the first three): +0.38%. One step at 5/8/11 frames: -0.72/-0.68/-0.53%; two steps (8, 2): -1.12%. Searched: "patent increasing quantization parameter for P frames near end of GOP fewer remaining frames reference propagation closed GOP encoder": closest are bit-budget rate control (US8223836, US12413738, US20120076202), not propagation. Being a frame-level QP schedule, it should also be compared against x264 with MB-tree on (CRF) in the final report. |
| Key-frame and B-frame QP offsets re-tuned with this branch's tools: key 3/5/6, B 3/5/6, combinations (both were 4) | Key 5 and B 5: -0.76% Y, -1.62% YUV, -3.96% SSIM with retrained tables; kept | Key 3: +1.29%; B 3: +1.02%. Key 6 with B 5: -0.63% Y, -2.33% YUV, -4.98% SSIM (luma slightly worse than key 5: an SSIM-leaning choice). Key 5 with B 6: -0.54%. |
| 13b. QP per key-frame macroblock (H.263 DQUANT syntax, -2..+2): -2 / -1 where the mean absolute difference to the same area of the rest of the GOP is below 1 / 2 (static, so many frames predict from it) | Untrained: -0.89% Y, -1.52% YUV, -3.53% SSIM (akiyo SSIM -7.6%); dropped | Thresholds 0.7-3 and +1 for changing areas all within 0.05%; a -3 level traded luma for SSIM. Dropped for patent risk: per-macroblock QP from a look-ahead is close to active patents (MB-tree / temporal importance maps: US11228766, US10757419, US11159801; look-ahead QP: EP2377321A1, US9774848), and no pre-2004 source was found for the macroblock-level rule (Ramchandran, Ortega and Vetterli 1994 is frame-level; Puri and Aravind 1991 could not be confirmed to cover it). Revisit only with a clear prior-art source. |
| From the bit breakdown (intra modes 6.8% of akiyo's bits at QP 28): the intra mode's model chosen by MPEG-4 Part 2's DC-gradient direction rule (3 contexts per frame type) | Untrained: -0.19%; retrained tables: -0.23% Y, -0.27% YUV, -0.18% SSIM; kept | The direction MPEG-4 derives predicts the RD-chosen mode only weakly. A context from the neighbours' modes was not tried (close to H.264's most probable mode). |
| 14. Time-domain lapped transform on key-frame luma (Tran, Liang and Tu 2003): prefilter 1/2 W diag(I, V) W across each internal 8x8 boundary, V = J C_II^T D_s C_IV J, post-filter with V^-1 in Q12 in place of deblocking | s = 1.0: +4.53% Y; s = 1.2: +2.08%; s = 1.6: +2.66% (akiyo +5%); dropped | Helps at low rates (akiyo QP 28: +0.12 dB at 1% more bytes) and hurts at high rates (QP 8: -0.65 dB). This shape has an inexact round trip (float prefilter, Q12 post-filter), RD decisions in the prefiltered domain and Annex I prediction not adapted to lapping; an exactly invertible lifting version with weighted RD would be a different, larger design. |
| After the queue, bit breakdown 1 (luma AC levels 30-60% of the bits): Z and G of luma AC levels also by how many of the two previous levels in zigzag order are nonzero (C. Tu and T. D. Tran, IEEE TIP 11(11), 2002) | Untrained: -0.32%; retrained tables: -0.36% Y, -0.25% YUV, -0.87% SSIM; kept | First of the "three in a row under 0.5%" of the done rule. |
| Bit breakdown 2 (band flags 9-12% of the bits): band flag also by how many of the left and above blocks code the same band (Tu and Tran 2002, space-frequency context) | Untrained: -0.44%; retrained tables: -0.73% Y, -0.55% YUV, -0.90% SSIM; kept | Over 0.5%: the count of three restarts. |
| Bit breakdown 3 (chroma up to 13% of the bits): chroma coded flag also by whether the same plane's left block is coded | Untrained: -0.11%; retrained tables: -0.21% Y, -0.15% YUV, -0.17% SSIM; kept | First of three under 0.5% after the restart. |
| Bit breakdown 4 (MV deltas up to 19% of foreman's bits): MV-delta category decisions also by the prediction agreement class | Untrained: -0.07%; retrained tables: -0.05% Y, -0.06% YUV, 0.00% SSIM (foreman +0.12%); dropped | Second of three under 0.5%. The agreement class already sets the zero and "above one" models. |
| Bit breakdown 5 (macroblock mode up to 11% of foreman's bits): mode model also by how many of the left and above macroblocks are SPLIT | Untrained: -0.24%; retrained tables: -0.04% Y, -0.12% YUV, -0.18% SSIM; dropped | Third under 0.5% in a row: the done rule's experiments are complete. |
