# Encoder techniques

This guide describes how `sinterenc` makes its decisions: the motion vectors, the modes, the quantized levels and the QP of each frame. None of these methods is part of the bitstream. A decoder doesn't know or care how the encoder chose, so the same methods can be used in an encoder for another format, for example H.264.

For each technique the guide gives what it does, where it is in the code, what it saved, and what carries over to other formats. Most of these methods are well known: x264, for example, has its own versions of trellis quantization and RD mode decisions. The value here is a record of what worked in this codec, with measured numbers.

Savings are luma BD-rates (and YUV where it differs), the mean of foreman, mobile and akiyo (Xiph CIF, 60 frames, 2 B-frames), each against the step before it. [EXPERIMENTS.md](EXPERIMENTS.md) and the git history have the full record. Each method has a source from before 2004 or a recorded patent search. These records are a map, not a patent clearance.

All functions named below are in `libsinter/codec/encoder.c` unless another file is given.

## 1. Decide by trial encoding

Every decision compares candidates by the same cost:

    J = D + lambda * R

- **D** is the squared error of the candidate's reconstruction against the source.
- **R** is the bits the candidate would take: mode flags, vector deltas and coded residual.
- **lambda** is `RD_LAMBDA * step^2` with `RD_LAMBDA = 0.14`, where step is the quantizer step. Because the DCT is orthonormal, squared error in the coefficient domain equals squared error in pixels. Lambda 0.7x, 1.4x and 2x of this were all worse.

R comes from the entropy models in their current state: `cdf_bits` gives `log2(total) - log2(count)` of a model. `code_luma` has an estimate mode (no symbol sink) that runs the full coding of a block, levels included, and adds up these costs without emitting anything or changing the models. So a candidate's cost is close to what it would really take.

**Mode decision.** `decide_block` tries, for each 8x8 block, skip, inter coding at each vector candidate (section 2), and intra coding in each mode, and keeps the lowest J. `eval_mb16` tries the 16x16 candidates: skip, and inter with one vector for the whole macroblock. `code_macroblock` then compares split (the four best 8x8 choices) against the 16x16 ones. Trial encoding of the 8x8 modes saved 7.7% when it replaced three thresholds in the first round.

**Pitfalls:**
- **Provisional state.** The second block of a macroblock is decided with the first block's chosen reconstruction, vectors and DC already written in. Its neighbours' predictions and contexts depend on them. Emission later overwrites this state.
- **The 16x16 candidate's own blocks.** Inside a 16x16 candidate, a block's left or above neighbour may be another block of the same candidate. Its contexts (DC level, coded bands) and its overlapped prediction must use the candidate's own values for those blocks, not the provisional split decision (`inter16_cost`, and the `in16` flag of `pred_block`).
- **Frozen models.** Estimates use the models as they stand before the block. Emission then adapts them symbol by symbol, so later blocks see slightly different costs. This is accurate enough, and much simpler than tracking adaptation inside an estimate.

**Carries over:** RD mode decision is standard practice (Sullivan and Wiegand, IEEE Signal Processing Magazine, 1998). The pitfalls apply to any encoder with neighbour-dependent prediction or contexts.

## 2. Motion vectors

The search finds a shortlist; trial encoding chooses among the vectors; then the winner is refined. Together this saved **6.0%** (YUV 5.5%), the largest single gain of the second round.

**Search** (`search_mv`, `core/prediction/motion_est.c`): a full integer search over +/-16 pixels around zero, plus the predicted vector, by SAD + `sqrt(lambda)` x the exp-Golomb length of the vector delta. Then a half-sample and a quarter-sample step, trying the 8 neighbours each time. Using the entropy models' real vector costs in this search instead of exp-Golomb lengths changed nothing (+0.06%): the search only makes a shortlist.

**Candidates** (`mv_candidates`). Each is coded in full by trial encoding (`inter8_cost`, `inter16_cost`):
- the searched vector;
- the predicted vector, which costs no vector delta;
- the skip motion (the left block's motion, or the global motion);
- the vectors of the blocks to the left, above and above right;
- for an 8x8 block, the macroblock's 16x16 searched vector; for the 16x16 candidate, the vectors its four 8x8 blocks chose;
- in B-frames, each candidate forward, backward and bi-directionally, plus the best forward with the best backward (`mv_cand_best`);
- in P-frames with several references, the candidates again on each other reference, each with its own search.

Duplicates are skipped (`mv_cand_dup`).

How much each step added:

| Step | Luma |
|------|------|
| The predicted vector alone | -1.99% |
| With the skip motion and the left and above vectors | -2.51% |
| With the above-right vector and the other block size's vectors | -3.74% |
| With the bi-directional pair of the best one-way vectors | -0.17% more |
| With refinement (below) | -2.16% more |

**Refinement** (end of `decide_block` and `eval_mb16`): the 8 quarter-sample neighbours of each vector the best inter candidate uses, each costed by trial encoding, with the other direction's vector held. For bi-directional prediction this is S.-W. Wu and A. Gersho's joint estimation of forward and backward vectors (IEEE TIP 3(5), 1994). A second refinement pass saved only 0.2% more, for 30% more encoder time.

**Why it works.** The SAD search picks the vector that best matches the pixels, but the best vector to code is the one with the lowest total cost. A neighbour's vector or the predictor often codes in fewer bits, and its residual is nearly as small. Only trial encoding sees that trade.

**Carries over:** directly. An H.264 encoder can use the same candidates (predictor, neighbours, co-located vector, other partition sizes) costed with CABAC or CAVLC bit estimates.

## 3. RD-optimized quantization

`code_luma` (with `rd_candidates`) chooses each level by J instead of rounding it. It saved **1.5%** (YUV 2.0%):

- **Per coefficient**, the candidates are the rounded level, one step toward zero, and zero. The DC level is chosen the same way. This is the thresholding of K. Ramchandran and M. Vetterli ("Rate-distortion optimal fast thresholding with complete JPEG/MPEG decoder compatibility", IEEE TIP 3(5), 1994). It saved 1.25% on its own.
- **Per band.** Each band of the zigzag scan is then coded or zeroed. A band flag's model depends on whether the previous band is coded, so the choices for all 16 bands are made at once, as a shortest path through two states per band (Viterbi). This follows J. Wen, M. Luttrell and J. Villasenor, "Trellis-based R-D optimal quantization in H.263+" (IEEE TIP 9(8), 2000). The DC level's zero or nonzero choice is the path's first stage.
- **Rounding.** The first candidate is rounded with an offset of step/6, the same offset the decoder adds back. Changing the rounding offset together with the decoder's reconstruction offset (to step/8, /10 or /14) made results worse (+0.9% to +3.0%): the step/6 rounding acts as a useful dead zone.
- **One special case:** when every earlier level of a band is zero, the last level is known to be nonzero (its zero decision isn't coded), so zero isn't a candidate for it.

**Chroma** (`chroma_rd`) uses the same choice with its own lambda, `CHROMA_LAMBDA` = 0.25 x the luma lambda. A sweep, measured before retraining the starting tables:

| Chroma lambda | Luma | YUV |
|---------------|------|-----|
| 0 (plain rounding) | -2.28% | -2.31% |
| 0.25 x luma | -3.53% | -2.61% |
| 0.5 x luma | -4.92% | -2.48% |
| 1.0 x luma | -6.53% | -0.47% |

Above 0.25 the luma gain is bought with chroma quality (section 6).

**Carries over:** yes; x264's trellis is based on the same 2000 paper. In H.264 the dependency between levels comes from the significance map and the last-coefficient flag rather than band flags, but the method is the same. Later RDOQ designs from the HEVC era have patents of their own: follow the 1994 and 2000 papers.

## 4. Chroma in the mode decisions

`chroma_cost` adds each candidate's two 4x4 chroma blocks to its cost, costed as the chroma pass will code them (`chroma_rd`):
- skip: the squared error of the chroma prediction;
- inter: the RD-coded residual;
- intra: pixels minus 128, with the DC prediction estimated from the source, since chroma is coded after all luma.

It saved 0.3% on luma and 0.6% on YUV. A chroma weight of 1.0 (plain squared error) was best:
- 0.5 and 0.75 traded YUV for luma;
- 2.0 and 4.0 made the modes follow chroma at luma's expense (+7% at 4.0).

**A lesson in order.** Tried before chroma had a coded flag, this lost SSIM. At that time every non-skipped block paid 16 chroma symbols, so chroma's cost pushed blocks toward skip and lost luma detail. Once chroma had a coded flag, the same decision helped.

**Carries over:** yes. x264 also counts chroma in its RD decisions.

## 5. QP across frames

Each frame header carries its QP, so the encoder can choose it per frame without any format change.

- **Key frames 5 steps finer, B-frames 5 coarser** (`KEY_QP_OFFSET`, `B_QP_OFFSET`). Every frame of a GOP predicts from the key frame, so its quality pays off many times. Nothing predicts from a B-frame. Both offsets were 4 in the first round, when the finer key frame alone saved 6.4%. Re-tuned for the second round's tools, 5 and 5 saved another 0.8% (YUV 1.6%, SSIM 4.0%). Key 6 gave 1% more on SSIM but lost 0.1% on luma.
- **Coarser late anchors** (`gop_qp_steps`). An anchor's errors spread to the frames after it in its closed GOP, so anchors near the end of a GOP are coded coarser: +1 QP with 8 or fewer frames left, +2 with 2 or fewer, +3 for the last frame. This is dependent quantization (K. Ramchandran, A. Ortega and M. Vetterli, "Bit allocation for dependent quantization with applications to multiresolution and MPEG video coders", IEEE TIP 3(5), 1994). It saved 1.25% (YUV 1.4%, SSIM 2.0%) and needs the caller's key-frame interval (`SntEncoderCtx.keyint`). Finer early anchors were worse (+0.4%). A patent search found only bit-budget rate control by remaining frames (EXPERIMENTS.md, item 13a).

**Not kept: a QP per macroblock by look-ahead.** Key-frame macroblocks whose area stays the same through the GOP were coded finer. It saved 0.9% on luma and 3.5% on SSIM, but it is close to active patents on look-ahead adaptive quantization (MB-tree and similar), and no pre-2004 source was found for it.

**Carries over:** x264's `ipratio` and `pbratio` do what the key-frame and B-frame offsets do. The late-anchor steps apply to any closed-GOP encode. Comparisons against x264 must match its mode: at a fixed QP, x264 turns MB-tree and adaptive quantization off.

## 6. Global motion estimation

From the first round (`core/prediction/global_motion.c`):
- One translation per frame and reference, by Lucas-Kanade (1981): Gauss-Newton steps on a three-level image pyramid.
- Pixels whose residual is over 2.5x the mean are left out, so moving foreground doesn't pull the estimate.
- The encoder then picks the quarter-sample vector near the estimate that minimizes the sum over 8x8 blocks of min(SSE, 30 lambda): the vector that lets the most blocks skip, not the best average match.

With the format's global-motion skip, it saved 2.1%. After the second round only 0.3–1% of skipped blocks still use it; the rest repeat the left block's motion.

**Carries over:** to formats with global motion (MPEG-4 Part 2's GMC). H.264 has no global motion.

## 7. Speed

The second round's encoder kept within its budget of 650 ms per CIF frame (single-threaded) through three changes, each found by profiling first:

| Change | Output | Slowest frame |
|--------|--------|---------------|
| SAD summed row by row in the integer search, stopping once the cost can't win | Identical | 628 -> 484 ms, with the next change |
| `log2` of counts from a table in `cdf_bits` | Within 0.04% | (together with the previous change) |
| The last 16 plain 8x8 predictions cached by block and motion (`cached_pred8`); OBMC's neighbour predictions are the same for every candidate of a block | Identical | 623 -> 558 ms |

The encoder now takes 250–600 ms per CIF frame. Trial encoding is the main cost; a real-time encoder would need fewer candidates and SIMD.

## 8. Measuring

Most of what this round learned about method:

- **One harness.** `tools/eval/run.sh` builds a commit or the working tree and runs the tuning points; `tools/eval/bd.py` gives the BD-rates. The harness was checked against `sinterenc`, `sinterdec` and ffmpeg (same sizes, PSNR and SSIM).
- **Luma and YUV together.** A parameter that moves bits between luma and chroma (lambda, chroma lambda, chroma QP) shows a luma gain with a YUV loss, or the reverse. That's a trade. Keep a change only if both improve.
- **Refresh the baseline after every commit.** Once a baseline was a commit behind, and OBMC was first reported at -3.8% instead of its real -2.0%.
- **Like with like for new models.** When a change adds models, compare two builds without trained tables, or retrain the tables for the new build. Gains often grow with retrained tables (more contexts need trained starting points) and sometimes shrink.
- **Check everything else too:**
  - encoder and decoder reconstructions match at every test point;
  - stills checked by eye for every format change;
  - held-out clips measured only at the end;
  - tools that adapt the QP also compared against x264 in CRF mode.
- **Measure the upper bound first.** Before building a large tool, count what it could save. Affine global motion and 32x32 skip blocks were dropped this way without being built: under 1% of skipped blocks use the global motion, and the skip runs already cost 0.15 bits per 2x2 macroblock group.

## Sources

- K. Ramchandran and M. Vetterli, "Rate-distortion optimal fast thresholding with complete JPEG/MPEG decoder compatibility", IEEE TIP 3(5), 1994
- K. Ramchandran, A. Ortega and M. Vetterli, "Bit allocation for dependent quantization with applications to multiresolution and MPEG video coders", IEEE TIP 3(5), 1994
- S.-W. Wu and A. Gersho, "Joint estimation of forward and backward motion vectors for interpolative prediction of video", IEEE TIP 3(5), 1994
- G. J. Sullivan and T. Wiegand, "Rate-distortion optimization for video compression", IEEE Signal Processing Magazine, 1998
- J. Wen, M. Luttrell and J. Villasenor, "Trellis-based R-D optimal quantization in H.263+", IEEE TIP 9(8), 2000
- B. D. Lucas and T. Kanade, "An iterative image registration technique with an application to stereo vision", IJCAI 1981
