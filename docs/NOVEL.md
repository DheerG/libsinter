# Novel tools

Phase 3 of the clean-format work tried tools of Sinter's own design. Each entry gives the idea, why it should help, the search for the same idea in patents and papers, and the measured result. A new tool is kept only if the search finds no patent or application with the same idea. The searches are checks for obvious collisions, not patent clearance.

BD-rates are luma (and YUV), mean of foreman, mobile and akiyo (Xiph CIF, 60 frames, 2 B-frames), against the codec just before the change.

## 1. A decoder that learns while decoding (dropped after the search)

**Idea.** Code each binary decision with several context models at once (coarse and fast, fine and slow) and combine their probabilities by logistic mixing, with the mixing weights learned online from the decoded data (as in the PAQ compressors). A new context could then be added without diluting the statistics: the mixer learns how much to trust it.

**Why it should help.** Richer contexts lost before because each one relearns slowly after a key frame; mixing lets a fast general model cover for slow specific ones.

**Search.**
- Query: "patent video coding entropy context mixing logistic mixing weights updated online probability estimation transform coefficients". Found Google's [EP3677027A1](https://patents.google.com/patent/EP3677027A1/en) (granted as EP3677027B1) and [US10887595B2](https://patents.google.com/patent/US10887595B2/en), "Mixing of probabilities for entropy coding in video compression" / "Using multiple probability models for entropy coding in video compression" (priority 2017-08-29, active): two or more probability models for the tokens of quantized transform coefficients, combined into a mixed probability with weights updated online from the code length each model would have produced. Also [US10735736](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/10735736), "Selective mixing for entropy coding in video compression", and [US11259053B2](https://patents.google.com/patent/US11259053) (machine-learned probability refinement).

**Result.** Dropped before implementation: adaptive mixing of several models' probabilities for coefficient coding is the same idea, although these patents mix linearly rather than in the logistic domain. The prior-art route would be Context Tree Weighting (F. Willems, Y. Shtarkov and T. Tjalkens, IEEE Trans. Information Theory 41(3), 1995), which weights context models by their accumulated code lengths and predates these patents by 22 years. But that weighting is essentially the patents' mechanism, so building it would follow the letter of the source rule and not its purpose. It is left for a patent professional to judge.

## 2. Contexts from the prediction (kept)

**Idea.** Before an inter block's residual is decoded, the decoder already has the block's motion-compensated prediction. Take its DCT and compare each coefficient's magnitude with the quantizer step (classes: below 1 step, below 4 steps, above), and each band's summed magnitude (below 2 steps, below 8, above). These classes select the models for the band's coded flag and for each level's zero and "above one" decisions (`levels.h`). Intra blocks use class 0. Nothing extra is sent.

**Why it should help.** A residual is roughly the prediction's misalignment times its gradient, so it tends to be nonzero at the frequencies where the prediction has energy, and zero where the prediction is flat.

**Search.**
- "patent entropy coding residual coefficients context derived from motion compensated prediction block energy frequency": closest [US7974343B2](https://patents.google.com/patent/US7974343/zh-cn) (code tables chosen by macroblock type), [US20160295242A1](https://patents.google.com/patent/US20160295242), [US20100046626A1](https://patents.google.com/patent/US20100046626A1/en). None derives contexts from the prediction.
- "patent selecting context for coding transform coefficients of residual based on prediction block characteristics ... variance gradient": closest [US6680974B1](https://patents.google.com/patent/US6680974) (contexts from previously reconstructed neighbouring coefficients; expired), [US8891630B2](https://ipverse.greyb.com/patents/US8891630B2) (partition sets).
- '"prediction signal" OR "predicted block" used to derive context for entropy coding of residual transform coefficients': closest [US10397577](https://patents.justia.com/patent/10397577) (contexts by prediction mode and transform type, not the prediction's content), [WO2019029951A1](https://patents.google.com/patent/WO2019029951A1/en), [US9621894B2](https://patents.google.com/patent/US9621894).

No patent or application found that selects contexts from the content or spectrum of the prediction block.

**Result.** Kept. -0.9% (-0.7%): foreman -0.9%, mobile -1.2%, akiyo -0.5%. Lowering or raising the class thresholds by 2x was no better.

## 3. Per-clip starting tables (dropped after the search)

**Idea.** A first encoder pass collects the clip's symbol statistics; the sequence or key-frame header then sends starting probability tables as small changes from the trained defaults.

**Search.**
- "patent video encoder transmit initial probability tables context initialization in bitstream header per sequence first pass statistics arithmetic coding": closest [US20260122287A1](https://patents.google.com/patent/US20260122287A1/en) (frame context initialization from averaged tile contexts), [US8344917B2](https://patents.google.com/patent/US8344917B2/en) and [US20170250709A1](https://patents.google.com/patent/US20170250709) (initialization from QP), [WO2016196307A1](https://patents.google.com/patent/WO2016196307A1/en).
- "VP8 coefficient probability updates frame header": VP8 ([RFC 6386](https://datatracker.ietf.org/doc/rfc6386/), 2010) sends exactly this: default token probabilities plus per-frame updates in the frame header, which persist until the next key frame.

**Result.** Dropped: sending probability changes from defaults in the header is VP8's design, and VP8 is an excluded source. Choosing one of a few predefined tables instead would be H.264's `cabac_init_idc`, also excluded. JPEG's two-pass optimized tables (1992) would be the pre-2004 form, but sending complete tables for hundreds of adaptive models costs more than the trained starting tables leave to gain.

## 4. Prediction agreement for mode flags (kept)

Chosen from a measured breakdown of where the bits go: at QP 28 the macroblock mode and the per-block inter and skip flags took 6-16% of the bits (foreman 15.8%, akiyo 8%), coded with models that only knew the history of those flags.

**Idea.** Before a macroblock's mode (or a block's inter and skip flags) is read, the decoder can already form two predictions from the references alone. In P-frames: the skip prediction (the motion a skipped block would use) and the prediction at the median motion vector predictor. In B-frames: the forward prediction at the forward predictor and the backward prediction at the backward predictor. Their mean absolute difference, against the quantizer step (below step/8, below step/2, above), selects the mode models. Nothing extra is sent.

**Why it should help.** When the two predictions agree, the motion around the block is consistent (in B-frames, consistent in time from both sides), so skipping is likely to be good enough; when they disagree, the block probably has its own motion and needs coding.

**Search.**
- "patent skip flag context derived from difference between two motion compensated predictions decoder side similarity forward backward prediction entropy coding mode": closest [US 11,310,522](https://patents.justia.com/patent/11310522) (ByteDance, decoder-side motion vector refinement: the SAD between the two bi-predictions decides whether to skip the refinement, not a context for coding), [US9420300B2](https://patents.google.com/patent/US9420300), [US9451270B2](https://patents.google.com/patent/US9451270B2/en), [JP2006518564A](https://patents.google.com/patent/JP2006518564A/en).
- "entropy coding mode decision probability estimated from agreement of candidate predictions computed at decoder context macroblock skip": closest [US8520732](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/8520732) and [US8428133](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/8428133) (contexts from whether neighbours use one or two directions), [US10630978](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/10630978) (decoder reorders intra modes by redundancy), [US9967594B2](https://patents.google.com/patent/US9967594).

No patent or application found that selects the models for mode or skip flags from the agreement of two decoder-side predictions.

**Result.** Kept. -0.6% (-0.5%): foreman -0.6%, mobile -0.2%, akiyo -0.9%. Thresholds twice as coarse (step/4, step) were worse (-0.5%). The decoder forms two 8x8 predictions per flag decision; it still decodes a CIF frame in 3-6 ms.
