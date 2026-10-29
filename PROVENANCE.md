# Provenance of Sinter's format tools

A format tool is anything the decoder must run to reproduce the output: syntax, entropy coding, prediction, transforms, filters and reconstruction. Each one must either come from a publication dated before 2004, implemented as published with small listed changes, or be new to Sinter, with a record of the search for the same idea in patents and papers. Nothing comes from H.264/AVC, its drafts, or any later standard or codec (VP8 onward, Daala, AV1).

Encoder-only methods (motion search, mode decisions, rate control, quantizer choices) don't affect the format and aren't listed.

This file records where each tool comes from. It is not a patent clearance.

| Tool | Source | Date | Changes from the source |
|------|--------|------|-------------------------|
| Range coder | G. N. N. Martin, "Range encoding: an algorithm for removing redundancy from a digitised message", Video & Data Recording Conference, Southampton | 1979 | Bytes leave the top of a 32-bit range. A carry is resolved by holding back bytes that it could still change (the byte analogue of Witten, Neal & Cleary's bits_to_follow). The final value is chosen to end in zero bytes, which are not written. |
| Adaptive symbol models | I. H. Witten, R. M. Neal and J. G. Cleary, "Arithmetic coding for data compression", Communications of the ACM 30(6) | 1987 | A count per symbol, incremented on use, with all counts halved when their total passes a limit, as published. The starting count, increment and limit are set per alphabet (the paper uses 1, 1 and 2^14 - 1): 1, 128 and 65000 for coefficient-like alphabets, and 96, 96 and 8192 for mode flags. |
| Luma interpolation (quarter-pel motion compensation) | M. Unser, A. Aldroubi and M. Eden, "Fast B-spline transforms for continuous image representation and interpolation", IEEE Trans. PAMI 13(3); and "B-spline signal processing, Part II", IEEE Trans. Signal Processing 41(2) | 1991, 1993 | Cubic B-spline with the published recursive prefilter (pole sqrt(3) - 2, mirror boundaries), in fixed point: constants in Q16, coefficients in Q8. The frame is first extended by repeating its edge samples 32 pixels on each side. Whole-sample positions copy the reference, which equals the exact spline there. |
| Skipped blocks | ITU-T H.263, "Video coding for low bit rate communication" (COD = 1: a macroblock not coded, copied from the reference at zero motion) | 1996 | Applies to 8x8 blocks and 16x16 macroblocks. In B-frames a skipped block is the rounded average of both references at zero motion. A skipped block counts as zero motion for its neighbours' motion vector prediction, as in H.263. |

`libsinter/core/entropy/rans.c` is no longer part of the format or the build.
