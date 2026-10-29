# Provenance of Sinter's format tools

A format tool is anything the decoder must run to reproduce the output: syntax, entropy coding, prediction, transforms, filters and reconstruction. Each one must either come from a publication dated before 2004, implemented as published with small listed changes, or be new to Sinter, with a record of the search for the same idea in patents and papers. Nothing comes from H.264/AVC, its drafts, or any later standard or codec (VP8 onward, Daala, AV1).

Encoder-only methods (motion search, mode decisions, rate control, quantizer choices) don't affect the format and aren't listed.

This file records where each tool comes from. It is not a patent clearance.

| Tool | Source | Date | Changes from the source |
|------|--------|------|-------------------------|
| Range coder | G. N. N. Martin, "Range encoding: an algorithm for removing redundancy from a digitised message", Video & Data Recording Conference, Southampton | 1979 | Bytes leave the top of a 32-bit range. A carry is resolved by holding back bytes that it could still change (the byte analogue of Witten, Neal & Cleary's bits_to_follow). The final value is chosen to end in zero bytes, which are not written. |
| Adaptive symbol models | I. H. Witten, R. M. Neal and J. G. Cleary, "Arithmetic coding for data compression", Communications of the ACM 30(6) | 1987 | A count per symbol, incremented on use, with all counts halved when their total passes a limit, as published. The starting count, increment and limit are set per alphabet (the paper uses 1, 1 and 2^14 - 1): 1, 128 and 65000 for coefficient-like alphabets, and 96, 96 and 8192 for mode flags. |

`libsinter/core/entropy/rans.c` is no longer part of the format or the build.
