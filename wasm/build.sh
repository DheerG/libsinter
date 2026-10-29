#!/bin/bash
# Build libsinter decoder as WebAssembly using Emscripten
#
# Prerequisites: install Emscripten SDK (https://emscripten.org/docs/getting_started/)
#   git clone https://github.com/emscripten-core/emsdk.git
#   cd emsdk && ./emsdk install latest && ./emsdk activate latest
#   source ./emsdk_env.sh
#
# Then run this script from the repo root:
#   bash wasm/build.sh
#
# Output: wasm/sinter_decoder.js + wasm/sinter_decoder.wasm

set -e

SRCDIR="libsinter"
OUTDIR="wasm"

# All decoder-side source files (no encoder needed for playback)
DECODER_SRC="
  $SRCDIR/sinter.c
  $SRCDIR/codec/decoder.c
  $SRCDIR/codec/encoder.c
  $SRCDIR/codec/rdo.c
  $SRCDIR/codec/rate_control.c
  $SRCDIR/core/transform/lapping.c
  $SRCDIR/core/transform/dct.c
  $SRCDIR/core/quantize/pvq.c
  $SRCDIR/core/quantize/qp_table.c
  $SRCDIR/core/entropy/rangecoder.c
  $SRCDIR/core/bitstream/bitstream.c
  $SRCDIR/core/prediction/intra.c
  $SRCDIR/core/prediction/obmc.c
  $SRCDIR/core/prediction/motion_est.c
  $SRCDIR/core/prediction/mc.c
  $SRCDIR/core/filter/deblock.c
  $SRCDIR/core/prediction/intra_spatial.c
  $SRCDIR/core/partition/quadtree.c
  $SRCDIR/core/grain/grain_synth.c
  $SRCDIR/core/tile/tile.c
  $SRCDIR/perceptual/csf_weights.c
  $SRCDIR/perceptual/activity_mask.c
  $SRCDIR/perceptual/ptd.c
  $SRCDIR/perceptual/grain_analysis.c
  $SRCDIR/perceptual/butteraugli_target.c
"

echo "Building libsinter WASM decoder..."

emcc -O2 -std=c11 -I$SRCDIR \
  $DECODER_SRC \
  -s WASM=1 \
  -s EXPORTED_FUNCTIONS='["_snt_decoder_create","_snt_decoder_decode","_snt_decoder_destroy","_malloc","_free"]' \
  -s EXPORTED_RUNTIME_METHODS='["ccall","cwrap","getValue","setValue"]' \
  -s ALLOW_MEMORY_GROWTH=1 \
  -s INITIAL_MEMORY=16777216 \
  -o $OUTDIR/sinter_decoder.js

echo "Built: $OUTDIR/sinter_decoder.js + $OUTDIR/sinter_decoder.wasm"
echo "Size: $(wc -c < $OUTDIR/sinter_decoder.wasm) bytes (WASM)"
