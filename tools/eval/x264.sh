#!/bin/sh
# Run x264 on the clips and measure its output with the harness's
# measures. Lines in build/eval/NAME.txt as run.sh writes them (the
# frame-type, mismatch and timing fields are 0).
#
# usage: tools/eval/x264.sh NAME [crf]
#   default: -preset medium -tune psnr -g 30 -qp Q (constant QP; x264 turns
#   MB-tree and adaptive quantization off), Q in QPS (default 16 20 ... 44)
#   crf: -crf Q instead, with MB-tree on (for tools that adapt the QP)
#   environment: CLIPS (default "foreman mobile akiyo bus news"), FRAMES (150)
set -e
name=$1; mode=${2:-qp}
root=$(cd "$(dirname "$0")/../.." && pwd)
out=$root/build/eval
mkdir -p "$out/$name.tmp"
CLIPS=${CLIPS:-"foreman mobile akiyo bus news"}
QPS=${QPS:-"16 20 24 28 32 36 40 44"}
FRAMES=${FRAMES:-150}
obj=$out/x264harness.obj
make -s -C "$root" BUILD="$obj" CFLAGS="-std=c11 -O2 -I libsinter" "$obj/libsinter.a" >/dev/null
gcc -std=c11 -O2 -I "$root/libsinter" "$root/tools/eval/harness.c" -L"$obj" -lsinter -lm -o "$obj/harness"
rm -f "$out/$name.txt"
for c in $CLIPS; do for q in $QPS; do echo "$c $q"; done; done |
  xargs -P "${JOBS:-5}" -L 1 sh -c '
    t='"$out/$name.tmp"'/${0}_$1
    ffmpeg -loglevel error -y -i '"$root"'/clips/$0.y4m -frames:v '"$FRAMES"' -c:v libx264 \
      -preset medium -tune psnr -g 30 -'"$mode"' $1 -f h264 $t.h264
    ffmpeg -loglevel error -y -i $t.h264 -pix_fmt yuv420p -f yuv4mpegpipe $t.y4m
    set -- $0 $1 $(stat -f %z $t.h264) $('"$obj"'/harness -compare '"$root"'/clips/$0.y4m $t.y4m '"$FRAMES"')
    rm -f $t.y4m
    echo "$1 $2 '"$FRAMES"' $3 0 0 0 $4 $5 $6 $7 0 0 0"' >> "$out/$name.txt"
sort -k1,1 -k2n "$out/$name.txt" -o "$out/$name.txt"
rm -rf "$out/$name.tmp"
