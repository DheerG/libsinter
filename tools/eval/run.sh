#!/bin/sh
# Build the working tree (or git REF) and run the tuning points through
# the harness. One line per point in build/eval/NAME.txt: the clip, then
# the harness's fields (harness.c).
#
# usage: tools/eval/run.sh NAME [REF] [CFLAGS]
#   REF: a commit to measure instead of the working tree ("." = working tree)
#   environment: CLIPS (default "foreman mobile akiyo"), QPS (default
#   "8 12 16 20 24 28"), FRAMES (default 60), BITS=1 (stats build; bit
#   breakdowns in build/eval/NAME.bits/, see bits.py), JOBS (default 10)
set -e
name=$1; ref=${2:-.}; flags=$3
[ -n "$name" ] || { echo "usage: $0 NAME [REF] [CFLAGS]"; exit 1; }
root=$(cd "$(dirname "$0")/../.." && pwd)
out=$root/build/eval
mkdir -p "$out"
CLIPS=${CLIPS:-"foreman mobile akiyo"}
QPS=${QPS:-"8 12 16 20 24 28"}
FRAMES=${FRAMES:-60}
JOBS=${JOBS:-10}
[ "$BITS" = 1 ] && flags="$flags -DSNT_TRAIN_STATS"

src=$root
if [ "$ref" != . ]; then
  src=$out/$name.src
  rm -rf "$src" && mkdir -p "$src"
  git -C "$root" archive "$ref" | tar -x -C "$src"
fi
obj=$out/$name.obj
make -s -C "$src" BUILD="$obj" CFLAGS="-std=c11 -O2 -I libsinter $flags" "$obj/libsinter.a" >/dev/null
gcc -std=c11 -O2 -I "$src/libsinter" $flags "$root/tools/eval/harness.c" -L"$obj" -lsinter -lm -o "$obj/harness"

for c in $CLIPS; do
  [ -f "$root/clips/$c.y4m" ] || { echo "missing clips/$c.y4m (tools/eval/fetch_clips.sh)"; exit 1; }
done
rm -f "$out/$name.txt"
[ "$BITS" = 1 ] && rm -rf "$out/$name.bits" && mkdir -p "$out/$name.bits"
for c in $CLIPS; do for q in $QPS; do echo "$c $q"; done; done |
  xargs -P "$JOBS" -L 1 sh -c '
    extra=""; [ "'"$BITS"'" = 1 ] && extra="-bits '"$out/$name.bits"'/${0}_q$1.txt"
    echo "$0 $('"$obj"'/harness '"$root"'/clips/$0.y4m $1 '"$FRAMES"' $extra)"' >> "$out/$name.txt"
sort -k1,1 -k2n "$out/$name.txt" -o "$out/$name.txt"
awk '{m += $12} END {print NR " points, encoder/decoder mismatches: " m}' "$out/$name.txt"
