#!/bin/sh
# Regenerate libsinter/codec/model_tables.h: collect symbol counts from
# 60 frames of each training clip at seven QPs, in parallel, then write
# the tables (README.md gives the steps one by one).
#
# usage: tools/train/train.sh [CLIPDIR]
#   CLIPDIR (default clips/train) holds the ten training clips as
#   NAME_cif.y4m. Never put foreman, mobile, akiyo, bus or news there.
set -e
root=$(cd "$(dirname "$0")/../.." && pwd)
dir=${1:-$root/clips/train}
out=$root/build/train
flags="-std=c11 -O2 -I libsinter -DSNT_TRAIN_STATS -DSNT_NO_TRAINED_TABLES"
make -s -C "$root" BUILD="$out" CFLAGS="$flags" "$out/libsinter.a" >/dev/null
gcc $flags -I "$root/libsinter" "$root/tools/train/collect_stats.c" -L"$out" -lsinter -lm -o "$out/collect_stats"
rm -f "$out"/*.stats
for clip in "$dir"/*_cif.y4m; do
  case $(basename "$clip") in foreman*|mobile*|akiyo*|bus*|news*) echo "not a training clip: $clip"; exit 1;; esac
  echo "$clip"
done | xargs -P "${JOBS:-10}" -n 1 sh -c '"$0" "$1/$(basename "$2" .y4m).stats" 60 8,12,16,20,24,28,32 "$2" 2>/dev/null' \
  "$out/collect_stats" "$out"
python3 "$root/tools/train/make_tables.py" "$out"/*.stats > "$root/libsinter/codec/model_tables.h"
echo "wrote libsinter/codec/model_tables.h from $(ls "$out"/*.stats | wc -l | tr -d ' ') clips"
