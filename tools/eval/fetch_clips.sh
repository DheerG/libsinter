#!/bin/sh
# Download the Xiph CIF test clips into clips/, first 150 frames.
set -e
root=$(cd "$(dirname "$0")/../.." && pwd)
mkdir -p "$root/clips"
for c in ${CLIPS:-foreman mobile akiyo bus news}; do
  [ -f "$root/clips/$c.y4m" ] && continue
  curl -sfL "https://media.xiph.org/video/derf/y4m/${c}_cif.y4m" -o "$root/clips/$c.full.y4m"
  ffmpeg -loglevel error -y -i "$root/clips/$c.full.y4m" -frames:v 150 -f yuv4mpegpipe "$root/clips/$c.y4m"
  rm "$root/clips/$c.full.y4m"
done
