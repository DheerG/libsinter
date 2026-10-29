# Training the starting tables

`libsinter/codec/model_tables.h` holds the probabilities each entropy model starts from when its set restarts at a key frame. It is generated from encodes of Xiph CIF clips that are not used to evaluate the codec. Never train on foreman, mobile, akiyo, bus or news.

```bash
# 1. A library that counts symbols, and the collector
make BUILD=build_train CFLAGS="-std=c11 -O2 -I libsinter -DSNT_TRAIN_STATS -DSNT_NO_TRAINED_TABLES" build_train/libsinter.a
gcc -std=c11 -O2 -I libsinter -DSNT_TRAIN_STATS -DSNT_NO_TRAINED_TABLES tools/train/collect_stats.c \
    -Lbuild_train -lsinter -lm -o build_train/collect_stats

# 2. Counts from 60 frames of each training clip at seven QPs (one file per clip)
for c in coastguard container hall_monitor mother_daughter paris silent tempete flower waterfall highway; do
  curl -sO https://media.xiph.org/video/derf/y4m/${c}_cif.y4m
  build_train/collect_stats ${c}_cif.stats 60 8,12,16,20,24,28,32 ${c}_cif.y4m
done

# 3. The tables
python3 tools/train/make_tables.py *_cif.stats > libsinter/codec/model_tables.h
```

The training build starts every model uniform (`SNT_NO_TRAINED_TABLES`), so the counts don't depend on the previous tables. Regenerate the tables whenever the model set in `cdf_set.h` changes (a static assertion in the header catches a mismatch in the number of models) or the encoder's decisions change much.
