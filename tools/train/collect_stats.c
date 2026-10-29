/*
 * Collect symbol statistics for the trained starting tables.
 *
 * Encodes each Y4M clip at each QP, in the same GOP structure as
 * sinterenc (key frame every KEYINT frames, closed GOPs, BFRAMES B-frames
 * between anchors), with a library built with -DSNT_TRAIN_STATS, then
 * writes the symbol counts of every model: one line per model set, QP
 * class and model, with 256 counts.
 *
 * usage: collect_stats OUT.txt FRAMES QP[,QP...] CLIP.y4m...
 *
 * Train only on clips that are not used to evaluate the codec.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "codec/encoder.h"
#include "core/sinter_types.h"

#define KEYINT  30
#define BFRAMES 2

static int parse_hdr(FILE *f, int *w, int *h) {
    char line[512];
    int pos = 0, c;
    while ((c = fgetc(f)) != '\n' && c != EOF && pos < 511) line[pos++] = (char)c;
    line[pos] = 0;
    if (strncmp(line, "YUV4MPEG2", 9)) return 0;
    for (char *p = line + 9; *p; p++) {
        if (p[-1] == ' ' && *p == 'W') *w = atoi(p + 1);
        if (p[-1] == ' ' && *p == 'H') *h = atoi(p + 1);
    }
    return *w > 0 && *h > 0;
}

static int read_frame(FILE *f, uint8_t *buf, size_t fs) {
    int c = fgetc(f);
    if (c != 'F') return 0;
    while ((c = fgetc(f)) != '\n' && c != EOF);
    return fread(buf, 1, fs, f) == fs;
}

static void free_ctx(SntEncoderCtx *c) {
    free(c->ref_y); free(c->ref_u); free(c->ref_v); free(c->cdf_state);
    free(c->ref_old_y); free(c->ref_old_u); free(c->ref_old_v); free(c->cdf_state_b);
    free(c->last_y); free(c->last_u); free(c->last_v);
    free(c->ref_coef); free(c->ref_old_coef);
    free(c->key_y); free(c->key_u); free(c->key_v); free(c->key_coef);
    free(c->anchor_mv); free(c->anchor_dir);
}

static void encode_clip(const char *path, int max_frames, int qp) {
    FILE *f = fopen(path, "rb");
    int w = 0, h = 0;
    if (!f || !parse_hdr(f, &w, &h)) { fprintf(stderr, "bad y4m %s\n", path); exit(1); }
    size_t ls = (size_t)w * h, cs = ls / 4, fs = ls + 2 * cs;
    uint8_t *src = malloc(fs * max_frames);
    int n = 0;
    while (n < max_frames && read_frame(f, src + fs * n, fs)) n++;
    fclose(f);

    /* Coding order: per GOP the key frame, then each anchor followed by
     * the B-frames before it */
    int *order = malloc(sizeof(int) * n), *bidi = malloc(sizeof(int) * n), no = 0;
    for (int g = 0; g < n; g += KEYINT) {
        int len = n - g < KEYINT ? n - g : KEYINT;
        order[no] = g; bidi[no++] = 0;
        int prev = 0;
        while (prev < len - 1) {
            int next = prev + BFRAMES + 1 < len - 1 ? prev + BFRAMES + 1 : len - 1;
            order[no] = g + next; bidi[no++] = 0;
            for (int b = prev + 1; b < next; b++) { order[no] = g + b; bidi[no++] = 1; }
            prev = next;
        }
    }

    SntEncoderCtx ctx;
    memset(&ctx, 0, sizeof ctx);
    ctx.header.base_qp = (uint8_t)qp;
    ctx.header.num_tile_cols = 1;
    ctx.header.num_tile_rows = 1;
    ctx.header.chroma_qp_offset = 4;
    ctx.keyint = KEYINT;
    uint8_t *bs = malloc(ls * 4);
    for (int i = 0; i < no; i++) {
        int d = order[i];
        uint8_t *y = src + fs * d, *u = y + ls, *v = u + cs;
        if (d % KEYINT == 0) { ctx.num_refs = 0; ctx.cdf_state_valid = 0; }
        ctx.code_bidi = bidi[i];
        ctx.header.frame_number = (uint32_t)d;
        uint32_t len = 0;
        if (snte_encode_frame(&ctx, y, u, v, w, h, w, w / 2, bs, (uint32_t)(ls * 4), &len)) {
            fprintf(stderr, "encode failed\n");
            exit(1);
        }
    }
    free_ctx(&ctx);
    free(bs); free(src); free(order); free(bidi);
}

int main(int argc, char **argv) {
    if (argc < 5) {
        fprintf(stderr, "usage: %s OUT.txt FRAMES QP[,QP...] CLIP.y4m...\n", argv[0]);
        return 1;
    }
    int frames = atoi(argv[2]);
    for (int c = 4; c < argc; c++) {
        char qps[256];
        snprintf(qps, sizeof qps, "%s", argv[3]);
        for (char *t = strtok(qps, ","); t; t = strtok(NULL, ",")) {
            fprintf(stderr, "%s QP %s\n", argv[c], t);
            encode_clip(argv[c], frames, atoi(t));
        }
    }
    FILE *out = fopen(argv[1], "w");
    if (!out) return 1;
    snte_train_dump(out);
    fclose(out);
    return 0;
}
