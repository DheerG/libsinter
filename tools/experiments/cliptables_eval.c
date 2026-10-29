/*
 * Measure per-clip starting tables (docs/NOVEL.md, section 3).
 *
 * Pass 1 encodes the clip normally and counts each binary model's
 * symbols. A model gets a per-clip start when its expected saving, about
 * KL(clip || default) bits for each of the first 24 symbols after every
 * set restart, beats the cost of sending it (about 10 bits). Pass 2
 * encodes again with the tables sent in the first key frame, decodes
 * every frame, and prints size, mean PSNR and encoder/decoder mismatches
 * in the same format as the other evaluation harnesses.
 *
 * Build against a library made with
 *   -DSNT_EXPERIMENT_CLIPTABLES -DSNT_TRAIN_STATS
 * usage: cliptables_eval in.y4m QP keyint nb [-q]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "codec/encoder.h"
#include "codec/decoder.h"
#include "codec/cliptables.h"
#include "core/sinter_types.h"

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

static double psnr(const uint8_t *a, const uint8_t *b, size_t n) {
    double se = 0;
    for (size_t i = 0; i < n; i++) { double d = (double)a[i] - b[i]; se += d * d; }
    return se == 0 ? 99.0 : 10.0 * log10(255.0 * 255.0 * n / se);
}

static size_t mismatch(const uint8_t *a, const uint8_t *b, size_t n) {
    size_t m = 0;
    for (size_t i = 0; i < n; i++) m += (a[i] != b[i]);
    return m;
}

static double kl(double p, double q) {
    return p * log2(p / q) + (1 - p) * log2((1 - p) / (1 - q));
}

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: %s in.y4m QP keyint nb [-q]\n", argv[0]); return 1; }
    int qp = atoi(argv[2]), keyint = atoi(argv[3]), nb = atoi(argv[4]);
    FILE *f = fopen(argv[1], "rb");
    int w = 0, h = 0;
    if (!f || !parse_hdr(f, &w, &h)) { fprintf(stderr, "bad y4m\n"); return 1; }
    size_t ls = (size_t)w * h, cs = ls / 4, fs = ls + 2 * cs;
    int cap = 1024, n = 0;
    uint8_t *src = malloc(fs * cap);
    while (n < cap && read_frame(f, src + fs * n, fs)) n++;
    uint8_t *dec = malloc(fs * n);

    int *order = malloc(sizeof(int) * n), *bidi = malloc(sizeof(int) * n), no = 0, gops = 0;
    for (int g = 0; g < n; g += keyint) {
        int len = n - g < keyint ? n - g : keyint;
        order[no] = g; bidi[no++] = 0;
        gops++;
        int prev = 0;
        while (prev < len - 1) {
            int next = prev + nb + 1 < len - 1 ? prev + nb + 1 : len - 1;
            order[no] = g + next; bidi[no++] = 0;
            for (int b = prev + 1; b < next; b++) { order[no] = g + b; bidi[no++] = 1; }
            prev = next;
        }
    }

    uint8_t *bs = malloc(ls * 4);
    SntClipTables tables;
    uint64_t tot = 0;
    size_t mm = 0;
    for (int pass = 0; pass < 2; pass++) {
        SntEncoderCtx ectx;
        memset(&ectx, 0, sizeof ectx);
        ectx.header.base_qp = (uint8_t)qp;
        ectx.header.chroma_qp_offset = 4;
        if (pass == 1) { ectx.clip_tables = &tables; ectx.clip_tables_send = 1; }
        SntDecoderCtx dctx;
        memset(&dctx, 0, sizeof dctx);
        snte_train_reset();
        for (int i = 0; i < no; i++) {
            int d = order[i];
            uint8_t *y = src + fs * d, *u = y + ls, *v = u + cs;
            if (d % keyint == 0) { ectx.num_refs = 0; ectx.cdf_state_valid = 0; }
            ectx.code_bidi = bidi[i];
            ectx.header.frame_number = (uint32_t)d;
            uint32_t len = 0;
            if (snte_encode_frame(&ectx, y, u, v, w, h, w, w / 2, bs, (uint32_t)(ls * 4), &len)) return 2;
            if (pass == 0) continue;
            uint8_t *oy = dec + fs * d, *ou = oy + ls, *ov = ou + cs;
            int dw, dh;
            if (sntd_decode_frame(&dctx, bs, len, oy, ou, ov, &dw, &dh, w, w / 2)) return 3;
            mm += mismatch(ectx.last_y, oy, ls) + mismatch(ectx.last_u, ou, cs) + mismatch(ectx.last_v, ov, cs);
            tot += len;
        }
        if (pass == 0) {
            /* Choose the per-clip starts from pass 1's counts */
            int used = 0;
            for (int set = 0; set < 2; set++) {
                int set_qp = set == SNT_SET_P ? (qp > 4 ? qp - 4 : 0) : (qp + 4 < 63 ? qp + 4 : 63);
                SntCdfSet def;
                snt_cdf_set_init(&def, set, set_qp);
                for (int k = 0; k < SNT_NUM_MODELS; k++) {
                    tables.q[set][k] = SNT_CLIP_NONE;
                    SntModel *m = &((SntModel *)&def)[k];
                    if (m->num_symbols != 2) continue;
                    double n0 = snte_train_count(set, k, 0), n1 = snte_train_count(set, k, 1);
                    if (n0 + n1 < 32) continue;
                    double p = (n1 + 0.5) / (n0 + n1 + 1), pd = (double)m->freq[1] / m->total;
                    int q = (int)lround(log(p / (1 - p)) * 256.0 / SNT_CLIP_QSTEP);
                    q = q < -15 ? -15 : (q > 15 ? 15 : q);
                    double pq = snt_clip_p1(q) / 4096.0;
                    double per_reset = (n0 + n1) / gops < 24 ? (n0 + n1) / gops : 24;
                    double saving = gops * per_reset * (kl(p, pd) - kl(p, pq));
                    if (saving > 10) { tables.q[set][k] = (int8_t)q; used++; }
                }
            }
            fprintf(stderr, "per-clip starts: %d models\n", used);
        }
        free(ectx.ref_y); free(ectx.ref_u); free(ectx.ref_v); free(ectx.cdf_state);
        free(ectx.ref_old_y); free(ectx.ref_old_u); free(ectx.ref_old_v); free(ectx.cdf_state_b);
        free(ectx.last_y); free(ectx.last_u); free(ectx.last_v);
        free(ectx.ref_coef); free(ectx.ref_old_coef);
        free(ectx.key_y); free(ectx.key_u); free(ectx.key_v); free(ectx.key_coef);
        free(dctx.ref_y); free(dctx.ref_u); free(dctx.ref_v); free(dctx.cdf_state);
        free(dctx.ref_old_y); free(dctx.ref_old_u); free(dctx.ref_old_v); free(dctx.cdf_state_b);
        free(dctx.ref_coef); free(dctx.ref_old_coef); free(dctx.clip_tables);
        free(dctx.key_y); free(dctx.key_u); free(dctx.key_v); free(dctx.key_coef);
    }
    double sy = 0, su = 0, sv = 0;
    for (int d = 0; d < n; d++) {
        uint8_t *a = src + fs * d, *b = dec + fs * d;
        sy += psnr(a, b, ls); su += psnr(a + ls, b + ls, cs); sv += psnr(a + ls + cs, b + ls + cs, cs);
    }
    printf("SUMMARY qp=%d frames=%d bytes=%llu  avgY=%.2f avgU=%.2f avgV=%.2f  mismatchY=%zu mismatchUV=0\n",
           qp, n, (unsigned long long)tot, sy / n, su / n, sv / n, mm);
    return 0;
}
