/*
 * Evaluation harness: encode a Y4M clip through the internal API in
 * sinterenc's GOP structure (key frame every KEYINT frames, closed GOPs,
 * BFRAMES B-frames between anchors), decode every frame, and report the
 * size, quality and encoder/decoder agreement on one line:
 *
 *   QP FRAMES BYTES KEY_BYTES P_BYTES B_BYTES Y U V SSIM MISMATCH ENC_MS DEC_MS
 *
 * BYTES is the size of the .sntr file sinterenc would write (container
 * header and frame lengths included). Y, U and V are the mean of per-frame
 * PSNR; SSIM is the mean of per-frame luma SSIM, computed as ffmpeg does
 * (8x8 windows on a 4-pixel grid). MISMATCH counts samples where the
 * encoder's reconstruction differs from the decoder's output. ENC_MS and
 * DEC_MS are milliseconds per frame.
 *
 * usage: harness CLIP.y4m QP FRAMES [-k KEYINT] [-b BFRAMES] [-bits FILE]
 *        harness -compare SOURCE.y4m DECODED.y4m FRAMES
 *
 * -bits (library built with -DSNT_TRAIN_STATS): write the ideal bits of
 * every model, with the model groups of SntCdfSet, for bits.py.
 * -compare: the same quality measures for another codec's output, printed
 * as "Y U V SSIM".
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <time.h>
#include "codec/encoder.h"
#include "codec/decoder.h"
#include "codec/cdf_set.h"
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

/* Up to max_frames frames of a Y4M clip, as consecutive I420 frames */
static uint8_t *read_clip(const char *path, int max_frames, int *w, int *h, int *n) {
    FILE *f = fopen(path, "rb");
    if (!f || !parse_hdr(f, w, h)) { fprintf(stderr, "bad y4m %s\n", path); exit(1); }
    size_t fs = (size_t)*w * *h * 3 / 2;
    uint8_t *buf = malloc(fs * max_frames);
    *n = 0;
    while (*n < max_frames && read_frame(f, buf + fs * *n, fs)) (*n)++;
    fclose(f);
    return buf;
}

static double psnr(const uint8_t *a, const uint8_t *b, size_t n) {
    double se = 0;
    for (size_t i = 0; i < n; i++) { double d = (double)a[i] - b[i]; se += d * d; }
    return se == 0 ? 99.0 : 10.0 * log10(255.0 * 255.0 * n / se);
}

/* Luma SSIM as ffmpeg's ssim filter computes it: sums over 4x4 blocks,
 * combined into 8x8 windows at every 4-pixel step */
static double ssim(const uint8_t *a, const uint8_t *b, int w, int h) {
    int bw = w / 4, bh = h / 4;
    int (*s)[4] = malloc(sizeof(int[4]) * bw * bh);
    for (int by = 0; by < bh; by++)
        for (int bx = 0; bx < bw; bx++) {
            int s1 = 0, s2 = 0, ss = 0, s12 = 0;
            for (int y = 0; y < 4; y++)
                for (int x = 0; x < 4; x++) {
                    int p = a[(by * 4 + y) * w + bx * 4 + x], q = b[(by * 4 + y) * w + bx * 4 + x];
                    s1 += p; s2 += q; ss += p * p + q * q; s12 += p * q;
                }
            int *o = s[by * bw + bx];
            o[0] = s1; o[1] = s2; o[2] = ss; o[3] = s12;
        }
    const double c1 = .01 * .01 * 255 * 255 * 64, c2 = .03 * .03 * 255 * 255 * 64 * 63;
    double sum = 0;
    for (int by = 0; by + 1 < bh; by++)
        for (int bx = 0; bx + 1 < bw; bx++) {
            double t[4] = {0, 0, 0, 0};
            for (int k = 0; k < 4; k++)
                for (int j = 0; j < 4; j++)
                    t[j] += s[(by + (k >> 1)) * bw + bx + (k & 1)][j];
            double vars = t[2] * 64 - t[0] * t[0] - t[1] * t[1];
            double covar = t[3] * 64 - t[0] * t[1];
            sum += (2 * t[0] * t[1] + c1) * (2 * covar + c2)
                 / ((t[0] * t[0] + t[1] * t[1] + c1) * (vars + c2));
        }
    free(s);
    return sum / ((bw - 1) * (bh - 1));
}

static size_t mismatch(const uint8_t *a, const uint8_t *b, size_t n) {
    size_t m = 0;
    for (size_t i = 0; i < n; i++) m += a[i] != b[i];
    return m;
}

static double now_ms(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec * 1e3 + t.tv_nsec / 1e6;
}

static int compare(const char *src_path, const char *dec_path, int max_frames) {
    int w, h, n, w2, h2, n2;
    uint8_t *src = read_clip(src_path, max_frames, &w, &h, &n);
    uint8_t *dec = read_clip(dec_path, max_frames, &w2, &h2, &n2);
    if (w != w2 || h != h2 || n != n2) { fprintf(stderr, "clips differ\n"); return 1; }
    size_t ls = (size_t)w * h, cs = ls / 4, fs = ls + 2 * cs;
    double q[4] = {0, 0, 0, 0};
    for (int i = 0; i < n; i++) {
        const uint8_t *a = src + fs * i, *b = dec + fs * i;
        q[0] += psnr(a, b, ls);
        q[1] += psnr(a + ls, b + ls, cs);
        q[2] += psnr(a + ls + cs, b + ls + cs, cs);
        q[3] += ssim(a, b, w, h);
    }
    printf("%.4f %.4f %.4f %.6f\n", q[0] / n, q[1] / n, q[2] / n, q[3] / n);
    return 0;
}

int main(int argc, char **argv) {
    if (argc == 5 && !strcmp(argv[1], "-compare")) return compare(argv[2], argv[3], atoi(argv[4]));
    if (argc < 4) {
        fprintf(stderr, "usage: %s CLIP.y4m QP FRAMES [-k KEYINT] [-b BFRAMES] [-bits FILE]\n"
                        "       %s -compare SOURCE.y4m DECODED.y4m FRAMES\n", argv[0], argv[0]);
        return 1;
    }
    int qp = atoi(argv[2]), keyint = 30, bframes = 2;
    const char *bits_path = NULL;
    for (int i = 4; i + 1 < argc; i += 2) {
        if (!strcmp(argv[i], "-k")) keyint = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "-b")) bframes = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "-bits")) bits_path = argv[i + 1];
    }
    int w, h, n;
    uint8_t *src = read_clip(argv[1], atoi(argv[3]), &w, &h, &n);
    size_t ls = (size_t)w * h, cs = ls / 4, fs = ls + 2 * cs;

    /* Coding order, as sinterenc: per GOP the key frame, then each anchor
     * followed by the B-frames before it */
    int *order = malloc(sizeof(int) * n), *bidi = malloc(sizeof(int) * n), no = 0;
    for (int g = 0; g < n; g += keyint) {
        int len = n - g < keyint ? n - g : keyint;
        order[no] = g; bidi[no++] = 0;
        int prev = 0;
        while (prev < len - 1) {
            int next = prev + bframes + 1 < len - 1 ? prev + bframes + 1 : len - 1;
            order[no] = g + next; bidi[no++] = 0;
            for (int b = prev + 1; b < next; b++) { order[no] = g + b; bidi[no++] = 1; }
            prev = next;
        }
    }

    SntEncoderCtx ectx;
    memset(&ectx, 0, sizeof ectx);
    ectx.header.base_qp = (uint8_t)qp;
    ectx.header.dual_rans = 1;
    ectx.header.num_tile_cols = 1;
    ectx.header.num_tile_rows = 1;
    ectx.header.chroma_qp_offset = 4;
    ectx.psy_cfg.alpha = 0.5f;
    ectx.psy_cfg.qp_clamp = 6;
    SntDecoderCtx dctx;
    memset(&dctx, 0, sizeof dctx);

    uint8_t *bs = malloc(ls * 4), *out = malloc(fs);
    double q[4] = {0, 0, 0, 0}, enc_ms = 0, dec_ms = 0;
    uint64_t bytes = 20 + 4ull * n, type_bytes[3] = {0, 0, 0};
    size_t mm = 0;
    for (int i = 0; i < n; i++) {
        int d = order[i];
        const uint8_t *y = src + fs * d, *u = y + ls, *v = u + cs;
        if (d % keyint == 0) { ectx.num_refs = 0; ectx.cdf_state_valid = 0; }
        ectx.code_bidi = bidi[i];
        ectx.header.frame_number = (uint32_t)d;
        uint32_t len = 0;
        double t0 = now_ms();
        if (snte_encode_frame(&ectx, y, u, v, w, h, w, w / 2, bs, (uint32_t)(ls * 4), &len)) {
            fprintf(stderr, "encode failed at frame %d\n", d);
            return 2;
        }
        double t1 = now_ms();
        memset(out + ls, 128, 2 * cs);
        int dw, dh;
        if (sntd_decode_frame(&dctx, bs, len, out, out + ls, out + ls + cs, &dw, &dh, w, w / 2)) {
            fprintf(stderr, "decode failed at frame %d\n", d);
            return 3;
        }
        dec_ms += now_ms() - t1;
        enc_ms += t1 - t0;
        bytes += len;
        type_bytes[ectx.header.frame_type == SNT_FRAME_KEY ? 0 : (bidi[i] ? 2 : 1)] += len;
        mm += mismatch(ectx.last_y, out, ls) + mismatch(ectx.last_u, out + ls, cs)
            + mismatch(ectx.last_v, out + ls + cs, cs);
        q[0] += psnr(y, out, ls);
        q[1] += psnr(u, out + ls, cs);
        q[2] += psnr(v, out + ls + cs, cs);
        q[3] += ssim(y, out, w, h);
    }
    printf("%d %d %llu %llu %llu %llu %.4f %.4f %.4f %.6f %zu %.1f %.2f\n", qp, n,
           (unsigned long long)bytes, (unsigned long long)type_bytes[0],
           (unsigned long long)type_bytes[1], (unsigned long long)type_bytes[2],
           q[0] / n, q[1] / n, q[2] / n, q[3] / n, mm, enc_ms / n, dec_ms / n);

    if (bits_path) {
#ifdef SNT_TRAIN_STATS
        FILE *bf = fopen(bits_path, "w");
        if (!bf) return 1;
        /* Model groups of SntCdfSet: name, first model, number of models */
#define G(f) fprintf(bf, "group %s %zu %zu\n", #f, offsetof(SntCdfSet, f) / sizeof(SntModel), \
                     sizeof(((SntCdfSet *)0)->f) / sizeof(SntModel))
        G(dc); G(acz); G(acs); G(acg); G(acx); G(acm); G(band); G(imode); G(mbmode); G(bdir);
        G(inter); G(iskip); G(mvx); G(mvy); G(cdc); G(cac);
#undef G
        fprintf(bf, "frametypes %llu %llu %llu\n", (unsigned long long)type_bytes[0],
                (unsigned long long)type_bytes[1], (unsigned long long)type_bytes[2]);
        snte_bits_dump(bf);
        fclose(bf);
#else
        fprintf(stderr, "-bits needs a library built with -DSNT_TRAIN_STATS\n");
        return 1;
#endif
    }
    return 0;
}
