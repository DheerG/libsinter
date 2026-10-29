/*
 * sinterenc — Sinter Codec Reference Encoder CLI
 *
 * Usage:
 *   sinterenc -i input.y4m -o output.sntr [-q QP] [-k KEYINT] [-b BFRAMES]
 *   sinterenc -i input.yuv -w WIDTH -h HEIGHT -o output.sntr [-q QP] [-n FRAMES] [-k KEYINT]
 *   cat input.y4m | sinterenc -o output.sntr [-q QP]
 *
 * Reads Y4M (auto-detected) or raw YUV 4:2:0, encodes, writes .sntr container.
 *
 * GOP structure: a key frame every KEYINT frames (closed GOPs), and up to
 * BFRAMES B-frames (default 2) between anchor frames. The last frame of
 * each GOP is an anchor, so no B-frame predicts across a key frame.
 *
 * .sntr container format:
 *   [4] magic "SNTR"
 *   [4] version (5)
 *   [4] width
 *   [4] height
 *   [4] num_frames
 *   Per frame, in decode order (each frame header carries its display
 *   number; decoders reorder):
 *     [4] frame_size (bytes)
 *     [frame_size] encoded frame data
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "codec/encoder.h"
#include "core/sinter_types.h"

static void write32(FILE *f, uint32_t v) {
    fputc(v & 0xFF, f); fputc((v >> 8) & 0xFF, f);
    fputc((v >> 16) & 0xFF, f); fputc((v >> 24) & 0xFF, f);
}

/* Parse Y4M header. Returns 1 on success, 0 on failure.
 * Reads byte-by-byte to avoid text/binary mixing with fgets. */
static int parse_y4m_header(FILE *f, int *width, int *height, int *fps_num, int *fps_den) {
    char line[256];
    int pos = 0, c;
    while ((c = fgetc(f)) != '\n' && c != EOF && pos < 255)
        line[pos++] = (char)c;
    line[pos] = '\0';
    if (strncmp(line, "YUV4MPEG2", 9) != 0) return 0;

    *width = 0; *height = 0; *fps_num = 30; *fps_den = 1;

    char *p = line + 9;
    while (*p) {
        while (*p == ' ') p++;
        if (*p == 'W') { *width = atoi(p + 1); }
        else if (*p == 'H') { *height = atoi(p + 1); }
        else if (*p == 'F') {
            *fps_num = atoi(p + 1);
            char *colon = strchr(p, ':');
            if (colon) *fps_den = atoi(colon + 1);
        }
        while (*p && *p != ' ' && *p != '\n') p++;
    }
    return (*width > 0 && *height > 0);
}

/* Read one Y4M frame. Returns 1 on success, 0 on EOF/error.
 * Uses binary I/O only (no fgets) to avoid text/binary mixing issues. */
static int read_y4m_frame(FILE *f, uint8_t *y, uint8_t *u, uint8_t *v,
                           int width, int height) {
    /* Read "FRAME" tag byte by byte until newline */
    int c = fgetc(f);
    if (c == EOF) return 0;
    if (c != 'F') return 0;
    /* Skip rest of FRAME line (consume until newline) */
    while ((c = fgetc(f)) != '\n' && c != EOF);

    size_t luma_size = (size_t)width * height;
    size_t chroma_size = luma_size / 4;
    if (fread(y, 1, luma_size, f) != luma_size) return 0;
    if (fread(u, 1, chroma_size, f) != chroma_size) return 0;
    if (fread(v, 1, chroma_size, f) != chroma_size) return 0;
    return 1;
}

int main(int argc, char **argv) {
    const char *input_path = NULL, *output_path = NULL;
    int width = 0, height = 0;
    uint8_t qp = 20;
    int num_frames = 0; /* 0 = read until EOF */
    int keyint = 30;
    int bframes = 2;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-i") == 0 && i + 1 < argc) input_path = argv[++i];
        else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) output_path = argv[++i];
        else if (strcmp(argv[i], "-w") == 0 && i + 1 < argc) width = atoi(argv[++i]);
        else if (strcmp(argv[i], "-h") == 0 && i + 1 < argc) height = atoi(argv[++i]);
        else if (strcmp(argv[i], "-q") == 0 && i + 1 < argc) qp = (uint8_t)atoi(argv[++i]);
        else if (strcmp(argv[i], "-n") == 0 && i + 1 < argc) num_frames = atoi(argv[++i]);
        else if (strcmp(argv[i], "-k") == 0 && i + 1 < argc) keyint = atoi(argv[++i]);
        else if (strcmp(argv[i], "-b") == 0 && i + 1 < argc) bframes = atoi(argv[++i]);
    }

    if (!output_path) {
        fprintf(stderr, "Usage: sinterenc -i input.y4m -o output.sntr [-q QP] [-k KEYINT] [-b BFRAMES]\n");
        fprintf(stderr, "       sinterenc -i input.yuv -w W -h H -o output.sntr [-q QP] [-n N]\n");
        return 1;
    }

    FILE *fin = input_path ? fopen(input_path, "rb") : stdin;
    if (!fin) { fprintf(stderr, "Cannot open input\n"); return 1; }

    /* Try Y4M header */
    int is_y4m = 0, fps_num = 30, fps_den = 1;
    {
        int y4m_w = 0, y4m_h = 0;
        if (parse_y4m_header(fin, &y4m_w, &y4m_h, &fps_num, &fps_den)) {
            is_y4m = 1;
            if (width == 0) width = y4m_w;
            if (height == 0) height = y4m_h;
        } else {
            rewind(fin);
        }
    }

    if (keyint < 1) keyint = 1;
    if (bframes < 0) bframes = 0;

    if (width <= 0 || height <= 0) {
        fprintf(stderr, "Error: width/height not specified and no Y4M header found\n");
        return 1;
    }

    FILE *fout = fopen(output_path, "wb");
    if (!fout) { fprintf(stderr, "Cannot open output\n"); return 1; }

    size_t luma_size = (size_t)width * height;
    size_t chroma_size = luma_size / 4;
    size_t frame_size = luma_size + 2 * chroma_size;
    int chroma_w = width / 2;
    /* One GOP of input frames, encoded in coding order once read */
    uint8_t *gop = (uint8_t *)malloc(frame_size * (size_t)keyint);
    int *order = (int *)malloc(sizeof(int) * (size_t)keyint);
    int *is_b = (int *)malloc(sizeof(int) * (size_t)keyint);
    uint8_t *bitstream = (uint8_t *)malloc(luma_size * 4);

    SntEncoderCtx ctx;
    memset(&ctx, 0, sizeof(ctx));
    ctx.header.base_qp = qp;
    ctx.header.dual_rans = 1;
    ctx.header.num_tile_cols = 1;
    ctx.header.num_tile_rows = 1;
    ctx.header.chroma_qp_offset = 4;
    ctx.keyint = keyint;
    ctx.psy_cfg.alpha = 0.5f;
    ctx.psy_cfg.qp_clamp = 6;

    /* Write container header (placeholder num_frames, updated at end) */
    long header_pos = 0;
    fwrite("SNTR", 1, 4, fout);
    write32(fout, 5); /* version */
    write32(fout, (uint32_t)width);
    write32(fout, (uint32_t)height);
    header_pos = ftell(fout);
    write32(fout, 0); /* placeholder */

    uint64_t total_bytes = 0;
    int frames_encoded = 0;

    for (int gop_start = 0; num_frames == 0 || gop_start < num_frames; gop_start += keyint) {
        /* Read up to one GOP */
        int len = 0;
        while (len < keyint && (num_frames == 0 || gop_start + len < num_frames)) {
            uint8_t *y = gop + frame_size * len, *u = y + luma_size, *v = u + chroma_size;
            int got_frame;
            if (is_y4m) {
                got_frame = read_y4m_frame(fin, y, u, v, width, height);
            } else {
                got_frame = (fread(y, 1, luma_size, fin) == luma_size);
                if (got_frame) {
                    fread(u, 1, chroma_size, fin);
                    fread(v, 1, chroma_size, fin);
                }
            }
            if (!got_frame) break;
            len++;
        }
        if (len == 0) break;

        /* Coding order: the key frame, then each anchor followed by the
         * B-frames before it in display order */
        int n = 0, prev = 0;
        order[n] = 0; is_b[n++] = 0;
        while (prev < len - 1) {
            int next = prev + bframes + 1 < len - 1 ? prev + bframes + 1 : len - 1;
            order[n] = next; is_b[n++] = 0;
            for (int b = prev + 1; b < next; b++) { order[n] = b; is_b[n++] = 1; }
            prev = next;
        }

        for (int i = 0; i < n; i++) {
            int f = gop_start + order[i];
            uint8_t *y = gop + frame_size * order[i], *u = y + luma_size, *v = u + chroma_size;
            if (order[i] == 0) {
                ctx.num_refs = 0;  /* key frame */
                ctx.cdf_state_valid = 0;
            }
            ctx.code_bidi = is_b[i];
            ctx.header.frame_number = (uint32_t)f;

            uint32_t bs_len;
            int ret = snte_encode_frame(&ctx, y, u, v, width, height, width, chroma_w,
                                         bitstream, (uint32_t)(luma_size * 4), &bs_len);
            if (ret != 0) {
                fprintf(stderr, "Encode error at frame %d\n", f);
                goto done;
            }

            write32(fout, bs_len);
            fwrite(bitstream, 1, bs_len, fout);
            total_bytes += bs_len;
            frames_encoded++;

            const char *type = ctx.header.frame_type == SNT_FRAME_KEY ? "I" :
                               ctx.header.frame_type == SNT_FRAME_BIDI ? "B" : "P";
            fprintf(stderr, "Frame %d [%s]: %u bytes\n", f, type, bs_len);
        }
        if (len < keyint) break;
    }
done:

    /* Update num_frames in header */
    fseek(fout, header_pos, SEEK_SET);
    write32(fout, (uint32_t)frames_encoded);

    fprintf(stderr, "Encoded %d frames, %llu bytes total (%.1f KB, %.2f kbps at %d/%d fps)\n",
            frames_encoded, (unsigned long long)total_bytes,
            (double)total_bytes / 1024.0,
            (double)total_bytes * 8.0 * fps_num / fps_den / frames_encoded / 1000.0,
            fps_num, fps_den);

    if (fin != stdin) fclose(fin);
    fclose(fout);
    free(gop); free(order); free(is_b); free(bitstream);
    free(ctx.ref_y); free(ctx.ref_u); free(ctx.ref_v); free(ctx.cdf_state);
    free(ctx.ref_old_y); free(ctx.ref_old_u); free(ctx.ref_old_v); free(ctx.cdf_state_b);
    free(ctx.last_y); free(ctx.last_u); free(ctx.last_v);
    free(ctx.ref_coef); free(ctx.ref_old_coef);
    free(ctx.key_y); free(ctx.key_u); free(ctx.key_v); free(ctx.key_coef);
    free(ctx.anchor_mv); free(ctx.anchor_dir);
    return 0;
}
