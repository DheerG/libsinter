/*
 * sinterdec — Sinter Codec Reference Decoder CLI
 *
 * Usage:
 *   sinterdec -i input.sntr -o output.y4m
 *   sinterdec -i input.sntr -o output.yuv [--raw]
 *   sinterdec -i input.sntr -o - | ffplay -f yuv4mpeg -
 *
 * Reads .sntr bitstream, decodes, writes Y4M (default) or raw YUV 4:2:0.
 * Frames are stored in decode order; output is reordered to display order
 * by the frame number in each frame header.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "codec/decoder.h"
#include "core/sinter_types.h"

/* Frames decoded but not yet written. A B-frame is decoded after the
 * anchor that follows it in display order, so anchors wait here. */
#define MAX_PENDING 16
typedef struct { uint32_t num; uint8_t *buf; } Pending;

static void write_frame(FILE *fout, int raw_output, const uint8_t *buf, size_t frame_size) {
    if (!raw_output)
        fprintf(fout, "FRAME\n");
    fwrite(buf, 1, frame_size, fout);
}

static uint32_t read32(FILE *f) {
    uint8_t b[4];
    if (fread(b, 1, 4, f) != 4) return 0;
    return (uint32_t)b[0] | ((uint32_t)b[1] << 8) |
           ((uint32_t)b[2] << 16) | ((uint32_t)b[3] << 24);
}

int main(int argc, char **argv) {
    const char *input_path = NULL, *output_path = NULL;
    int raw_output = 0;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-i") == 0 && i + 1 < argc) input_path = argv[++i];
        else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) output_path = argv[++i];
        else if (strcmp(argv[i], "--raw") == 0) raw_output = 1;
    }

    if (!input_path || !output_path) {
        fprintf(stderr, "Usage: sinterdec -i input.sntr -o output.y4m\n");
        fprintf(stderr, "       sinterdec -i input.sntr -o output.yuv --raw\n");
        return 1;
    }

    FILE *fin = fopen(input_path, "rb");
    FILE *fout = (strcmp(output_path, "-") == 0) ? stdout : fopen(output_path, "wb");
    if (!fin || !fout) { fprintf(stderr, "Error opening files\n"); return 1; }

    /* Read container header */
    char magic[4];
    fread(magic, 1, 4, fin);
    if (memcmp(magic, "SNTR", 4) != 0) {
        fprintf(stderr, "Not a .sntr file\n");
        return 1;
    }
    uint32_t version = read32(fin);
    uint32_t width = read32(fin);
    uint32_t height = read32(fin);
    uint32_t num_frames = read32(fin);

    fprintf(stderr, "SNTR v%u: %ux%u, %u frames\n", version, width, height, num_frames);
    if (version != 4) {
        fprintf(stderr, "Unsupported .sntr version %u (this decoder reads version 4)\n", version);
        return 1;
    }

    /* Write Y4M header */
    if (!raw_output) {
        fprintf(fout, "YUV4MPEG2 W%u H%u F30:1 Ip C420\n", width, height);
    }

    size_t luma_size = (size_t)width * height;
    size_t chroma_size = luma_size / 4;
    int chroma_w = (int)width / 2;
    uint8_t *bitstream = (uint8_t *)malloc(luma_size * 4);
    uint8_t *output_y = (uint8_t *)malloc(luma_size);
    uint8_t *output_u = (uint8_t *)calloc(chroma_size, 1);
    uint8_t *output_v = (uint8_t *)calloc(chroma_size, 1);

    size_t frame_size = luma_size + 2 * chroma_size;
    Pending pend[MAX_PENDING];
    int npend = 0;
    uint32_t next_out = 0;

    SntDecoderCtx ctx;
    memset(&ctx, 0, sizeof(ctx));

    for (uint32_t f = 0; f < num_frames; f++) {
        uint32_t frame_len = read32(fin);
        if (frame_len == 0 || frame_len > luma_size * 4) {
            fprintf(stderr, "Invalid frame size at frame %u\n", f);
            break;
        }

        if (fread(bitstream, 1, frame_len, fin) != frame_len) {
            fprintf(stderr, "Short read at frame %u\n", f);
            break;
        }

        int out_w, out_h;
        memset(output_u, 128, chroma_size);
        memset(output_v, 128, chroma_size);

        int ret = sntd_decode_frame(&ctx, bitstream, frame_len, output_y,
                                     output_u, output_v, &out_w, &out_h,
                                     (int)width, chroma_w);
        if (ret != 0) {
            fprintf(stderr, "Decode error at frame %u\n", f);
            break;
        }

        fprintf(stderr, "Frame %u decoded: %dx%d (display %u)\n", f, out_w, out_h,
                ctx.header.frame_number);

        /* Queue the frame, then write every frame that is next in
         * display order. A full queue (malformed stream) writes its
         * lowest-numbered frame. */
        Pending *q = &pend[npend++];
        q->num = ctx.header.frame_number;
        q->buf = (uint8_t *)malloc(frame_size);
        memcpy(q->buf, output_y, luma_size);
        memcpy(q->buf + luma_size, output_u, chroma_size);
        memcpy(q->buf + luma_size + chroma_size, output_v, chroma_size);
        for (;;) {
            int k = -1;
            for (int j = 0; j < npend; j++)
                if (pend[j].num == next_out) k = j;
            if (k < 0 && npend == MAX_PENDING) {
                k = 0;
                for (int j = 1; j < npend; j++)
                    if (pend[j].num < pend[k].num) k = j;
            }
            if (k < 0) break;
            write_frame(fout, raw_output, pend[k].buf, frame_size);
            next_out = pend[k].num + 1;
            free(pend[k].buf);
            pend[k] = pend[--npend];
        }
    }

    /* Write what is left, in display order */
    while (npend > 0) {
        int k = 0;
        for (int j = 1; j < npend; j++)
            if (pend[j].num < pend[k].num) k = j;
        write_frame(fout, raw_output, pend[k].buf, frame_size);
        free(pend[k].buf);
        pend[k] = pend[--npend];
    }

    if (fin) fclose(fin);
    if (fout != stdout) fclose(fout);
    free(bitstream); free(output_y); free(output_u); free(output_v);
    free(ctx.ref_y); free(ctx.ref_u); free(ctx.ref_v); free(ctx.cdf_state);
    free(ctx.ref_old_y); free(ctx.ref_old_u); free(ctx.ref_old_v); free(ctx.cdf_state_b);
    free(ctx.ref_coef); free(ctx.ref_old_coef);
    return 0;
}
