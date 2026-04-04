/*
 * sinterenc — Sinter Codec Reference Encoder CLI
 *
 * Usage: sinterenc -i input.yuv -w WIDTH -h HEIGHT -q QP -o output.sntr
 *
 * Reads raw YUV 4:2:0 (luma only for v0.3), encodes using Sinter,
 * writes .sntr bitstream with simple IVF-style framing.
 *
 * .sntr container format:
 *   [4] magic "SNTR"
 *   [4] version (1)
 *   [4] width
 *   [4] height
 *   [4] num_frames
 *   Per frame:
 *     [4] frame_size (bytes)
 *     [frame_size] encoded frame data
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "codec/encoder.h"
#include "core/sinter_types.h"

static void write32(FILE *f, uint32_t v) {
    fputc(v & 0xFF, f);
    fputc((v >> 8) & 0xFF, f);
    fputc((v >> 16) & 0xFF, f);
    fputc((v >> 24) & 0xFF, f);
}

int main(int argc, char **argv) {
    const char *input_path = NULL, *output_path = NULL;
    int width = 0, height = 0;
    uint8_t qp = 28;
    int num_frames = 1;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-i") == 0 && i + 1 < argc) input_path = argv[++i];
        else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) output_path = argv[++i];
        else if (strcmp(argv[i], "-w") == 0 && i + 1 < argc) width = atoi(argv[++i]);
        else if (strcmp(argv[i], "-h") == 0 && i + 1 < argc) height = atoi(argv[++i]);
        else if (strcmp(argv[i], "-q") == 0 && i + 1 < argc) qp = (uint8_t)atoi(argv[++i]);
        else if (strcmp(argv[i], "-n") == 0 && i + 1 < argc) num_frames = atoi(argv[++i]);
    }

    if (!input_path || !output_path || width <= 0 || height <= 0) {
        fprintf(stderr, "Usage: sinterenc -i input.yuv -w WIDTH -h HEIGHT [-q QP] [-n FRAMES] -o output.sntr\n");
        return 1;
    }

    FILE *fin = fopen(input_path, "rb");
    FILE *fout = fopen(output_path, "wb");
    if (!fin || !fout) {
        fprintf(stderr, "Error opening files\n");
        return 1;
    }

    /* Write container header */
    fwrite("SNTR", 1, 4, fout);
    write32(fout, 1); /* version */
    write32(fout, (uint32_t)width);
    write32(fout, (uint32_t)height);
    write32(fout, (uint32_t)num_frames);

    /* Allocate buffers */
    size_t frame_size = (size_t)width * height;
    size_t chroma_size = frame_size / 4; /* per plane, 4:2:0 */
    int chroma_w = width / 2;
    uint8_t *input_y = (uint8_t *)malloc(frame_size);
    uint8_t *input_u = (uint8_t *)malloc(chroma_size);
    uint8_t *input_v = (uint8_t *)malloc(chroma_size);
    uint8_t *bitstream = (uint8_t *)malloc(frame_size * 4); /* generous */

    SntEncoderCtx ctx;
    memset(&ctx, 0, sizeof(ctx));
    ctx.header.base_qp = qp;
    ctx.header.dual_rans = 1;
    ctx.header.num_tile_cols = 1;
    ctx.header.num_tile_rows = 1;
    ctx.header.chroma_qp_offset = 4;
    ctx.psy_cfg.alpha = 0.5f;
    ctx.psy_cfg.qp_clamp = 6;

    uint64_t total_bytes = 0;
    for (int f = 0; f < num_frames; f++) {
        /* Read luma + chroma planes (YUV 4:2:0) */
        size_t read = fread(input_y, 1, frame_size, fin);
        if (read < frame_size) {
            fprintf(stderr, "Short read at frame %d\n", f);
            break;
        }
        fread(input_u, 1, chroma_size, fin);
        fread(input_v, 1, chroma_size, fin);

        ctx.header.frame_number = (uint32_t)f;

        uint32_t bs_len;
        int ret = snte_encode_frame(&ctx, input_y, input_u, input_v,
                                     width, height, width, chroma_w,
                                     bitstream, (uint32_t)(frame_size * 4), &bs_len);
        if (ret != 0) {
            fprintf(stderr, "Encode error at frame %d\n", f);
            break;
        }

        write32(fout, bs_len);
        fwrite(bitstream, 1, bs_len, fout);
        total_bytes += bs_len;

        fprintf(stderr, "Frame %d: %u bytes\n", f, bs_len);
    }

    fprintf(stderr, "Encoded %d frames, total %llu bytes (%.2f kbps at 30fps)\n",
            num_frames, (unsigned long long)total_bytes,
            (double)total_bytes * 8.0 * 30.0 / num_frames / 1000.0);

    fclose(fin);
    fclose(fout);
    free(input_y);
    free(bitstream);
    return 0;
}
