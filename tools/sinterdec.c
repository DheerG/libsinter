/*
 * sinterdec — Sinter Codec Reference Decoder CLI
 *
 * Usage: sinterdec -i input.sntr -o output.yuv
 *
 * Reads .sntr bitstream, decodes, writes raw YUV 4:2:0 (luma only for v0.3).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "codec/decoder.h"
#include "core/sinter_types.h"

static uint32_t read32(FILE *f) {
    uint8_t b[4];
    if (fread(b, 1, 4, f) != 4) return 0;
    return (uint32_t)b[0] | ((uint32_t)b[1] << 8) |
           ((uint32_t)b[2] << 16) | ((uint32_t)b[3] << 24);
}

int main(int argc, char **argv) {
    const char *input_path = NULL, *output_path = NULL;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-i") == 0 && i + 1 < argc) input_path = argv[++i];
        else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) output_path = argv[++i];
    }

    if (!input_path || !output_path) {
        fprintf(stderr, "Usage: sinterdec -i input.sntr -o output.yuv\n");
        return 1;
    }

    FILE *fin = fopen(input_path, "rb");
    FILE *fout = fopen(output_path, "wb");
    if (!fin || !fout) {
        fprintf(stderr, "Error opening files\n");
        return 1;
    }

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

    size_t frame_size = (size_t)width * height;
    uint8_t *bitstream = (uint8_t *)malloc(frame_size * 2);
    uint8_t *output_y = (uint8_t *)malloc(frame_size);
    uint8_t *chroma_uv = (uint8_t *)calloc(frame_size / 2, 1); /* blank chroma */

    SntDecoderCtx ctx;
    memset(&ctx, 0, sizeof(ctx));

    for (uint32_t f = 0; f < num_frames; f++) {
        uint32_t frame_len = read32(fin);
        if (frame_len == 0 || frame_len > frame_size * 2) {
            fprintf(stderr, "Invalid frame size at frame %u\n", f);
            break;
        }

        fread(bitstream, 1, frame_len, fin);

        int out_w, out_h;
        int ret = sntd_decode_frame(&ctx, bitstream, frame_len, output_y,
                                     NULL, NULL, &out_w, &out_h,
                                     (int)width, 0);
        if (ret != 0) {
            fprintf(stderr, "Decode error at frame %u\n", f);
            break;
        }

        /* Write YUV 4:2:0: luma + blank chroma */
        fwrite(output_y, 1, frame_size, fout);
        fwrite(chroma_uv, 1, frame_size / 2, fout);

        fprintf(stderr, "Frame %u decoded: %dx%d\n", f, out_w, out_h);
    }

    fclose(fin);
    fclose(fout);
    free(bitstream);
    free(output_y);
    free(chroma_uv);
    return 0;
}
