/*
 * libsinter — Public C API
 * Sinter Video Codec Reference Implementation
 *
 * This is the stable public API. Keep surface area minimal (<10 functions).
 * All internal types are opaque. Caller owns input/output buffers.
 *
 * Usage:
 *   SntEncoder *enc = snt_encoder_create(width, height, &opts);
 *   snt_encoder_encode(enc, y, u, v, stride_y, stride_uv, &pkt);
 *   // pkt.data contains encoded bitstream, pkt.size bytes
 *   snt_encoder_destroy(enc);
 *
 *   SntDecoder *dec = snt_decoder_create();
 *   snt_decoder_decode(dec, pkt.data, pkt.size, &frame);
 *   // frame.y/u/v contain decoded planes
 *   snt_decoder_destroy(dec);
 */
#ifndef SINTER_H
#define SINTER_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ── Version ──────────────────────────────────────────────────────── */
#define SINTER_VERSION_MAJOR 0
#define SINTER_VERSION_MINOR 3
#define SINTER_VERSION_PATCH 0

/* ── Opaque handles ───────────────────────────────────────────────── */
typedef struct SntEncoder SntEncoder;
typedef struct SntDecoder SntDecoder;

/* ── Encoder options ──────────────────────────────────────────────── */
typedef struct {
    uint8_t  qp;              /* base QP [0..63], default 28 */
    float    crf;             /* CRF target (0=disabled, uses fixed QP) */
    float    butteraugli;     /* Butteraugli target (0=disabled) */
    uint8_t  speed;           /* speed preset [0..10], default 6 */
    uint8_t  num_tile_cols;   /* tile columns, default 1 */
    uint8_t  num_tile_rows;   /* tile rows, default 1 */
    int      grain_enabled;   /* enable film grain estimation */
} SntEncoderOpts;

/* Initialize options to defaults. Always call this first. */
void snt_encoder_opts_default(SntEncoderOpts *opts);

/* ── Encoded packet ───────────────────────────────────────────────── */
typedef struct {
    uint8_t  *data;     /* encoded bitstream (owned by encoder, valid until next encode) */
    uint32_t  size;     /* bitstream size in bytes */
    uint32_t  frame_num;
    int       is_key;   /* 1 if keyframe */
} SntPacket;

/* ── Decoded frame ────────────────────────────────────────────────── */
typedef struct {
    uint8_t  *y;        /* luma plane (owned by decoder, valid until next decode) */
    uint8_t  *u;        /* chroma-U plane */
    uint8_t  *v;        /* chroma-V plane */
    int       width;
    int       height;
    int       stride_y;
    int       stride_uv;
} SntDecodedFrame;

/* ── Encoder API ──────────────────────────────────────────────────── */

/* Create an encoder for the given resolution. Returns NULL on failure. */
SntEncoder *snt_encoder_create(int width, int height, const SntEncoderOpts *opts);

/* Encode one YUV420 frame. u/v may be NULL for luma-only.
 * On success, pkt is filled with the encoded data. Returns 0 on success. */
int snt_encoder_encode(SntEncoder *enc,
                        const uint8_t *y, const uint8_t *u, const uint8_t *v,
                        int stride_y, int stride_uv,
                        SntPacket *pkt);

/* Destroy encoder and free all resources. */
void snt_encoder_destroy(SntEncoder *enc);

/* ── Decoder API ──────────────────────────────────────────────────── */

/* Create a decoder. Returns NULL on failure. */
SntDecoder *snt_decoder_create(void);

/* Decode one frame from bitstream data.
 * On success, frame is filled with decoded planes. Returns 0 on success. */
int snt_decoder_decode(SntDecoder *dec,
                        const uint8_t *data, uint32_t size,
                        SntDecodedFrame *frame);

/* Destroy decoder and free all resources. */
void snt_decoder_destroy(SntDecoder *dec);

#ifdef __cplusplus
}
#endif

#endif /* SINTER_H */
