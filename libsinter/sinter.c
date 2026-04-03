/*
 * libsinter — Public C API Implementation
 * Wraps the internal encoder/decoder with opaque handles and stable ABI.
 */

#include <stdlib.h>
#include <string.h>
#include "sinter.h"
#include "codec/encoder.h"
#include "codec/decoder.h"
#include "core/sinter_types.h"

/* ── Opaque encoder struct ────────────────────────────────────────── */
struct SntEncoder {
    SntEncoderCtx ctx;
    int width, height;
    uint8_t *bitstream_buf;
    uint32_t bitstream_cap;
    uint32_t frame_num;
};

/* ── Opaque decoder struct ────────────────────────────────────────── */
struct SntDecoder {
    SntDecoderCtx ctx;
    uint8_t *y_buf;
    uint8_t *u_buf;
    uint8_t *v_buf;
    int buf_w, buf_h;
};

/* ── Options defaults ─────────────────────────────────────────────── */
void snt_encoder_opts_default(SntEncoderOpts *opts) {
    memset(opts, 0, sizeof(*opts));
    opts->qp = 28;
    opts->crf = 0;
    opts->butteraugli = 0;
    opts->speed = 6;
    opts->num_tile_cols = 1;
    opts->num_tile_rows = 1;
    opts->grain_enabled = 0;
}

/* ── Encoder ──────────────────────────────────────────────────────── */

SntEncoder *snt_encoder_create(int width, int height, const SntEncoderOpts *opts) {
    SntEncoder *enc = (SntEncoder *)calloc(1, sizeof(SntEncoder));
    if (!enc) return NULL;

    enc->width = width;
    enc->height = height;
    enc->frame_num = 0;

    /* Generous bitstream buffer */
    enc->bitstream_cap = (uint32_t)(width * height * 3);
    enc->bitstream_buf = (uint8_t *)malloc(enc->bitstream_cap);
    if (!enc->bitstream_buf) { free(enc); return NULL; }

    /* Configure internal context */
    memset(&enc->ctx, 0, sizeof(enc->ctx));
    enc->ctx.header.base_qp = opts ? opts->qp : 28;
    enc->ctx.header.dual_rans = 1;
    enc->ctx.header.num_tile_cols = opts ? opts->num_tile_cols : 1;
    enc->ctx.header.num_tile_rows = opts ? opts->num_tile_rows : 1;
    enc->ctx.header.chroma_qp_offset = 4;
    enc->ctx.header.bit_depth = SNT_BIT_DEPTH_8;

    enc->ctx.psy_cfg.alpha = 0.5f;
    enc->ctx.psy_cfg.qp_clamp = 6;
    enc->ctx.psy_cfg.chroma_qp_offset = 4;
    enc->ctx.psy_cfg.speed_preset = opts ? opts->speed : 6;
    enc->ctx.psy_cfg.grain_enabled = opts ? opts->grain_enabled : 0;
    enc->ctx.psy_cfg.butteraugli_target = opts ? opts->butteraugli : 0;
    enc->ctx.psy_cfg.butteraugli_iters = 3;

    return enc;
}

int snt_encoder_encode(SntEncoder *enc,
                        const uint8_t *y, const uint8_t *u, const uint8_t *v,
                        int stride_y, int stride_uv,
                        SntPacket *pkt) {
    enc->ctx.header.frame_number = enc->frame_num;

    uint32_t bs_len;
    int ret = snte_encode_frame(&enc->ctx, y, u, v,
                                 enc->width, enc->height, stride_y, stride_uv,
                                 enc->bitstream_buf, enc->bitstream_cap, &bs_len);
    if (ret != 0) return ret;

    pkt->data = enc->bitstream_buf;
    pkt->size = bs_len;
    pkt->frame_num = enc->frame_num;
    pkt->is_key = (enc->ctx.header.frame_type == SNT_FRAME_KEY) ? 1 : 0;

    enc->frame_num++;
    return 0;
}

void snt_encoder_destroy(SntEncoder *enc) {
    if (!enc) return;
    free(enc->bitstream_buf);
    free(enc);
}

/* ── Decoder ──────────────────────────────────────────────────────── */

SntDecoder *snt_decoder_create(void) {
    SntDecoder *dec = (SntDecoder *)calloc(1, sizeof(SntDecoder));
    return dec;
}

static void decoder_ensure_buffers(SntDecoder *dec, int w, int h) {
    if (dec->buf_w >= w && dec->buf_h >= h) return;

    free(dec->y_buf); free(dec->u_buf); free(dec->v_buf);
    dec->buf_w = w; dec->buf_h = h;
    dec->y_buf = (uint8_t *)calloc((size_t)w * h, 1);
    dec->u_buf = (uint8_t *)calloc((size_t)(w/2) * (h/2), 1);
    dec->v_buf = (uint8_t *)calloc((size_t)(w/2) * (h/2), 1);
}

int snt_decoder_decode(SntDecoder *dec,
                        const uint8_t *data, uint32_t size,
                        SntDecodedFrame *frame) {
    /* Pre-allocate for expected size (will resize if needed) */
    decoder_ensure_buffers(dec, 4096, 2160);

    int out_w, out_h;
    int ret = sntd_decode_frame(&dec->ctx, data, size,
                                 dec->y_buf, dec->u_buf, dec->v_buf,
                                 &out_w, &out_h, out_w, out_w / 2);

    /* Resize if actual size differs from pre-allocation */
    if (ret == 0) {
        frame->y = dec->y_buf;
        frame->u = dec->u_buf;
        frame->v = dec->v_buf;
        frame->width = out_w;
        frame->height = out_h;
        frame->stride_y = out_w;
        frame->stride_uv = out_w / 2;
    }

    return ret;
}

void snt_decoder_destroy(SntDecoder *dec) {
    if (!dec) return;
    free(dec->y_buf);
    free(dec->u_buf);
    free(dec->v_buf);
    free(dec);
}
