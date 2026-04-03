/*
 * FFmpeg libsinter encoder wrapper
 *
 * This file would live at libavcodec/libsinter_enc.c in the FFmpeg tree.
 * It wraps libsinter's public C API for FFmpeg's AVCodec interface.
 *
 * Build: ffmpeg configure --enable-libsinter (requires pkg-config libsinter)
 * Usage: ffmpeg -i input.mp4 -c:v libsinter -qp 28 output.sntr
 *
 * Based on the pattern from libavcodec/libaomenc.c (AV1).
 */

/*
 * NOTE: This is a standalone reference file showing the FFmpeg integration
 * pattern. It cannot be compiled without FFmpeg headers. It serves as
 * documentation for how to write the actual FFmpeg patch.
 */

#ifdef HAVE_FFMPEG /* Only compile within FFmpeg build system */

#include "libavutil/opt.h"
#include "libavutil/imgutils.h"
#include "libavcodec/avcodec.h"
#include "libavcodec/internal.h"
#include "sinter.h"

typedef struct SinterEncContext {
    AVClass *class;
    SntEncoder *encoder;
    int qp;
    float crf;
    float butteraugli;
    int speed;
    int tile_cols;
    int tile_rows;
    int grain;
} SinterEncContext;

static av_cold int sinter_encode_init(AVCodecContext *avctx) {
    SinterEncContext *ctx = avctx->priv_data;

    SntEncoderOpts opts;
    snt_encoder_opts_default(&opts);
    opts.qp = ctx->qp;
    opts.crf = ctx->crf;
    opts.butteraugli = ctx->butteraugli;
    opts.speed = ctx->speed;
    opts.num_tile_cols = ctx->tile_cols;
    opts.num_tile_rows = ctx->tile_rows;
    opts.grain_enabled = ctx->grain;

    ctx->encoder = snt_encoder_create(avctx->width, avctx->height, &opts);
    if (!ctx->encoder)
        return AVERROR(ENOMEM);

    avctx->pix_fmt = AV_PIX_FMT_YUV420P;
    return 0;
}

static int sinter_encode_frame(AVCodecContext *avctx, AVPacket *pkt,
                                const AVFrame *frame, int *got_packet) {
    SinterEncContext *ctx = avctx->priv_data;
    SntPacket snt_pkt;

    int ret = snt_encoder_encode(ctx->encoder,
                                  frame->data[0], frame->data[1], frame->data[2],
                                  frame->linesize[0], frame->linesize[1],
                                  &snt_pkt);
    if (ret != 0)
        return AVERROR_EXTERNAL;

    ret = ff_get_encode_buffer(avctx, pkt, snt_pkt.size, 0);
    if (ret < 0) return ret;

    memcpy(pkt->data, snt_pkt.data, snt_pkt.size);
    if (snt_pkt.is_key)
        pkt->flags |= AV_PKT_FLAG_KEY;

    *got_packet = 1;
    return 0;
}

static av_cold int sinter_encode_close(AVCodecContext *avctx) {
    SinterEncContext *ctx = avctx->priv_data;
    snt_encoder_destroy(ctx->encoder);
    ctx->encoder = NULL;
    return 0;
}

#define OFFSET(x) offsetof(SinterEncContext, x)
#define VE AV_OPT_FLAG_VIDEO_PARAM | AV_OPT_FLAG_ENCODING_PARAM
static const AVOption options[] = {
    { "qp",          "Quantization parameter",   OFFSET(qp),          AV_OPT_TYPE_INT,   {.i64 = 28},  0, 63, VE },
    { "sinter-crf",  "Constant rate factor",     OFFSET(crf),         AV_OPT_TYPE_FLOAT, {.dbl = 0},   0, 63, VE },
    { "sinter-speed","Speed preset",             OFFSET(speed),       AV_OPT_TYPE_INT,   {.i64 = 6},   0, 10, VE },
    { "tile-columns","Number of tile columns",   OFFSET(tile_cols),   AV_OPT_TYPE_INT,   {.i64 = 1},   1, 16, VE },
    { "tile-rows",   "Number of tile rows",      OFFSET(tile_rows),   AV_OPT_TYPE_INT,   {.i64 = 1},   1, 16, VE },
    { "grain",       "Enable grain synthesis",   OFFSET(grain),       AV_OPT_TYPE_BOOL,  {.i64 = 0},   0,  1, VE },
    { "butteraugli", "Butteraugli quality target",OFFSET(butteraugli), AV_OPT_TYPE_FLOAT, {.dbl = 0},   0, 10, VE },
    { NULL },
};

static const AVClass sinter_enc_class = {
    .class_name = "libsinter",
    .item_name  = av_default_item_name,
    .option     = options,
    .version    = LIBAVUTIL_VERSION_INT,
};

AVCodec ff_libsinter_encoder = {
    .name           = "libsinter",
    .long_name      = "libsinter Sinter video encoder",
    .type           = AVMEDIA_TYPE_VIDEO,
    .id             = AV_CODEC_ID_SINTER, /* would be registered in codec_id.h */
    .priv_data_size = sizeof(SinterEncContext),
    .init           = sinter_encode_init,
    .encode2        = sinter_encode_frame,
    .close          = sinter_encode_close,
    .pix_fmts       = (const enum AVPixelFormat[]){ AV_PIX_FMT_YUV420P, AV_PIX_FMT_NONE },
    .priv_class     = &sinter_enc_class,
    .capabilities   = AV_CODEC_CAP_AUTO_THREADS,
    .wrapper_name   = "libsinter",
};

#endif /* HAVE_FFMPEG */
