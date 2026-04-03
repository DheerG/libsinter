/*
 * FFmpeg libsinter decoder wrapper
 *
 * This file would live at libavcodec/libsinter_dec.c in the FFmpeg tree.
 * Usage: ffmpeg -i input.sntr -c:v rawvideo output.y4m
 *
 * Based on the pattern from libavcodec/libaomdec.c (AV1 / dav1d).
 */

#ifdef HAVE_FFMPEG

#include "libavutil/imgutils.h"
#include "libavcodec/avcodec.h"
#include "libavcodec/internal.h"
#include "sinter.h"

typedef struct SinterDecContext {
    SntDecoder *decoder;
} SinterDecContext;

static av_cold int sinter_decode_init(AVCodecContext *avctx) {
    SinterDecContext *ctx = avctx->priv_data;
    ctx->decoder = snt_decoder_create();
    if (!ctx->decoder)
        return AVERROR(ENOMEM);
    avctx->pix_fmt = AV_PIX_FMT_YUV420P;
    return 0;
}

static int sinter_decode_frame(AVCodecContext *avctx, void *data,
                                int *got_frame, AVPacket *pkt) {
    SinterDecContext *ctx = avctx->priv_data;
    AVFrame *frame = data;
    SntDecodedFrame snt_frame;

    int ret = snt_decoder_decode(ctx->decoder, pkt->data, pkt->size, &snt_frame);
    if (ret != 0)
        return AVERROR_INVALIDDATA;

    avctx->width  = snt_frame.width;
    avctx->height = snt_frame.height;

    ret = ff_get_buffer(avctx, frame, 0);
    if (ret < 0) return ret;

    /* Copy decoded planes to AVFrame */
    av_image_copy_plane(frame->data[0], frame->linesize[0],
                        snt_frame.y, snt_frame.stride_y,
                        snt_frame.width, snt_frame.height);
    av_image_copy_plane(frame->data[1], frame->linesize[1],
                        snt_frame.u, snt_frame.stride_uv,
                        snt_frame.width / 2, snt_frame.height / 2);
    av_image_copy_plane(frame->data[2], frame->linesize[2],
                        snt_frame.v, snt_frame.stride_uv,
                        snt_frame.width / 2, snt_frame.height / 2);

    *got_frame = 1;
    return pkt->size;
}

static av_cold int sinter_decode_close(AVCodecContext *avctx) {
    SinterDecContext *ctx = avctx->priv_data;
    snt_decoder_destroy(ctx->decoder);
    ctx->decoder = NULL;
    return 0;
}

AVCodec ff_libsinter_decoder = {
    .name           = "libsinter",
    .long_name      = "libsinter Sinter video decoder",
    .type           = AVMEDIA_TYPE_VIDEO,
    .id             = AV_CODEC_ID_SINTER,
    .priv_data_size = sizeof(SinterDecContext),
    .init           = sinter_decode_init,
    .decode         = sinter_decode_frame,
    .close          = sinter_decode_close,
    .capabilities   = AV_CODEC_CAP_DR1,
    .wrapper_name   = "libsinter",
};

#endif /* HAVE_FFMPEG */
