/*
 * Sinter Codec — Common Types
 *
 * Shared data structures used across encoder and decoder.
 * int16_t for 8-bit content coefficients, int32_t for 10-bit+.
 * Never bare 'int' for codec data.
 */
#ifndef SNT_TYPES_H
#define SNT_TYPES_H

#include <stdint.h>
#include <stdbool.h>

/* ── Block sizes ──────────────────────────────────────────────────── */
#define SNT_MAX_SB_SIZE   128  /* maximum superblock dimension */
#define SNT_MIN_BLK_SIZE  4    /* minimum block dimension */
#define SNT_MAX_BLK_LOG2  7    /* log2(128) */
#define SNT_MIN_BLK_LOG2  2    /* log2(4) */

/* Number of supported block sizes: 4, 8, 16, 32, 64, 128 */
#define SNT_NUM_BLK_SIZES 6

/* ── Bit depth ────────────────────────────────────────────────────── */
#define SNT_BIT_DEPTH_8   8
#define SNT_BIT_DEPTH_10  10

/* ── Frame types ──────────────────────────────────────────────────── */
typedef enum {
    SNT_FRAME_KEY    = 0,  /* intra only */
    SNT_FRAME_INTER  = 1,  /* inter (P-frame) */
    SNT_FRAME_BIDI   = 2,  /* bi-directional (B-frame) */
} SntFrameType;

/* ── Transfer functions (HDR) ─────────────────────────────────────── */
typedef enum {
    SNT_TF_SDR    = 0,  /* BT.709 gamma */
    SNT_TF_PQ     = 1,  /* SMPTE ST 2084 (HDR10) */
    SNT_TF_HLG    = 2,  /* ARIB STD-B67 */
} SntTransferFunc;

/* ── Color primaries ──────────────────────────────────────────────── */
typedef enum {
    SNT_CP_BT709  = 0,
    SNT_CP_BT2020 = 1,
} SntColorPrimaries;

/* ── Transform modes ──────────────────────────────────────────────── */
typedef enum {
    SNT_TX_LAPPED = 0,  /* lapped transform (TDLT): pre/post filter + DCT */
    SNT_TX_BLOCK  = 1,  /* plain block DCT (identity pre/post filter) */
} SntTransformMode;

/* ── Prediction modes ─────────────────────────────────────────────── */
typedef enum {
    SNT_PRED_INTRA_DC   = 0,
    SNT_PRED_INTRA_FREQ = 1,  /* frequency-domain prediction */
    SNT_PRED_INTER      = 2,
} SntPredMode;

/* ── Motion vector ────────────────────────────────────────────────── */
typedef struct {
    int16_t x;  /* quarter-pel luma units; chroma reads the same value as 1/8-pel */
    int16_t y;
} SntMV;

/* ── Affine motion parameters (per superblock) ────────────────────── */
typedef struct {
    /* 6-parameter affine: x' = a*x + b*y + c, y' = d*x + e*y + f
     * Stored in Q16 fixed-point (16 fractional bits) */
    int32_t a, b, c;
    int32_t d, e, f;
} SntAffineParams;

/* ── Block info ───────────────────────────────────────────────────── */
typedef struct {
    uint16_t       x, y;          /* position in pixels */
    uint8_t        log2_w, log2_h;/* log2 of width/height */
    SntTransformMode tx_mode;
    SntPredMode    pred_mode;
    int8_t         qp;            /* quantization parameter */
    int8_t         qp_delta;      /* perceptual QP offset */
    SntMV          mv;            /* motion vector (inter only) */
    uint16_t       pvq_k;         /* PVQ pulse count */
    uint32_t       pvq_gain;      /* quantized PVQ gain */
    bool           skip;          /* zero-coefficient block */
} SntBlockInfo;

/* ── Tile info ────────────────────────────────────────────────────── */
typedef struct {
    uint16_t col, row;           /* tile grid position */
    uint16_t x, y;               /* pixel position */
    uint16_t width, height;      /* tile dimensions in pixels */
} SntTileInfo;

/* ── Frame buffer ─────────────────────────────────────────────────── */
/* Padded and aligned for SIMD. Padding = SNT_MAX_SB_SIZE on each edge. */
typedef struct {
    uint8_t  *y;     /* luma plane */
    uint8_t  *u;     /* chroma-U plane */
    uint8_t  *v;     /* chroma-V plane */
    int16_t  *y_coeff; /* luma coefficient buffer */
    int16_t  *u_coeff; /* chroma-U coefficient buffer */
    int16_t  *v_coeff; /* chroma-V coefficient buffer */
    uint32_t  width;
    uint32_t  height;
    uint32_t  stride_y;   /* luma stride (includes padding) */
    uint32_t  stride_uv;  /* chroma stride */
    uint32_t  pad;        /* padding on each edge */
} SntFrame;

/* ── Film grain parameters (per-frame) ────────────────────────────── */
typedef struct {
    bool     enabled;
    /* AR(3) luma coefficients */
    int16_t  ar_coeff_y[3];    /* Q8 fixed-point */
    /* AR(4) chroma coefficients (4th tap = co-located luma grain) */
    int16_t  ar_coeff_u[4];    /* Q8 */
    int16_t  ar_coeff_v[4];    /* Q8 */
    uint8_t  strength_y;       /* luma grain strength [0..255] */
    uint8_t  strength_u;       /* chroma-U grain strength */
    uint8_t  strength_v;       /* chroma-V grain strength */
    uint16_t seed;             /* derived from frame_number */
} SntGrainParams;

/* ── Frame header ─────────────────────────────────────────────────── */
typedef struct {
    uint32_t      width;
    uint32_t      height;
    SntFrameType  frame_type;
    uint8_t       base_qp;       /* [0..63] */
    uint8_t       num_tile_cols;
    uint8_t       num_tile_rows;
    bool          dual_rans;     /* normative: 1-bit flag, dual interleaved rANS */
    int8_t        chroma_qp_offset; /* picture-level chroma QP delta */
    uint8_t       bit_depth;     /* 8 or 10 */
    SntTransferFunc transfer_func; /* SDR, PQ, or HLG */
    SntColorPrimaries color_primaries; /* BT.709 or BT.2020 */
    uint8_t       temporal_layer; /* 0=base, 1-3=enhancement (for B-frames) */
    SntGrainParams grain;
    uint32_t      frame_number;
} SntFrameHeader;

/* ── Perceptual encoder config (non-normative) ────────────────────── */
typedef struct {
    float    alpha;              /* activity masking strength (default 0.5) */
    int8_t   qp_clamp;          /* max |delta_QP| (default 6) */
    int8_t   chroma_qp_offset;  /* default +4 to +6 */
    float    butteraugli_target; /* 0 = disabled */
    uint8_t  butteraugli_iters;  /* max iterations (default 3) */
    bool     grain_enabled;
    bool     ptd_enabled;        /* use PTD in RDO */
    uint8_t  speed_preset;       /* 0 (slowest) to 10 (fastest) */
} SntPerceptualConfig;

/* ── Encoder context ──────────────────────────────────────────────── */
typedef struct {
    SntFrameHeader     header;
    SntPerceptualConfig psy_cfg;
    SntFrame           *cur_frame;   /* current input frame */
    SntFrame           *rec_frame;   /* reconstructed (reference) */
    SntFrame           *ref_frames[4]; /* reference frame buffer */
    uint8_t            num_refs;
    /* Simple flat reference buffers for inter prediction (no padding) */
    uint8_t            *ref_y;
    uint8_t            *ref_u;
    uint8_t            *ref_v;
    int                 ref_w, ref_h;
    /* Cross-frame CDF state (opaque blob, managed by encoder) */
    void               *cdf_state;
    int                 cdf_state_valid;
} SntEncoderCtx;

/* ── Decoder context ──────────────────────────────────────────────── */
typedef struct {
    SntFrameHeader     header;
    SntFrame           *rec_frame;
    SntFrame           *ref_frames[4];
    uint8_t            num_refs;
    /* Simple flat reference buffers for inter prediction (no padding) */
    uint8_t            *ref_y;
    uint8_t            *ref_u;
    uint8_t            *ref_v;
    int                 ref_w, ref_h;
    /* Cross-frame CDF state (opaque blob, managed by decoder) */
    void               *cdf_state;
    int                 cdf_state_valid;
} SntDecoderCtx;

/* ── Frame buffer allocation/free ─────────────────────────────────── */
#include "../util/memory.h"

static inline SntFrame *snt_frame_alloc(uint32_t width, uint32_t height) {
    SntFrame *f = (SntFrame *)snt_calloc(1, sizeof(SntFrame));
    if (!f) return NULL;

    f->width  = width;
    f->height = height;
    f->pad    = SNT_MAX_SB_SIZE;

    uint32_t padded_w = width + 2 * f->pad;
    uint32_t padded_h = height + 2 * f->pad;

    /* Stride aligned to 64 bytes */
    f->stride_y  = (padded_w + 63) & ~63u;
    f->stride_uv = ((padded_w / 2) + 63) & ~63u;

    size_t y_size  = (size_t)f->stride_y * padded_h;
    size_t uv_size = (size_t)f->stride_uv * (padded_h / 2);

    f->y = (uint8_t *)snt_calloc(1, y_size);
    f->u = (uint8_t *)snt_calloc(1, uv_size);
    f->v = (uint8_t *)snt_calloc(1, uv_size);

    size_t coeff_y_size  = (size_t)f->stride_y * padded_h * sizeof(int16_t);
    size_t coeff_uv_size = (size_t)f->stride_uv * (padded_h / 2) * sizeof(int16_t);

    f->y_coeff = (int16_t *)snt_calloc(1, coeff_y_size);
    f->u_coeff = (int16_t *)snt_calloc(1, coeff_uv_size);
    f->v_coeff = (int16_t *)snt_calloc(1, coeff_uv_size);

    if (!f->y || !f->u || !f->v ||
        !f->y_coeff || !f->u_coeff || !f->v_coeff) {
        /* Partial alloc — caller should call snt_frame_free */
        return f;
    }
    return f;
}

static inline void snt_frame_free(SntFrame *f) {
    if (!f) return;
    snt_free(f->y);
    snt_free(f->u);
    snt_free(f->v);
    snt_free(f->y_coeff);
    snt_free(f->u_coeff);
    snt_free(f->v_coeff);
    snt_free(f);
}

/* Get pointer to pixel at (x, y) in luma plane, accounting for padding */
static inline uint8_t *snt_frame_pixel_y(SntFrame *f, uint32_t x, uint32_t y) {
    return f->y + (y + f->pad) * f->stride_y + (x + f->pad);
}

/* Get pointer to pixel at (x, y) in chroma plane */
static inline uint8_t *snt_frame_pixel_uv(uint8_t *plane, SntFrame *f,
                                            uint32_t x, uint32_t y) {
    uint32_t pad_uv = f->pad / 2;
    return plane + (y + pad_uv) * f->stride_uv + (x + pad_uv);
}

#endif /* SNT_TYPES_H */
