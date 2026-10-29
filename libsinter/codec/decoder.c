/*
 * Sinter Codec — Decoder with range coding + inter prediction
 * Mirrors encoder symbol order exactly.
 * No lapping post-filter: the encoder does not pre-filter (see encoder.c).
 */

#include <string.h>
#include <stdlib.h>
#include "decoder.h"
#include "cdf_set.h"
#ifdef SNT_EXPERIMENT_CLIPTABLES
#include "cliptables.h"
#endif
#ifdef SNT_EXPERIMENT_MIXING
#include "mixer.h"
/* Mixing state of the model set in use (one decoder per process) */
static SntMixState dec_mix[2];
static SntMixState *cur_mix;
#endif
#include "scan.h"
#include "../util/fixed_point.h"
#include "../core/transform/dct.h"
#include "../core/quantize/pvq.h"
#include "../core/quantize/qp_table.h"
#include "../core/bitstream/bitstream.h"
#include "../core/entropy/rangecoder.h"
#include "../core/prediction/motion_est.h"
#include "../core/prediction/mc.h"
#include "../core/filter/deblock.h"
#include "../core/prediction/intra_acdc.h"

#define BLK 8

/* Decode an MV delta of a block of size16 (0 8x8, 1 16x16) and add it to
 * *mv: x, then y, each as binary decisions (levels.h, MV models) */
static void dec_mvd(SntRangeDec *rd, SntCdfSet *cs, int size16, int ac, SntMV *mv) {
    int16_t d[2];
    for (int c = 0; c < 2; c++) {
        int xz = c && d[0] == 0;
        d[c] = 0;
        if (!snt_rc_decode_sym(rd, &cs->vz[c][size16][xz][ac])) continue;
        int neg = snt_rc_decode_sym(rd, &cs->vs[c]);
        int32_t mag = 1;
        if (snt_rc_decode_sym(rd, &cs->vg[c][size16][xz][ac])) {
            int k = 1;
            while (k < SNT_LEVEL_MAX_CAT && snt_rc_decode_sym(rd, &cs->vx[c][k - 1]))
                k++;
            uint32_t a1 = 1;
            for (int i = k - 2; i >= 0; i--)
                a1 = (a1 << 1) | (uint32_t)snt_rc_decode_sym(rd, &cs->vm[c][k]);
            mag = (int32_t)a1 + 1;
            if (mag > 32767) mag = 32767;  /* damaged stream */
        }
        d[c] = (int16_t)(neg ? -mag : mag);
    }
    mv->x = (int16_t)(mv->x + d[0]);
    mv->y = (int16_t)(mv->y + d[1]);
}

/* Prediction of a skipped block (bx, by): as the block to its left if
 * that block is inter-predicted (MPEG-1's rule for skipped macroblocks
 * in B-pictures), otherwise with the global motion */
static void skip_motion(const uint8_t *blk_dir, const SntMV *mvf_field, const SntMV *mvb_field,
                        int blk_cols, int bx, int by, int is_bidi, SntMV gmf, SntMV gmb,
                        int *dir, SntMV *mvf, SntMV *mvb) {
    int left = by * blk_cols + bx - 1;
    if (bx > 0 && blk_dir[left]) {
        *dir = blk_dir[left];
        *mvf = mvf_field[left];
        *mvb = mvb_field[left];
    } else {
        *dir = is_bidi ? SNT_DIR_BI : SNT_DIR_FWD;
        *mvf = gmf;
        *mvb = gmb;
    }
}

/* Prediction agreement of 8x8 block (bx, by) (cdf_set.h): P-frames compare
 * the skip prediction with the prediction at the forward predictor pf;
 * B-frames compare forward at pf with backward at pb. Zero for a
 * malformed stream that lacks a reference. */
static uint32_t agree_sad8(const SntSplinePlane *rf, const SntSplinePlane *rb,
                           const SntSplinePlane *r2, int bx, int by,
                           int is_bidi, int sdir, SntMV smf, SntMV smb, SntMV pf, SntMV pb) {
    SntMV zero = {0, 0};
    if (!rf->coef || (is_bidi && !rb->coef)) return 0;
    if (is_bidi)
        return snt_pred_sad8(rf, rb, r2, bx * BLK, by * BLK, SNT_DIR_FWD, pf, zero, SNT_DIR_BWD, zero, pb);
    return snt_pred_sad8(rf, rb, r2, bx * BLK, by * BLK, sdir, smf, smb, SNT_DIR_FWD, pf, zero);
}

/* Decode a luma AC level at zigzag index k from its binary decisions
 * (levels.h; the same decisions and models as snt_level_bins) */
#ifdef SNT_EXPERIMENT_MIXING
/* Z or G decision with the mixed probability (mixer.h) */
static int dec_mixed(SntRangeDec *rd, SntCdfSet *cs, int kind, int t, int k, int band, int pc,
                     int nzprev) {
    SntMixBin mb;
    snt_mix_prepare(cur_mix, cs, kind, t, k, band, pc, nzprev, &mb);
    int bit = snt_rc_decode_bit(rd, mb.p);
    snt_mix_update(&mb, bit);
    return bit;
}
#define DEC_Z(rd, cs, t, k, pc) dec_mixed(rd, cs, 0, t, k, band, pc, nzprev)
#define DEC_G(rd, cs, t, k, pc) dec_mixed(rd, cs, 1, t, k, band, pc, nzprev)
#else
#define DEC_Z(rd, cs, t, k, pc) snt_rc_decode_sym(rd, &(cs)->acz[t][k][pc])
#define DEC_G(rd, cs, t, k, pc) snt_rc_decode_sym(rd, &(cs)->acg[t][k][pc])
#endif

static int16_t dec_level(SntRangeDec *rd, SntCdfSet *cs, int t, int k, int pc, int z_known,
                         int band, int nzprev) {
    (void)band; (void)nzprev;
    if (!z_known && !DEC_Z(rd, cs, t, k, pc)) return 0;
    int neg = snt_rc_decode_sym(rd, &cs->acs[t]);
    if (!DEC_G(rd, cs, t, k, pc)) return (int16_t)(neg ? -1 : 1);
    int c = 1;
    while (c < SNT_LEVEL_MAX_CAT && snt_rc_decode_sym(rd, &cs->acx[t][k > SNT_LEVEL_KLOW][c - 1]))
        c++;
    uint32_t a1 = 1;
    for (int i = c - 2; i >= 0; i--)
        a1 = (a1 << 1) | (uint32_t)snt_rc_decode_sym(rd, &cs->acm[t][c]);
    int32_t mag = (int32_t)a1 + 1;
    if (mag > 32767) mag = 32767;  /* damaged stream */
    return (int16_t)(neg ? -mag : mag);
}

/* A luma DC level (levels.h, DC models) in DC context da (snt_dc_context) */
static int16_t dec_dclevel(SntRangeDec *rd, SntCdfSet *cs, int t, int da) {
    if (!snt_rc_decode_sym(rd, &cs->dz[t][da])) return 0;
    int neg = snt_rc_decode_sym(rd, &cs->ds[t][da]);
    if (!snt_rc_decode_sym(rd, &cs->dg[t][da][neg])) return (int16_t)(neg ? -1 : 1);
    int c = 1;
    while (c < SNT_LEVEL_MAX_CAT && snt_rc_decode_sym(rd, &cs->dx[t][c - 1]))
        c++;
    uint32_t a1 = 1;
    for (int k = c - 2; k >= 0; k--)
        a1 = (a1 << 1) | (uint32_t)snt_rc_decode_sym(rd, &cs->dm[t][c]);
    int32_t mag = (int32_t)a1 + 1;
    if (mag > 32767) mag = 32767;  /* damaged stream */
    return (int16_t)(neg ? -mag : mag);
}

/* A chroma level at zigzag index i (levels.h, chroma models) */
static int16_t dec_clevel(SntRangeDec *rd, SntCdfSet *cs, int t, int i, int z_known) {
    if (!z_known && !snt_rc_decode_sym(rd, &cs->cz[t][i])) return 0;
    int neg = snt_rc_decode_sym(rd, &cs->cs[t]);
    if (!snt_rc_decode_sym(rd, &cs->cg[t][i])) return (int16_t)(neg ? -1 : 1);
    int c = 1;
    while (c < SNT_LEVEL_MAX_CAT && snt_rc_decode_sym(rd, &cs->cx[t][i > 0][c - 1]))
        c++;
    uint32_t a1 = 1;
    for (int k = c - 2; k >= 0; k--)
        a1 = (a1 << 1) | (uint32_t)snt_rc_decode_sym(rd, &cs->cm[t][c]);
    int32_t mag = (int32_t)a1 + 1;
    if (mag > 32767) mag = 32767;  /* damaged stream */
    return (int16_t)(neg ? -mag : mag);
}

/* A P-frame inter block's reference flags (two or three references):
 * returns 0, SNT_DIR_REF1 or SNT_DIR_REF2. left_dir: the direction of the
 * block to the left (0 if intra or outside), the flags' context. */
static int dec_ref(SntRangeDec *rd, SntCdfSet *cs, int three, int left_dir) {
    int left = (left_dir & SNT_DIR_REF2) ? 2 : (left_dir & SNT_DIR_REF1) ? 1 : 0;
    if (!snt_rc_decode_sym(rd, &cs->ref1[left > 0])) return 0;
    if (three && snt_rc_decode_sym(rd, &cs->ref2[left == 2])) return SNT_DIR_REF2;
    return SNT_DIR_REF1;
}

/* Inter predictions (luma overlapped with the above and left blocks'
 * motion, snt_pred_luma8_obmc), or zeros for a malformed stream that lacks
 * a reference */
static void dec_pred_luma(const SntSplinePlane *rf, const SntSplinePlane *rb,
                          const SntSplinePlane *r2, int x0, int y0,
                          SntMotion cur, SntMotion above, SntMotion left, uint8_t *pred) {
    for (int k = 0; k < 3; k++) {
        int dir = (k == 0 ? cur : k == 1 ? above : left).dir;
        if (((dir & SNT_DIR_FWD) && !rf->coef) || ((dir & (SNT_DIR_BWD | SNT_DIR_REF1)) && !rb->coef) ||
            ((dir & SNT_DIR_REF2) && !r2->coef)) {
            memset(pred, 0, 64);
            return;
        }
    }
    snt_pred_luma8_obmc(rf, rb, r2, x0, y0, cur, above, left, pred);
}

static void dec_pred_chroma(const uint8_t *rf, const uint8_t *rb, const uint8_t *r2, int w, int h,
                            int x0, int y0, int dir, SntMV mvf, SntMV mvb, uint8_t *pred) {
    if (((dir & SNT_DIR_FWD) && !rf) || ((dir & (SNT_DIR_BWD | SNT_DIR_REF1)) && !rb) ||
        ((dir & SNT_DIR_REF2) && !r2)) {
        memset(pred, 128, 16);
        return;
    }
    snt_pred_chroma4(rf, rb, r2, w, h, x0, y0, dir, mvf, mvb, pred);
}

/*
 * Symbol consumption order (must match encoder exactly):
 *
 * FOR each 16x16 macroblock (raster):
 *   IF inter frame AND full macroblock: DEC(mbmode)
 *     IF INTER16: [B: DEC(bdir)] [P, two or three references: dec_ref], then an MV
 *                 delta per used direction
 *   FOR each 8x8 block in Z order (SKIP16/INTER16: no per-block mode symbols):
 *   IF inter frame (split macroblock):
 *     DEC(inter) → inter_flag
 *     IF inter_flag:
 *       DEC(iskip) → skip_flag
 *       IF skip: prediction with the global motion, continue (B: bi-directional)
 *       IF coded: [B: DEC(bdir)] [P, two or three references: dec_ref], then an MV
 *                 delta per used direction
 *   IF intra block: DEC(imode[is_inter_frame]) → mode (SNT_IPRED_*, intra_acdc.h)
 *   DC level's decisions (DC models, context: the left and above blocks' DC
 *     levels) → dc_level (residual for inter, absolute for intra)
 *   FOR band b=0..15: band[t][b][nz], then (if coded) ac[t][b/4] per level
 * FOR each chroma plane, each 4x4 block not skipped: cflag[ct][luma_nz][pctx],
 *   then (if coded) 16 levels in zigzag order (chroma models)
 * (contexts: see cdf_set.h and the encoder's symbol-order comment)
 */

int sntd_decode_frame(SntDecoderCtx *ctx, const uint8_t *bitstream,
                       uint32_t bs_len, uint8_t *output_y,
                       uint8_t *output_u, uint8_t *output_v,
                       int *out_width, int *out_height,
                       int stride_y, int stride_uv) {
    SntBitReader br;
    snt_br_init(&br, bitstream, bs_len);
    snt_read_frame_header(&br, &ctx->header);

    int width = (int)ctx->header.width;
    int height = (int)ctx->header.height;
    *out_width = width;
    *out_height = height;

    int pad_w = ((width + BLK - 1) / BLK) * BLK;
    int pad_h = ((height + BLK - 1) / BLK) * BLK;
    int blk_cols = pad_w / BLK;
    int blk_rows = pad_h / BLK;

    int is_inter = (ctx->header.frame_type != SNT_FRAME_KEY);
    int is_bidi = (ctx->header.frame_type == SNT_FRAME_BIDI);
    /* P-frames predict from the newest anchor; B-frames forward from the
     * older anchor and backward from the newest */
    int same_ref = ctx->ref_w == width && ctx->ref_h == height;
    SntSplinePlane sp_new = { ctx->ref_y, same_ref ? ctx->ref_coef : NULL, width, height };
    SntSplinePlane sp_old = { ctx->ref_old_y, same_ref ? ctx->ref_old_coef : NULL, width, height };
    SntSplinePlane sp_key = { ctx->key_y, same_ref ? ctx->key_coef : NULL, width, height };
    SntSplinePlane sp_none = { NULL, NULL, width, height };
    /* P-frames also predict from the anchor before the newest when it is in
     * the same GOP, as ref_b (SNT_DIR_REF1) */
    int has_ref1 = is_inter && !is_bidi && ctx->num_anchors >= 2 && sp_old.coef;
    const SntSplinePlane *ref_f = is_bidi ? &sp_old : &sp_new;
    /* From the third anchor on, also from the GOP's key frame, as ref_2
     * (SNT_DIR_REF2) */
    int has_ref2 = has_ref1 && ctx->num_anchors >= 3 && sp_key.coef;
    const SntSplinePlane *ref_b = is_bidi ? &sp_new : (has_ref1 ? &sp_old : &sp_none);
    const SntSplinePlane *ref_2 = has_ref2 ? &sp_key : &sp_none;
    uint16_t step = snt_qp_to_step(ctx->header.base_qp);

    /* Fixed fields: chroma present, global motion against each reference
     * (skipped blocks predict with it) */
    uint32_t chroma_present = snt_br_get_bits(&br, 1);
    SntMV gmf = {0, 0}, gmb = {0, 0};
    if (is_inter) {
        gmf.x = (int16_t)snt_br_get_sexpgolomb(&br);
        gmf.y = (int16_t)snt_br_get_sexpgolomb(&br);
    }
    if (is_bidi) {
        gmb.x = (int16_t)snt_br_get_sexpgolomb(&br);
        gmb.y = (int16_t)snt_br_get_sexpgolomb(&br);
    }
    snt_br_byte_align(&br);

    /* The range-coded bytes run to the end of the frame */
    uint32_t rc_start = br.byte_pos < bs_len ? br.byte_pos : bs_len;
    const uint8_t *rc_buf = bitstream + rc_start;
    uint32_t rc_bytes = bs_len - rc_start;

    /* ── Init range decoder ────────────────────────────────────────── */
    SntRangeDec rdec;
    snt_rc_dec_init(&rdec, rc_buf, rc_bytes);

    /* ── Entropy contexts: reset on key frames, carried across inter ──
     * B-frames keep their own set, carried from B-frame to B-frame. */
    if (!ctx->cdf_state)
        ctx->cdf_state = malloc(sizeof(SntCdfSet));
    if (!ctx->cdf_state_b)
        ctx->cdf_state_b = malloc(sizeof(SntCdfSet));
    if (!is_inter)
        ctx->cdf_state_b_valid = 0;
    SntCdfSet *cs = (SntCdfSet *)(is_bidi ? ctx->cdf_state_b : ctx->cdf_state);
    if (!is_inter || !(is_bidi ? ctx->cdf_state_b_valid : ctx->cdf_state_valid)) {
        snt_cdf_set_init(cs, is_bidi ? SNT_SET_B : SNT_SET_P, ctx->header.base_qp);
#ifdef SNT_EXPERIMENT_MIXING
        snt_mix_init(&dec_mix[is_bidi]);
#endif
#ifdef SNT_EXPERIMENT_CLIPTABLES
        if (ctx->clip_tables)
            snt_clip_apply(cs, (const SntClipTables *)ctx->clip_tables, is_bidi);
#endif
    }
#ifdef SNT_EXPERIMENT_CLIPTABLES
    /* Key frames start with a flag: per-clip tables follow */
    if (!is_inter && snt_rc_decode_bit(&rdec, 2048)) {
        if (!ctx->clip_tables) ctx->clip_tables = malloc(sizeof(SntClipTables));
        snt_clip_decode(&rdec, (SntClipTables *)ctx->clip_tables);
        snt_clip_apply(cs, (const SntClipTables *)ctx->clip_tables, SNT_SET_P);
    }
#endif
#ifdef SNT_EXPERIMENT_MIXING
    cur_mix = &dec_mix[is_bidi];
#endif

    #define DEC(m) ((uint16_t)snt_rc_decode_sym(&rdec, &(m)))

    /* ── Decode luma blocks ────────────────────────────────────────── */
    int16_t *recon = (int16_t *)calloc((size_t)pad_w * pad_h, sizeof(int16_t));
    int total_blocks_d = blk_cols * blk_rows;
    /* Per-block state for neighbours' prediction: coefficients of the
     * reconstruction (intra_acdc.h), and whether the block is coded yet */
    int16_t *acdc = (int16_t *)calloc((size_t)total_blocks_d * SNT_ACDC_STATE, sizeof(int16_t));
    uint8_t *coded = (uint8_t *)calloc(total_blocks_d, 1);
    uint8_t *luma_nz = (uint8_t *)calloc(total_blocks_d, 1);  /* any nonzero luma level? */
    int16_t *dc_lev = (int16_t *)calloc(total_blocks_d, sizeof(int16_t));  /* luma DC levels */

    /* MV field and block flags for chroma inter derivation */
    SntMV *mv_field = (SntMV *)calloc(blk_cols * blk_rows, sizeof(SntMV));   /* forward MVs */
    SntMV *mvb_field = (SntMV *)calloc(blk_cols * blk_rows, sizeof(SntMV));  /* backward MVs */
    uint8_t *blk_dir = (uint8_t *)calloc(blk_cols * blk_rows, 1);            /* SNT_DIR_*, 0 intra */
    int *blk_is_inter = (int *)calloc(blk_cols * blk_rows, sizeof(int));
    int *blk_is_skip = (int *)calloc(blk_cols * blk_rows, sizeof(int));

    #define NUM_BANDS 16
    static const int band_bounds[NUM_BANDS + 1] = {0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,59,63};

    /* Blocks are coded in 16x16 macroblocks in raster order, and the four
     * 8x8 blocks of a macroblock in Z order (top-left, top-right,
     * bottom-left, bottom-right). */
    int mb_cols = (blk_cols + 1) / 2, mb_rows = (blk_rows + 1) / 2;
    for (int mb = 0; mb < mb_cols * mb_rows; mb++) {
        /* Macroblock mode (inter frames, full macroblocks only). The 16x16
         * predictor comes from neighbours outside the macroblock. */
        int bx0 = (mb % mb_cols) * 2, by0 = (mb / mb_cols) * 2;
        int mb_mode = SNT_MB_SPLIT, mb_dir = SNT_DIR_FWD, mac = 0;
        SntMV mb_mvf = {0, 0}, mb_mvb = {0, 0};
        if (is_inter && bx0 + 1 < blk_cols && by0 + 1 < blk_rows) {
            mb_mvf = snt_mv_predict_at(mv_field, blk_cols, bx0, by0);
            if (is_bidi)
                mb_mvb = snt_mv_predict_at(mvb_field, blk_cols, bx0, by0);
            /* The mode's model: agreement of two predictions over the
             * macroblock, with its skip motion and predictors */
            int sdir;
            SntMV smf, smb;
            skip_motion(blk_dir, mv_field, mvb_field, blk_cols, bx0, by0, is_bidi,
                        gmf, gmb, &sdir, &smf, &smb);
            uint32_t sad = 0;
            for (int z = 0; z < 4; z++)
                sad += agree_sad8(ref_f, ref_b, ref_2, bx0 + (z & 1), by0 + (z >> 1), is_bidi,
                                  sdir, smf, smb, mb_mvf, mb_mvb);
            mac = snt_agree_class(sad, 256, step);
            mb_mode = DEC(cs->mbmode[snt_flag_context(mac, bx0 > 0 && blk_is_skip[by0 * blk_cols + bx0 - 1],
                                                      by0 > 0 && blk_is_skip[(by0 - 1) * blk_cols + bx0])]);
            /* SKIP16 predicts with the global motion: bi-directional in
             * B-frames */
            mb_dir = is_bidi ? SNT_DIR_BI : SNT_DIR_FWD;
            if (mb_mode == SNT_MB_SKIP16)
                skip_motion(blk_dir, mv_field, mvb_field, blk_cols, bx0, by0, is_bidi,
                            gmf, gmb, &mb_dir, &mb_mvf, &mb_mvb);
            if (mb_mode == SNT_MB_INTER16) {
                mb_dir = is_bidi ? DEC(cs->bdir) + 1 : SNT_DIR_FWD;
                if (has_ref1)
                    mb_dir |= dec_ref(&rdec, cs, has_ref2, bx0 > 0 ? blk_dir[by0 * blk_cols + bx0 - 1] : 0);
                /* An unused direction's MV is zero (H.263 Annex O) */
                if (!(mb_dir & SNT_DIR_FWD)) mb_mvf = (SntMV){0, 0};
                if (!(mb_dir & SNT_DIR_BWD)) mb_mvb = (SntMV){0, 0};
                if (mb_dir & SNT_DIR_FWD) dec_mvd(&rdec, cs, 1, mac, &mb_mvf);
                if (mb_dir & SNT_DIR_BWD) dec_mvd(&rdec, cs, 1, mac, &mb_mvb);
            }
        }
        for (int z = 0; z < 4; z++) {
            int bx = (mb % mb_cols) * 2 + (z & 1), by = (mb / mb_cols) * 2 + (z >> 1);
            if (bx >= blk_cols || by >= blk_rows) continue;
            int blk_idx = by * blk_cols + bx;
            int px0 = bx * BLK, py0 = by * BLK;
            int block_is_inter = 0;
            int block_is_skip = 0;
            int dir = 0;
            SntMV best_mv = {0, 0}, best_mvb = {0, 0};
            uint8_t mc_blk[BLK * BLK];  /* MC prediction (skip and inter coded blocks) */

            /* ── Inter mode decode ────────────────────────────── */
            if (is_inter) {
                if (mb_mode == SNT_MB_SPLIT) {
                    /* The flags' models: agreement of two predictions */
                    int sdir;
                    SntMV smf, smb, pf, pb = {0, 0};
                    skip_motion(blk_dir, mv_field, mvb_field, blk_cols, bx, by, is_bidi,
                                gmf, gmb, &sdir, &smf, &smb);
                    pf = snt_mv_predict_at(mv_field, blk_cols, bx, by);
                    if (is_bidi)
                        pb = snt_mv_predict_at(mvb_field, blk_cols, bx, by);
                    int ac = snt_agree_class(agree_sad8(ref_f, ref_b, ref_2, bx, by, is_bidi,
                                                        sdir, smf, smb, pf, pb), 64, step);
                    int fc = snt_flag_context(ac, bx > 0 && blk_is_skip[blk_idx - 1],
                                              by > 0 && blk_is_skip[blk_idx - blk_cols]);
                    block_is_inter = DEC(cs->inter[fc]) == 1;
                    if (block_is_inter) {
                        /* A skip block predicts with the global motion; a
                         * coded block's MVs are coded against the median of
                         * its neighbours' (an unused direction's MV is zero) */
                        block_is_skip = DEC(cs->iskip[fc]) == 1;
                        dir = is_bidi ? SNT_DIR_BI : SNT_DIR_FWD;
                        if (block_is_skip)
                            skip_motion(blk_dir, mv_field, mvb_field, blk_cols, bx, by, is_bidi,
                                        gmf, gmb, &dir, &best_mv, &best_mvb);
                        if (!block_is_skip) {
                            best_mv = snt_mv_predict_at(mv_field, blk_cols, bx, by);
                            if (is_bidi)
                                best_mvb = snt_mv_predict_at(mvb_field, blk_cols, bx, by);
                            if (is_bidi)
                                dir = DEC(cs->bdir) + 1;
                            else
                                dir = SNT_DIR_FWD;
                            if (has_ref1)
                                dir |= dec_ref(&rdec, cs, has_ref2, bx > 0 ? blk_dir[blk_idx - 1] : 0);
                            if (!(dir & SNT_DIR_FWD)) best_mv = (SntMV){0, 0};
                            if (!(dir & SNT_DIR_BWD)) best_mvb = (SntMV){0, 0};
                            if (dir & SNT_DIR_FWD) dec_mvd(&rdec, cs, 0, ac, &best_mv);
                            if (dir & SNT_DIR_BWD) dec_mvd(&rdec, cs, 0, ac, &best_mvb);
                        }
                    }
                } else {
                    /* 16x16 skip or inter: every block takes the MB's MVs */
                    block_is_inter = 1;
                    block_is_skip = (mb_mode == SNT_MB_SKIP16);
                    dir = mb_dir;
                    best_mv = mb_mvf;
                    best_mvb = mb_mvb;
                }
                blk_is_inter[blk_idx] = block_is_inter;
                blk_is_skip[blk_idx] = block_is_skip;
                if (block_is_inter) {
                    mv_field[blk_idx] = best_mv;
                    mvb_field[blk_idx] = best_mvb;
                    blk_dir[blk_idx] = (uint8_t)dir;
                    /* Overlapped with the motion of the inter blocks above
                     * and to the left (an intra or missing neighbour counts
                     * as this block's motion) */
                    SntMotion cur = { dir, best_mv, best_mvb }, am = cur, lm = cur;
                    int ai = blk_idx - blk_cols, li = blk_idx - 1;
                    if (by > 0 && blk_dir[ai]) am = (SntMotion){ blk_dir[ai], mv_field[ai], mvb_field[ai] };
                    if (bx > 0 && blk_dir[li]) lm = (SntMotion){ blk_dir[li], mv_field[li], mvb_field[li] };
                    dec_pred_luma(ref_f, ref_b, ref_2, px0, py0, cur, am, lm, mc_blk);
                }

                if (block_is_skip) {
                    /* Copy the prediction to recon */
                    int16_t skip_blk[64];
                    for (int dy = 0; dy < BLK; dy++)
                        for (int dx = 0; dx < BLK; dx++) {
                            int py = py0 + dy, px_v = px0 + dx;
                            skip_blk[dy * BLK + dx] = mc_blk[dy * BLK + dx];
                            if (py < height && px_v < width)
                                recon[py * pad_w + px_v] = mc_blk[dy * BLK + dx];
                        }
                    /* Neighbours' intra prediction reads this block's
                     * coefficients */
                    int16_t c[64];
                    snt_dct_fwd_2d_ref(skip_blk, BLK, BLK, c);
                    snt_acdc_store(c, &acdc[blk_idx * SNT_ACDC_STATE]);
                    coded[blk_idx] = 1;
                    continue;
                }
            }

            /* Intra: coefficient prediction from the neighbours in the
             * block's mode; inter: the MC prediction is in mc_blk */
            int16_t td_pred[64];
            memset(td_pred, 0, sizeof(td_pred));
            if (!block_is_inter) {
                int imode = DEC(cs->imode[is_inter]);
                if (imode >= SNT_NUM_IPRED) imode = SNT_IPRED_DC;  /* damaged stream */
                snt_acdc_predict(imode,
                                 bx > 0 ? &acdc[(blk_idx - 1) * SNT_ACDC_STATE] : NULL,
                                 by > 0 ? &acdc[(blk_idx - blk_cols) * SNT_ACDC_STATE] : NULL,
                                 td_pred);
            }
            int has_pred_blk = block_is_inter && !block_is_skip;

            /* ── Decode DC + bands (same for intra and inter coded) ── */
            int16_t recon_coeffs[64];

            /* DC prediction: intra blocks from their mode's prediction;
             * inter residuals are zero-mean */
            int16_t dc_pred = has_pred_blk ? 0 : td_pred[0];

            int16_t dc_level = dec_dclevel(&rdec, cs, block_is_inter,
                                           snt_dc_context(bx > 0 ? dc_lev[blk_idx - 1] : 0,
                                                          by > 0 ? dc_lev[blk_idx - blk_cols] : 0));
            dc_lev[blk_idx] = dc_level;

            int16_t dc_recon_val = pvq_dequantize_dc(dc_level, step) + dc_pred;
            recon_coeffs[0] = dc_recon_val;

            /* Bands: a coded flag (context: is the previous band, for
             * band 0 the DC level, nonzero?), then the levels */
            /* Contexts from the prediction (inter blocks; levels.h) */
            uint8_t kcls[64] = {0}, bcls[NUM_BANDS] = {0};
            if (has_pred_blk) {
                int16_t p16[64], pdct[64];
                for (int i = 0; i < 64; i++) p16[i] = mc_blk[i];
                snt_dct_fwd_2d_ref(p16, BLK, BLK, pdct);
                snt_pred_classes(pdct, step, snt_zigzag8, band_bounds, NUM_BANDS, kcls, bcls);
            }
            int prev_nz = dc_level != 0;
            luma_nz[blk_idx] = dc_level != 0;
            int16_t lv[64] = {0};  /* levels so far, zigzag order */
            lv[0] = dc_level;
            for (int b = 0; b < NUM_BANDS; b++) {
                int start = band_bounds[b];
                int n = band_bounds[b + 1] - start;
                int coded = DEC(cs->band[block_is_inter][b][prev_nz][bcls[b]]) == 1;
                int any = 0;
                for (int i = 0; i < n; i++) {
                    int16_t level = 0;
                    if (coded) {
                        /* The last level is known nonzero if the others are 0 */
                        int k = 1 + start + i;
                        level = dec_level(&rdec, cs, block_is_inter, k, kcls[k], i == n - 1 && !any,
                                          b, (lv[k - 1] != 0) + (k >= 2 && lv[k - 2] != 0));
                        lv[k] = level;
                        any |= level != 0;
                    }
                    recon_coeffs[snt_zigzag8[1 + start + i]] = pvq_dequantize_dc(level, step);
                }
                prev_nz = coded;
                luma_nz[blk_idx] |= coded;
            }

            /* Add back the intra AC prediction (first row or column) */
            if (!block_is_inter)
                for (int i = 1; i < 64; i++)
                    recon_coeffs[i] += td_pred[i];

            /* IDCT */
            int16_t recon_blk[64];
            snt_dct_inv_2d_ref(recon_coeffs, BLK, recon_blk, BLK);

            /* Write to recon: add the prediction block (inter or spatial intra) */
            for (int dy = 0; dy < BLK; dy++)
                for (int dx = 0; dx < BLK; dx++) {
                    int py = py0 + dy, px_v = px0 + dx;
                    int16_t val = recon_blk[dy * BLK + dx];
                    if (has_pred_blk)
                        val += mc_blk[dy * BLK + dx];
                    recon[py * pad_w + px_v] = val;
                    recon_blk[dy * BLK + dx] = val;
                }

            /* Neighbours' intra prediction reads this block's coefficients */
            int16_t c[64];
            snt_dct_fwd_2d_ref(recon_blk, BLK, BLK, c);
            snt_acdc_store(c, &acdc[blk_idx * SNT_ACDC_STATE]);
            coded[blk_idx] = 1;
        }
    }

    /* Copy to output */
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
            output_y[y * stride_y + x] = snt_clip_pixel((int32_t)recon[y * pad_w + x]);
    snt_deblock_luma(output_y, width, height, stride_y, step,
                         is_inter ? blk_is_skip : NULL, blk_cols);

    /* Decode chroma with inter prediction support */
    if (chroma_present && output_u && output_v) {
        int cw = width / 2, ch_h = height / 2;
        int cblk = 4;
        int cpad_w = ((cw + cblk - 1) / cblk) * cblk;
        int cpad_h = ((ch_h + cblk - 1) / cblk) * cblk;
        int cblk_cols = cpad_w / cblk;
        /* Reduced chroma QP offset for inter frames (matches encoder) */
        int cqp_off = is_inter ? 1 : ctx->header.chroma_qp_offset;
        int cqp = ctx->header.base_qp + cqp_off;
        if (cqp < 0) cqp = 0; if (cqp > 63) cqp = 63;
        uint16_t cstep = snt_qp_to_step((uint8_t)cqp);

        uint8_t *chroma_out[2] = {output_u, output_v};
        uint8_t *u_coded = (uint8_t *)calloc((size_t)cblk_cols * (cpad_h / cblk), 1);
        const uint8_t *chroma_ref[2] = {is_bidi ? ctx->ref_old_u : ctx->ref_u,
                                         is_bidi ? ctx->ref_old_v : ctx->ref_v};
        const uint8_t *chroma_ref_b[2] = {is_bidi ? ctx->ref_u : (has_ref1 ? ctx->ref_old_u : NULL),
                                           is_bidi ? ctx->ref_v : (has_ref1 ? ctx->ref_old_v : NULL)};
        const uint8_t *chroma_ref_2[2] = {has_ref2 ? ctx->key_u : NULL, has_ref2 ? ctx->key_v : NULL};
        for (int plane = 0; plane < 2; plane++) {
            const uint8_t *cr = chroma_ref[plane];
            const uint8_t *crb = chroma_ref_b[plane];
            int16_t c_dc_left = 0;
            int16_t *c_dc_above = (int16_t *)calloc(cblk_cols, sizeof(int16_t));

            for (int cby = 0; cby < cpad_h; cby += cblk) {
                c_dc_left = 0;
                for (int cbx = 0; cbx < cpad_w; cbx += cblk) {
                    int ci = cbx / cblk;
                    int cj = cby / cblk;

                    /* Check co-located luma block's inter/skip status */
                    int luma_blk_idx = (ci < blk_cols && cj < blk_rows) ?
                        cj * blk_cols + ci : 0;
                    int c_is_skip = is_inter && blk_is_skip[luma_blk_idx];
                    int c_is_inter = is_inter && blk_is_inter[luma_blk_idx];

                    /* Chroma skip: MC prediction at the luma block's MV */
                    if (c_is_skip && cr) {
                        uint8_t c_skip[16];
                        dec_pred_chroma(cr, crb, chroma_ref_2[plane], cw, ch_h, cbx, cby, blk_dir[luma_blk_idx],
                                        mv_field[luma_blk_idx], mvb_field[luma_blk_idx], c_skip);
                        for (int dy = 0; dy < 4; dy++)
                            for (int dx = 0; dx < 4; dx++) {
                                int py = cby + dy, px = cbx + dx;
                                if (py < ch_h && px < cw)
                                    chroma_out[plane][py*stride_uv+px] = c_skip[dy*4+dx];
                            }
                        c_dc_left = snt_pred_dc4c(c_skip);
                        c_dc_above[ci] = c_dc_left;
                        continue;
                    }

                    int c_bt = (c_is_inter && cr) ? 1 : 0;  /* chroma context: inter-predicted? */
                    uint8_t c_mc[16];  /* chroma MC prediction: luma MVs as 1/8-pel */
                    if (c_bt)
                        dec_pred_chroma(cr, crb, chroma_ref_2[plane], cw, ch_h, cbx, cby, blk_dir[luma_blk_idx],
                                        mv_field[luma_blk_idx], mvb_field[luma_blk_idx], c_mc);
                    /* DC prediction: intra blocks from the neighbours'
                     * pixel-domain DC; inter residuals are zero-mean */
                    int16_t c_dc_pred = 0;
                    if (c_bt)
                        c_dc_pred = 0;
                    else if (cbx > 0 && cby > 0)
                        c_dc_pred = (c_dc_left + c_dc_above[ci] + 1) >> 1;
                    else if (cbx > 0) c_dc_pred = c_dc_left;
                    else if (cby > 0) c_dc_pred = c_dc_above[ci];

                    /* Coded flag, then the levels in zigzag order */
                    int bi = cj * cblk_cols + ci;
                    int pctx = plane == 0 ? 0 : 1 + u_coded[bi];
                    int cc = DEC(cs->cflag[c_bt][luma_nz[luma_blk_idx]][pctx]);
                    if (plane == 0) u_coded[bi] = (uint8_t)cc;
                    int16_t ccoeffs[16] = {0};
                    if (cc) {
                        int any = 0;
                        for (int i = 0; i < 16; i++) {
                            int16_t level = dec_clevel(&rdec, cs, c_bt, i, i == 15 && !any);
                            ccoeffs[snt_zigzag4[i]] = pvq_dequantize_dc(level, cstep);
                            any |= level != 0;
                        }
                    }
                    ccoeffs[0] += c_dc_pred;
                    int16_t cdc_recon = ccoeffs[0];
                    /* Neighbours see the pixel-domain DC (inter: add the prediction's) */
                    c_dc_left = c_bt ? cdc_recon + snt_pred_dc4c(c_mc) : cdc_recon;
                    c_dc_above[ci] = c_dc_left;

                    int16_t cblock[16];
                    snt_dct_inv_2d_ref(ccoeffs, 4, cblock, 4);

                    /* Output: add MC reference for inter, add 128 for intra */
                    for (int dy = 0; dy < 4; dy++)
                        for (int dx = 0; dx < 4; dx++) {
                            int py = cby + dy, px = cbx + dx;
                            if (py < ch_h && px < cw) {
                                int32_t val = (int32_t)cblock[dy*4+dx];
                                if (c_is_inter && cr) {
                                    val += c_mc[dy*4+dx];
                                } else {
                                    val += 128;
                                }
                                chroma_out[plane][py*stride_uv+px] = snt_clip_pixel(val);
                            }
                        }
                }
            }
            free(c_dc_above);
        }
        free(u_coded);
        snt_deblock_chroma(output_u, cw, ch_h, stride_uv, cstep, is_inter ? blk_is_skip : NULL, blk_cols);
        snt_deblock_chroma(output_v, cw, ch_h, stride_uv, cstep, is_inter ? blk_is_skip : NULL, blk_cols);
    }

    /* Models now carry into the next frame of the same kind */
    if (is_bidi)
        ctx->cdf_state_b_valid = 1;
    else
        ctx->cdf_state_valid = 1;

    /* ── Anchors (key and P-frames) become the newest reference; the
     *    previous newest becomes the older one. B-frames are not stored. ── */
    if (!is_bidi) {
        int cw_dec = width / 2, ch_dec = height / 2;
        free(ctx->ref_old_y); free(ctx->ref_old_u); free(ctx->ref_old_v);
        free(ctx->ref_old_coef);
        int same_size = ctx->ref_w == width && ctx->ref_h == height;
        ctx->ref_old_y = same_size ? ctx->ref_y : NULL;
        ctx->ref_old_u = same_size ? ctx->ref_u : NULL;
        ctx->ref_old_v = same_size ? ctx->ref_v : NULL;
        ctx->ref_old_coef = same_size ? ctx->ref_coef : NULL;
        if (!same_size) {
            free(ctx->ref_y); free(ctx->ref_u); free(ctx->ref_v); free(ctx->ref_coef);
        }
        ctx->ref_y = (uint8_t *)malloc((size_t)width * height);
        ctx->ref_u = (uint8_t *)calloc((size_t)cw_dec * ch_dec, 1);
        ctx->ref_v = (uint8_t *)calloc((size_t)cw_dec * ch_dec, 1);
        ctx->ref_w = width;
        ctx->ref_h = height;
        for (int y = 0; y < height; y++)
            for (int x = 0; x < width; x++)
                ctx->ref_y[y * width + x] = output_y[y * stride_y + x];
        ctx->ref_coef = snt_spline_coefs(ctx->ref_y, width, height);
        if (output_u && output_v) {
            for (int y = 0; y < ch_dec; y++)
                for (int x = 0; x < cw_dec; x++) {
                    ctx->ref_u[y * cw_dec + x] = output_u[y * stride_uv + x];
                    ctx->ref_v[y * cw_dec + x] = output_v[y * stride_uv + x];
                }
        }
        ctx->num_refs = 1;
        ctx->num_anchors = is_inter ? ctx->num_anchors + 1 : 1;
        /* A key frame is also kept as the GOP's long-term reference */
        if (!is_inter) {
            free(ctx->key_y); free(ctx->key_u); free(ctx->key_v); free(ctx->key_coef);
            ctx->key_y = (uint8_t *)malloc((size_t)width * height);
            ctx->key_u = (uint8_t *)malloc((size_t)cw_dec * ch_dec);
            ctx->key_v = (uint8_t *)malloc((size_t)cw_dec * ch_dec);
            memcpy(ctx->key_y, ctx->ref_y, (size_t)width * height);
            memcpy(ctx->key_u, ctx->ref_u, (size_t)cw_dec * ch_dec);
            memcpy(ctx->key_v, ctx->ref_v, (size_t)cw_dec * ch_dec);
            ctx->key_coef = snt_spline_coefs(ctx->key_y, width, height);
        }
    }

    free(recon);
    free(acdc);
    free(coded);
    free(luma_nz);
    free(dc_lev);
    free(mv_field);
    free(mvb_field);
    free(blk_dir);
    free(blk_is_inter);
    free(blk_is_skip);
    return 0;
}
