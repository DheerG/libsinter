/*
 * Sinter Codec — Bitstream Reader/Writer
 * Implements bitstream syntax elements.
 *
 * Frame header layout:
 *   [32] width  [32] height  [1] frame_type  [6] base_qp
 *   [4] num_tile_cols  [4] num_tile_rows  [1] dual_rans
 *   [8] chroma_qp_offset (signed)  [32] frame_number
 *   [1] grain_enabled  [grain params if enabled]
 *
 * Per-block QP deltas are rANS-coded (13 symbols for [-6..+6]).
 * Byte-align after each entropy-coded segment (Gary's requirement).
 */

#include <string.h>
#include "bitstream.h"

/* ── Writer ───────────────────────────────────────────────────────── */

void snt_bw_init(SntBitWriter *bw, uint8_t *buf, uint32_t capacity) {
    bw->buf = buf;
    bw->capacity = capacity;
    bw->byte_pos = 0;
    bw->bit_pos = 0;
    bw->bits_written = 0;
    memset(buf, 0, capacity);
}

void snt_bw_put_bits(SntBitWriter *bw, uint32_t val, int nbits) {
    for (int i = nbits - 1; i >= 0; i--) {
        if (bw->byte_pos >= bw->capacity) return;
        uint8_t bit = (val >> i) & 1;
        bw->buf[bw->byte_pos] |= (bit << (7 - bw->bit_pos));
        bw->bit_pos++;
        bw->bits_written++;
        if (bw->bit_pos == 8) {
            bw->bit_pos = 0;
            bw->byte_pos++;
        }
    }
}

void snt_bw_put_byte(SntBitWriter *bw, uint8_t val) {
    snt_bw_byte_align(bw);
    if (bw->byte_pos < bw->capacity) {
        bw->buf[bw->byte_pos++] = val;
        bw->bits_written += 8;
    }
}

void snt_bw_byte_align(SntBitWriter *bw) {
    if (bw->bit_pos > 0) {
        bw->bits_written += (8 - bw->bit_pos);
        bw->bit_pos = 0;
        bw->byte_pos++;
    }
}

uint32_t snt_bw_bytes_written(const SntBitWriter *bw) {
    return bw->byte_pos + (bw->bit_pos > 0 ? 1 : 0);
}

/* ── Reader ───────────────────────────────────────────────────────── */

void snt_br_init(SntBitReader *br, const uint8_t *buf, uint32_t len) {
    br->buf = buf;
    br->buf_len = len;
    br->byte_pos = 0;
    br->bit_pos = 0;
    br->bits_read = 0;
}

uint32_t snt_br_get_bits(SntBitReader *br, int nbits) {
    uint32_t val = 0;
    for (int i = 0; i < nbits; i++) {
        if (br->byte_pos >= br->buf_len) return val;
        uint8_t bit = (br->buf[br->byte_pos] >> (7 - br->bit_pos)) & 1;
        val = (val << 1) | bit;
        br->bit_pos++;
        br->bits_read++;
        if (br->bit_pos == 8) {
            br->bit_pos = 0;
            br->byte_pos++;
        }
    }
    return val;
}

uint8_t snt_br_get_byte(SntBitReader *br) {
    snt_br_byte_align(br);
    if (br->byte_pos < br->buf_len) {
        br->bits_read += 8;
        return br->buf[br->byte_pos++];
    }
    return 0;
}

void snt_br_byte_align(SntBitReader *br) {
    if (br->bit_pos > 0) {
        br->bits_read += (8 - br->bit_pos);
        br->bit_pos = 0;
        br->byte_pos++;
    }
}

int snt_br_eof(const SntBitReader *br) {
    return br->byte_pos >= br->buf_len;
}

/* ── Frame header ─────────────────────────────────────────────────── */

void snt_write_frame_header(SntBitWriter *bw, const SntFrameHeader *hdr) {
    snt_bw_put_bits(bw, hdr->width, 32);
    snt_bw_put_bits(bw, hdr->height, 32);
    snt_bw_put_bits(bw, hdr->frame_type, 1);
    snt_bw_put_bits(bw, hdr->base_qp, 6);
    snt_bw_put_bits(bw, hdr->num_tile_cols, 4);
    snt_bw_put_bits(bw, hdr->num_tile_rows, 4);
    snt_bw_put_bits(bw, hdr->dual_rans ? 1 : 0, 1);
    snt_bw_put_bits(bw, (uint32_t)(uint8_t)hdr->chroma_qp_offset, 8);
    snt_bw_put_bits(bw, hdr->frame_number, 32);
    snt_bw_put_bits(bw, hdr->grain.enabled ? 1 : 0, 1);
    if (hdr->grain.enabled) {
        snt_write_grain_params(bw, &hdr->grain);
    }
    snt_bw_byte_align(bw);
}

void snt_read_frame_header(SntBitReader *br, SntFrameHeader *hdr) {
    hdr->width = snt_br_get_bits(br, 32);
    hdr->height = snt_br_get_bits(br, 32);
    hdr->frame_type = (SntFrameType)snt_br_get_bits(br, 1);
    hdr->base_qp = (uint8_t)snt_br_get_bits(br, 6);
    hdr->num_tile_cols = (uint8_t)snt_br_get_bits(br, 4);
    hdr->num_tile_rows = (uint8_t)snt_br_get_bits(br, 4);
    hdr->dual_rans = snt_br_get_bits(br, 1) != 0;
    hdr->chroma_qp_offset = (int8_t)(uint8_t)snt_br_get_bits(br, 8);
    hdr->frame_number = snt_br_get_bits(br, 32);
    hdr->grain.enabled = snt_br_get_bits(br, 1) != 0;
    if (hdr->grain.enabled) {
        snt_read_grain_params(br, &hdr->grain);
    }
    snt_br_byte_align(br);
}

/* ── Grain parameters ─────────────────────────────────────────────── */

void snt_write_grain_params(SntBitWriter *bw, const SntGrainParams *g) {
    /* AR(3) luma coefficients: 3 x 16-bit signed */
    for (int i = 0; i < 3; i++)
        snt_bw_put_bits(bw, (uint32_t)(uint16_t)g->ar_coeff_y[i], 16);
    /* AR(4) chroma-U: 4 x 16-bit signed */
    for (int i = 0; i < 4; i++)
        snt_bw_put_bits(bw, (uint32_t)(uint16_t)g->ar_coeff_u[i], 16);
    /* AR(4) chroma-V: 4 x 16-bit signed */
    for (int i = 0; i < 4; i++)
        snt_bw_put_bits(bw, (uint32_t)(uint16_t)g->ar_coeff_v[i], 16);
    /* Strengths: 3 x 8-bit */
    snt_bw_put_bits(bw, g->strength_y, 8);
    snt_bw_put_bits(bw, g->strength_u, 8);
    snt_bw_put_bits(bw, g->strength_v, 8);
    /* Seed: 16-bit */
    snt_bw_put_bits(bw, g->seed, 16);
}

void snt_read_grain_params(SntBitReader *br, SntGrainParams *g) {
    for (int i = 0; i < 3; i++)
        g->ar_coeff_y[i] = (int16_t)(uint16_t)snt_br_get_bits(br, 16);
    for (int i = 0; i < 4; i++)
        g->ar_coeff_u[i] = (int16_t)(uint16_t)snt_br_get_bits(br, 16);
    for (int i = 0; i < 4; i++)
        g->ar_coeff_v[i] = (int16_t)(uint16_t)snt_br_get_bits(br, 16);
    g->strength_y = (uint8_t)snt_br_get_bits(br, 8);
    g->strength_u = (uint8_t)snt_br_get_bits(br, 8);
    g->strength_v = (uint8_t)snt_br_get_bits(br, 8);
    g->seed = (uint16_t)snt_br_get_bits(br, 16);
}
