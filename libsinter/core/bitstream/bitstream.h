/*
 * Sinter Codec — Bitstream Reader/Writer
 * Implements bitstream syntax for frame headers, tile headers, and block data.
 * Explicit methods per syntax element (Gary's requirement).
 */
#ifndef SNT_BITSTREAM_H
#define SNT_BITSTREAM_H

#include <stdint.h>
#include "../sinter_types.h"

/* ── Bitstream Writer ─────────────────────────────────────────────── */
typedef struct {
    uint8_t  *buf;
    uint32_t  capacity;
    uint32_t  byte_pos;
    uint8_t   bit_pos;   /* bits written in current byte (0-7) */
    uint32_t  bits_written;
} SntBitWriter;

void snt_bw_init(SntBitWriter *bw, uint8_t *buf, uint32_t capacity);
void snt_bw_put_bits(SntBitWriter *bw, uint32_t val, int nbits);
void snt_bw_put_byte(SntBitWriter *bw, uint8_t val);
void snt_bw_byte_align(SntBitWriter *bw);
uint32_t snt_bw_bytes_written(const SntBitWriter *bw);

/* ── Bitstream Reader ─────────────────────────────────────────────── */
typedef struct {
    const uint8_t *buf;
    uint32_t  buf_len;
    uint32_t  byte_pos;
    uint8_t   bit_pos;
    uint32_t  bits_read;
} SntBitReader;

void snt_br_init(SntBitReader *br, const uint8_t *buf, uint32_t len);
uint32_t snt_br_get_bits(SntBitReader *br, int nbits);
uint8_t  snt_br_get_byte(SntBitReader *br);
void snt_br_byte_align(SntBitReader *br);
int snt_br_eof(const SntBitReader *br);

/* ── Exp-Golomb variable-length coding (unsigned, order 0) ────────── */
/* Small values use few bits: 0→1bit, 1-2→3bits, 3-6→5bits, etc.
 * Patent-free: prior art in H.261 (1988). */
void snt_bw_put_expgolomb(SntBitWriter *bw, uint32_t val);
uint32_t snt_br_get_expgolomb(SntBitReader *br);

/* Signed exp-Golomb: interleave positive/negative → unsigned */
void snt_bw_put_sexpgolomb(SntBitWriter *bw, int32_t val);
int32_t snt_br_get_sexpgolomb(SntBitReader *br);

/* ── Frame header syntax ──────────────────────────────────────────── */
void snt_write_frame_header(SntBitWriter *bw, const SntFrameHeader *hdr);
void snt_read_frame_header(SntBitReader *br, SntFrameHeader *hdr);

/* ── Grain parameters syntax ──────────────────────────────────────── */
void snt_write_grain_params(SntBitWriter *bw, const SntGrainParams *grain);
void snt_read_grain_params(SntBitReader *br, SntGrainParams *grain);

#endif /* SNT_BITSTREAM_H */
