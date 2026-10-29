/*
 * Sinter Codec — Range coder with adaptive frequency-count models (normative)
 *
 * Range coding: G. N. N. Martin, "Range encoding: an algorithm for removing
 * redundancy from a digitised message", Video & Data Recording Conference,
 * Southampton, 1979.
 * Adaptive models: I. H. Witten, R. M. Neal and J. G. Cleary, "Arithmetic
 * coding for data compression", Communications of the ACM 30(6), 1987:
 * a frequency count per symbol, incremented on each use, with all counts
 * halved when their total reaches a limit.
 *
 * Changes from the sources (see PROVENANCE.md): byte-wise renormalisation
 * with a 32-bit range; the carry is resolved by holding back bytes that a
 * carry could still change, the byte analogue of Witten, Neal & Cleary's
 * bits_to_follow; the increment and the halving limit are set per model
 * (Witten, Neal & Cleary used 1 and 2^14 - 1).
 */
#ifndef SNT_RANGECODER_H
#define SNT_RANGECODER_H

#include <stdint.h>

#define SNT_RC_MAX_SYMBOLS 256
#define SNT_RC_MAX_TOTAL   (1u << 16)  /* keeps range / total >= 2^8 */

/* Model parameters: starting count, increment per symbol, halving limit.
 * Coefficient-like alphabets (256 symbols) start nearly empty and learn
 * fast; mode flags start with a prior of one observation per symbol
 * (the starting count equals the increment, as in Witten, Neal & Cleary).
 * Chosen by BD-rate on the Xiph CIF tuning clips. */
#define SNT_RC_INIT_COEFF   1
#define SNT_RC_INC_COEFF    128
#define SNT_RC_LIMIT_COEFF  65000
#define SNT_RC_INIT_MODE    96
#define SNT_RC_INC_MODE     96
#define SNT_RC_LIMIT_MODE   8192

/* Adaptive frequency-count model */
typedef struct {
    uint32_t total;
    uint16_t num_symbols;
    uint16_t inc;
    uint32_t limit;
    uint16_t freq[SNT_RC_MAX_SYMBOLS];
} SntModel;

/* Every symbol starts with count init. */
void snt_model_init(SntModel *m, int num_symbols, int init, int inc, uint32_t limit);

/* Count one occurrence of sym; halve all counts (keeping them >= 1) once
 * the total exceeds the limit. */
void snt_model_update(SntModel *m, int sym);

/* Cumulative count of the symbols below sym */
static inline uint32_t snt_model_cum(const SntModel *m, int sym) {
    uint32_t c = 0;
    for (int i = 0; i < sym; i++) c += m->freq[i];
    return c;
}

/* ── Encoder ──────────────────────────────────────────────────────── */
typedef struct {
    uint8_t  *buf;
    uint32_t  cap, pos;
    uint64_t  low;      /* 33 bits: bit 32 is a pending carry */
    uint32_t  range;
    uint8_t   cache;    /* last byte not yet written (a carry may change it) */
    uint32_t  pending;  /* bytes held back: cache, then pending - 1 0xFF bytes */
} SntRangeEnc;

void snt_rc_enc_init(SntRangeEnc *e, uint8_t *buf, uint32_t cap);
void snt_rc_encode(SntRangeEnc *e, uint32_t cum, uint32_t freq, uint32_t total);
/* Encode sym with model m, then update m */
void snt_rc_encode_sym(SntRangeEnc *e, SntModel *m, int sym);
/* Write the final bytes; returns the number of bytes in buf */
uint32_t snt_rc_enc_finish(SntRangeEnc *e);

/* ── Decoder ──────────────────────────────────────────────────────── */
typedef struct {
    const uint8_t *buf;
    uint32_t len, pos;
    uint32_t code;      /* offset of the coded value from the range's low end */
    uint32_t range;
} SntRangeDec;

void snt_rc_dec_init(SntRangeDec *d, const uint8_t *buf, uint32_t len);
/* Decode one symbol with model m, then update m */
int snt_rc_decode_sym(SntRangeDec *d, SntModel *m);

/* Binary decision with an explicit probability p1 (of a 1) in 1/4096ths,
 * 1..4095: encode bit, or decode and return it */
void snt_rc_encode_bit(SntRangeEnc *e, int bit, uint32_t p1);
int snt_rc_decode_bit(SntRangeDec *d, uint32_t p1);

#endif /* SNT_RANGECODER_H */
