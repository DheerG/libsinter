/*
 * Sinter Codec — rANS Entropy Coder
 * See rans.c for implementation details.
 */
#ifndef SNT_RANS_H
#define SNT_RANS_H

#include <stdint.h>

/* ── Constants ────────────────────────────────────────────────────── */
#define RANS_PROB_BITS   15         /* CDF resolution: 15-bit */
#define RANS_PROB_TOTAL  (1 << RANS_PROB_BITS)  /* 32768 */
#define RANS_STATE_BITS  32
#define RANS_RENORM_BITS 16         /* renormalize when state < (1<<16) */
#define RANS_RENORM_LO   (1u << RANS_RENORM_BITS)  /* 65536 */

/* Probability adaptation rates (PARA) */
#define RANS_ADAPT_RATE_COEFF  5   /* for transform coefficients */
#define RANS_ADAPT_RATE_MODE   4   /* for mode decisions */

/* Maximum symbols per CDF */
#define RANS_MAX_SYMBOLS  256

/* ── rANS state ───────────────────────────────────────────────────── */
typedef uint32_t RansState;

/* ── CDF (Cumulative Distribution Function) ───────────────────────── */
typedef struct {
    uint16_t cdf[RANS_MAX_SYMBOLS + 1]; /* cdf[i] = cumulative prob of symbols 0..i-1 */
    uint16_t num_symbols;
    uint8_t  adapt_rate;                /* PARA adaptation rate */
} RansCDF;

/* ── Bitstream buffer for rANS ────────────────────────────────────── */
typedef struct {
    uint8_t  *buf;       /* buffer start */
    uint32_t  capacity;  /* buffer size in bytes */
    uint32_t  pos;       /* current read/write position */
} RansBitBuf;

/* ── Encoder (writes in reverse, LIFO) ────────────────────────────── */
typedef struct {
    RansState state0;    /* primary rANS state (for even-indexed symbols) */
    RansState state1;    /* secondary rANS state (for odd-indexed symbols) */
    int       dual;      /* 1 = dual interleaved, 0 = single state */
    /* Reverse buffer: encoder writes from end to start */
    uint16_t *rev_buf;   /* 16-bit output words, written in reverse */
    uint32_t  rev_cap;   /* capacity in 16-bit words */
    uint32_t  rev_pos;   /* next write position (counts down from rev_cap-1) */
} RansEncoder;

/* ── Decoder (reads forward, FIFO) ────────────────────────────────── */
typedef struct {
    RansState state0;
    RansState state1;
    int       dual;
    uint32_t  sym_count; /* number of symbols decoded so far */
    const uint8_t *buf;  /* input bitstream */
    uint32_t  buf_len;
    uint32_t  pos;       /* current read position */
} RansDecoder;

/* ── CDF operations ───────────────────────────────────────────────── */

/* Initialize CDF to uniform distribution for num_symbols symbols. */
void rans_cdf_init_uniform(RansCDF *cdf, uint16_t num_symbols, uint8_t adapt_rate);

/* Adapt CDF after observing a symbol (PARA exponential moving average).
 * cdf->cdf[i] += ((target - cdf[i]) + (1 << (rate-1))) >> rate */
void rans_cdf_adapt(RansCDF *cdf, uint16_t symbol);

/* Get probability range for a symbol: [cdf_lo, cdf_lo + freq) */
static inline void rans_cdf_symbol_range(const RansCDF *cdf, uint16_t symbol,
                                          uint16_t *cdf_lo, uint16_t *freq) {
    *cdf_lo = cdf->cdf[symbol];
    *freq   = cdf->cdf[symbol + 1] - cdf->cdf[symbol];
}

/* Find symbol from cumulative value (for decoding) */
uint16_t rans_cdf_lookup(const RansCDF *cdf, uint16_t cum_val);

/* ── Encoder operations ───────────────────────────────────────────── */

/* Initialize encoder. dual=1 for dual interleaved rANS. */
void rans_enc_init(RansEncoder *enc, uint16_t *rev_buf, uint32_t rev_cap, int dual);

/* Encode a symbol. Must be called in REVERSE order relative to decode.
 * fwd_idx: the forward-order index of this symbol (0-based).
 * In dual mode, even indices use state0, odd use state1. */
void rans_enc_put(RansEncoder *enc, uint16_t cdf_lo, uint16_t freq, uint32_t fwd_idx);

/* Flush final state to buffer. Call after all symbols are encoded. */
void rans_enc_flush(RansEncoder *enc);

/* Get encoded data: returns pointer and length of the encoded bitstream. */
const uint16_t *rans_enc_data(const RansEncoder *enc, uint32_t *len_words);

/* ── Decoder operations ───────────────────────────────────────────── */

/* Initialize decoder from bitstream. */
void rans_dec_init(RansDecoder *dec, const uint8_t *buf, uint32_t buf_len, int dual);

/* Get cumulative frequency from current state (for CDF lookup). */
uint16_t rans_dec_get(const RansDecoder *dec);

/* Advance decoder state after identifying the symbol. */
void rans_dec_advance(RansDecoder *dec, uint16_t cdf_lo, uint16_t freq);

/* Convenience: decode one symbol given a CDF.
 * IMPORTANT: This function does NOT adapt the CDF. Callers MUST call
 * rans_cdf_adapt(cdf, symbol) after each call, or CDF states will
 * diverge between encoder and decoder. */
uint16_t rans_dec_symbol(RansDecoder *dec, RansCDF *cdf);

#endif /* SNT_RANS_H */
