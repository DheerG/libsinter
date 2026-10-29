/*
 * Sinter Codec — Range coder with adaptive frequency-count models (normative)
 * Sources and changes: see rangecoder.h.
 *
 * The coder keeps an interval [low, low + range) inside [0, 2^32), scaled
 * by 2^8 each time a byte leaves the top. Coding a symbol with cumulative
 * count cum and count freq out of total narrows it to
 *   low += (range / total) * cum,  range = (range / total) * freq,
 * then bytes are shifted out while range < 2^24.
 */

#include <string.h>
#include "rangecoder.h"

#define RC_TOP (1u << 24)

/* ── Models ───────────────────────────────────────────────────────── */

void snt_model_init(SntModel *m, int num_symbols, int init, int inc, uint32_t limit) {
    m->num_symbols = (uint16_t)num_symbols;
    m->inc = (uint16_t)inc;
    m->limit = limit;
    for (int i = 0; i < num_symbols; i++) m->freq[i] = (uint16_t)init;
    m->total = (uint32_t)(num_symbols * init);
}

void snt_model_update(SntModel *m, int sym) {
    m->freq[sym] += m->inc;
    m->total += m->inc;
    if (m->total > m->limit) {
        m->total = 0;
        for (int i = 0; i < m->num_symbols; i++) {
            m->freq[i] = (uint16_t)((m->freq[i] + 1) >> 1);
            m->total += m->freq[i];
        }
    }
}

/* ── Encoder ──────────────────────────────────────────────────────── */

void snt_rc_enc_init(SntRangeEnc *e, uint8_t *buf, uint32_t cap) {
    e->buf = buf;
    e->cap = cap;
    e->pos = 0;
    e->low = 0;
    e->range = 0xFFFFFFFFu;
    e->cache = 0;
    e->pending = 0;
}

static void put_byte(SntRangeEnc *e, uint8_t b) {
    if (e->pos < e->cap) e->buf[e->pos] = b;
    e->pos++;
}

/* Move the top byte of low out. A byte other than 0xFF, or a carry,
 * settles every held byte; a 0xFF could still turn into 0x00 with a
 * carry, so it is held back too. */
static void shift_low(SntRangeEnc *e) {
    uint32_t carry = (uint32_t)(e->low >> 32);
    uint8_t top = (uint8_t)(e->low >> 24);
    if (top != 0xFF || carry) {
        if (e->pending) {
            put_byte(e, (uint8_t)(e->cache + carry));
            for (uint32_t k = 1; k < e->pending; k++)
                put_byte(e, (uint8_t)(0xFF + carry));
        }
        e->cache = top;
        e->pending = 1;
    } else if (e->pending) {
        e->pending++;
    } else {
        e->cache = top;
        e->pending = 1;
    }
    e->low = (e->low & 0x00FFFFFFu) << 8;
}

void snt_rc_encode(SntRangeEnc *e, uint32_t cum, uint32_t freq, uint32_t total) {
    uint32_t r = e->range / total;
    e->low += (uint64_t)r * cum;
    e->range = r * freq;
    while (e->range < RC_TOP) {
        e->range <<= 8;
        shift_low(e);
    }
}

void snt_rc_encode_sym(SntRangeEnc *e, SntModel *m, int sym) {
    snt_rc_encode(e, snt_model_cum(m, sym), m->freq[sym], m->total);
    snt_model_update(m, sym);
}

uint32_t snt_rc_enc_finish(SntRangeEnc *e) {
    /* Any value in [low, low + range) decodes correctly, and the decoder
     * reads zeros past the end: take the one ending in three zero bytes
     * (range >= 2^24, so one exists), write its top byte, then settle the
     * held bytes. Trailing zeros are implicit, so drop them. */
    e->low = (e->low + 0x00FFFFFFu) & ~(uint64_t)0x00FFFFFFu;
    shift_low(e);
    shift_low(e);
    uint32_t n = e->pos < e->cap ? e->pos : e->cap;
    while (n > 0 && e->buf[n - 1] == 0) n--;
    return n;
}

/* ── Decoder ──────────────────────────────────────────────────────── */

static uint8_t get_byte(SntRangeDec *d) {
    return d->pos < d->len ? d->buf[d->pos++] : 0;
}

void snt_rc_dec_init(SntRangeDec *d, const uint8_t *buf, uint32_t len) {
    d->buf = buf;
    d->len = len;
    d->pos = 0;
    d->range = 0xFFFFFFFFu;
    d->code = 0;
    for (int i = 0; i < 4; i++)
        d->code = (d->code << 8) | get_byte(d);
}

int snt_rc_decode_sym(SntRangeDec *d, SntModel *m) {
    uint32_t r = d->range / m->total;
    uint32_t v = d->code / r;
    if (v >= m->total) v = m->total - 1;  /* only in a damaged stream */
    int sym = 0;
    uint32_t cum = 0;
    while (sym < m->num_symbols - 1 && v >= cum + m->freq[sym])
        cum += m->freq[sym++];
    d->code -= r * cum;
    d->range = r * m->freq[sym];
    while (d->range < RC_TOP) {
        d->code = (d->code << 8) | get_byte(d);
        d->range <<= 8;
    }
    snt_model_update(m, sym);
    return sym;
}

void snt_rc_encode_bit(SntRangeEnc *e, int bit, uint32_t p1) {
    if (bit) snt_rc_encode(e, 4096 - p1, p1, 4096);
    else     snt_rc_encode(e, 0, 4096 - p1, 4096);
}

int snt_rc_decode_bit(SntRangeDec *d, uint32_t p1) {
    uint32_t r = d->range / 4096;
    uint32_t v = d->code / r;
    int bit = v >= 4096 - p1;
    uint32_t cum = bit ? 4096 - p1 : 0, freq = bit ? p1 : 4096 - p1;
    d->code -= r * cum;
    d->range = r * freq;
    while (d->range < RC_TOP) {
        d->code = (d->code << 8) | get_byte(d);
        d->range <<= 8;
    }
    return bit;
}
