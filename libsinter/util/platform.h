/*
 * Sinter Codec — Platform Detection
 *
 * SIMD detection and endianness. All DSP functions have _ref (scalar)
 * versions first; SIMD variants are added later.
 */
#ifndef SNT_PLATFORM_H
#define SNT_PLATFORM_H

#include <stdint.h>

/* Endianness detection */
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
  #define SNT_BIG_ENDIAN 1
#else
  #define SNT_BIG_ENDIAN 0
#endif

/* Compiler hints */
#ifdef __GNUC__
  #define SNT_INLINE    static inline __attribute__((always_inline))
  #define SNT_UNLIKELY(x) __builtin_expect(!!(x), 0)
  #define SNT_LIKELY(x)   __builtin_expect(!!(x), 1)
#else
  #define SNT_INLINE    static inline
  #define SNT_UNLIKELY(x) (x)
  #define SNT_LIKELY(x)   (x)
#endif

#endif /* SNT_PLATFORM_H */
