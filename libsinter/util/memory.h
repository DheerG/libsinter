/*
 * Sinter Codec — Memory Utilities
 *
 * Aligned allocation for SIMD-friendly frame buffers.
 * 64-byte alignment for AVX-512 compatibility.
 */
#ifndef SNT_MEMORY_H
#define SNT_MEMORY_H

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define SNT_ALIGN 64 /* bytes — AVX-512 friendly */

/* Aligned allocation. Returns NULL on failure. */
static inline void *snt_alloc(size_t size) {
    void *ptr = NULL;
#if defined(_WIN32)
    ptr = _aligned_malloc(size, SNT_ALIGN);
#else
    if (posix_memalign(&ptr, SNT_ALIGN, size) != 0)
        return NULL;
#endif
    return ptr;
}

/* Aligned allocation, zero-initialized */
static inline void *snt_calloc(size_t count, size_t size) {
    size_t total = count * size;
    void *ptr = snt_alloc(total);
    if (ptr) memset(ptr, 0, total);
    return ptr;
}

static inline void snt_free(void *ptr) {
#if defined(_WIN32)
    _aligned_free(ptr);
#else
    free(ptr);
#endif
}

#endif /* SNT_MEMORY_H */
