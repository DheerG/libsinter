/*
 * Sinter Codec — Quad-Tree Block Partitioning
 * Implements Section 1.1: Block Partitioning (quad-tree, 32x32 down to 4x4)
 *
 * The encoder decides per-superblock whether to split into 4 sub-blocks.
 * Split decisions are signaled as 1 bit per node in the quad-tree.
 * The decoder mirrors the split logic to reconstruct the partition.
 */
#ifndef SNT_QUADTREE_H
#define SNT_QUADTREE_H

#include <stdint.h>
#include <stdbool.h>

/* Supported block sizes: 32, 16, 8, 4 */
#define SNT_QT_MAX_SIZE  32
#define SNT_QT_MIN_SIZE  4
#define SNT_QT_MAX_DEPTH 3  /* 32->16->8->4 = 3 splits */

/* Quad-tree node: represents a block that may be split into 4 children */
typedef struct SntQTNode {
    uint16_t x, y;         /* top-left position in pixels */
    uint8_t  log2_size;    /* log2 of block width/height */
    bool     split;        /* true if split into 4 children */
    uint8_t  tx_mode;      /* 0=lapped, 1=block-DCT (per-block decision) */
} SntQTNode;

/* Callback for processing a leaf block (not split further).
 * Called by the quad-tree traversal with the block's position and size. */
typedef void (*SntQTLeafFunc)(int x, int y, int size, int tx_mode, void *ctx);

/* Traverse a quad-tree, calling leaf_func for each leaf block.
 * Reads split decisions from the bitstream (decoder) or signals them (encoder). */
void snt_qt_traverse(int x, int y, int size, int frame_w, int frame_h,
                      const uint8_t *split_flags, int *flag_idx,
                      SntQTLeafFunc leaf_func, void *ctx);

/* Encoder: decide quad-tree splits using RDO.
 * Returns the number of split flags written. */
int snte_qt_decide_splits(const int16_t *pre_filtered, int stride,
                           int frame_w, int frame_h,
                           uint8_t base_qp, float lambda,
                           uint8_t *split_flags, uint8_t *tx_modes);

#endif /* SNT_QUADTREE_H */
