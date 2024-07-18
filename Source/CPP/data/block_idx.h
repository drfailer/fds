#ifndef LOOP_FDS_BLOCK_IDX_H
#define LOOP_FDS_BLOCK_IDX_H

#include "block_idx_types.h"
#include <cstddef>
#include <memory>

template <BlockIdxTypes BlockType> struct BlockIdx {
  BlockIdx(size_t i, size_t j, size_t k, size_t iEnd, size_t jEnd, size_t kEnd,
           size_t nbBlocks, size_t ibar, size_t jbar, size_t kbar)
      : i(i), j(j), k(k), iEnd(iEnd), jEnd(jEnd), kEnd(kEnd),
        nbBlocks(nbBlocks), ibar(ibar), jbar(jbar), kbar(kbar) {}

  template <BlockIdxTypes Other>
  BlockIdx(BlockIdx<Other> const &other)
      : BlockIdx(other.i, other.j, other.k, other.iEnd, other.jEnd, other.kEnd,
                 other.nbBlocks, other.ibar, other.jbar, other.kbar) {}

  template <BlockIdxTypes Other>
  BlockIdx(std::shared_ptr<BlockIdx<Other>> const &other)
      : BlockIdx(other.i, other.j, other.k, other.iEnd, other.jEnd, other.kEnd,
                 other.nbBlocks, other.ibar, other.jbar, other.kbar) {}

  template <BlockIdxTypes Other>
  BlockIdx(std::shared_ptr<BlockIdx<Other>> &&other) noexcept
      : BlockIdx(other->i, other->j, other->k, other->iEnd, other->jEnd,
                 other->kEnd, other->nbBlocks, other->ibar, other->jbar,
                 other->kbar) {}

  size_t i = 0;
  size_t j = 0;
  size_t k = 0;
  size_t iEnd = 0;
  size_t jEnd = 0;
  size_t kEnd = 0;
  size_t nbBlocks = 0;
  size_t ibar = 0;
  size_t jbar = 0;
  size_t kbar = 0;
};

#endif // LOOP_FDS_BLOCK_IDX_H
