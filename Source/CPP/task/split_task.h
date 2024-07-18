#ifndef LOOP_FDS_SPLIT_TASK_H
#define LOOP_FDS_SPLIT_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/block_idx.h"
#include "../data/parameters.h"

#define STInNb 1
#define STIn Parameters
#define STOut BlockIdx<X>, BlockIdx<Y>, BlockIdx<Z>

class SplitTask : public hh::AbstractTask<STInNb, STIn, STOut> {
 public:
  SplitTask() :
          hh::AbstractTask<STInNb, STIn, STOut>("Split Task", 1) {}

  void execute(std::shared_ptr<Parameters> params) override {
    size_t ibar = params->ibp1 - 1;
    size_t jbar = params->jbp1 - 1;
    size_t kbar = params->kbp1 - 1;

    // we add 1 because the loops in fortran have an end condition equivalent to <= and not <
    for (size_t kBlock = 0; kBlock < params->nbBlocksK; ++kBlock) {
      for (size_t jBlock = 0; jBlock < params->nbBlocksJ; ++jBlock) {
        for (size_t iBlock = 0; iBlock < params->nbBlocksI; ++iBlock) {
          size_t i = iBlock * params->iBlockSize;
          size_t j = jBlock * params->jBlockSize;
          size_t k = kBlock * params->kBlockSize;
          size_t iEnd = std::min(ibar, i + params->iBlockSize - 1);
          size_t jEnd = std::min(jbar, j + params->jBlockSize - 1);
          size_t kEnd = std::min(kbar, k + params->kBlockSize - 1);
          size_t nbBlocks = params->nbBlocksI * params->nbBlocksJ * params->nbBlocksK;

          // we start at 1 by default, so we have less measures to make in the
          // rest of the code
          i = (i == 0) ? 1 : i;
          j = (j == 0) ? 1 : j;
          k = (k == 0) ? 1 : k;

          this->addResult(std::make_shared<BlockIdx<X>>(i, j, k, iEnd, jEnd, kEnd, nbBlocks, ibar, jbar, kbar));
          this->addResult(std::make_shared<BlockIdx<Y>>(i, j, k, iEnd, jEnd, kEnd, nbBlocks, ibar, jbar, kbar));
          this->addResult(std::make_shared<BlockIdx<Z>>(i, j, k, iEnd, jEnd, kEnd, nbBlocks, ibar, jbar, kbar));
        }
      }
    }
  }

  std::shared_ptr<hh::AbstractTask<STInNb, STIn, STOut>>
  copy() override {
    return std::make_shared<SplitTask>();
  }
};


#endif //LOOP_FDS_SPLIT_TASK_H
