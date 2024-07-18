#ifndef LOOP_FDS_COMPUTE_Y_FLUY_H
#define LOOP_FDS_COMPUTE_Y_FLUY_H
#include "hedgehog/hedgehog.h"
#include "../../data/block_idx.h"
#include "../../fortran/cxx/loops.h"

#define CYFTInNb 1
#define CYFTIn BlockIdx<Y>
#define CYFTOut BlockIdx<Y>

class ComputeYFluxTask : public hh::AbstractTask<CYFTInNb, CYFTIn, CYFTOut> {
 public:
  explicit ComputeYFluxTask(size_t nbThreads):
          hh::AbstractTask<CYFTInNb, CYFTIn, CYFTOut>("Compute Y Task", nbThreads) {}

  void execute(std::shared_ptr<BlockIdx<Y>> block) override {
    size_t j = block->j == 1 ? 0 : block->j;
    computeYFlux(&block->i, &j, &block->k,
                 &block->iEnd,
                 &block->jEnd,
                 &block->kEnd);
    this->addResult(block);
  }

  std::shared_ptr<hh::AbstractTask<CYFTInNb, CYFTIn, CYFTOut>>
  copy() override {
    return std::make_shared<ComputeYFluxTask>(this->numberThreads());
  }
};

#endif //LOOP_FDS_COMPUTE_Y_FLUY_H
