#ifndef LOOP_FDS_COMPUTE_X_FLUX_TASK_H
#define LOOP_FDS_COMPUTE_X_FLUX_TASK_H
#include "hedgehog/hedgehog.h"
#include "../../data/block_idx.h"
#include "../../fortran/cxx/loops.h"

#define CXFTInNb 1
#define CXFTIn BlockIdx<X>
#define CXFTOut BlockIdx<X>

class ComputeXFluxTask : public hh::AbstractTask<CXFTInNb, CXFTIn, CXFTOut> {
 public:
  explicit ComputeXFluxTask(size_t nbThreads):
   hh::AbstractTask<CXFTInNb, CXFTIn, CXFTOut>("Compute X Flux Task", nbThreads) {}

  void execute(std::shared_ptr<BlockIdx<X>> block) override {
    size_t i = block->i == 1 ? 0 : block->i;
    computeXFlux(&i, &block->j, &block->k,
                 &block->iEnd,
                 &block->jEnd,
                 &block->kEnd);
    this->addResult(block);
  }

  std::shared_ptr<hh::AbstractTask<CXFTInNb, CXFTIn, CXFTOut>>
  copy() override {
    return std::make_shared<ComputeXFluxTask>(this->numberThreads());
  }
};

#endif //LOOP_FDS_COMPUTE_X_FLUX_TASK_H
