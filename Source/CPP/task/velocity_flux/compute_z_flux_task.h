#ifndef LOOP_FDS_COMPUTE_Z_FLUZ_H
#define LOOP_FDS_COMPUTE_Z_FLUZ_H
#include "hedgehog/hedgehog.h"
#include "../../data/block_idx.h"
#include "../../fortran/cxx/loops.h"

#define CZFTInNb 1
#define CZFTIn BlockIdx<Z>
#define CZFTOut BlockIdx<Z>

class ComputeZFluxTask : public hh::AbstractTask<CZFTInNb, CZFTIn, CZFTOut> {
 public:
  explicit ComputeZFluxTask(size_t nbThreads):
          hh::AbstractTask<CZFTInNb, CZFTIn, CZFTOut>("Compute Z Task", nbThreads) {}

  void execute(std::shared_ptr<BlockIdx<Z>> block) override {
    size_t k = block->k == 1 ? 0 : block->k;
    computeZFlux(&block->i, &block->j, &k,
                 &block->iEnd,
                 &block->jEnd,
                 &block->kEnd);
    this->addResult(block);
  }

  std::shared_ptr<hh::AbstractTask<CZFTInNb, CZFTIn, CZFTOut>>
  copy() override {
    return std::make_shared<ComputeZFluxTask>(this->numberThreads());
  }
};

#endif //LOOP_FDS_COMPUTE_Z_FLUZ_H
