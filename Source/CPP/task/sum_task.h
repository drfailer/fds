#ifndef LOOP_FDS_SUM_TASK_H
#define LOOP_FDS_SUM_TASK_H
#include <hedgehog/hedgehog.h>
#include <memory>
#include <utility>
#include "../data/block_idx.h"
#include "../../fortran/cxx/loops.h"

#define SumTInNb 3
#define SumTIn BlockIdx<X>, BlockIdx<Y>, BlockIdx<Z>
#define SumTOut                     \
  std::pair<BlockIdx<X>, double>, \
  std::pair<BlockIdx<Y>, double>, \
  std::pair<BlockIdx<Z>, double>

// note: the calculus may be done in fortran

class SumTask : public hh::AbstractTask<SumTInNb, SumTIn, SumTOut> {
 public:
  explicit SumTask(size_t nbThreads):
   hh::AbstractTask<SumTInNb, SumTIn, SumTOut>("Sum Task", nbThreads) {}

  void execute(std::shared_ptr<BlockIdx<X>> block) override {
    auto result = computePartialSum(block);
    this->addResult(std::make_shared<std::pair<BlockIdx<X>, double>>(*block, result));
  }

  void execute(std::shared_ptr<BlockIdx<Y>> block) override {
    auto result = computePartialSum(block);
    this->addResult(std::make_shared<std::pair<BlockIdx<Y>, double>>(*block, result));
  }

  void execute(std::shared_ptr<BlockIdx<Z>> block) override {
    auto result = computePartialSum(block);
    this->addResult(std::make_shared<std::pair<BlockIdx<Z>, double>>(*block, result));
  }

  std::shared_ptr<hh::AbstractTask<SumTInNb, SumTIn, SumTOut>>
  copy() override {
    return std::make_shared<SumTask>(this->numberThreads());
  }

  template <BlockIdxTypes BT>
    double computePartialSum(std::shared_ptr<BlockIdx<BT>> block) {
      double result = 0;
      size_t direction = BT;
      size_t i = block->i;
      size_t j = block->j;
      size_t k = block->k;
      size_t iEnd = block->iEnd;
      size_t jEnd = block->jEnd;
      size_t kEnd = block->kEnd;

      partialSum(&direction, &i, &j, &k, &iEnd, &jEnd, &kEnd, &result);
      return result;
    }

};

#endif //LOOP_FDS_SUM_TASK_H
