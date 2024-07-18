#ifndef LOOP_FDS_RESULT_STATE_H
#define LOOP_FDS_RESULT_STATE_H
#include <hedgehog/hedgehog.h>
#include <utility>
#include <tuple>
#include "../data/block_idx.h"

#define RSInNb 3
#define RSIn                      \
  std::pair<BlockIdx<X>, double>, \
  std::pair<BlockIdx<Y>, double>, \
  std::pair<BlockIdx<Z>, double>
#define RSOut std::tuple<double, double, double>

class ResultState : public hh::AbstractState<RSInNb, RSIn, RSOut > {
 public:
  explicit ResultState(size_t ibar, size_t jbar, size_t kbar)
    : hh::AbstractState<RSInNb, RSIn, RSOut >(),
    total_(ibar * jbar * kbar),
    results_(0, 0, 0)
  {}

  void execute(std::shared_ptr<std::pair<BlockIdx<X>, double>> partialSum) override {
    tryInit(partialSum->first);
    std::get<X>(results_) += partialSum->second;
    --nbBlocksX_;
    trySendResult();
  }

  void execute(std::shared_ptr<std::pair<BlockIdx<Y>, double>> partialSum) override {
    tryInit(partialSum->first);
    std::get<Y>(results_) += partialSum->second;
    --nbBlocksY_;
    trySendResult();
  }

  void execute(std::shared_ptr<std::pair<BlockIdx<Z>, double>> partialSum) override {
    tryInit(partialSum->first);
    std::get<Z>(results_) += partialSum->second;
    --nbBlocksZ_;
    trySendResult();
  }

  template <BlockIdxTypes BT>
    void tryInit(BlockIdx<BT> const &block) {
      if (!nbBlockSet_) {
        nbBlocksX_ = block.nbBlocks;
        nbBlocksY_ = block.nbBlocks;
        nbBlocksZ_ = block.nbBlocks;
        nbBlockSet_ = true;
      }
    }

  void trySendResult() {
    if (isDone()) {
      this->addResult(computeResult());
    }
  }

  [[nodiscard]] bool isDone() const {
    return nbBlockSet_ && !nbBlocksX_ && !nbBlocksY_ && !nbBlocksZ_;
  }

  [[nodiscard]] std::shared_ptr<std::tuple<double, double, double>>
  computeResult() const {
    return std::make_shared<std::tuple<double, double, double>>(
        std::get<X>(results_) / total_,
        std::get<Y>(results_) / total_,
        std::get<Z>(results_) / total_
        );
  }

 private:
  size_t total_ = 0;
  std::tuple<double, double, double> results_;
  size_t nbBlocksX_ = 0;
  size_t nbBlocksY_ = 0;
  size_t nbBlocksZ_ = 0;
  bool nbBlockSet_ = false;
};

#endif //LOOP_FDS_RESULT_STATE_H
