#ifndef SIGNAL_FILTER_H
#define SIGNAL_FILTER_H
#include <hedgehog/hedgehog.h>
#include "../data/signal.h"

#define SFTaskInNb 1
#define SFTaskIn Signal<Sig>
#define SFTaskOut SFTaskIn

template <Sigs Sig>
class SignalFilter : public hh::AbstractTask<SFTaskInNb, SFTaskIn, SFTaskOut> {
  public:
  SignalFilter():
    hh::AbstractTask<SFTaskInNb, SFTaskIn, SFTaskOut>("Signal Filter", 1) {}

  void execute(std::shared_ptr<Signal<Sigs::Stop>> signal) override {
    this->addResult(signal);
  }

  std::shared_ptr<hh::AbstractTask<SFTaskInNb, SFTaskIn, SFTaskOut>>
  copy() override {
    return std::make_shared<SignalFilter<Sig>>();
  }
};

template <Sigs ...Sig>
class SignalsFilter : public SignalFilter<Sig>... {};

#endif
