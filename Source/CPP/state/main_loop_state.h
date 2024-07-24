#ifndef MAIN_LOOP_STATE_H
#define MAIN_LOOP_STATE_H
#include "../../Fortran/include/fds.h"
#include "../data/parameters.h"
#include "../data/signal.h"
#include <hedgehog/hedgehog.h>

class TimeStepState;

#define MainLoopStateInNb 1
#define MainLoopStateIn Parameters<ParameterIds::None>
#define MainLoopStateOut                                                       \
  Parameters<ParameterIds::None>, bool, Signal<Sigs::Stop>

class MainLoopState
    : public hh::AbstractState<MainLoopStateInNb, MainLoopStateIn,
                               MainLoopStateOut> {
public:
  MainLoopState()
      : hh::AbstractState<MainLoopStateInNb, MainLoopStateIn,
                          MainLoopStateOut>() {}

  void
  execute(std::shared_ptr<Parameters<ParameterIds::None>> parameters) override {
    if (stopMainLoop()) {
      INFO("HH STOP");
      stop_ = true;
      this->addResult(std::make_shared<bool>(true)); // todo: bool is temporary
      // send stop signal to the TimeStepState
      this->addResult(std::make_shared<Signal<Sigs::Stop>>());
    } else {
      this->addResult(parameters);
    }
  }

  [[nodiscard]] bool isDone() { return stop_; }

private:
  bool stop_ = false;
};

#endif
