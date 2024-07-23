#ifndef MAIN_LOOP_STATE_H
#define MAIN_LOOP_STATE_H
#include "../../Fortran/include/fds.h"
#include "../data/parameters.h"
#include <hedgehog/hedgehog.h>

#define MainLoopStateInNb 1
#define MainLoopStateIn Parameters
#define MainLoopStateOut Parameters, bool

class MainLoopState
    : public hh::AbstractState<MainLoopStateInNb, MainLoopStateIn,
                               MainLoopStateOut> {
public:
  MainLoopState()
      : hh::AbstractState<MainLoopStateInNb, MainLoopStateIn,
                          MainLoopStateOut>() {}

  void execute(std::shared_ptr<Parameters> parameters) override {
    if (stopMainLoop()) {
      std::cout << "HH STOP" << std::endl;
      stop_ = true;
      this->addResult(std::make_shared<bool>(true)); // todo: bool is temporary
    } else {
      this->addResult(parameters);
    }
  }

  [[nodiscard]] bool isDone() { return stop_; }

private:
  bool stop_ = false;
};

#endif
