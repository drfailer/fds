#ifndef TIME_STEP_STATE_MANAGER_H
#define TIME_STEP_STATE_MANAGER_H
#include "time_step_state.h"
#include <hedgehog/hedgehog.h>

class TimeStepStateManager
    : public hh::StateManager<TSStateInNb, TSStateIn, TSStateOut> {
public:
  TimeStepStateManager(std::shared_ptr<TimeStepState> const &state)
      : hh::StateManager<TSStateInNb, TSStateIn, TSStateOut>(
            state, "Time State Manager") {}

  [[nodiscard]] bool canTerminate() const override {
    this->state()->lock();
    auto ret =
        std::dynamic_pointer_cast<TimeStepState>(this->state())->isDone();
    this->state()->unlock();
    return ret;
  }
};

#endif
