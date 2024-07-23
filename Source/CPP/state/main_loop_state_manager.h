#ifndef MAIN_LOOP_STATE_MANAGER_H
#define MAIN_LOOP_STATE_MANAGER_H
#include "main_loop_state.h"
#include <hedgehog/hedgehog.h>

class MainLoopStateManager
    : public hh::StateManager<MainLoopStateInNb, MainLoopStateIn,
                              MainLoopStateOut> {
public:
  MainLoopStateManager(std::shared_ptr<MainLoopState> const &state)
      : hh::StateManager<MainLoopStateInNb, MainLoopStateIn, MainLoopStateOut>(
            state, "Main Loop State Manager") {}

  [[nodiscard]] bool canTerminate() const override {
    this->state()->lock();
    auto ret =
        std::dynamic_pointer_cast<MainLoopState>(this->state())->isDone();
    this->state()->unlock();
    return ret;
  }
};

#endif
