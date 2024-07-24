#ifndef TIME_STEP_STATE_H
#define TIME_STEP_STATE_H
#include "../../../Fortran/include/fds.h"
#include "../../data/parameters.h"
#include <hedgehog/hedgehog.h>
#include <iostream>

#define TSStateInNb 2
#define TSStateIn                                                              \
  Parameters<ParameterIds::InstabilityCheck>,                                  \
      Parameters<ParameterIds::TimeStepReduced>
#define TSStateOut                                                             \
  Parameters<ParameterIds::None>, Parameters<ParameterIds::UpdateTS>,          \
      Parameters<ParameterIds::Finished>

class TimeStepState
    : public hh::AbstractState<TSStateInNb, TSStateIn, TSStateOut> {
public:
  TimeStepState() : hh::AbstractState<TSStateInNb, TSStateIn, TSStateOut>() {}

  void execute(std::shared_ptr<Parameters<ParameterIds::InstabilityCheck>>
                   parameters) override {
    if (verifyInstability() == 1) {
      ERROR("instability error");
      this->addResult(std::make_shared<Parameters<ParameterIds::Finished>>(
          std::move(parameters)));
    } else {
      this->addResult(std::make_shared<Parameters<ParameterIds::UpdateTS>>(
          std::move(parameters)));
    }
  }

  void execute(std::shared_ptr<Parameters<ParameterIds::TimeStepReduced>>
                   parameters) override {
    if (timeStepReduced() == 0) {
      INFO_GRP("timeStepReduced", TIME_STEP_LOOP_INFO_GRP);
      this->addResult(std::make_shared<Parameters<ParameterIds::Finished>>(
          std::move(parameters)));
    } else {
      INFO_GRP("continue time step loop", TIME_STEP_LOOP_INFO_GRP);
      this->addResult(std::make_shared<Parameters<ParameterIds::None>>(
          std::move(parameters)));
    }
  }

  [[nodiscard]] bool isDone() { return !loop_; }

private:
  bool loop_ = true;
};

#endif
