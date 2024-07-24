#ifndef TIME_STEP_STATE_H
#define TIME_STEP_STATE_H
#include "../../../Fortran/include/fds.h"
#include "../../data/parameters.h"
#include "../../data/signal.h"
#include <hedgehog/hedgehog.h>

#define TSStateInNb 3
#define TSStateIn                                                              \
  Parameters<ParameterIds::InstabilityCheck>,                                  \
      Parameters<ParameterIds::TimeStepReduced>, Signal<Sigs::Stop>
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

  void execute(std::shared_ptr<Signal<Sigs::Stop>>) override {
    isDone_ = true;
  }

  [[nodiscard]] bool isDone() { return isDone_; }

private:
  bool isDone_ = false;
};

#endif
