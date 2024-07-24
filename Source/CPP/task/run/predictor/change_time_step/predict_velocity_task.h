#ifndef PREDICT_VELOCITY_TASK_H
#define PREDICT_VELOCITY_TASK_H
#include "../../../../../Fortran/include/fds.h"
#include "../../../../data/parameters.h"
#include <hedgehog/hedgehog.h>

#define PredictVelocityTaskInNb 1
#define PredictVelocityTaskIn Parameters<ParameterIds::None>
#define PredictVelocityTaskOut Parameters<ParameterIds::InstabilityCheck>

class PredictVelocityTask
    : public hh::AbstractTask<PredictVelocityTaskInNb, PredictVelocityTaskIn,
                              PredictVelocityTaskOut> {
public:
  PredictVelocityTask(size_t nbThreads)
      : hh::AbstractTask<PredictVelocityTaskInNb, PredictVelocityTaskIn,
                         PredictVelocityTaskOut>("Predict Velocity Task",
                                                 nbThreads) {}

  void
  execute(std::shared_ptr<Parameters<ParameterIds::None>> parameters) override {
    INFO_GRP("predictVelocity", TIME_STEP_LOOP_INFO_GRP);
    predictVelocity();
    this->addResult(
        std::make_shared<Parameters<ParameterIds::InstabilityCheck>>(
            parameters));
  }

  std::shared_ptr<hh::AbstractTask<
      PredictVelocityTaskInNb, PredictVelocityTaskIn, PredictVelocityTaskOut>>
  copy() override {
    return std::make_shared<PredictVelocityTask>(this->numberThreads());
  }
};

#endif
