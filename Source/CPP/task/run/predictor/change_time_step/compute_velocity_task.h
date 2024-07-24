#ifndef COMPUTE_VELOCITY_TASK_H
#define COMPUTE_VELOCITY_TASK_H
#include "../../../../../Fortran/include/fds.h"
#include "../../../../data/parameters.h"
#include <hedgehog/hedgehog.h>

#define ComputeVelocityTaskInNb 1
#define ComputeVelocityTaskIn Parameters<ParameterIds::None>
#define ComputeVelocityTaskOut Parameters<ParameterIds::None>

class ComputeVelocityTask
    : public hh::AbstractTask<ComputeVelocityTaskInNb, ComputeVelocityTaskIn,
                              ComputeVelocityTaskOut> {
public:
  ComputeVelocityTask(size_t nbThreads)
      : hh::AbstractTask<ComputeVelocityTaskInNb, ComputeVelocityTaskIn,
                         ComputeVelocityTaskOut>("Compute Velocity Task",
                                                 nbThreads) {}

  void
  execute(std::shared_ptr<Parameters<ParameterIds::None>> parameters) override {
    INFO("computeVelocity");
    computeVelocity();
    this->addResult(parameters);
  }

  std::shared_ptr<hh::AbstractTask<
      ComputeVelocityTaskInNb, ComputeVelocityTaskIn, ComputeVelocityTaskOut>>
  copy() override {
    return std::make_shared<ComputeVelocityTask>(this->numberThreads());
  }
};

#endif
