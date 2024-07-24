#ifndef COMPUTE_WALL_BC_TASK_H
#define COMPUTE_WALL_BC_TASK_H
#include "../../../../../Fortran/include/fds.h"
#include "../../../../data/parameters.h"
#include <hedgehog/hedgehog.h>

#define ComputeWallBCTaskInNb 1
#define ComputeWallBCTaskIn Parameters<ParameterIds::None>
#define ComputeWallBCTaskOut Parameters<ParameterIds::None>

class ComputeWallBCTask
    : public hh::AbstractTask<ComputeWallBCTaskInNb, ComputeWallBCTaskIn,
                              ComputeWallBCTaskOut> {
public:
  ComputeWallBCTask(size_t nbThreads)
      : hh::AbstractTask<ComputeWallBCTaskInNb, ComputeWallBCTaskIn,
                         ComputeWallBCTaskOut>("Compute Wall BC Task",
                                               nbThreads) {}

  void
  execute(std::shared_ptr<Parameters<ParameterIds::None>> parameters) override {
    INFO("computeWallBC");
    computeWallBC();
    this->addResult(parameters);
  }

  std::shared_ptr<hh::AbstractTask<ComputeWallBCTaskInNb, ComputeWallBCTaskIn,
                                   ComputeWallBCTaskOut>>
  copy() override {
    return std::make_shared<ComputeWallBCTask>(this->numberThreads());
  }
};

#endif
