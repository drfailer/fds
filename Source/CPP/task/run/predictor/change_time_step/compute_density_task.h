#ifndef COMPUTE_DENSITY_TASK_H
#define COMPUTE_DENSITY_TASK_H
#include "../../../../../Fortran/include/fds.h"
#include "../../../../data/parameters.h"
#include <hedgehog/hedgehog.h>

#define ComputeDensityTaskInNb 1
#define ComputeDensityTaskIn Parameters<ParameterIds::None>
#define ComputeDensityTaskOut Parameters<ParameterIds::None>

class ComputeDensityTask
    : public hh::AbstractTask<ComputeDensityTaskInNb, ComputeDensityTaskIn,
                              ComputeDensityTaskOut> {
public:
  ComputeDensityTask(size_t nbThreads)
      : hh::AbstractTask<ComputeDensityTaskInNb, ComputeDensityTaskIn,
                         ComputeDensityTaskOut>("Compute Density Task",
                                                nbThreads) {}

  void
  execute(std::shared_ptr<Parameters<ParameterIds::None>> parameters) override {
    INFO("computeDensity");
    computeDensity();
    this->addResult(parameters);
  }

  std::shared_ptr<hh::AbstractTask<ComputeDensityTaskInNb, ComputeDensityTaskIn,
                                   ComputeDensityTaskOut>>
  copy() override {
    return std::make_shared<ComputeDensityTask>(this->numberThreads());
  }
};

#endif
