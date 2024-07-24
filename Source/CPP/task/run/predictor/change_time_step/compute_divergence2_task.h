#ifndef COMPUTE_DIVERGENCE2_TASK_H
#define COMPUTE_DIVERGENCE2_TASK_H
#include "../../../../../Fortran/include/fds.h"
#include "../../../../data/parameters.h"
#include <hedgehog/hedgehog.h>

#define ComputeDivergence2TaskInNb 1
#define ComputeDivergence2TaskIn Parameters<ParameterIds::None>
#define ComputeDivergence2TaskOut Parameters<ParameterIds::None>

class ComputeDivergence2Task
    : public hh::AbstractTask<ComputeDivergence2TaskInNb,
                              ComputeDivergence2TaskIn,
                              ComputeDivergence2TaskOut> {
public:
  ComputeDivergence2Task(size_t nbThreads)
      : hh::AbstractTask<ComputeDivergence2TaskInNb, ComputeDivergence2TaskIn,
                         ComputeDivergence2TaskOut>("Compute Divergence 2 Task",
                                                    nbThreads) {}

  void
  execute(std::shared_ptr<Parameters<ParameterIds::None>> parameters) override {
    INFO("computeDivergencePhase2");
    computeDivergencePhase2();
    this->addResult(parameters);
  }

  std::shared_ptr<
      hh::AbstractTask<ComputeDivergence2TaskInNb, ComputeDivergence2TaskIn,
                       ComputeDivergence2TaskOut>>
  copy() override {
    return std::make_shared<ComputeDivergence2Task>(this->numberThreads());
  }
};

#endif
