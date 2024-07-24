#ifndef SOLVE_PRESSURE_TASK_H
#define SOLVE_PRESSURE_TASK_H
#include "../../../../../Fortran/include/fds.h"
#include "../../../../data/parameters.h"
#include <hedgehog/hedgehog.h>

#define SolvePressureTaskInNb 1
#define SolvePressureTaskIn Parameters<ParameterIds::None>
#define SolvePressureTaskOut Parameters<ParameterIds::None>

class SolvePressureTask
    : public hh::AbstractTask<SolvePressureTaskInNb, SolvePressureTaskIn,
                              SolvePressureTaskOut> {
public:
  SolvePressureTask(size_t nbThreads)
      : hh::AbstractTask<SolvePressureTaskInNb, SolvePressureTaskIn,
                         SolvePressureTaskOut>("Solve Pressure Task",
                                               nbThreads) {}

  void
  execute(std::shared_ptr<Parameters<ParameterIds::None>> parameters) override {
    INFO("solvePressure");
    solvePressure();
    this->addResult(parameters);
  }

  std::shared_ptr<hh::AbstractTask<SolvePressureTaskInNb, SolvePressureTaskIn,
                                   SolvePressureTaskOut>>
  copy() override {
    return std::make_shared<SolvePressureTask>(this->numberThreads());
  }
};

#endif
