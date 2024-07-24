#ifndef HVAC_SOLVER_TASK_H
#define HVAC_SOLVER_TASK_H
#include "../../../../../Fortran/include/fds.h"
#include "../../../../data/parameters.h"
#include <hedgehog/hedgehog.h>

#define HVACSolverTaskInNb 1
#define HVACSolverTaskIn Parameters<ParameterIds::None>
#define HVACSolverTaskOut Parameters<ParameterIds::None>

class HVACSolverTask
    : public hh::AbstractTask<HVACSolverTaskInNb, HVACSolverTaskIn,
                              HVACSolverTaskOut> {
public:
  HVACSolverTask(size_t nbThreads)
      : hh::AbstractTask<HVACSolverTaskInNb, HVACSolverTaskIn,
                         HVACSolverTaskOut>("HVAC Solver Task", nbThreads) {}

  void
  execute(std::shared_ptr<Parameters<ParameterIds::None>> parameters) override {
    INFO_GRP("hvacSolver", TIME_STEP_LOOP_INFO_GRP);
    hvacSolver();
    this->addResult(parameters);
  }

  std::shared_ptr<
      hh::AbstractTask<HVACSolverTaskInNb, HVACSolverTaskIn, HVACSolverTaskOut>>
  copy() override {
    return std::make_shared<HVACSolverTask>(this->numberThreads());
  }
};

#endif
