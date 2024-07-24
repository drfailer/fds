#ifndef RUN_CORRECTOR_TASK_H
#define RUN_CORRECTOR_TASK_H
#include "../../../Fortran/include/fds.h"
#include "../../data/parameters.h"
#include <hedgehog/hedgehog.h>

#define RunCorrectorTaskInNb 1
#define RunCorrectorTaskIn Parameters<ParameterIds::None>
#define RunCorrectorTaskOut Parameters<ParameterIds::None>

class RunCorrectorTask
    : public hh::AbstractTask<RunCorrectorTaskInNb, RunCorrectorTaskIn,
                              RunCorrectorTaskOut> {
public:
  RunCorrectorTask(size_t nbThreads)
      : hh::AbstractTask<RunCorrectorTaskInNb, RunCorrectorTaskIn,
                         RunCorrectorTaskOut>("Run Corrector Taks", nbThreads) {
  }

  void
  execute(std::shared_ptr<Parameters<ParameterIds::None>> parameters) override {
    INFO("run corrector");
    runCorrector();
    this->addResult(parameters);
  }

  std::shared_ptr<hh::AbstractTask<RunCorrectorTaskInNb, RunCorrectorTaskIn,
                                   RunCorrectorTaskOut>>
  copy() override {
    return std::make_shared<RunCorrectorTask>(this->numberThreads());
  }
};

#endif
