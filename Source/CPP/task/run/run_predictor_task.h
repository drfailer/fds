#ifndef RUN_PREDICTOR_TASK_H
#define RUN_PREDICTOR_TASK_H
#include "../../../Fortran/include/fds.h"
#include "../../data/parameters.h"
#include <hedgehog/hedgehog.h>

#define RunPredictorTaskInNb 1
#define RunPredictorTaskIn Parameters
#define RunPredictorTaskOut Parameters

class RunPredictorTask
    : public hh::AbstractTask<RunPredictorTaskInNb, RunPredictorTaskIn,
                       RunPredictorTaskOut> {
public:
  RunPredictorTask(size_t nbThreads)
      : hh::AbstractTask<RunPredictorTaskInNb, RunPredictorTaskIn,
                         RunPredictorTaskOut>("Run Predictor Taks", nbThreads) {
  }

  void execute(std::shared_ptr<Parameters> parameters) override {
    runPredictor();
    this->addResult(parameters);
  }

  std::shared_ptr<hh::AbstractTask<RunPredictorTaskInNb, RunPredictorTaskIn,
                                   RunPredictorTaskOut>>
  copy() override {
    return std::make_shared<RunPredictorTask>(this->numberThreads());
  }
};

#endif
