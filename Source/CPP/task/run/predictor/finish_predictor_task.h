#ifndef FINISH_PREDICTOR_TASK_H
#define FINISH_PREDICTOR_TASK_H
#include "../../../../Fortran/include/fds.h"
#include "../../../data/parameters.h"
#include <hedgehog/hedgehog.h>

#define FinishPredictorTaskInNb 1
#define FinishPredictorTaskIn Parameters<ParameterIds::Finished>
#define FinishPredictorTaskOut Parameters<ParameterIds::None>

class FinishPredictorTask
    : public hh::AbstractTask<FinishPredictorTaskInNb, FinishPredictorTaskIn,
                              FinishPredictorTaskOut> {
public:
  FinishPredictorTask(size_t nbThreads)
      : hh::AbstractTask<FinishPredictorTaskInNb, FinishPredictorTaskIn,
                         FinishPredictorTaskOut>("Finish Predictor Task",
                                                 nbThreads) {}

  void execute(
      std::shared_ptr<Parameters<ParameterIds::Finished>> parameters) override {
    INFO_GRP("finish predictor", PREDICTOR_GRP);
    finishPredictor();
    this->addResult(
        std::make_shared<Parameters<ParameterIds::None>>(parameters));
  }

  std::shared_ptr<hh::AbstractTask<
      FinishPredictorTaskInNb, FinishPredictorTaskIn, FinishPredictorTaskOut>>
  copy() override {
    return std::make_shared<FinishPredictorTask>(this->numberThreads());
  }
};

#endif
