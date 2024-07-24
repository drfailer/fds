#ifndef START_PREDICTOR_TASK_H
#define START_PREDICTOR_TASK_H
#include "../../../../Fortran/include/fds.h"
#include "../../../data/parameters.h"
#include <hedgehog/hedgehog.h>

#define StartPredictorTaskInNb 1
#define StartPredictorTaskIn Parameters<ParameterIds::Start>
#define StartPredictorTaskOut Parameters<ParameterIds::None>

class StartPredictorTask
    : public hh::AbstractTask<StartPredictorTaskInNb, StartPredictorTaskIn,
                              StartPredictorTaskOut> {
public:
  StartPredictorTask(size_t nbThreads)
      : hh::AbstractTask<StartPredictorTaskInNb, StartPredictorTaskIn,
                         StartPredictorTaskOut>("Start Predictor Task",
                                                nbThreads) {}

  void execute(
      std::shared_ptr<Parameters<ParameterIds::Start>> parameters) override {
    INFO_GRP("start predictor", PREDICTOR_GRP);
    startPredictor();
    this->addResult(
        std::make_shared<Parameters<ParameterIds::None>>(parameters));
  }

  std::shared_ptr<hh::AbstractTask<StartPredictorTaskInNb, StartPredictorTaskIn,
                                   StartPredictorTaskOut>>
  copy() override {
    return std::make_shared<StartPredictorTask>(this->numberThreads());
  }
};

#endif
