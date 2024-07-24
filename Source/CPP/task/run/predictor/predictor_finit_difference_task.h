#ifndef PREDICTOR_FINIT_DIFFERENCE_TASK_H
#define PREDICTOR_FINIT_DIFFERENCE_TASK_H
#include "../../../../Fortran/include/fds.h"
#include "../../../data/parameters.h"
#include <hedgehog/hedgehog.h>

#define PredictorFinitDifferenceTaskInNb 1
#define PredictorFinitDifferenceTaskIn Parameters<ParameterIds::None>
#define PredictorFinitDifferenceTaskOut Parameters<ParameterIds::None>

class PredictorFinitDifferenceTask
    : public hh::AbstractTask<PredictorFinitDifferenceTaskInNb,
                              PredictorFinitDifferenceTaskIn,
                              PredictorFinitDifferenceTaskOut> {
public:
  PredictorFinitDifferenceTask(size_t nbThreads)
      : hh::AbstractTask<PredictorFinitDifferenceTaskInNb,
                         PredictorFinitDifferenceTaskIn,
                         PredictorFinitDifferenceTaskOut>(
            "Predictor Finit Difference Task", nbThreads) {}

  void
  execute(std::shared_ptr<Parameters<ParameterIds::None>> parameters) override {
    predictorFinitDifference();
    this->addResult(parameters);
  }

  std::shared_ptr<hh::AbstractTask<PredictorFinitDifferenceTaskInNb,
                                   PredictorFinitDifferenceTaskIn,
                                   PredictorFinitDifferenceTaskOut>>
  copy() override {
    return std::make_shared<PredictorFinitDifferenceTask>(
        this->numberThreads());
  }
};

#endif
