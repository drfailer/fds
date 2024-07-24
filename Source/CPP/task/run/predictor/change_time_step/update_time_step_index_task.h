#ifndef UPDATE_TIME_STEP_INDEX_TASK_H
#define UPDATE_TIME_STEP_INDEX_TASK_H
#include "../../../../../Fortran/include/fds.h"
#include "../../../../data/parameters.h"
#include <hedgehog/hedgehog.h>

#define UpdateTimeStepIndexTaskInNb 1
#define UpdateTimeStepIndexTaskIn Parameters<ParameterIds::UpdateTS>
#define UpdateTimeStepIndexTaskOut Parameters<ParameterIds::TimeStepReduced>

class UpdateTimeStepIndexTask
    : public hh::AbstractTask<UpdateTimeStepIndexTaskInNb,
                              UpdateTimeStepIndexTaskIn,
                              UpdateTimeStepIndexTaskOut> {
public:
  UpdateTimeStepIndexTask(size_t nbThreads)
      : hh::AbstractTask<UpdateTimeStepIndexTaskInNb, UpdateTimeStepIndexTaskIn,
                         UpdateTimeStepIndexTaskOut>(
            "Update Time Step Index Task", nbThreads) {}

  void execute(
      std::shared_ptr<Parameters<ParameterIds::UpdateTS>> parameters) override {
    INFO("updateTimeStepIndex");
    updateTimeStepIndex();
    this->addResult(std::make_shared<Parameters<ParameterIds::TimeStepReduced>>(
        parameters));
  }

  std::shared_ptr<
      hh::AbstractTask<UpdateTimeStepIndexTaskInNb, UpdateTimeStepIndexTaskIn,
                       UpdateTimeStepIndexTaskOut>>
  copy() override {
    return std::make_shared<UpdateTimeStepIndexTask>(this->numberThreads());
  }
};

#endif
