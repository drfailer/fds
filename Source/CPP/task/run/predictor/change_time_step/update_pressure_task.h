#ifndef UPDATE_PRESSURE_TASK_H
#define UPDATE_PRESSURE_TASK_H
#include "../../../../../Fortran/include/fds.h"
#include "../../../../data/parameters.h"
#include <hedgehog/hedgehog.h>

#define UpdatePressureTaskInNb 1
#define UpdatePressureTaskIn Parameters<ParameterIds::None>
#define UpdatePressureTaskOut Parameters<ParameterIds::None>

class UpdatePressureTask
    : public hh::AbstractTask<UpdatePressureTaskInNb, UpdatePressureTaskIn,
                              UpdatePressureTaskOut> {
public:
  UpdatePressureTask(size_t nbThreads)
      : hh::AbstractTask<UpdatePressureTaskInNb, UpdatePressureTaskIn,
                         UpdatePressureTaskOut>("Update Pressure Task",
                                                nbThreads) {}

  void
  execute(std::shared_ptr<Parameters<ParameterIds::None>> parameters) override {
    INFO("updateGlobalPressure");
    updateGlobalPressure();
    this->addResult(parameters);
  }

  std::shared_ptr<hh::AbstractTask<UpdatePressureTaskInNb, UpdatePressureTaskIn,
                                   UpdatePressureTaskOut>>
  copy() override {
    return std::make_shared<UpdatePressureTask>(this->numberThreads());
  }
};

#endif
