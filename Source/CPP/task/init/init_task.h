#ifndef INIT_TASK_H
#define INIT_TASK_H
#include "../../../Fortran/include/fds.h"
#include "../../data/parameters.h"
#include "hedgehog/src/api/task/abstract_task.h"
#include <hedgehog/hedgehog.h>

#define InitTaskInNb 1
#define InitTaskIn Parameters<ParameterIds::None>
#define InitTaskOut Parameters<ParameterIds::None>

class InitTask
    : public hh::AbstractTask<InitTaskInNb, InitTaskIn, InitTaskOut> {
public:
  InitTask(size_t nbThreads)
      : hh::AbstractTask<InitTaskInNb, InitTaskIn, InitTaskOut>("Init Task",
                                                                nbThreads) {}

  void
  execute(std::shared_ptr<Parameters<ParameterIds::None>> parameters) override {
    INFO("start");
    start(parameters->fnInput.data());
    this->addResult(parameters);
  }

  std::shared_ptr<hh::AbstractTask<InitTaskInNb, InitTaskIn, InitTaskOut>>
  copy() override {
    return std::make_shared<InitTask>(this->numberThreads());
  }
};

#endif
