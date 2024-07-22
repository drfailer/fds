#ifndef RUN_TASK_H
#define RUN_TASK_H
#include "../../../Fortran/include/fds.h"
#include "../../data/parameters.h"
#include "hedgehog/src/api/task/abstract_task.h"
#include <hedgehog/hedgehog.h>

#define RunTaskInNb 1
#define RunTaskIn Parameters
#define RunTaskOut Parameters

class RunTask
    : public hh::AbstractTask<RunTaskInNb, RunTaskIn, RunTaskOut> {
public:
  RunTask(size_t nbThreads)
      : hh::AbstractTask<RunTaskInNb, RunTaskIn, RunTaskOut>("Run Task",
                                                                nbThreads) {}

  void execute(std::shared_ptr<Parameters> parameters) override {
    run();
    this->addResult(parameters);
  }

  std::shared_ptr<hh::AbstractTask<RunTaskInNb, RunTaskIn, RunTaskOut>>
  copy() override {
    return std::make_shared<RunTask>(this->numberThreads());
  }
};

#endif
