#ifndef FINISH_TASK_H
#define FINISH_TASK_H
#include "../../../Fortran/include/fds.h"
#include "../../data/parameters.h"
#include "hedgehog/src/api/task/abstract_task.h"
#include <hedgehog/hedgehog.h>

#define FinishTaskInNb 1
#define FinishTaskIn Parameters
#define FinishTaskOut Parameters

class FinishTask
    : public hh::AbstractTask<FinishTaskInNb, FinishTaskIn, FinishTaskOut> {
public:
  FinishTask(size_t nbThreads)
      : hh::AbstractTask<FinishTaskInNb, FinishTaskIn, FinishTaskOut>(
            "Finish Task", nbThreads) {}

  void execute(std::shared_ptr<Parameters> fileName) override {
    finish();
    this->addResult(fileName);
  }

  std::shared_ptr<hh::AbstractTask<FinishTaskInNb, FinishTaskIn, FinishTaskOut>>
  copy() override {
    return std::make_shared<FinishTask>(this->numberThreads());
  }
};

#endif
