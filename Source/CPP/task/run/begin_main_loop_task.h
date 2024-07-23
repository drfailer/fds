#ifndef BEGIN_MAIN_LOOP_TASK_H
#define BEGIN_MAIN_LOOP_TASK_H
#include "../../../Fortran/include/fds.h"
#include "../../data/parameters.h"
#include <hedgehog/hedgehog.h>

#define BeginMainLoopTaskInNb 1
#define BeginMainLoopTaskIn Parameters
#define BeginMainLoopTaskOut Parameters

class BeginMainLoopTask
    : public hh::AbstractTask<BeginMainLoopTaskInNb, BeginMainLoopTaskIn,
                              BeginMainLoopTaskOut> {
public:
  BeginMainLoopTask(size_t nbThreads)
      : hh::AbstractTask<BeginMainLoopTaskInNb, BeginMainLoopTaskIn,
                         BeginMainLoopTaskOut>("Begin Main Loop Task",
                                               nbThreads) {}

  void execute(std::shared_ptr<Parameters> parameters) override {
    mainLoopBegin();
    this->addResult(parameters);
  }

  std::shared_ptr<hh::AbstractTask<BeginMainLoopTaskInNb, BeginMainLoopTaskIn,
                                   BeginMainLoopTaskOut>>
  copy() override {
    return std::make_shared<BeginMainLoopTask>(this->numberThreads());
  }
};

#endif
