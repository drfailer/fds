#ifndef RUN_GRAPH_H
#define RUN_GRAPH_H
#include "../task/run/run_task.h"
#include <hedgehog/hedgehog.h>

#define RunGraphInNb 1
#define RunGraphIn Parameters
#define RunGraphOut Parameters

class RunGraph : public hh::Graph<RunGraphInNb, RunGraphIn, RunGraphOut> {
public:
  RunGraph()
      : hh::Graph<RunGraphInNb, RunGraphIn, RunGraphOut>("Run Graph") {
    auto runTask = std::make_shared<RunTask>(1);

    this->inputs(runTask);
    this->outputs(runTask);
  }
};

#endif
