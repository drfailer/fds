#ifndef FINISH_GRAPH_H
#define FINISH_GRAPH_H
#include "../task/finish/finish_task.h"
#include <hedgehog/hedgehog.h>

#define FinishGraphInNb 1
#define FinishGraphIn Parameters
#define FinishGraphOut Parameters

class FinishGraph : public hh::Graph<FinishGraphInNb, FinishGraphIn, FinishGraphOut> {
public:
  FinishGraph()
      : hh::Graph<FinishGraphInNb, FinishGraphIn, FinishGraphOut>("Finish Graph") {
    auto finishTask = std::make_shared<FinishTask>(1);

    this->inputs(finishTask);
    this->outputs(finishTask);
  }
};

#endif
