#ifndef INIT_GRAPH_H
#define INIT_GRAPH_H
#include "../task/init/init_task.h"
#include <hedgehog/hedgehog.h>

#define InitGraphInNb 1
#define InitGraphIn Parameters<ParameterIds::None>
#define InitGraphOut Parameters<ParameterIds::None>

class InitGraph : public hh::Graph<InitGraphInNb, InitGraphIn, InitGraphOut> {
public:
  InitGraph()
      : hh::Graph<InitGraphInNb, InitGraphIn, InitGraphOut>("Init Graph") {
    auto initTask = std::make_shared<InitTask>(1);

    this->inputs(initTask);
    this->outputs(initTask);
  }
};

#endif
