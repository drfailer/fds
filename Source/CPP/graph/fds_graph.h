#ifndef FDS_GRAPH_H
#define FDS_GRAPH_H
#include "finish_graph.h"
#include "init_graph.h"
#include "run_graph.h"
#include <hedgehog/hedgehog.h>

#define FDSGraphInNb 1
#define FDSGraphIn Parameters<ParameterIds::None>
#define FDSGraphOut bool

class FDSGraph : public hh::Graph<FDSGraphInNb, FDSGraphIn, FDSGraphOut> {
public:
  FDSGraph() : hh::Graph<FDSGraphInNb, FDSGraphIn, FDSGraphOut>("FDS Graph") {
    auto initGraph = std::make_shared<InitGraph>();
    auto runGraph = std::make_shared<RunGraph>();
    auto finishGraph = std::make_shared<FinishGraph>();

    this->inputs(initGraph);
    this->edges(initGraph, runGraph);
    this->edges(runGraph, finishGraph);
    this->outputs(finishGraph);
  }
};

#endif
