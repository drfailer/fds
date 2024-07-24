#ifndef RUN_GRAPH_H
#define RUN_GRAPH_H
#include "../state/main_loop_state.h"
#include "../state/main_loop_state_manager.h"
#include "../task/run/begin_main_loop_task.h"
#include "../task/run/dump_output_files_task.h"
#include "../task/run/run_corrector_task.h"
#include "predictor_graph.h"
#include <hedgehog/hedgehog.h>

#define RunGraphInNb 1
#define RunGraphIn Parameters<ParameterIds::None>
#define RunGraphOut bool

class RunGraph : public hh::Graph<RunGraphInNb, RunGraphIn, RunGraphOut> {
public:
  RunGraph() : hh::Graph<RunGraphInNb, RunGraphIn, RunGraphOut>("Run Graph") {
    // tasks
    auto beginMainLoopTask = std::make_shared<BeginMainLoopTask>(1);
    auto predictorGraph = std::make_shared<PredictorGraph>();
    auto runCorrectorTask = std::make_shared<RunCorrectorTask>(1);
    auto dumpOutputFilesTask = std::make_shared<DumpOutputFilesTask>(1);
    // states
    auto mainLoopStateManager = std::make_shared<MainLoopStateManager>(
        std::make_shared<MainLoopState>());

    this->inputs(beginMainLoopTask);

    this->edges(beginMainLoopTask, predictorGraph);
    this->edges(predictorGraph, runCorrectorTask);
    this->edges(runCorrectorTask, dumpOutputFilesTask);
    this->edges(dumpOutputFilesTask, mainLoopStateManager);
    this->edges(mainLoopStateManager, beginMainLoopTask);

    // WARN: this can cause a deadlock if these 2 state have more that 1 type in
    // common
    this->edges(mainLoopStateManager, predictorGraph);

    this->outputs(mainLoopStateManager);
  }
};

#endif
