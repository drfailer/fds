#ifndef RUN_GRAPH_H
#define RUN_GRAPH_H
#include "../state/main_loop_state.h"
#include "../state/main_loop_state_manager.h"
#include "../task/run/begin_main_loop_task.h"
#include "../task/run/dump_output_files_task.h"
#include "../task/run/run_corrector_task.h"
#include "../task/run/run_predictor_task.h"
#include <hedgehog/hedgehog.h>

#define RunGraphInNb 1
#define RunGraphIn Parameters
#define RunGraphOut bool

class RunGraph : public hh::Graph<RunGraphInNb, RunGraphIn, RunGraphOut> {
public:
  RunGraph() : hh::Graph<RunGraphInNb, RunGraphIn, RunGraphOut>("Run Graph") {
    // tasks
    auto beginMainLoopTask = std::make_shared<BeginMainLoopTask>(1);
    auto runPredictorTask = std::make_shared<RunPredictorTask>(1);
    auto runCorrectorTask = std::make_shared<RunCorrectorTask>(1);
    auto dumpOutputFilesTask = std::make_shared<DumpOutputFilesTask>(1);
    // states
    auto mainLoopStateManager = std::make_shared<MainLoopStateManager>(
        std::make_shared<MainLoopState>());

    this->inputs(beginMainLoopTask);

    this->edges(beginMainLoopTask, runPredictorTask);
    this->edges(runPredictorTask, runCorrectorTask);
    this->edges(runCorrectorTask, dumpOutputFilesTask);
    this->edges(dumpOutputFilesTask, mainLoopStateManager);
    this->edges(mainLoopStateManager, beginMainLoopTask);

    this->outputs(mainLoopStateManager);
  }
};

#endif
