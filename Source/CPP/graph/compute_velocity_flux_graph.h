#ifndef LOOP_FDS_COMPUTE_FLUX_GRAPH_H
#define LOOP_FDS_COMPUTE_FLUX_GRAPH_H

#include <hedgehog/hedgehog.h>
#include "../data/parameters.h"
#include "../task/split_task.h"
#include "../task/sum_task.h"
#include "../task/velocity_flux/compute_x_flux_task.h"
#include "../task/velocity_flux/compute_y_flux_task.h"
#include "../task/velocity_flux/compute_z_flux_task.h"
#include "../state/result_state.h"

class ComputeFluxGraph : public hh::Graph<1, Parameters, std::tuple<double, double, double>> {
 public:
  ComputeFluxGraph(size_t nbThreadsComputeTasks,
                   size_t nbThreadsSumTask,
                   size_t ibar, size_t jbar, size_t kbar)
    : hh::Graph<1, Parameters, std::tuple<double, double, double>>("Compute Flux Graph") {
    auto splitTask = std::make_shared<SplitTask>();
    auto computeXTask = std::make_shared<ComputeXFluxTask>(nbThreadsComputeTasks);
    auto computeYTask = std::make_shared<ComputeYFluxTask>(nbThreadsComputeTasks);
    auto computeZTask = std::make_shared<ComputeZFluxTask>(nbThreadsComputeTasks);
    auto sumTask = std::make_shared<SumTask>(nbThreadsSumTask);
    auto resultState = std::make_shared<ResultState>(ibar, jbar, kbar);
    auto resultStateManager = std::make_shared<hh::StateManager<RSInNb, RSIn, RSOut>>(resultState, "Result State Manager");

    this->inputs(splitTask);

    this->edges(splitTask, computeXTask);
    this->edges(splitTask, computeYTask);
    this->edges(splitTask, computeZTask);

    this->edges(computeXTask, sumTask);
    this->edges(computeYTask, sumTask);
    this->edges(computeZTask, sumTask);

    this->edges(sumTask, resultStateManager);

    this->outputs(resultStateManager);
  }
};

#endif //LOOP_FDS_COMPUTE_FLUX_GRAPH_H
