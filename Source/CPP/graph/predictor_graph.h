#ifndef PREDICTOR_GRAPH_H
#define PREDICTOR_GRAPH_H
#include "../data/parameters.h"
#include "../state/predictor/time_step_state.h"
#include "../state/predictor/time_step_state_manager.h"
#include "../task/run/predictor/change_time_step/compute_density_task.h"
#include "../task/run/predictor/change_time_step/compute_divergence2_task.h"
#include "../task/run/predictor/change_time_step/compute_velocity_task.h"
#include "../task/run/predictor/change_time_step/compute_wall_bc_task.h"
#include "../task/run/predictor/change_time_step/exchange_mesh_divergence_task.h"
#include "../task/run/predictor/change_time_step/exchange_values_task.h"
#include "../task/run/predictor/change_time_step/hvac_solver_task.h"
#include "../task/run/predictor/change_time_step/predict_velocity_task.h"
#include "../task/run/predictor/change_time_step/solve_pressure_task.h"
#include "../task/run/predictor/change_time_step/update_pressure_task.h"
#include "../task/run/predictor/change_time_step/update_time_step_index_task.h"
#include "../task/run/predictor/finish_predictor_task.h"
#include "../task/run/predictor/start_predictor_task.h"
#include "../task/signal_filter.h"
#include "../data/signal.h"
#include <hedgehog/hedgehog.h>

#define PredictorGraphInNb 2
#define PredictorGraphIn Parameters<ParameterIds::None>, Signal<Sigs::Stop>
#define PredictorGraphOut Parameters<ParameterIds::None>

class PredictorGraph : public hh::Graph<PredictorGraphInNb, PredictorGraphIn,
                                        PredictorGraphOut> {
public:
  PredictorGraph()
      : hh::Graph<PredictorGraphInNb, PredictorGraphIn, PredictorGraphOut>(
            "Predictor Graph") {
    // tasks
    auto startPredictorTask = std::make_shared<StartPredictorTask>(1);
    auto exchangeValuesTask = std::make_shared<ExchangeValuesTask>(1); // MPI
    auto computeDensityTask = std::make_shared<ComputeDensityTask>(1);
    auto computeVelocityTask = std::make_shared<ComputeVelocityTask>(1);
    auto hvacSolverTask = std::make_shared<HVACSolverTask>(1);
    auto exchangeMeshDivergenceTask =
        std::make_shared<ExchangeMeshDivergenceTask>(1); // MPI
    auto computeWallBCTask = std::make_shared<ComputeWallBCTask>(1);
    auto updateGlobalPressureTask = std::make_shared<UpdatePressureTask>(1);
    auto computeDivergence2Task = std::make_shared<ComputeDivergence2Task>(1);
    auto solvePressureTask = std::make_shared<SolvePressureTask>(1);
    auto predictVelocityTask = std::make_shared<PredictVelocityTask>(1);
    auto updateTimeStepIndexTask = std::make_shared<UpdateTimeStepIndexTask>(1);
    auto finishPredictorTask = std::make_shared<FinishPredictorTask>(1);
    // state managers
    auto timeStepStateManager = std::make_shared<TimeStepStateManager>(
        std::make_shared<TimeStepState>());
    // signals
    auto signalFilter = std::make_shared<SignalFilter<Sigs::Stop>>();

    // catch signals
    this->inputs(signalFilter);
    this->edges(signalFilter, timeStepStateManager);

    // the graph
    this->inputs(startPredictorTask);

    this->edges(startPredictorTask, computeDensityTask);
    /* this->edges(computeDensityTask, exchangeValuesTask); // MPI */
    this->edges(computeDensityTask, computeVelocityTask);
    this->edges(computeVelocityTask, hvacSolverTask);
    this->edges(hvacSolverTask, computeWallBCTask);
    /* this->edges(computeWallBCTask, exchangeMeshDivergenceTask); // MPI */
    this->edges(computeWallBCTask, updateGlobalPressureTask);
    this->edges(updateGlobalPressureTask, computeDivergence2Task);
    this->edges(computeDivergence2Task, solvePressureTask);
    this->edges(solvePressureTask, predictVelocityTask);
    this->edges(predictVelocityTask, timeStepStateManager); // instability ?
    this->edges(timeStepStateManager, updateTimeStepIndexTask);
    this->edges(updateTimeStepIndexTask, timeStepStateManager); // time step
                                                                // reduced ?

    this->edges(timeStepStateManager, computeDensityTask);  // loop
    this->edges(timeStepStateManager, finishPredictorTask); // end

    this->outputs(finishPredictorTask);
  }
};

#endif
