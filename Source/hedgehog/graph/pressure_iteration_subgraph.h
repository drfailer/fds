#ifndef PRESSURE_ITERATION_SUBGRAPH_H
#define PRESSURE_ITERATION_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/termination_data.h"
#include "../task/pressure_parallel_task.h"
#include "../state/pressure_convergence_state.h"
#include "../tool/thread_budget.h"

/// Build the pressure iteration sub-graph.
///
/// Shared between predictor and corrector phases. Phase routing is handled
/// internally via MeshData::phase (0=predictor, 1=corrector).
///
/// Pipeline:
///   Baroclinic → PreSolveExch → [external exchange] → SolvePhase → Solve
///   Solve → PostSolveExch → [external exchange] → VelErrorPhase → VelError
///   VelError → Pressure → Convergence → cycle or exit
///
/// Exchange I/O (PreSolveExch/PostSolveExch ↔ SolvePhase/VelErrorPhase)
/// is exposed as subgraph inputs/outputs, wired to an ExchangeGraph at the
/// fds_graph level.
inline auto buildPressureIterationSubgraph(int nmeshes,
                                            const ThreadBudget &budget,
                                            int presFlag = 0) {
    using SubGraphType = hh::Graph<5,
        MeshData<MeshState::PredictorPressure>,
        MeshData<MeshState::CorrectorPressure>,
        MeshData<MeshState::SolvePhase>,
        MeshData<MeshState::VelErrorPhase>,
        TerminationData,
        MeshData<MeshState::PredictorPressure>,
        MeshData<MeshState::CorrectorPressure>,
        MeshData<MeshState::PreSolveExch>,
        MeshData<MeshState::PostSolveExch>>;
    auto subgraph = std::make_shared<SubGraphType>("PressureIteration");

    // --- Packed parallel task (threads from budget) ---
    auto pressureParallelTask = std::make_shared<PressureParallelTask>(
        budget.pressureParallel, presFlag);

    // --- Convergence barrier (convergence check only) ---
    auto convergenceSM = std::make_shared<PressureConvergenceManager>(
        std::make_shared<PressureConvergenceState>(nmeshes),
        "PressureConvergence");

    // Entry: all mesh data types → pressureParallelTask
    subgraph->template input<MeshData<MeshState::PredictorPressure>>(pressureParallelTask);
    subgraph->template input<MeshData<MeshState::CorrectorPressure>>(pressureParallelTask);
    subgraph->template input<MeshData<MeshState::SolvePhase>>(pressureParallelTask);
    subgraph->template input<MeshData<MeshState::VelErrorPhase>>(pressureParallelTask);
    subgraph->template input<TerminationData>(convergenceSM);
    subgraph->template input<TerminationData>(pressureParallelTask);

    // VelError output (Pressure) → Convergence
    subgraph->edges(pressureParallelTask, convergenceSM);

    // Cycle: MeshData<Pressure> → back to pressureParallelTask (baroclinic phase)
    subgraph->template edge<MeshData<MeshState::Pressure>>(convergenceSM, pressureParallelTask);

    // Exchange outputs: PreSolveExch + PostSolveExch → subgraph output
    subgraph->template output<MeshData<MeshState::PreSolveExch>>(pressureParallelTask);
    subgraph->template output<MeshData<MeshState::PostSolveExch>>(pressureParallelTask);

    // Exit: PredPressure + CorrPressure → subgraph output
    subgraph->outputs(convergenceSM);

    return subgraph;
}

#endif // PRESSURE_ITERATION_SUBGRAPH_H
