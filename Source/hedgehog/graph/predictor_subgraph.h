#ifndef PREDICTOR_SUBGRAPH_H
#define PREDICTOR_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../data/pred_fork_data.h"
#include "../state/collector_state.h"
#include "../state/pred_step1_state.h"
#include "../state/div_setup_state.h"
#include "../state/pred_fork_state.h"
#include "../task/barrier_tasks.h"
#include "../task/pred_step1_kernel_task.h"
#include "../task/mass_fd_kernel_task.h"
#include "../task/density_pred_kernel_task.h"
#include "../task/div_setup_kernel_task.h"
#include "../task/pred_wall_div_kernel_task.h"
#include "../task/pred_fork_tasks.h"
#include "compute_viscosity_block_subgraph.h"
#include "velocity_flux_block_subgraph.h"
#include "../task/divergence_part2_kernel_task.h"
#include "divergence_part2_block_subgraph.h"
#include "velocity_predictor_block_subgraph.h"
#include "change_timestep_subgraph.h"
#include "velocity_bc_subgraph.h"
#include "wallbc_subgraph.h"
#include "wallbc_block_subgraph.h"
#include "pressure_iteration_subgraph.h"
#include "density_block_subgraph.h"
#include "pred_fork_vflux_subgraph.h"
#include "pred_fork_div_subgraph.h"

/// Build the Predictor sub-graph.
///
/// Implements the full predictor phase of the FDS time-stepping loop.
///
/// Non-CC_IBM (pipelined):
///   PredStep1 -> MESH_EXCHANGE(1) -> HVAC+InitDiv -> DIV_P1_prefork -> Fork
///     Branch A: VFLUX -> PART_MOM
///     Branch B: WallBC -> DIV_P1_early (WORK_BRANCH=2)
///   -> Join -> DIV_P1_late (WORK_BRANCH=2) -> DivExchange -> DivP2 ->
///   PressureIteration -> VelocityPredictor -> ChangeTimeStep ->
///   MESH_EXCHANGE(3) -> PredFinal -> PhaseTransition
///
/// CC_IBM (sequential):
///   PredStep1 -> MESH_EXCHANGE(1) -> VFLUX -> HVAC+InitDiv -> WallBC ->
///   PredWallDiv(PMOM+DIV_P1) -> DivExchange -> ... (unchanged)
///
/// @param nmeshes Number of meshes
/// @param tEnd Simulation end time (passed to ChangeTimeStep sub-graph)
/// @param kernelThreads Number of threads for parallel kernel tasks
/// @param termSignal Shared termination signal for pressure iteration sub-graph
/// @return Shared pointer to the constructed sub-graph
inline auto buildPredictorSubgraph(int nmeshes, double tEnd, size_t kernelThreads,
                                    size_t blockThreads, int numBlocks,
                                    std::shared_ptr<TerminationSignal> termSignal) {
    auto subgraph = std::make_shared<hh::Graph<1, MeshData, MeshData>>("Predictor");

    // --- Kernel sub-graph components ---

    // PredStep1: sequential INSERT_ALL_PARTICLES -> parallel kernels
    auto predStep1OrchSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
        std::make_shared<PredStep1Orchestrator>(nmeshes), "PredStep1Orch");
    auto predStep1KernelTask = std::make_shared<PredStep1KernelTask>(kernelThreads);

    // Viscosity block decomposition: if non-DEARDORFF/DYNSMAG/CC_IBM, use K-block parallel
    bool canBlockVisc = fds_compute_viscosity_can_block_decompose() != 0;

    // Density block decomposition: if non-CC_IBM and non-MMS, use K-block parallel
    bool canBlockDensity = fds_density_can_block_decompose() != 0 && numBlocks > 1;

    // DensityPred: parallel DENSITY_KERNEL (mesh-level fallback)
    auto densPredKernelTask = std::make_shared<DensityPredKernelTask>(kernelThreads);

    // Feature flags
    bool ccIBM = fds_is_cc_ibm() != 0;
    bool canBlockFlux = fds_velocity_flux_can_block_decompose(1) != 0;
    bool canBlockWallBC = fds_wall_bc_can_block_decompose() != 0;

    // PredDivPart2: parallel DIVERGENCE_PART_2_KERNEL
    bool canBlockDivP2 = fds_divergence_part_2_can_block_decompose() != 0 && numBlocks > 1;
    auto predDivP2KernelTask = std::make_shared<DivergencePart2KernelTask>(kernelThreads);

    // VelocityPredictor: block-decomposed kernel
    auto velPredSubgraph = buildVelocityPredictorBlockSubgraph(
        blockThreads, numBlocks);

    // PredFinal sub-graph (Pattern B, outputs BarrierData)
    auto predFinalSubgraph = buildPredFinalSubgraph(nmeshes, kernelThreads, blockThreads, numBlocks);

    // ChangeTimeStep sub-graph (CFL retry loop)
    auto changeTimeStepSubgraph = buildChangeTimeStepSubgraph(tEnd, nmeshes, kernelThreads);

    // --- Barrier tasks ---

    auto collector1SM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "Collector(1)");
    auto meshExchange1 = std::make_shared<MeshExchangeTask>(1, /*ccDensity=*/ccIBM);

    auto predHvacCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "PredHvacCollector");
    auto hvacInitDivTask = std::make_shared<HvacInitDivTask>(1);

    auto predDivCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "PredDivCollector");
    auto predDivExchangeTask = std::make_shared<DivergenceExchangeTask>(/*corrector=*/false);

    auto predPressureCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "PredPressureCollector");

    auto changeTimeStepCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "ChangeTimeStepCollector");

    auto collector3SM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "Collector(3)");
    auto meshExchange3 = std::make_shared<MeshExchangeTask>(3, /*ccDensity=*/false, /*ccEndStep=*/ccIBM);

    auto phaseTransTask = std::make_shared<PhaseTransitionTask>();

    // --- Wire the sub-graph ---

    subgraph->inputs(predStep1OrchSM);

    // PredStep1: orchestrator (INSERT_ALL_PARTICLES) -> viscosity -> mass_fd -> Density
    if (canBlockVisc) {
        auto predViscBlockSubgraph = buildComputeViscosityBlockSubgraph(
            nmeshes, blockThreads, numBlocks);
        auto predMassFDKernelTask = std::make_shared<MassFDKernelTask>(kernelThreads);
        subgraph->edges(predStep1OrchSM, predViscBlockSubgraph);
        subgraph->edges(predViscBlockSubgraph, predMassFDKernelTask);
        if (canBlockDensity) {
            auto predDensityBlockSubgraph = buildDensityBlockSubgraph(
                nmeshes, blockThreads, numBlocks);
            subgraph->edges(predMassFDKernelTask, predDensityBlockSubgraph);
            subgraph->edges(predDensityBlockSubgraph, collector1SM);
        } else {
            subgraph->edges(predMassFDKernelTask, densPredKernelTask);
            subgraph->edges(densPredKernelTask, collector1SM);
        }
    } else {
        subgraph->edges(predStep1OrchSM, predStep1KernelTask);
        if (canBlockDensity) {
            auto predDensityBlockSubgraph = buildDensityBlockSubgraph(
                nmeshes, blockThreads, numBlocks);
            subgraph->edges(predStep1KernelTask, predDensityBlockSubgraph);
            subgraph->edges(predDensityBlockSubgraph, collector1SM);
        } else {
            subgraph->edges(predStep1KernelTask, densPredKernelTask);
            subgraph->edges(densPredKernelTask, collector1SM);
        }
    }

    // MESH_EXCHANGE(1)
    subgraph->edges(collector1SM, meshExchange1);

    // --- Predictor middle section: Fork (non-CC_IBM) or Sequential (CC_IBM) ---

    if (!ccIBM) {
        // Pipelined path: HVAC+InitDiv → DIV_P1_prefork → Fork → Join → DIV_P1_late
        // Reorder: HVAC collects directly from MeshExchange(1) (safe: HVAC ⊥ VFLUX)
        subgraph->edges(meshExchange1, predHvacCollectorSM);
        subgraph->edges(predHvacCollectorSM, hvacInitDivTask);

        auto divP1PreforkTask = std::make_shared<DivP1PreforkTask>(kernelThreads);
        subgraph->edges(hvacInitDivTask, divP1PreforkTask);

        // Fork: (VFLUX + PART_MOM) || (WallBC + DIV_P1_early)
        auto predForkSM = std::make_shared<hh::StateManager<
            1, MeshData, PredForkVFluxWork, PredForkDivWork>>(
            std::make_shared<PredForkState>(), "PredFork");
        auto predForkVFluxSG = buildPredForkVFluxSubgraph(
            nmeshes, kernelThreads, blockThreads, numBlocks, canBlockFlux);
        auto predForkDivSG = buildPredForkDivSubgraph(
            nmeshes, kernelThreads, blockThreads, numBlocks, canBlockWallBC);
        auto predJoinSM = std::make_shared<hh::StateManager<
            2, PredForkVFluxResult, PredForkDivResult, MeshData>>(
            std::make_shared<PredJoinState>(), "PredJoin");

        subgraph->edges(divP1PreforkTask, predForkSM);
        // Branch A: VFLUX + PART_MOM
        subgraph->edges(predForkSM, predForkVFluxSG);
        subgraph->edges(predForkVFluxSG, predJoinSM);
        // Branch B: WallBC + DIV_P1_early (WORK_BRANCH=2)
        subgraph->edges(predForkSM, predForkDivSG);
        subgraph->edges(predForkDivSG, predJoinSM);

        // After join: DIV_P1_late (WORK_BRANCH=2, copies RTRM to WORK1 for DIV_P2)
        auto divP1LateTask = std::make_shared<DivP1LateTask>(kernelThreads);
        subgraph->edges(predJoinSM, divP1LateTask);
        subgraph->edges(divP1LateTask, predDivCollectorSM);
    } else {
        // CC_IBM sequential path: VFLUX → HVAC+InitDiv → WallBC → PredWallDiv
        if (canBlockFlux) {
            auto predDivSetupBlockSubgraph = buildVelocityFluxBlockSubgraph(
                nmeshes, blockThreads, numBlocks);
            subgraph->edges(meshExchange1, predDivSetupBlockSubgraph);
            subgraph->edges(predDivSetupBlockSubgraph, predHvacCollectorSM);
        } else {
            auto predDivSetupOrchSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
                std::make_shared<PredDivSetupOrchestrator>(nmeshes), "PredDivSetupOrch");
            auto predDivSetupKernelTask = std::make_shared<DivSetupKernelTask>(kernelThreads);
            subgraph->edges(meshExchange1, predDivSetupOrchSM);
            subgraph->edges(predDivSetupOrchSM, predDivSetupKernelTask);
            subgraph->edges(predDivSetupKernelTask, predHvacCollectorSM);
        }
        subgraph->edges(predHvacCollectorSM, hvacInitDivTask);

        auto predWallBCSubgraph = canBlockWallBC
            ? buildWallBCBlockSubgraph(nmeshes, blockThreads, numBlocks)
            : buildWallBCSubgraph(nmeshes, kernelThreads);
        subgraph->edges(hvacInitDivTask, predWallBCSubgraph);

        auto predWallDivKernelTask = std::make_shared<PredWallDivKernelTask>(kernelThreads);
        subgraph->edges(predWallBCSubgraph, predWallDivKernelTask);
        subgraph->edges(predWallDivKernelTask, predDivCollectorSM);
    }

    // --- Common downstream: DivExchange → DivP2 → Pressure → VelPred → ... ---

    subgraph->edges(predDivCollectorSM, predDivExchangeTask);

    // PredDivPart2 -> Pressure (block-decomposed or mesh-level)
    if (canBlockDivP2) {
        auto predDivP2BlockSubgraph = buildDivergencePart2BlockSubgraph(
            nmeshes, blockThreads, numBlocks);
        subgraph->edges(predDivExchangeTask, predDivP2BlockSubgraph);
        subgraph->edges(predDivP2BlockSubgraph, predPressureCollectorSM);
    } else {
        subgraph->edges(predDivExchangeTask, predDivP2KernelTask);
        subgraph->edges(predDivP2KernelTask, predPressureCollectorSM);
    }

    // Pressure iteration: parallel sub-graph or sequential fallback
    bool useParallelPressure = fds_use_pressure_subgraph() != 0;
    if (useParallelPressure) {
        auto predPressureSubgraph = buildPressureIterationSubgraph(
            tEnd, nmeshes, kernelThreads, /*predictor=*/true, termSignal,
            fds_get_pres_flag());
        subgraph->edges(predPressureCollectorSM, predPressureSubgraph);
        subgraph->edges(predPressureSubgraph, velPredSubgraph);
    } else {
        auto predPressureTask = std::make_shared<PressureIterationTask>(/*predictor=*/true);
        subgraph->edges(predPressureCollectorSM, predPressureTask);
        subgraph->edges(predPressureTask, velPredSubgraph);
    }

    // VelocityPredictor -> ChangeTimeStep
    subgraph->edges(velPredSubgraph, changeTimeStepCollectorSM);
    subgraph->edges(changeTimeStepCollectorSM, changeTimeStepSubgraph);
    subgraph->edges(changeTimeStepSubgraph, collector3SM);
    subgraph->edges(collector3SM, meshExchange3);

    // PredFinal (outputs BarrierData) -> PhaseTransition (no collector needed)
    subgraph->edges(meshExchange3, predFinalSubgraph);
    subgraph->edges(predFinalSubgraph, phaseTransTask);

    subgraph->outputs(phaseTransTask);

    return subgraph;
}

#endif // PREDICTOR_SUBGRAPH_H
