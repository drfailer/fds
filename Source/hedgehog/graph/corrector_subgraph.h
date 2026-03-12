#ifndef CORRECTOR_SUBGRAPH_H
#define CORRECTOR_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../state/collector_state.h"
#include "../state/div_setup_state.h"
#include "../state/corr_particle_state.h"
#include "../state/velocity_corrector_state.h"
#include "../task/barrier_tasks.h"
#include "../task/corr_step1_kernel_task.h"
#include "../task/div_setup_kernel_task.h"
#include "../task/corr_condens_kernel_task.h"
#include "../task/corr_particle_kernel_task.h"
#include "../task/corr_div_part1_kernel_task.h"
#include "../task/divergence_part2_kernel_task.h"
#include "../task/velocity_corrector_kernel_task.h"
#include "wallbc_subgraph.h"
#include "velocity_bc_subgraph.h"
#include "corr_radiation_subgraph.h"

/// Build the Corrector sub-graph.
///
/// Implements the full corrector phase of the FDS time-stepping loop:
///   CorrStep1 -> MESH_EXCHANGE(4) -> CorrDivSetup -> Combustion -> HVAC ->
///   CorrCondens -> CorrParticle -> MESH_EXCHANGE(7) -> WallBC ->
///   MESH_EXCHANGE(6) -> CorrRadiation -> MESH_EXCHANGE(2) -> InitDivIntegrals ->
///   CorrDivPart1 -> DivergenceExchange -> CorrDivPart2 -> PressureIteration ->
///   VelocityCorrector -> MESH_EXCHANGE(6) -> CorrFinal
///
/// @param nmeshes Number of meshes
/// @param kernelThreads Number of threads for parallel kernel tasks
/// @return Shared pointer to the constructed sub-graph
inline auto buildCorrectorSubgraph(int nmeshes, size_t kernelThreads) {
    auto subgraph = std::make_shared<hh::Graph<1, MeshData, MeshData>>("Corrector");

    // --- Kernel tasks (MeshData -> MeshData, no orchestrator/collector needed) ---

    auto corrStep1KernelTask = std::make_shared<CorrStep1KernelTask>(kernelThreads);
    auto corrCondensKernelTask = std::make_shared<CorrCondensKernelTask>(kernelThreads);
    auto corrDivP1KernelTask = std::make_shared<CorrDivPart1KernelTask>(kernelThreads);
    auto corrDivP2KernelTask = std::make_shared<DivergencePart2KernelTask>(kernelThreads);

    // --- Sub-graphs with orchestrators (sequential pre-processing required) ---

    // CorrDivSetup: parallel kernel (+ sequential CC_VELOCITY_BC if CC_IBM)
    bool ccIBM = fds_is_cc_ibm() != 0;
    auto corrDivSetupKernelTask = std::make_shared<DivSetupKernelTask>(kernelThreads);

    // CorrParticle: sequential MASS_ENERGY + MOVE -> parallel MOMENTUM kernel
    auto corrParticleOrchSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
        std::make_shared<CorrParticleOrchestrator>(nmeshes), "CorrParticleOrch");
    auto corrParticleKernelTask = std::make_shared<CorrParticleKernelTask>(kernelThreads);

    // VelocityCorrector: parallel kernel (+ CC_PROJECT_VELOCITY orch/collector if CC_IBM)
    auto velCorrKernelTask = std::make_shared<VelocityCorrectorKernelTask>(kernelThreads);

    // --- Named sub-graphs ---

    auto wallBCSubgraph = buildWallBCSubgraph(nmeshes, kernelThreads);
    auto corrRadiationSubgraph = buildCorrRadiationSubgraph(nmeshes, kernelThreads);
    auto corrFinalSubgraph = buildCorrFinalSubgraph(nmeshes, kernelThreads);

    // --- Barrier tasks ---

    auto collector4SM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "Collector(4)");
    auto meshExchange4 = std::make_shared<MeshExchangeTask>(4);

    auto combustionCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "CombustionCollector");
    auto combustionTask = std::make_shared<CombustionTask>();

    auto corrHvacCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "CorrHvacCollector");
    auto corrHvacTask = std::make_shared<HvacTask>(1);

    auto collector7SM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "Collector(7)");
    auto meshExchange7 = std::make_shared<MeshExchangeTask>(7);

    auto collector6aSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "Collector(6a)");
    auto meshExchange6a = std::make_shared<MeshExchangeTask>(6);

    auto collector2SM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "Collector(2)");
    auto meshExchange2 = std::make_shared<MeshExchangeTask>(2);

    auto corrInitDivCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "CorrInitDivCollector");
    auto corrInitDivTask = std::make_shared<InitDivIntegralsTask>();

    auto corrDivCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "CorrDivCollector");
    auto corrDivExchangeTask = std::make_shared<DivergenceExchangeTask>(/*corrector=*/true);

    auto corrPressureCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "CorrPressureCollector");
    auto corrPressureTask = std::make_shared<PressureIterationTask>(/*predictor=*/false);

    auto collector6bSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "Collector(6b)");
    auto meshExchange6b = std::make_shared<MeshExchangeTask>(6);

    // --- Wire the sub-graph ---

    subgraph->inputs(corrStep1KernelTask);

    // CorrStep1 -> MESH_EXCHANGE(4)
    subgraph->edges(corrStep1KernelTask, collector4SM);
    subgraph->edges(collector4SM, meshExchange4);

    // CorrDivSetup: parallel kernel (with optional CC_VELOCITY_BC orchestrator if CC_IBM)
    if (ccIBM) {
        auto corrDivSetupOrchSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
            std::make_shared<CorrDivSetupOrchestrator>(nmeshes), "CorrDivSetupOrch");
        subgraph->edges(meshExchange4, corrDivSetupOrchSM);
        subgraph->edges(corrDivSetupOrchSM, corrDivSetupKernelTask);
    } else {
        subgraph->edges(meshExchange4, corrDivSetupKernelTask);
    }
    subgraph->edges(corrDivSetupKernelTask, combustionCollectorSM);
    subgraph->edges(combustionCollectorSM, combustionTask);
    subgraph->edges(combustionTask, corrHvacCollectorSM);
    subgraph->edges(corrHvacCollectorSM, corrHvacTask);

    // CorrCondens -> CorrParticle: orchestrator (particle ops) -> parallel kernel
    subgraph->edges(corrHvacTask, corrCondensKernelTask);
    subgraph->edges(corrCondensKernelTask, corrParticleOrchSM);
    subgraph->edges(corrParticleOrchSM, corrParticleKernelTask);
    subgraph->edges(corrParticleKernelTask, collector7SM);
    subgraph->edges(collector7SM, meshExchange7);

    // WallBC sub-graph
    subgraph->edges(meshExchange7, wallBCSubgraph);
    subgraph->edges(wallBCSubgraph, collector6aSM);
    subgraph->edges(collector6aSM, meshExchange6a);

    // CorrRadiation sub-graph
    subgraph->edges(meshExchange6a, corrRadiationSubgraph);
    subgraph->edges(corrRadiationSubgraph, collector2SM);
    subgraph->edges(collector2SM, meshExchange2);
    subgraph->edges(meshExchange2, corrInitDivCollectorSM);
    subgraph->edges(corrInitDivCollectorSM, corrInitDivTask);

    // CorrDivPart1 -> DivExchange
    subgraph->edges(corrInitDivTask, corrDivP1KernelTask);
    subgraph->edges(corrDivP1KernelTask, corrDivCollectorSM);
    subgraph->edges(corrDivCollectorSM, corrDivExchangeTask);

    // CorrDivPart2 -> Pressure
    subgraph->edges(corrDivExchangeTask, corrDivP2KernelTask);
    subgraph->edges(corrDivP2KernelTask, corrPressureCollectorSM);
    subgraph->edges(corrPressureCollectorSM, corrPressureTask);

    // VelocityCorrector: parallel kernel (+ CC_PROJECT_VELOCITY orch/collector if CC_IBM)
    if (ccIBM) {
        auto velCorrCCOrchSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
            std::make_shared<VelocityCorrectorCCOrchestrator>(nmeshes), "VelCorrCCOrch");
        auto velCorrCCCollSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
            std::make_shared<VelocityCorrectorCCCollector>(nmeshes), "VelCorrCCCollector");
        subgraph->edges(corrPressureTask, velCorrCCOrchSM);
        subgraph->edges(velCorrCCOrchSM, velCorrKernelTask);
        subgraph->edges(velCorrKernelTask, velCorrCCCollSM);
        subgraph->edges(velCorrCCCollSM, collector6bSM);
    } else {
        subgraph->edges(corrPressureTask, velCorrKernelTask);
        subgraph->edges(velCorrKernelTask, collector6bSM);
    }
    subgraph->edges(collector6bSM, meshExchange6b);

    // CorrFinal sub-graph
    subgraph->edges(meshExchange6b, corrFinalSubgraph);

    subgraph->outputs(corrFinalSubgraph);

    return subgraph;
}

#endif // CORRECTOR_SUBGRAPH_H
