#ifndef CORRECTOR_SUBGRAPH_H
#define CORRECTOR_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../data/corr_step1_data.h"
#include "../data/div_setup_data.h"
#include "../data/corr_condens_data.h"
#include "../data/corr_particle_data.h"
#include "../data/corr_div_part1_data.h"
#include "../data/divergence_part2_data.h"
#include "../data/velocity_corrector_data.h"
#include "../state/collector_state.h"
#include "../state/mesh_barrier_state.h"
#include "../state/corr_step1_state.h"
#include "../state/div_setup_state.h"
#include "../state/corr_condens_state.h"
#include "../state/corr_particle_state.h"
#include "../state/corr_div_part1_state.h"
#include "../state/divergence_part2_state.h"
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

    // --- Kernel sub-graph components ---

    // CorrStep1: parallel VISCOSITY + MASS_FD + DENSITY kernels
    auto corrStep1OrchSM = std::make_shared<hh::StateManager<1, MeshData, CorrStep1Work>>(
        std::make_shared<CorrStep1Orchestrator>(nmeshes), "CorrStep1Orch");
    auto corrStep1KernelTask = std::make_shared<CorrStep1KernelTask>(kernelThreads);
    auto corrStep1CollectorSM = std::make_shared<hh::StateManager<1, CorrStep1Work, MeshData>>(
        std::make_shared<CorrStep1Collector>(nmeshes), "CorrStep1Collector");

    // CorrDivSetup: sequential VISCOSITY_BC + AGGLOMERATION + parallel VELOCITY_FLUX_KERNEL
    auto corrDivSetupOrchSM = std::make_shared<hh::StateManager<1, MeshData, DivSetupWork>>(
        std::make_shared<CorrDivSetupOrchestrator>(nmeshes), "CorrDivSetupOrch");
    auto corrDivSetupKernelTask = std::make_shared<DivSetupKernelTask>(kernelThreads);
    auto corrDivSetupCollectorSM = std::make_shared<hh::StateManager<1, DivSetupWork, MeshData>>(
        std::make_shared<DivSetupCollector>(nmeshes), "CorrDivSetupCollector");

    // CorrCondens: parallel CONDENSATION_EVAPORATION_KERNEL
    auto corrCondensOrchSM = std::make_shared<hh::StateManager<1, MeshData, CorrCondensWork>>(
        std::make_shared<CorrCondensOrchestrator>(nmeshes), "CorrCondensOrch");
    auto corrCondensKernelTask = std::make_shared<CorrCondensKernelTask>(kernelThreads);
    auto corrCondensCollectorSM = std::make_shared<hh::StateManager<1, CorrCondensWork, MeshData>>(
        std::make_shared<CorrCondensCollector>(nmeshes), "CorrCondensCollector");
    auto corrCondensBarrierSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
        std::make_shared<PassthroughBarrierState>(nmeshes), "CorrCondensBarrier");

    // CorrParticle: sequential MASS_ENERGY + MOVE + parallel MOMENTUM kernel
    auto corrParticleOrchSM = std::make_shared<hh::StateManager<1, MeshData, CorrParticleWork>>(
        std::make_shared<CorrParticleOrchestrator>(nmeshes), "CorrParticleOrch");
    auto corrParticleKernelTask = std::make_shared<CorrParticleKernelTask>(kernelThreads);
    auto corrParticleCollectorSM = std::make_shared<hh::StateManager<1, CorrParticleWork, MeshData>>(
        std::make_shared<CorrParticleCollector>(nmeshes), "CorrParticleCollector");

    // WallBC sub-graph (Pattern B)
    auto wallBCSubgraph = buildWallBCSubgraph(nmeshes, kernelThreads);

    // CorrRadiation sub-graph (Pattern A with global accumulators)
    auto corrRadiationSubgraph = buildCorrRadiationSubgraph(nmeshes, kernelThreads);

    // CorrDivPart1: sequential COMBUSTION_BC + parallel DIVERGENCE_PART_1_KERNEL
    auto corrDivP1OrchSM = std::make_shared<hh::StateManager<1, MeshData, CorrDivPart1Work>>(
        std::make_shared<CorrDivPart1Orchestrator>(nmeshes), "CorrDivP1Orch");
    auto corrDivP1KernelTask = std::make_shared<CorrDivPart1KernelTask>(kernelThreads);
    auto corrDivP1CollectorSM = std::make_shared<hh::StateManager<1, CorrDivPart1Work, MeshData>>(
        std::make_shared<CorrDivPart1Collector>(nmeshes), "CorrDivP1Collector");

    // CorrDivPart2: parallel DIVERGENCE_PART_2_KERNEL
    auto corrDivP2OrchSM = std::make_shared<hh::StateManager<1, MeshData, DivergencePart2Work>>(
        std::make_shared<DivergencePart2Orchestrator>(nmeshes), "CorrDivP2Orch");
    auto corrDivP2KernelTask = std::make_shared<DivergencePart2KernelTask>(kernelThreads);
    auto corrDivP2CollectorSM = std::make_shared<hh::StateManager<1, DivergencePart2Work, MeshData>>(
        std::make_shared<DivergencePart2Collector>(nmeshes), "CorrDivP2Collector");

    // VelocityCorrector: parallel VELOCITY_CORRECTOR_KERNEL
    auto velCorrOrchSM = std::make_shared<hh::StateManager<1, MeshData, VelocityCorrectorWork>>(
        std::make_shared<VelocityCorrectorOrchestrator>(nmeshes), "VelCorrOrch");
    auto velCorrKernelTask = std::make_shared<VelocityCorrectorKernelTask>(kernelThreads);
    auto velCorrCollectorSM = std::make_shared<hh::StateManager<1, VelocityCorrectorWork, MeshData>>(
        std::make_shared<VelocityCorrectorCollector>(nmeshes), "VelCorrCollector");

    // CorrFinal sub-graph (Pattern B)
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

    subgraph->inputs(corrStep1OrchSM);

    // CorrStep1 sub-graph
    subgraph->edges(corrStep1OrchSM, corrStep1KernelTask);
    subgraph->edges(corrStep1KernelTask, corrStep1CollectorSM);
    subgraph->edges(corrStep1CollectorSM, collector4SM);
    subgraph->edges(collector4SM, meshExchange4);

    // CorrDivSetup sub-graph
    subgraph->edges(meshExchange4, corrDivSetupOrchSM);
    subgraph->edges(corrDivSetupOrchSM, corrDivSetupKernelTask);
    subgraph->edges(corrDivSetupKernelTask, corrDivSetupCollectorSM);
    subgraph->edges(corrDivSetupCollectorSM, combustionCollectorSM);
    subgraph->edges(combustionCollectorSM, combustionTask);
    subgraph->edges(combustionTask, corrHvacCollectorSM);
    subgraph->edges(corrHvacCollectorSM, corrHvacTask);

    // CorrCondens sub-graph
    subgraph->edges(corrHvacTask, corrCondensOrchSM);
    subgraph->edges(corrCondensOrchSM, corrCondensKernelTask);
    subgraph->edges(corrCondensKernelTask, corrCondensCollectorSM);
    subgraph->edges(corrCondensCollectorSM, corrCondensBarrierSM);

    // CorrParticle sub-graph
    subgraph->edges(corrCondensBarrierSM, corrParticleOrchSM);
    subgraph->edges(corrParticleOrchSM, corrParticleKernelTask);
    subgraph->edges(corrParticleKernelTask, corrParticleCollectorSM);
    subgraph->edges(corrParticleCollectorSM, collector7SM);
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

    // CorrDivPart1 sub-graph
    subgraph->edges(corrInitDivTask, corrDivP1OrchSM);
    subgraph->edges(corrDivP1OrchSM, corrDivP1KernelTask);
    subgraph->edges(corrDivP1KernelTask, corrDivP1CollectorSM);
    subgraph->edges(corrDivP1CollectorSM, corrDivCollectorSM);
    subgraph->edges(corrDivCollectorSM, corrDivExchangeTask);

    // CorrDivPart2 sub-graph
    subgraph->edges(corrDivExchangeTask, corrDivP2OrchSM);
    subgraph->edges(corrDivP2OrchSM, corrDivP2KernelTask);
    subgraph->edges(corrDivP2KernelTask, corrDivP2CollectorSM);
    subgraph->edges(corrDivP2CollectorSM, corrPressureCollectorSM);
    subgraph->edges(corrPressureCollectorSM, corrPressureTask);

    // VelocityCorrector sub-graph
    subgraph->edges(corrPressureTask, velCorrOrchSM);
    subgraph->edges(velCorrOrchSM, velCorrKernelTask);
    subgraph->edges(velCorrKernelTask, velCorrCollectorSM);
    subgraph->edges(velCorrCollectorSM, collector6bSM);
    subgraph->edges(collector6bSM, meshExchange6b);

    // CorrFinal sub-graph
    subgraph->edges(meshExchange6b, corrFinalSubgraph);

    subgraph->outputs(corrFinalSubgraph);

    return subgraph;
}

#endif // CORRECTOR_SUBGRAPH_H
