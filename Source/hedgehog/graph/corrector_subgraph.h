#ifndef CORRECTOR_SUBGRAPH_H
#define CORRECTOR_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../state/collector_state.h"
#include "../state/barrier_state.h"
#include "../state/div_setup_state.h"
#include "../state/fork_join_state.h"
#include "../state/pipeline_fork2_state.h"
#include "../task/barrier_tasks.h"
#include "../task/corr_step1_kernel_task.h"
#include "../task/mass_fd_kernel_task.h"
#include "../task/div_setup_kernel_task.h"
#include "../task/combustion_kernel_task.h"
#include "../task/pipeline_fork1_tasks.h"
#include "pipeline_fork1_vflux_subgraph.h"
#include "../task/corr_condens_kernel_task.h"
#include "../task/particle_mass_energy_kernel_task.h"
#include "../task/particle_momentum_kernel_task.h"
#include "../task/corr_div_part1_kernel_task.h"
#include "../task/divergence_part2_kernel_task.h"
#include "../task/velocity_corrector_kernel_task.h"
#include "wallbc_subgraph.h"
#include "velocity_bc_subgraph.h"
#include "corr_radiation_subgraph.h"
#include "pressure_iteration_subgraph.h"
#include "../task/pipeline_fork2_tasks.h"

/// Build the Corrector sub-graph.
inline auto buildCorrectorSubgraph(int nmeshes, double tEnd, size_t kernelThreads,
                                    std::shared_ptr<TerminationSignal> termSignal) {
    auto subgraph = std::make_shared<hh::Graph<1, MeshData, BarrierData>>("Corrector");

    size_t meshThreads = static_cast<size_t>(nmeshes);

    // --- Kernel tasks ---

    auto corrStep1KernelTask = std::make_shared<CorrStep1KernelTask>(meshThreads);
    auto corrCondensKernelTask = std::make_shared<CorrCondensKernelTask>(meshThreads);
    auto corrDivP1KernelTask = std::make_shared<CorrDivPart1KernelTask>(meshThreads);
    auto corrDivP2KernelTask = std::make_shared<DivergencePart2KernelTask>(meshThreads);
    auto velCorrKernelTask = std::make_shared<VelocityCorrectorKernelTask>(meshThreads);

    bool ccIBM = fds_is_cc_ibm() != 0;

    auto particleMassEnergyKernelTask = std::make_shared<ParticleMassEnergyKernelTask>(meshThreads);
    auto partMomKernelTask = std::make_shared<ParticleMomentumKernelTask>(meshThreads);

    // --- Sub-graphs ---

    auto wallBCSubgraph = buildWallBCSubgraph(nmeshes, meshThreads);
    auto corrRadiationSubgraph = buildCorrRadiationSubgraph(nmeshes, meshThreads);
    auto corrFinalSubgraph = buildCorrFinalSubgraph(nmeshes, meshThreads);

    // --- Barrier states ---

    auto meshExchange4SM = makeBarrierSM(nmeshes, "MeshExchange(4)",
        "CC_DENSITY\\nMESH_EXCHANGE(4)",
        [ccIBM](auto& meshes) {
            if (ccIBM) { fds_cc_density(meshes[0]->t, meshes[0]->dt); }
            fds_mesh_exchange(4);
        });

    auto sootHvacSM = makeBarrierSM(nmeshes, "Soot+Hvac",
        "SOOT_OXIDATION_LOOP\\nHVAC_CALC",
        [](auto& meshes) {
            fds_soot_oxidation_loop(meshes[0]->dt);
            fds_hvac_calc(meshes[0]->t, meshes[0]->dt, 1);
        });

    auto removeMoveSM = makeBarrierSM(nmeshes, "RemoveMove",
        "REMOVE_PARTICLES\\nMOVE_PARTICLES",
        [](auto& meshes) {
            for (auto &md : meshes) {
                fds_remove_particles(md->t, md->nm);
                fds_move_particles(md->t, md->dt, md->nm);
            }
        });

    auto meshExchange7SM = makeBarrierSM(nmeshes, "MeshExchange(7)",
        "MESH_EXCHANGE(7)",
        [](auto& meshes) { fds_mesh_exchange(7); });

    auto meshExchange6aSM = makeBarrierSM(nmeshes, "MeshExchange(6a)",
        "MESH_EXCHANGE(6)",
        [](auto& meshes) { fds_mesh_exchange(6); });

    auto meshExchange2 = std::make_shared<MeshExchangeTask>(2, false, false, ccIBM);

    auto corrDivExchangeSM = makeBarrierSM(nmeshes, "CorrDivExchange",
        "EXCHANGE_DIVERGENCE_INFO\\nRTE_SOURCE_CORRECTION\\nGLOBAL_MATRIX_REASSIGN",
        [](auto& meshes) {
            fds_exchange_divergence_info();
            fds_rte_source_correction();
            fds_global_matrix_reassign(0);
        });

    bool useParallelPressure = fds_use_pressure_subgraph() != 0;

    auto meshExchange6bSM = makeBarrierSM(nmeshes, "MeshExchange(6b)",
        "CC_END_STEP\\nMESH_EXCHANGE(6)",
        [ccIBM](auto& meshes) {
            if (ccIBM) { fds_cc_end_step(meshes[0]->t, meshes[0]->dt, 0); }
            fds_mesh_exchange(6);
        });

    // --- Wire the sub-graph ---

    // CorrStep1 -> MeshExchange(4)
    subgraph->inputs(corrStep1KernelTask);
    subgraph->edges(corrStep1KernelTask, meshExchange4SM);

    // --- Fork 1: VFLUX || COMBUSTION ---

    auto fork1VFluxSubgraph = buildFork1VFluxSubgraph(nmeshes, ccIBM);
    auto fork1CombTask = std::make_shared<Fork1CombKernelTask>(meshThreads);
    auto join1SM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
        std::make_shared<ForkJoinState>(2), "Join1");

    subgraph->edges(meshExchange4SM, fork1VFluxSubgraph);
    subgraph->edges(meshExchange4SM, fork1CombTask);
    subgraph->edges(fork1VFluxSubgraph, join1SM);
    subgraph->edges(fork1CombTask, join1SM);

    // After join: Soot+HVAC barrier
    subgraph->edges(join1SM, sootHvacSM);

    // CorrCondens -> Particle pipeline
    subgraph->edges(sootHvacSM, corrCondensKernelTask);
    subgraph->edges(corrCondensKernelTask, particleMassEnergyKernelTask);
    subgraph->edges(particleMassEnergyKernelTask, removeMoveSM);
    subgraph->edges(removeMoveSM, partMomKernelTask);
    subgraph->edges(partMomKernelTask, meshExchange7SM);

    // WallBC sub-graph
    subgraph->edges(meshExchange7SM, wallBCSubgraph);
    subgraph->edges(wallBCSubgraph, meshExchange6aSM);

    // --- Fork 2: RADIATION || DIV_P1 (or sequential for CC_IBM) ---
    if (ccIBM) {
        subgraph->edges(meshExchange6aSM, corrRadiationSubgraph);
        subgraph->edges(corrRadiationSubgraph, meshExchange2);
        subgraph->edges(meshExchange2, corrDivP1KernelTask);
        subgraph->edges(corrDivP1KernelTask, corrDivExchangeSM);
    } else {
        auto fork2SM = std::make_shared<hh::StateManager<
            1, MeshData, MeshData>>(
            std::make_shared<PipelineFork2State>(nmeshes), "Fork2");

        auto fork2DivP1Task = std::make_shared<Fork2DivP1KernelTask>(meshThreads);
        auto fork2DivP1CollSM = std::make_shared<hh::StateManager<
            1, MeshData, BarrierData>>(
            std::make_shared<Fork2DivP1CollectorState>(nmeshes),
            "Fork2DivP1Collector");
        auto join2SM = std::make_shared<hh::StateManager<
            1, BarrierData, BarrierData>>(
            std::make_shared<BarrierJoinState>(2), "Join2");
        auto qrAddTask = std::make_shared<DivP1QRAdditionTask>(meshThreads);

        subgraph->edges(meshExchange6aSM, fork2SM);
        subgraph->edges(fork2SM, corrRadiationSubgraph);
        subgraph->edges(corrRadiationSubgraph, join2SM);
        subgraph->edges(fork2SM, fork2DivP1Task);
        subgraph->edges(fork2DivP1Task, fork2DivP1CollSM);
        subgraph->edges(fork2DivP1CollSM, join2SM);
        subgraph->edges(join2SM, meshExchange2);
        subgraph->edges(meshExchange2, qrAddTask);
        subgraph->edges(qrAddTask, corrDivExchangeSM);
    }

    // --- Common downstream: DivExchange -> DivP2 -> Pressure -> VelCorr -> ... ---

    if (useParallelPressure) {
        auto corrPressureCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
            std::make_shared<CollectorState>(nmeshes), "CorrPressureCollector");
        auto corrPressureSubgraph = buildPressureIterationSubgraph(
            tEnd, nmeshes, meshThreads, false, termSignal,
            fds_get_pres_flag());
        subgraph->edges(corrDivExchangeSM, corrDivP2KernelTask);
        subgraph->edges(corrDivP2KernelTask, corrPressureCollectorSM);
        subgraph->edges(corrPressureCollectorSM, corrPressureSubgraph);
        subgraph->edges(corrPressureSubgraph, velCorrKernelTask);
    } else {
        auto corrPressureSM = makeBarrierSM(nmeshes, "CorrPressure",
            "PRESSURE_ITERATION",
            [](auto& meshes) {
                fds_pressure_iteration(meshes[0]->t, meshes[0]->dt);
            });
        subgraph->edges(corrDivExchangeSM, corrDivP2KernelTask);
        subgraph->edges(corrDivP2KernelTask, corrPressureSM);
        subgraph->edges(corrPressureSM, velCorrKernelTask);
    }

    subgraph->edges(velCorrKernelTask, meshExchange6bSM);

    // CorrFinal sub-graph
    subgraph->edges(meshExchange6bSM, corrFinalSubgraph);

    subgraph->outputs(corrFinalSubgraph);

    return subgraph;
}

#endif // CORRECTOR_SUBGRAPH_H
