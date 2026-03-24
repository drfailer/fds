#ifndef CORRECTOR_SUBGRAPH_H
#define CORRECTOR_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <service/comm_service.hpp>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../state/barrier_state.h"
#include "../state/div_setup_state.h"
#include "../state/fork_join_state.h"
#include "../task/barrier_tasks.h"
#include "../task/corr_step1_kernel_task.h"
#include "../task/mass_fd_kernel_task.h"
#include "../task/div_setup_kernel_task.h"
#include "../task/combustion_kernel_task.h"
#include "../task/pipeline_fork1_tasks.h"
#include "pipeline_fork1_vflux_subgraph.h"
#include "../task/corr_div_part1_kernel_task.h"
#include "../task/divergence_part2_kernel_task.h"
#include "../task/velocity_corrector_kernel_task.h"
#include "wallbc_subgraph.h"
#include "velocity_bc_subgraph.h"
#include "corr_radiation_subgraph.h"
#include "pressure_iteration_subgraph.h"
#include "../task/pipeline_fork2_tasks.h"

/// Build the Corrector sub-graph.
///
/// Optimizations applied:
///   - Opt 2: Condensation+PartME merged into Join1+Soot+Hvac barrier;
///            RemoveMove + PartMom + MeshExch(7) merged into single barrier
///   - Opt 3: MeshExchange(2) + CorrDivExchange merged (non-CC_IBM)
///   - Opt 4: MeshExchange(6b) removed — absorbed into CorrFinalOrchestrator
inline auto buildCorrectorSubgraph(int nmeshes, double tEnd, size_t kernelThreads,
                                    std::shared_ptr<TerminationSignal> termSignal,
                                    hh::comm::CommService *commService = nullptr) {
    auto subgraph = std::make_shared<hh::Graph<1, MeshData, BarrierData>>("Corrector");

    size_t meshThreads = static_cast<size_t>(nmeshes);

    // --- Kernel tasks ---

    auto corrStep1KernelTask = std::make_shared<CorrStep1KernelTask>(meshThreads);
    auto corrDivP2KernelTask = std::make_shared<DivergencePart2KernelTask>(meshThreads);
    auto velCorrKernelTask = std::make_shared<VelocityCorrectorKernelTask>(meshThreads);

    bool ccIBM = fds_is_cc_ibm() != 0;

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

    // Merged: Fork1 Join (2*N tokens) + Soot+HVAC + Condensation+PartME (Opt 2)
    auto sootHvacCondSM = makeBarrierSM(nmeshes, "Join1+Soot+Hvac+Condens",
        "SOOT_OXIDATION_LOOP\\nHVAC_CALC\\nCONDENSATION\\nPART_MASS_ENERGY",
        [](auto& meshes) {
            fds_soot_oxidation_loop(meshes[0]->dt);
            fds_hvac_calc(meshes[0]->t, meshes[0]->dt, 1);
            for (auto &md : meshes) {
                fds_condensation_kernel(md->nm, md->dt);
                fds_particle_mass_energy_kernel(md->nm, md->t, md->dt);
            }
        },
        2 * nmeshes);

    // Merged: Remove+Move + ParticleMom + MeshExchange(7) (Opt 2)
    auto removeMovePartMomMeshExch7SM = makeBarrierSM(nmeshes, "RemoveMove+PartMom+MeshExch7",
        "REMOVE_PARTICLES\\nMOVE_PARTICLES\\nPARTICLE_MOMENTUM\\nMESH_EXCHANGE(7)",
        [](auto& meshes) {
            for (auto &md : meshes) {
                fds_remove_particles(md->t, md->nm);
                fds_move_particles(md->t, md->dt, md->nm);
                fds_particle_momentum_kernel(md->nm, md->dt);
            }
            fds_mesh_exchange(7);
        });

    bool useParallelPressure = fds_use_pressure_subgraph() != 0;

    // --- Wire the sub-graph ---

    // CorrStep1 -> MeshExchange(4)
    subgraph->inputs(corrStep1KernelTask);
    subgraph->edges(corrStep1KernelTask, meshExchange4SM);

    // --- Fork 1: VFLUX || COMBUSTION → merged Join+Soot+HVAC+Condens barrier ---

    auto fork1VFluxSubgraph = buildFork1VFluxSubgraph(nmeshes, ccIBM);
    auto fork1CombTask = std::make_shared<Fork1CombKernelTask>(meshThreads);

    subgraph->edges(meshExchange4SM, fork1VFluxSubgraph);
    subgraph->edges(meshExchange4SM, fork1CombTask);

    // Merged: Join(2*N) + Soot+HVAC+Condens barrier (Opt 2)
    subgraph->edges(fork1VFluxSubgraph, sootHvacCondSM);
    subgraph->edges(fork1CombTask, sootHvacCondSM);

    // Merged: Remove+Move+PartMom+MeshExch7 (Opt 2) → WallBC
    subgraph->edges(sootHvacCondSM, removeMovePartMomMeshExch7SM);
    subgraph->edges(removeMovePartMomMeshExch7SM, wallBCSubgraph);

    // --- Fork 2: RADIATION || DIV_P1 (or sequential for CC_IBM) ---
    if (ccIBM) {
        auto meshExchange6aSM = makeBarrierSM(nmeshes, "MeshExchange(6a)",
            "MESH_EXCHANGE(6)",
            [](auto& meshes) { fds_mesh_exchange(6); });

        auto meshExchange2 = std::make_shared<MeshExchangeTask>(2, false, false, true);
        auto corrDivP1KernelTask = std::make_shared<CorrDivPart1KernelTask>(meshThreads);

        auto corrDivExchangeSM = makeBarrierSM(nmeshes, "CorrDivExchange",
            "EXCH_DIV_INFO\\nRTE_SOURCE_CORR\\nGLOBAL_MATRIX_REASSIGN",
            [](auto& meshes) {
                fds_exchange_divergence_info();
                fds_rte_source_correction();
                fds_global_matrix_reassign(0);
            });

        subgraph->edges(wallBCSubgraph, meshExchange6aSM);
        subgraph->edges(meshExchange6aSM, corrRadiationSubgraph);
        subgraph->edges(corrRadiationSubgraph, meshExchange2);
        subgraph->edges(meshExchange2, corrDivP1KernelTask);
        subgraph->edges(corrDivP1KernelTask, corrDivExchangeSM);
        subgraph->edges(corrDivExchangeSM, corrDivP2KernelTask);
    } else {
        // Merged: MeshExchange(6a) + InitDivIntegrals
        auto meshExch6aInitDivSM = makeBarrierSM(nmeshes, "MeshExch6a+InitDiv",
            "MESH_EXCHANGE(6)\\nINIT_DIV_INTEGRALS",
            [](auto& meshes) {
                fds_mesh_exchange(6);
                fds_initialize_divergence_integrals();
            });
        subgraph->edges(wallBCSubgraph, meshExch6aInitDivSM);

        auto fork2DivP1Task = std::make_shared<Fork2DivP1KernelTask>(meshThreads);
        auto fork2DivP1CollSM = std::make_shared<hh::StateManager<
            1, MeshData, BarrierData>>(
            std::make_shared<Fork2DivP1CollectorState>(nmeshes),
            "Fork2DivP1Collector");
        auto join2SM = std::make_shared<hh::StateManager<
            1, BarrierData, BarrierData>>(
            std::make_shared<BarrierJoinState>(2), "Join2");

        // Merged: MeshExch(2) + QR_ADD + DivExchange (Opt 3)
        auto corrMeshExch2DivExchTask = std::make_shared<CorrMeshExch2DivExchangeTask>();

        // Multicast to both branches (Hedgehog routes by type)
        subgraph->edges(meshExch6aInitDivSM, corrRadiationSubgraph);
        subgraph->edges(meshExch6aInitDivSM, fork2DivP1Task);
        subgraph->edges(corrRadiationSubgraph, join2SM);
        subgraph->edges(fork2DivP1Task, fork2DivP1CollSM);
        subgraph->edges(fork2DivP1CollSM, join2SM);
        // Merged: MeshExch(2) + QR + DivExch replaces join2→meshExch2→corrDivExch
        subgraph->edges(join2SM, corrMeshExch2DivExchTask);
        subgraph->edges(corrMeshExch2DivExchTask, corrDivP2KernelTask);
    }

    // --- Common downstream: DivP2 -> Pressure -> VelCorr -> CorrFinal ---

    if (useParallelPressure) {
        auto corrPressureSubgraph = buildPressureIterationSubgraph(
            tEnd, nmeshes, meshThreads, false, termSignal,
            commService, fds_get_pres_flag());
        subgraph->edges(corrDivP2KernelTask, corrPressureSubgraph);
        subgraph->edges(corrPressureSubgraph, velCorrKernelTask);
    } else {
        auto corrPressureSM = makeBarrierSM(nmeshes, "CorrPressure",
            "PRESSURE_ITERATION",
            [](auto& meshes) {
                fds_pressure_iteration(meshes[0]->t, meshes[0]->dt);
            });
        subgraph->edges(corrDivP2KernelTask, corrPressureSM);
        subgraph->edges(corrPressureSM, velCorrKernelTask);
    }

    // CorrFinal sub-graph (MeshExch6b merged into CorrFinalOrchestrator — Opt 4)
    subgraph->edges(velCorrKernelTask, corrFinalSubgraph);

    subgraph->outputs(corrFinalSubgraph);

    return subgraph;
}

#endif // CORRECTOR_SUBGRAPH_H
