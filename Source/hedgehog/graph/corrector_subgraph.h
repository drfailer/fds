#ifndef CORRECTOR_SUBGRAPH_H
#define CORRECTOR_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <service/comm_service.hpp>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../data/termination_data.h"
#include "../state/barrier_state.h"
#include "../state/div_setup_state.h"
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
#include "../task/wallbc_kernel_task.h"
#include "velocity_bc_subgraph.h"
#include "corr_radiation_subgraph.h"
#include "pressure_iteration_subgraph.h"
#include "../task/pipeline_fork2_tasks.h"
#include "../tool/thread_budget.h"

/// Build the Corrector sub-graph.
///
/// WallBC subgraph inlined: orchestrator+collector merged into adjacent barriers.
///   - Group A: Join1+Soot+Hvac+Condens+RemoveMove+PartMom+MeshExch7+WallBCOrch (2N→N)
///   - Group B: WallBCFinalize+ResetWallCounter+MeshExch6a+InitDiv (N→N)
///   - Group C: Join2+MeshExch2+DivExch (non-CC_IBM, 2N→N via makeBarrierSM)
inline auto buildCorrectorSubgraph(int nmeshes, const ThreadBudget &budget,
                                    std::shared_ptr<MeshDependencyGraph> depGraph = nullptr,
                                    hh::comm::CommService *commService = nullptr) {
    auto subgraph = std::make_shared<hh::Graph<2, MeshData, TerminationData, BarrierData>>("Corrector");

    // --- Kernel tasks (threads from budget) ---

    auto corrStep1KernelTask = std::make_shared<CorrStep1KernelTask>(budget.corrStep1);
    auto wallBCKernelTask = std::make_shared<WallBCKernelTask>(budget.corrWallBC);
    auto corrDivP2KernelTask = std::make_shared<DivergencePart2KernelTask>(budget.corrDivPart2);
    auto velCorrKernelTask = std::make_shared<VelocityCorrectorKernelTask>(budget.velCorrector);

    bool ccIBM = fds_is_cc_ibm() != 0;

    // --- Sub-graphs ---

    auto corrRadiationSubgraph = buildCorrRadiationSubgraph(nmeshes, budget.corrFork2Radiation);
    auto corrFinalSubgraph = buildCorrFinalSubgraph(nmeshes, budget.corrFinalVelBC);

    // --- Barrier states ---

    auto meshExchange4SM = makeBarrierSM(nmeshes, "MeshExchange(4)",
        "CC_DENSITY\\nMESH_EXCHANGE(4)",
        [ccIBM](auto& meshes) {
            if (ccIBM) { fds_cc_density(meshes[0]->t, meshes[0]->dt); }
            fds_mesh_exchange(4);
        });

    bool useParallelPressure = fds_use_pressure_subgraph() != 0;

    // --- Group A: Join1+Soot+Hvac+Condens+RemoveMove+PartMom+MeshExch7+WallBCOrch ---
    // Collects 2N tokens from fork1, then runs sequential operations and
    // computes WallBC global state (dt_bc, call_ht_1d) before parallel kernel.
    auto groupASM = makeBarrierSM(nmeshes, "Join1+RemoveMove+WallBCOrch",
        "SOOT+HVAC\\nCOND+PARTME\\nREMOVE+MOVE+PARTMOM\\nMESH_EXCHANGE(7)\\nWALLBC_ORCH",
        [](auto& meshes) {
            fds_soot_oxidation_loop(meshes[0]->dt);
            fds_hvac_calc(meshes[0]->t, meshes[0]->dt, 1);
            for (auto &md : meshes) {
                fds_condensation_kernel(md->nm, md->dt);
                fds_particle_mass_energy_kernel(md->nm, md->t, md->dt);
            }
            for (auto &md : meshes) {
                fds_remove_particles(md->t, md->nm);
                fds_move_particles(md->t, md->dt, md->nm);
                fds_particle_momentum_kernel(md->nm, md->dt);
            }
            fds_mesh_exchange(7);
            // WallBC global state (corrector phase)
            double dt_bc = fds_compute_wall_bc_dt_bc(meshes[0]->t);
            fds_increment_wall_counter();
            int call_ht_1d = fds_check_call_ht_1d();
            if (call_ht_1d) {
                fds_update_bc_clock(meshes[0]->t);
            }
            for (auto &md : meshes) {
                md->dt_bc = dt_bc;
                md->call_ht_1d = call_ht_1d;
            }
        },
        2 * nmeshes);

    // --- Wire the sub-graph ---

    // CorrStep1 -> MeshExchange(4)
    subgraph->inputs(corrStep1KernelTask);
    subgraph->edges(corrStep1KernelTask, meshExchange4SM);

    // --- Fork 1: VFLUX || COMBUSTION → Group A barrier ---

    auto fork1VFluxSubgraph = buildFork1VFluxSubgraph(nmeshes, ccIBM, budget.corrFork1DivSetup);
    auto fork1CombTask = std::make_shared<Fork1CombKernelTask>(budget.corrFork1Comb);

    subgraph->edges(meshExchange4SM, fork1VFluxSubgraph);
    subgraph->edges(meshExchange4SM, fork1CombTask);

    // Group A collects 2N tokens from fork1 branches
    subgraph->edges(fork1VFluxSubgraph, groupASM);
    subgraph->edges(fork1CombTask, groupASM);

    // Group A → WallBCKernel (parallel, inlined)
    subgraph->edges(groupASM, wallBCKernelTask);

    // --- Fork 2: RADIATION || DIV_P1 (or sequential for CC_IBM) ---
    if (ccIBM) {
        // Group B (CC_IBM): WallBCFinalize + MeshExch(6) — no InitDiv
        auto groupBSM = makeBarrierSM(nmeshes, "WallBCFin+MeshExch6a",
            "WALLBC_FINALIZE\\nRESET_WALL_COUNTER\\nMESH_EXCHANGE(6)",
            [](auto& meshes) {
                for (auto &md : meshes) {
                    fds_wall_bc_finalize(md->nm, md->t, md->dt_bc, md->call_ht_1d);
                }
                fds_reset_wall_counter();
                fds_mesh_exchange(6);
            });

        // MeshExch(2) + InitDiv barrier (radiation subgraph now emits MeshData)
        auto meshExch2SM = makeBarrierSM(nmeshes, "MeshExchange(2)",
            "MESH_EXCHANGE(2)\\nINIT_DIV_INTEGRALS",
            [](auto& meshes) {
                if (fds_exchange_radiation()) { fds_mesh_exchange(2); }
                fds_initialize_divergence_integrals();
            });

        auto corrDivP1KernelTask = std::make_shared<CorrDivPart1KernelTask>(budget.standalone(2));

        auto corrDivExchangeSM = makeBarrierSM(nmeshes, "CorrDivExchange",
            "EXCH_DIV_INFO\\nRTE_SOURCE_CORR\\nGLOBAL_MATRIX_REASSIGN\\nPRES_INIT+INCR",
            [useParallelPressure](auto& meshes) {
                fds_exchange_divergence_info();
                fds_rte_source_correction();
                fds_global_matrix_reassign(0);
                if (useParallelPressure) {
                    fds_pressure_iteration_init();
                    fds_pressure_iteration_increment();
                }
            });

        subgraph->edges(wallBCKernelTask, groupBSM);
        subgraph->edges(groupBSM, corrRadiationSubgraph);
        subgraph->edges(corrRadiationSubgraph, meshExch2SM);
        subgraph->edges(meshExch2SM, corrDivP1KernelTask);
        subgraph->edges(corrDivP1KernelTask, corrDivExchangeSM);
        subgraph->edges(corrDivExchangeSM, corrDivP2KernelTask);
    } else {
        // Group B (non-CC_IBM): WallBCFinalize + MeshExch(6) + InitDiv
        auto groupBSM = makeBarrierSM(nmeshes, "WallBCFin+MeshExch6a+InitDiv",
            "WALLBC_FINALIZE\\nRESET_WALL_COUNTER\\nMESH_EXCHANGE(6)\\nINIT_DIV_INTEGRALS",
            [](auto& meshes) {
                for (auto &md : meshes) {
                    fds_wall_bc_finalize(md->nm, md->t, md->dt_bc, md->call_ht_1d);
                }
                fds_reset_wall_counter();
                fds_mesh_exchange(6);
                fds_initialize_divergence_integrals();
            });

        subgraph->edges(wallBCKernelTask, groupBSM);

        auto fork2DivP1Task = std::make_shared<Fork2DivP1KernelTask>(budget.corrFork2DivP1);

        // Group C: Join2 + MeshExch(2) + QR + DivExchange (2N→N barrier)
        auto groupCSM = makeBarrierSM(nmeshes, "Join2+MeshExch2+DivExch",
            "MESH_EXCHANGE(2)\\nQR_ADD\\nEXCH_DIV_INFO\\nRTE_SOURCE_CORR\\nGLOBAL_MATRIX_REASSIGN",
            [useParallelPressure](auto& meshes) {
                if (fds_exchange_radiation()) { fds_mesh_exchange(2); }
                for (auto &md : meshes) {
                    fds_divergence_part_1_add_qr_b(md->nm);
                }
                fds_exchange_divergence_info();
                fds_rte_source_correction();
                fds_global_matrix_reassign(0);
                if (useParallelPressure) {
                    fds_pressure_iteration_init();
                    fds_pressure_iteration_increment();
                }
            },
            2 * nmeshes);

        // Multicast to both branches (Hedgehog routes by type)
        subgraph->edges(groupBSM, corrRadiationSubgraph);
        subgraph->edges(groupBSM, fork2DivP1Task);
        subgraph->edges(corrRadiationSubgraph, groupCSM);
        subgraph->edges(fork2DivP1Task, groupCSM);
        // Group C emits MeshData directly to DivP2
        subgraph->edges(groupCSM, corrDivP2KernelTask);
    }

    // --- Common downstream: DivP2 -> Pressure -> VelCorr -> CorrFinal ---

    if (useParallelPressure) {
        auto corrPressureSubgraph = buildPressureIterationSubgraph(
            nmeshes, budget, false,
            depGraph, commService, fds_get_pres_flag());
        subgraph->input<TerminationData>(corrPressureSubgraph);
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
        // Sink for TerminationData when parallel pressure is not used
        auto termSinkSM = std::make_shared<hh::StateManager<1, TerminationData, TerminationData>>(
            std::make_shared<TerminationDataSink>(), "TermDataSink");
        subgraph->input<TerminationData>(termSinkSM);
    }

    // CorrFinal sub-graph (MeshExch6b merged into CorrFinalOrchestrator — Opt 4)
    subgraph->edges(velCorrKernelTask, corrFinalSubgraph);

    subgraph->outputs(corrFinalSubgraph);

    return subgraph;
}

#endif // CORRECTOR_SUBGRAPH_H
