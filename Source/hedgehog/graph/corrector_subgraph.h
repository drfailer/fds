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
#include "../task/wallbc_kernel_task.h"
#include "../task/particle_ops_kernel_task.h"
#include "../task/qr_add_copy_kernel_task.h"
#include "velocity_bc_subgraph.h"
#include "corr_radiation_subgraph.h"
#include "../task/pipeline_fork2_tasks.h"
#include "../tool/thread_budget.h"

/// Build the Corrector sub-graph.
///
/// Barrier splits applied:
///   - Group A: Split into pre-barrier (soot+hvac) + ParticleOpsKernel + post-barrier (exchange+WallBC)
///   - Group B: WallBCFinalize stays in barrier (uses POINT_TO_MESH, not thread-safe)
///   - Group C: Split into pre-barrier (MeshExch2) + QRAddCopyKernel + post-barrier (DivExch)
template<MeshState PressureTag = MeshState::Default>
inline auto buildCorrectorSubgraphImpl(int nmeshes, const ThreadBudget &budget,
                                    std::shared_ptr<MeshDependencyGraph> depGraph = nullptr,
                                    hh::comm::CommService *commService = nullptr) {
    auto subgraph = std::make_shared<hh::Graph<2,
        MeshData<>, MeshData<MeshState::CorrectorPressure>,
        MeshData<>, MeshData<MeshState::CorrectorPressure>, BarrierData>>("Corrector");

    // --- Kernel tasks (threads from budget) ---

    auto corrStep1KernelTask = std::make_shared<CorrStep1KernelTask>(budget.corrStep1);
    auto wallBCKernelTask = std::make_shared<WallBCKernelTask>(budget.corrWallBC);
    auto corrDivP2KernelTask = std::make_shared<DivergencePart2KernelTask<PressureTag>>(budget.corrDivPart2);
    auto velCorrKernelTask = std::make_shared<VelocityCorrectorKernelTask<PressureTag>>(budget.velCorrector);

    bool ccIBM = fds_is_cc_ibm() != 0;

    // --- Sub-graphs ---

    auto corrRadiationSubgraph = buildCorrRadiationSubgraph(nmeshes, budget.corrFork2Radiation);
    auto corrFinalSubgraph = buildCorrFinalSubgraph(nmeshes, budget.corrFinalVelBC);

    // --- Barrier states ---

    auto meshExchange4SM = makeBarrierSM(nmeshes, "MeshExchange(4)",
        "MESH_EXCHANGE(4)",
        [](auto& meshes) {
            fds_mesh_exchange(4);
        });

    constexpr bool useParallelPressure = (PressureTag != MeshState::Default);

    // --- Group A split: pre-barrier + ParticleOpsKernel + post-barrier ---
    //
    // Pre-barrier (2N→N): joins fork1 branches, runs soot+hvac (global)
    auto groupAPreSM = makeBarrierSM(nmeshes, "Join1+Soot+Hvac",
        "SOOT+HVAC",
        [](auto& meshes) {
            fds_soot_oxidation_loop(meshes[0]->dt);
            fds_hvac_calc(meshes[0]->t, meshes[0]->dt, 1);
        },
        2 * nmeshes);

    // Extracted: condensation + particle mass/energy (parallel per-mesh)
    auto particleOpsKernelTask = std::make_shared<ParticleOpsKernelTask>(budget.corrParticleOps);

    // Post-barrier (N→N): exchange + WallBC orch (particle ops moved to kernel)
    auto groupAPostSM = makeBarrierSM(nmeshes, "MeshExch7+WallBCOrch",
        "MESH_EXCHANGE(7)\\nWALLBC_ORCH",
        [](auto& meshes) {
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
        });

    // --- Wire the sub-graph ---

    // CorrStep1 -> MeshExchange(4)
    subgraph->inputs(corrStep1KernelTask);
    subgraph->edges(corrStep1KernelTask, meshExchange4SM);

    // --- Fork 1: VFLUX || COMBUSTION → Group A split ---

    auto fork1VFluxSubgraph = buildFork1VFluxSubgraph(nmeshes, ccIBM, budget.corrFork1DivSetup);
    auto fork1CombTask = std::make_shared<Fork1CombKernelTask>(budget.corrFork1Comb);

    subgraph->edges(meshExchange4SM, fork1VFluxSubgraph);
    subgraph->edges(meshExchange4SM, fork1CombTask);

    // Group A pre-barrier collects 2N tokens from fork1 branches
    subgraph->edges(fork1VFluxSubgraph, groupAPreSM);
    subgraph->edges(fork1CombTask, groupAPreSM);

    // Group A: pre-barrier → ParticleOpsKernel → post-barrier → WallBCKernel
    subgraph->edges(groupAPreSM, particleOpsKernelTask);
    subgraph->edges(particleOpsKernelTask, groupAPostSM);
    subgraph->edges(groupAPostSM, wallBCKernelTask);

    // --- Fork 2: RADIATION || DIV_P1 (or sequential for CC_IBM) ---
    if (ccIBM) {
        // Group B (CC_IBM): WallBC kernel includes finalize → MeshExch(6)
        // RESET_WALL_COUNTER moved to CorrFinalOrch barrier
        auto groupBPostSM = makeBarrierSM(nmeshes, "MeshExch6a",
            "MESH_EXCHANGE(6)",
            [](auto& meshes) {
                fds_mesh_exchange(6);
            });

        // MeshExch(2) + InitDiv barrier
        auto meshExch2SM = makeBarrierSM(nmeshes, "MeshExchange(2)",
            "MESH_EXCHANGE(2)\\nINIT_DIV_INTEGRALS",
            [](auto& meshes) {
                if (fds_exchange_radiation()) { fds_mesh_exchange(2); }
                fds_initialize_divergence_integrals();
            });

        auto corrDivP1KernelTask = std::make_shared<CorrDivPart1KernelTask>(budget.standalone(2));

        // CC_IBM: per-mesh loop kept for GET_LINKED_VELOCITIES (cross-mesh writes)
        // rte_source_correction moved to CorrFinalOrch barrier
        auto corrDivExchangeSM = makeBarrierSM(nmeshes, "CorrDivExchange",
            "EXCH_DIV_INFO\\nZONE_OPS\\nGLOBAL_MATRIX_REASSIGN\\nPRES_INIT+INCR",
            [useParallelPressure](auto& meshes) {
                fds_exchange_divergence_info();
                // Zone ops + GET_LINKED_VELOCITIES (CC_IBM needs per-mesh for cross-mesh writes)
                for (auto &md : meshes) {
                    fds_divergence_part_2_preprocessing(md->nm, md->dt);
                }
                fds_global_matrix_reassign(0);
                if (useParallelPressure) {
                    fds_pressure_iteration_init();
                    fds_pressure_iteration_increment();
                }
            });

        subgraph->edges(wallBCKernelTask, groupBPostSM);
        subgraph->edges(groupBPostSM, corrRadiationSubgraph);
        subgraph->edges(corrRadiationSubgraph, meshExch2SM);
        subgraph->edges(meshExch2SM, corrDivP1KernelTask);
        subgraph->edges(corrDivP1KernelTask, corrDivExchangeSM);
        subgraph->edges(corrDivExchangeSM, corrDivP2KernelTask);
    } else {
        // Group B (non-CC_IBM): WallBC kernel includes finalize → MeshExch(6) + InitDiv
        // RESET_WALL_COUNTER moved to CorrFinalOrch barrier
        auto groupBPostSM = makeBarrierSM(nmeshes, "MeshExch6a+InitDiv",
            "MESH_EXCHANGE(6)\\nINIT_DIV_INTEGRALS",
            [](auto& meshes) {
                fds_mesh_exchange(6);
                fds_initialize_divergence_integrals();
            });

        subgraph->edges(wallBCKernelTask, groupBPostSM);

        auto fork2DivP1Task = std::make_shared<Fork2DivP1KernelTask>(budget.corrFork2DivP1);

        // Split Group C: pre-barrier(MeshExch2) → QRAddCopyKernel → post-barrier(DivExch)

        // Pre-barrier (2N→N): joins fork2, exchanges radiation data
        auto groupCPreSM = makeBarrierSM(nmeshes, "Join2+MeshExch2",
            "MESH_EXCHANGE(2)",
            [](auto& meshes) {
                if (fds_exchange_radiation()) { fds_mesh_exchange(2); }
            },
            2 * nmeshes);

        // Extracted: QR addition + WORK1 copy (parallel per-mesh)
        auto qrAddCopyKernelTask = std::make_shared<QRAddCopyKernelTask>(budget.corrQRAddCopy);

        // Post-barrier (N→N): divergence exchange + zone ops
        // rte_source_correction moved to CorrFinalOrch barrier
        // R_PBAR moved to block kernel (per-mesh, thread-safe)
        auto groupCPostSM = makeBarrierSM(nmeshes, "DivExch+ZoneOps",
            "EXCH_DIV_INFO\\nZONE_OPS\\nGLOBAL_MATRIX_REASSIGN\\nPRES_INIT+INCR",
            [useParallelPressure](auto& meshes) {
                fds_exchange_divergence_info();
                for (auto &md : meshes) {
                    fds_divergence_part_2_preprocessing(md->nm, md->dt);
                }
                fds_global_matrix_reassign(0);
                if (useParallelPressure) {
                    fds_pressure_iteration_init();
                    fds_pressure_iteration_increment();
                }
            });

        // Multicast to both branches (Hedgehog routes by type)
        subgraph->edges(groupBPostSM, corrRadiationSubgraph);
        subgraph->edges(groupBPostSM, fork2DivP1Task);
        // Group C split: pre → QR kernel → post
        subgraph->edges(corrRadiationSubgraph, groupCPreSM);
        subgraph->edges(fork2DivP1Task, groupCPreSM);
        subgraph->edges(groupCPreSM, qrAddCopyKernelTask);
        subgraph->edges(qrAddCopyKernelTask, groupCPostSM);
        // Post-barrier emits MeshData<> to DivP2
        subgraph->edges(groupCPostSM, corrDivP2KernelTask);
    }

    // --- Common downstream: DivP2 -> Pressure -> VelCorr -> CorrFinal ---

    if constexpr (useParallelPressure) {
        // Pressure handled externally via shared pressure subgraph.
        // DivP2<PressureTag> exits subgraph, VelCorr<PressureTag> receives from outside.
        subgraph->outputs(corrDivP2KernelTask);
        subgraph->template input<MeshData<PressureTag>>(velCorrKernelTask);
    } else {
        auto corrPressureSM = makeBarrierSM(nmeshes, "CorrPressure",
            "PRESSURE_ITERATION",
            [](auto& meshes) {
                fds_pressure_iteration(meshes[0]->t, meshes[0]->dt);
            });
        subgraph->edges(corrDivP2KernelTask, corrPressureSM);
        subgraph->edges(corrPressureSM, velCorrKernelTask);
    }

    // CorrFinal sub-graph (MeshExch6b in CorrFinalOrch task)
    subgraph->edges(velCorrKernelTask, corrFinalSubgraph);

    subgraph->outputs(corrFinalSubgraph);

    return subgraph;
}

/// Dispatch wrapper: selects the correct template instantiation at runtime.
inline auto buildCorrectorSubgraph(int nmeshes, const ThreadBudget &budget,
                                    std::shared_ptr<MeshDependencyGraph> depGraph = nullptr,
                                    hh::comm::CommService *commService = nullptr) {
    if (fds_use_pressure_subgraph()) {
        return buildCorrectorSubgraphImpl<MeshState::CorrectorPressure>(
            nmeshes, budget, depGraph, commService);
    }
    return buildCorrectorSubgraphImpl<MeshState::Default>(
        nmeshes, budget, depGraph, commService);
}

#endif // CORRECTOR_SUBGRAPH_H
