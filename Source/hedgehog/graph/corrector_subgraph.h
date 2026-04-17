#ifndef CORRECTOR_SUBGRAPH_H
#define CORRECTOR_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <service/comm_service.hpp>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../data/termination_data.h"
#include "../state/barrier_state.h"
#include "../state/fork_join_state.h"
#include "../task/barrier_tasks.h"
#include "../task/corr_step1_kernel_task.h"
#include "../task/mass_fd_kernel_task.h"
#include "../task/div_setup_kernel_task.h"
#include "../task/combustion_kernel_task.h"
#include "../task/pipeline_fork1_tasks.h"
#include "../task/corr_div_part1_kernel_task.h"
#include "../task/divergence_part2_kernel_task.h"
#include "../task/velocity_corrector_kernel_task.h"
#include "../task/wallbc_kernel_task.h"
#include "../task/particle_ops_kernel_task.h"
#include "../task/corr_div_parallel_task.h"
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
    auto subgraph = std::make_shared<hh::Graph<3,
        MeshData<>, MeshData<MeshState::CorrectorPressure>, TerminationData,
        MeshData<>, MeshData<MeshState::CorrectorPressure>, BarrierData>>("Corrector");

    // --- Kernel tasks (threads from budget) ---

    auto corrStep1KernelTask = std::make_shared<CorrStep1KernelTask>(budget.corrStep1);
    auto wallBCKernelTask = std::make_shared<WallBCKernelTask>(budget.corrWallBC);
    auto velCorrKernelTask = std::make_shared<VelocityCorrectorKernelTask<PressureTag>>(budget.velCorrector);

    bool ccIBM = fds_is_cc_ibm() != 0;
    bool ht3d = fds_is_ht3d() != 0;

    // --- Sub-graphs ---

    auto corrRadiationSubgraph = buildCorrRadiationSubgraph(nmeshes, budget.corrFork2Radiation);

    // --- Barrier states ---

    auto meshExchange4SM = makeBarrierSM(nmeshes, "MeshExchange(4)",
        "MESH_EXCHANGE(4)",
        [](auto& meshes) {
            fds_mesh_exchange(4);
        });

    constexpr bool useParallelPressure = (PressureTag != MeshState::Default);

    // --- Group A: ForkJoin(2N→N) → Fork(HVAC || ParticleOps) → ForkJoin ---
    //
    // Soot oxidation merged into Fork1CombKernelTask (per-mesh, thread-safe).
    // HVAC_CALC is truly global (network solve) — independent of ParticleOps.
    // HVAC modifies global DUCT/DUCTNODE arrays; ParticleOps modifies per-mesh
    // particle data. No data dependency → run in parallel.
    auto fork1JoinTask = std::make_shared<ForkJoinTask>(nmeshes, 2, 1, "Fork1Join");

    // HVAC: global network solve. Collects N MeshData, emits 1 BarrierData.
    auto hvacBarrier = makeBarrierCollectToOne(nmeshes, "HvacCalc",
        "HVAC_CALC",
        [](auto& meshes) {
            fds_hvac_calc(meshes[0]->t, meshes[0]->dt, 1);
        });

    // Extracted: condensation + particle mass/energy (parallel per-mesh)
    auto particleOpsKernelTask = std::make_shared<ParticleOpsKernelTask>(budget.corrParticleOps);

    // MeshExch(7) on ParticleOps branch — only needs particle data, not HVAC.
    // Collects N MeshData, emits 1 BarrierData (carries meshes for downstream).
    auto meshExch7Barrier = makeBarrierCollectToOne(nmeshes, "MeshExchange(7)",
        "MESH_EXCHANGE(7)",
        [](auto& meshes) {
            fds_mesh_exchange(7);
        });

    // Join 2 BarrierData (HVAC + MeshExch7), scatter N MeshData
    auto hvacPartJoinTask = std::make_shared<BarrierJoinTask>(2, "HvacPartJoin");

    // --- Wire the sub-graph ---

    // CorrStep1 -> MeshExchange(4)
    subgraph->inputs(corrStep1KernelTask);
    subgraph->edges(corrStep1KernelTask, meshExchange4SM);

    // --- Fork 1: VFLUX || COMBUSTION → Group A split ---

    auto fork1DivSetupTask = std::make_shared<DivSetupKernelTask>(budget.corrFork1DivSetup);
    auto fork1CombTask = std::make_shared<Fork1CombKernelTask>(budget.corrFork1Comb);

    subgraph->edges(meshExchange4SM, fork1DivSetupTask);
    subgraph->edges(meshExchange4SM, fork1CombTask);

    // ForkJoin: collects 2N MeshData (DivSetup + Combustion+Soot), emits N
    subgraph->edges(fork1DivSetupTask, fork1JoinTask);
    subgraph->edges(fork1CombTask, fork1JoinTask);

    // Fork: HVAC(→BarrierData) || (ParticleOps → MeshExch7 →BarrierData) → Join → WallBC
    subgraph->edges(fork1JoinTask, hvacBarrier);
    subgraph->edges(fork1JoinTask, particleOpsKernelTask);
    subgraph->edges(particleOpsKernelTask, meshExch7Barrier);
    subgraph->edges(hvacBarrier, hvacPartJoinTask);
    subgraph->edges(meshExch7Barrier, hvacPartJoinTask);
    subgraph->edges(hvacPartJoinTask, wallBCKernelTask);

    // --- Fork 2: RADIATION || DIV_P1 (or sequential for CC_IBM) ---
    if (ccIBM) {
        // Group B (CC_IBM): Exchange(6) barrier for back wall data (HT3D only).
        // Radiation forked from WallBC directly — doesn't need Exchange(6) data.
        // RESET_WALL_COUNTER moved to CorrFinalOrch barrier
        auto groupBPostSM = makeBarrierSM(nmeshes, "MeshExch6a",
            "MESH_EXCHANGE(6) [HT3D]",
            [ht3d](auto& meshes) {
                if (ht3d && meshes[0]->call_ht_1d) { fds_mesh_exchange(6); }
            });

        // Join (N MeshData + 1 BarrierData): Radiation(BarrierData) + Exchange(6)(MeshData)
        auto meshExch2SM = makeDualInputBarrier(nmeshes, "Join+MeshExchange(2)",
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
                    // Pre-loop: link cut-face velocity fluxes before pressure iterations
                    for (auto &md : meshes) {
                        fds_get_linked_fv(md->nm, 0); // DO_BAROCLINIC=FALSE
                    }
                }
            });

        // Fork: Radiation starts immediately from WallBC (no Exchange(6) dependency).
        //        Exchange(6) runs in parallel. Both join at meshExch2SM.
        subgraph->edges(wallBCKernelTask, corrRadiationSubgraph);
        subgraph->edges(wallBCKernelTask, groupBPostSM);
        subgraph->edges(corrRadiationSubgraph, meshExch2SM);
        subgraph->edges(groupBPostSM, meshExch2SM);
        subgraph->edges(meshExch2SM, corrDivP1KernelTask);
        subgraph->edges(corrDivP1KernelTask, corrDivExchangeSM);
        // CC_IBM keeps standalone DivPart2 kernel (not packed)
        auto corrDivP2KernelTask = std::make_shared<DivergencePart2KernelTask<PressureTag>>(budget.corrDivPart2);
        subgraph->edges(corrDivExchangeSM, corrDivP2KernelTask);

        // TerminationData sink (no cycle to break in CC_IBM)
        auto termSink = std::make_shared<TerminationDataSink>();
        subgraph->template input<TerminationData>(termSink);

        // Downstream: DivP2 → Pressure → VelCorr
        if constexpr (useParallelPressure) {
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
    } else {
        // Group B (non-CC_IBM): Exchange(6) [HT3D only] + InitDiv barrier.
        // Radiation forked from WallBC directly — doesn't need Exchange(6) data.
        // DivP1 waits for this barrier (needs InitDiv zeroed arrays).
        // Exchange(6) conditional: at this point velocities haven't been corrected yet,
        // so only new data since Exchange(3) is back wall info needed for HT3D.
        // RESET_WALL_COUNTER moved to CorrFinalOrch barrier
        auto groupBPostSM = makeBarrierSM(nmeshes, "MeshExch6a+InitDiv",
            "MESH_EXCHANGE(6) [HT3D]\\nINIT_DIV_INTEGRALS",
            [ht3d](auto& meshes) {
                if (ht3d && meshes[0]->call_ht_1d) { fds_mesh_exchange(6); }
                fds_initialize_divergence_integrals();
            });

        auto fork2DivP1Task = std::make_shared<Fork2DivP1KernelTask>(budget.corrFork2DivP1);

        // Divergence pipeline: packed parallel task + 2 retagging barriers.
        // The 3 parallel kernels (QRAddCopy, DivP2Pre, DivPart2) share one
        // thread pool via CorrDivParallelTask. Sequential barriers are separate.
        //
        // Pipeline: Join2+MeshExch2 → Parallel(QRAddCopy) → DivExchange barrier
        //   → Parallel(DivP2Pre) → GlobalMatrix barrier → Parallel(DivPart2)
        //   → downstream

        // Pre-barrier (N MeshData + 1 BarrierData): joins fork2 + radiation
        auto groupCPreSM = makeDualInputBarrier(nmeshes, "Join2+MeshExch2",
            "MESH_EXCHANGE(2)",
            [](auto& meshes) {
                if (fds_exchange_radiation()) { fds_mesh_exchange(2); }
            });

        auto corrDivParallelTask = std::make_shared<CorrDivParallelTask<PressureTag>>(
            budget.corrDivParallel);

        auto corrDivExchBarrier = makeRetaggingBarrier<MeshState::DivExch, MeshState::DivP2Pre>(
            nmeshes, "DivExchange",
            "EXCH_DIV_INFO",
            [](auto&) {
                fds_exchange_divergence_info();
            });

        auto corrGlobalMatBarrier = makeRetaggingBarrier<MeshState::GlobalMat, MeshState::DivPart2>(
            nmeshes, "GlobalMatrix+PressureInit",
            "GLOBAL_MATRIX_REASSIGN\\nPRES_INIT+INCR",
            [useParallelPressure](auto&) {
                fds_global_matrix_reassign(0);
                if (useParallelPressure) {
                    fds_pressure_iteration_init();
                    fds_pressure_iteration_increment();
                }
            });

        // Fork: Radiation starts immediately from WallBC (no Exchange(6) dependency).
        //        DivP1 waits for Exchange(6) barrier.
        subgraph->edges(wallBCKernelTask, corrRadiationSubgraph);
        subgraph->edges(wallBCKernelTask, groupBPostSM);
        subgraph->edges(groupBPostSM, fork2DivP1Task);

        // Group C: Join2 → Parallel ↔ DivExchange ↔ Parallel ↔ GlobalMatrix ↔ Parallel
        subgraph->edges(corrRadiationSubgraph, groupCPreSM);
        subgraph->edges(fork2DivP1Task, groupCPreSM);
        subgraph->edges(groupCPreSM, corrDivParallelTask);
        subgraph->edges(corrDivParallelTask, corrDivExchBarrier);
        subgraph->edges(corrDivExchBarrier, corrDivParallelTask);
        subgraph->edges(corrDivParallelTask, corrGlobalMatBarrier);
        subgraph->edges(corrGlobalMatBarrier, corrDivParallelTask);

        // TerminationData breaks structural cycle at shutdown
        subgraph->template input<TerminationData>(corrDivParallelTask);

        // Downstream: final output (MeshData<PressureTag>) → Pressure → VelCorr
        if constexpr (useParallelPressure) {
            subgraph->template output<MeshData<PressureTag>>(corrDivParallelTask);
            subgraph->template input<MeshData<PressureTag>>(velCorrKernelTask);
        } else {
            auto corrPressureSM = makeBarrierSM(nmeshes, "CorrPressure",
                "PRESSURE_ITERATION",
                [](auto& meshes) {
                    fds_pressure_iteration(meshes[0]->t, meshes[0]->dt);
                });
            subgraph->edges(corrDivParallelTask, corrPressureSM);
            subgraph->edges(corrPressureSM, velCorrKernelTask);
        }
    }

    // --- CorrFinal (inlined): Orch → Fork(VelBCEdges || RTE) → Dump ---

    auto corrFinalOrchTask = std::make_shared<CorrFinalOrchTask>(nmeshes, ccIBM);

    auto corrFinalVelBCTask = std::make_shared<VelocityBCEdgesTask>(
        budget.corrFinalVelBC, /*applyToEstimated=*/0, /*doIBEdges=*/1, /*isCorrFinal=*/true);

    auto rteChainTask = makeBarrierChainTask("RTESourceCorr",
        "RTE_SOURCE_CORR",
        [](auto& meshes) {
            fds_rte_source_correction();
        });

    auto corrFinalDumpTask = std::make_shared<CorrFinalDumpTask>(nmeshes);

    subgraph->edges(velCorrKernelTask, corrFinalOrchTask);
    // Fork: VelocityBCEdges (MeshData) || RTE_SOURCE_CORRECTION (BarrierData)
    subgraph->edges(corrFinalOrchTask, corrFinalVelBCTask);
    subgraph->edges(corrFinalOrchTask, rteChainTask);
    // Join+dump
    subgraph->edges(corrFinalVelBCTask, corrFinalDumpTask);
    subgraph->edges(rteChainTask, corrFinalDumpTask);

    subgraph->outputs(corrFinalDumpTask);

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
