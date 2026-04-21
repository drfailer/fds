#ifndef CORRECTOR_SUBGRAPH_H
#define CORRECTOR_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../data/termination_data.h"
#include "../state/barrier_state.h"
#include "../state/fork_join_state.h"
#include "../task/barrier_tasks.h"
#include "../task/mass_fd_kernel_task.h"
#include "../task/div_setup_kernel_task.h"
#include "../task/combustion_kernel_task.h"
#include "../task/corr_divsetup_comb_part_task.h"
#include "../task/corr_div_part1_kernel_task.h"
#include "../task/divergence_part2_kernel_task.h"
#include "../task/corr_final_kernel_task.h"
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
inline auto buildCorrectorSubgraphImpl(int nmeshes, const ThreadBudget &budget) {
    auto subgraph = std::make_shared<hh::Graph<3,
        MeshData<>, MeshData<MeshState::CorrectorPressure>, TerminationData,
        MeshData<>, MeshData<MeshState::CorrectorPressure>, BarrierData>>("Corrector");

    // --- Kernel tasks (threads from budget) ---

    auto corrFinalKernelTask = std::make_shared<CorrFinalKernelTask>(budget.corrFinal);

    bool ccIBM = fds_is_cc_ibm() != 0;
    bool ht3d = fds_is_ht3d() != 0;

    // --- Sub-graphs ---

    auto corrRadiationSubgraph = buildCorrRadiationSubgraph<MeshState::PostWallBC>(
        nmeshes, budget.corrFork2Radiation);

    // --- Barrier states ---

    // MeshExch4: terminable barrier in cycle with CorrDivSetupCombPart
    auto meshExchange4SM = makeTerminableRetaggingBarrier<
        MeshState::PostCorrStep1, MeshState::PostCorrStep1>(
        nmeshes, "MeshExchange(4)",
        "MESH_EXCHANGE(4)",
        [](auto&) {
            fds_mesh_exchange(4);
        });

    constexpr bool useParallelPressure = (PressureTag != MeshState::Default);

    // --- CorrStep1 + DivSetup||Comb + ParticleOps → MeshExch7 → HvacCalc ---
    //
    // CorrDivSetupCombPartTask merges CorrStep1 + Fork1 (DivSetup||Comb via
    // AsyncWorker) + ParticleOps. Phase 1 (CorrStep1) → MeshExch4 → Phase 2:
    //   MeshData<> → HvacCalc barrier (emitted before ParticleOps)
    //   MeshData<PostParticleOps> → MeshExch7 barrier (emitted after ParticleOps)
    // HvacCalc collects N MeshData + 1 BarrierData (from MeshExch7), runs
    // HVAC_CALC, and emits N MeshData only when both sources are done.
    auto corrDivSetupCombPartTask = std::make_shared<CorrDivSetupCombPartTask>(
        budget.corrDivSetupCombPart);

    // MeshExch(7) — collects N MeshData<PostParticleOps>, emits 1 BarrierData.
    auto meshExch7Barrier = makeBarrierCollectToOne<MeshState::PostParticleOps>(
        nmeshes, "MeshExchange(7)",
        "MESH_EXCHANGE(7)",
        [](auto& meshes) {
            fds_mesh_exchange(7);
        });

    // HvacCalc: collects N MeshData + 1 BarrierData (MeshExch7), runs HVAC_CALC,
    // emits N MeshData<PostHvac> when both HVAC and MeshExch7 are complete.
    // Terminable: in cycle with CorrDivSetupCombPart (Phase 2 → HVAC → Phase 3).
    auto hvacBarrier = makeTerminableEagerDualInputBarrier<MeshState::PostHvac>(
        nmeshes, "HvacCalc",
        "HVAC_CALC",
        [](auto& meshes) {
            fds_hvac_calc(meshes[0]->t, meshes[0]->dt, 1);
        });

    // --- Wire the sub-graph ---

    // Subgraph input → CorrDivSetupCombPart (Phase 1: CorrStep1)
    subgraph->inputs(corrDivSetupCombPartTask);

    // CorrDivSetupCombPart ↔ MeshExchange(4) cycle
    subgraph->edges(corrDivSetupCombPartTask, meshExchange4SM);
    subgraph->edges(meshExchange4SM, corrDivSetupCombPartTask);

    // TerminationData breaks CorrDivSetupCombPart ↔ MeshExch4 cycle
    subgraph->template input<TerminationData>(meshExchange4SM);

    // MeshData<> output → HvacCalc barrier (emitted before ParticleOps)
    subgraph->edges(corrDivSetupCombPartTask, hvacBarrier);
    // MeshData<PostParticleOps> output → MeshExch7 barrier (emitted after ParticleOps)
    subgraph->edges(corrDivSetupCombPartTask, meshExch7Barrier);
    // MeshExch7 (BarrierData) → HvacCalc (waits for both before emitting)
    subgraph->edges(meshExch7Barrier, hvacBarrier);
    // HvacCalc (MeshData<PostHvac>) → CorrDivSetupCombPart Phase 3 (WallBC)
    subgraph->edges(hvacBarrier, corrDivSetupCombPartTask);

    // TerminationData breaks CorrDivSetupCombPart ↔ HvacCalc cycle
    subgraph->template input<TerminationData>(hvacBarrier);

    // --- Fork 2: RADIATION || DIV_P1 (or sequential for CC_IBM) ---
    if (ccIBM) {
        // Group B (CC_IBM): Exchange(6) barrier for back wall data (HT3D only).
        // Radiation forked from CorrDivSetupCombPart — doesn't need Exchange(6) data.
        // RESET_WALL_COUNTER moved to CorrFinalOrch barrier
        auto groupBPostSM = makeRetaggingBarrier<MeshState::PostWallBC, MeshState::Default>(
            nmeshes, "MeshExch6a",
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

        // Fork: Radiation starts immediately (no Exchange(6) dependency).
        //        Exchange(6) runs in parallel. Both join at meshExch2SM.
        subgraph->edges(corrDivSetupCombPartTask, corrRadiationSubgraph);
        subgraph->edges(corrDivSetupCombPartTask, groupBPostSM);
        subgraph->edges(corrRadiationSubgraph, meshExch2SM);
        subgraph->edges(groupBPostSM, meshExch2SM);
        subgraph->edges(meshExch2SM, corrDivP1KernelTask);
        subgraph->edges(corrDivP1KernelTask, corrDivExchangeSM);
        // CC_IBM keeps standalone DivPart2 kernel (not packed)
        auto corrDivP2KernelTask = std::make_shared<DivergencePart2KernelTask<PressureTag>>(budget.corrDivPart2);
        subgraph->edges(corrDivExchangeSM, corrDivP2KernelTask);

        // Downstream: DivP2 → Pressure → CorrFinalKernel
        if constexpr (useParallelPressure) {
            subgraph->outputs(corrDivP2KernelTask);
            subgraph->template input<MeshData<MeshState::CorrectorPressure>>(corrFinalKernelTask);
        } else {
            auto corrPressureSM = makeRetaggingBarrier<MeshState::Default, MeshState::CorrectorPressure>(
                nmeshes, "CorrPressure",
                "PRESSURE_ITERATION",
                [](auto& meshes) {
                    fds_pressure_iteration(meshes[0]->t, meshes[0]->dt);
                });
            subgraph->edges(corrDivP2KernelTask, corrPressureSM);
            subgraph->edges(corrPressureSM, corrFinalKernelTask);
        }
    } else {
        // Group B (non-CC_IBM): Exchange(6) [HT3D only] + InitDiv barrier.
        // Radiation forked from CorrDivSetupCombPart — doesn't need Exchange(6) data.
        // DivP1 waits for this barrier (needs InitDiv zeroed arrays).
        // Exchange(6) conditional: at this point velocities haven't been corrected yet,
        // so only new data since Exchange(3) is back wall info needed for HT3D.
        // RESET_WALL_COUNTER moved to CorrFinalOrch barrier
        auto groupBPostSM = makeRetaggingBarrier<MeshState::PostWallBC, MeshState::Default>(
            nmeshes, "MeshExch6a+InitDiv",
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

        auto corrDivExchBarrier = makeTerminableRetaggingBarrier<MeshState::DivExch, MeshState::DivP2Pre>(
            nmeshes, "DivExchange",
            "EXCH_DIV_INFO",
            [](auto&) {
                fds_exchange_divergence_info();
            });

        auto corrGlobalMatBarrier = makeTerminableRetaggingBarrier<MeshState::GlobalMat, MeshState::DivPart2>(
            nmeshes, "GlobalMatrix+PressureInit",
            "GLOBAL_MATRIX_REASSIGN\\nPRES_INIT+INCR",
            [useParallelPressure](auto&) {
                fds_global_matrix_reassign(0);
                if (useParallelPressure) {
                    fds_pressure_iteration_init();
                    fds_pressure_iteration_increment();
                }
            });

        // Fork: Radiation starts immediately (no Exchange(6) dependency).
        //        DivP1 waits for Exchange(6) barrier.
        subgraph->edges(corrDivSetupCombPartTask, corrRadiationSubgraph);
        subgraph->edges(corrDivSetupCombPartTask, groupBPostSM);
        subgraph->edges(groupBPostSM, fork2DivP1Task);

        // Group C: Join2 → Parallel ↔ DivExchange ↔ Parallel ↔ GlobalMatrix ↔ Parallel
        subgraph->edges(corrRadiationSubgraph, groupCPreSM);
        subgraph->edges(fork2DivP1Task, groupCPreSM);
        subgraph->edges(groupCPreSM, corrDivParallelTask);
        subgraph->edges(corrDivParallelTask, corrDivExchBarrier);
        subgraph->edges(corrDivExchBarrier, corrDivParallelTask);
        subgraph->edges(corrDivParallelTask, corrGlobalMatBarrier);
        subgraph->edges(corrGlobalMatBarrier, corrDivParallelTask);

        // TerminationData breaks structural cycles at shutdown
        subgraph->template input<TerminationData>(corrDivExchBarrier);
        subgraph->template input<TerminationData>(corrGlobalMatBarrier);

        // Downstream: final output (MeshData<PressureTag>) → Pressure → CorrFinalKernel
        if constexpr (useParallelPressure) {
            subgraph->template output<MeshData<PressureTag>>(corrDivParallelTask);
            subgraph->template input<MeshData<MeshState::CorrectorPressure>>(corrFinalKernelTask);
        } else {
            auto corrPressureSM = makeRetaggingBarrier<MeshState::Default, MeshState::CorrectorPressure>(
                nmeshes, "CorrPressure",
                "PRESSURE_ITERATION",
                [](auto& meshes) {
                    fds_pressure_iteration(meshes[0]->t, meshes[0]->dt);
                });
            subgraph->edges(corrDivParallelTask, corrPressureSM);
            subgraph->edges(corrPressureSM, corrFinalKernelTask);
        }
    }

    // --- CorrFinal: CorrFinalKernel ↔ CorrFinalOrch → CorrFinalDump ---
    //
    // CorrFinalKernelTask has 3 phases via different input types:
    //   Phase 1: MeshData<CorrectorPressure> → VelCorr → MeshData<PostVelCorr> → Orch
    //   Phase 2: MeshData<> (from Orch) → VelBCEdges → MeshData<> → Dump
    //   Phase 3: BarrierData (from Orch) → RTE → BarrierData → Dump

    auto corrFinalOrchTask = std::make_shared<CorrFinalOrchTask>(nmeshes, ccIBM);
    auto corrFinalDumpTask = std::make_shared<CorrFinalDumpTask>(nmeshes);

    // CorrFinalKernel ↔ CorrFinalOrch cycle
    subgraph->edges(corrFinalKernelTask, corrFinalOrchTask);  // MeshData<PostVelCorr>
    subgraph->edges(corrFinalOrchTask, corrFinalKernelTask);  // MeshData<> + BarrierData
    // CorrFinalKernel → CorrFinalDump
    subgraph->edges(corrFinalKernelTask, corrFinalDumpTask);  // MeshData<> + BarrierData

    // TerminationData breaks CorrFinalKernel ↔ CorrFinalOrch cycle
    subgraph->template input<TerminationData>(corrFinalOrchTask);

    subgraph->outputs(corrFinalDumpTask);

    return subgraph;
}

/// Dispatch wrapper: selects the correct template instantiation at runtime.
inline auto buildCorrectorSubgraph(int nmeshes, const ThreadBudget &budget) {
    if (fds_use_pressure_subgraph()) {
        return buildCorrectorSubgraphImpl<MeshState::CorrectorPressure>(
            nmeshes, budget);
    }
    return buildCorrectorSubgraphImpl<MeshState::Default>(
        nmeshes, budget);
}

#endif // CORRECTOR_SUBGRAPH_H
