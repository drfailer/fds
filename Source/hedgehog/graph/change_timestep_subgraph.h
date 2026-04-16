#ifndef CHANGE_TIMESTEP_SUBGRAPH_H
#define CHANGE_TIMESTEP_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/barrier_data.h"
#include "../data/change_timestep_data.h"
#include "../data/mesh_data.h"
#include "../data/termination_data.h"
#include "../task/change_timestep_tasks.h"
#include "../task/divergence_part2_kernel_task.h"
#include "../task/velocity_predictor_kernel_task.h"
#include "../state/barrier_state.h"
#include "../state/change_timestep_state.h"
#include "../tool/thread_budget.h"

/// Build the time step retry sub-graph.
///
/// Architecture: parallel DivP2 and VelPred via barriers + kernel tasks.
///
///   RetryPreKernel → RetryMomDivKernel → RetryDivExchSM → RetryDivP2Kernel
///     → RetryPressureSM → RetryVelPredKernel → RetryCheckSM
///                  ↘ (bypass when done=true) ↗
///
/// RetryDivExchSM: barrier (exchange_divergence_info + div2 preprocessing)
/// RetryDivP2Kernel: parallel block kernel per mesh
/// RetryPressureSM: barrier (pressure_iteration + init_change_time_step)
/// RetryVelPredKernel: parallel velocity predictor kernel per mesh
/// RetryCheckSM: stop_check + retry check, on exit CC_END_STEP + MESH_EXCHANGE(3)
inline auto buildChangeTimeStepSubgraph(int nmeshes, const ThreadBudget &budget, bool ccIBM) {
    using SubGraphType = hh::Graph<2, MeshData<>, TerminationData, MeshData<>>;
    auto subgraph = std::make_shared<SubGraphType>("ChangeTimeStepSubgraph");

    auto retryPreKernel = std::make_shared<RetryPreKernelTask>(nmeshes);
    auto retryMomDivKernel = std::make_shared<RetryMomentumDivKernelTask>(
        budget.retryMomDiv);

    // --- Post-kernel: barriers + parallel kernels (was sequential in old RetryLoopState) ---

    auto retryDivExchSM = makeBarrierSM(nmeshes, "RetryDivExch",
        "EXCH_DIV_INFO\\nDIV2_PREPROC",
        [](auto& meshes) {
            fds_exchange_divergence_info();
            for (auto &md : meshes) {
                fds_divergence_part_2_preprocessing(md->nm, md->dt);
            }
        });

    auto retryDivP2Kernel = std::make_shared<DivergencePart2KernelTask<>>(
        budget.standalone(1));

    auto retryPressureSM = makeBarrierSM(nmeshes, "RetryPressure",
        "PRESSURE_ITERATION\\nINIT_CHANGE_TIME_STEP",
        [](auto& meshes) {
            fds_pressure_iteration(meshes[0]->t, meshes[0]->dt);
            fds_init_change_time_step(meshes[0]->dt);
        });

    auto retryVelPredKernel = std::make_shared<VelocityPredictorKernelTask<>>(
        budget.standalone(1));

    auto retryCheckSM = std::make_shared<RetryCheckStateManager>(
        std::make_shared<RetryCheckState>(nmeshes, ccIBM), "RetryCheck");

    // --- Wire the sub-graph ---

    // Entry point: MeshData<> from VelocityPredictor → RetryPreKernel (collects N)
    subgraph->input<MeshData<>>(retryPreKernel);
    subgraph->input<TerminationData>(retryCheckSM);

    // RetryPreKernel → kernel: MeshData<> scatter for parallel processing
    subgraph->template edge<MeshData<>>(retryPreKernel, retryMomDivKernel);

    // RetryPreKernel → check: RetrySequenceData ONLY (bypass path)
    subgraph->template edge<RetrySequenceData>(retryPreKernel, retryCheckSM);

    // Kernel → barriers → parallel kernels → check
    subgraph->edges(retryMomDivKernel, retryDivExchSM);
    subgraph->edges(retryDivExchSM, retryDivP2Kernel);
    subgraph->edges(retryDivP2Kernel, retryPressureSM);
    subgraph->edges(retryPressureSM, retryVelPredKernel);
    subgraph->edges(retryVelPredKernel, retryCheckSM);

    // Cycle: check → pre-kernel: RetrySequenceData ONLY
    subgraph->template edge<RetrySequenceData>(retryCheckSM, retryPreKernel);

    // Exit: RetryCheckState emits MeshData<> → subgraph output
    subgraph->outputs(retryCheckSM);

    return subgraph;
}

#endif // CHANGE_TIMESTEP_SUBGRAPH_H
