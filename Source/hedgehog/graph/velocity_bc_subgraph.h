#ifndef VELOCITY_BC_SUBGRAPH_H
#define VELOCITY_BC_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <chrono>
#include <iomanip>
#include <memory>
#include <sstream>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../task/velocity_bc_edges_task.h"
#include "../task/barrier_tasks.h"
#include "../state/barrier_state.h"
#include "../state/collector_state.h"
#include "../state/fork_join_state.h"

/// Build the PredFinal sub-graph.
///
/// PhaseTransition converted from state (PredFinalCollector) to task:
///   CollectorState (N MeshData<> → 1 BarrierData) + PhaseTransitionTask (BarrierData → N MeshData<>)
///
///   VelocityBCEdgesTask → CollectorState → PhaseTransitionTask → output
inline auto buildPredFinalSubgraph(int nmeshes, size_t kernelThreads) {
    using SubGraphType = hh::Graph<1, MeshData<>, MeshData<>>;
    auto subgraph = std::make_shared<SubGraphType>("PredFinal");

    auto kernelTask = std::make_shared<VelocityBCEdgesTask>(kernelThreads, /*applyToEstimated=*/1);
    auto collectorTask = std::make_shared<CollectorTask>(nmeshes, "PredFinalCollector");
    auto phaseTransTask = std::make_shared<PhaseTransitionTask>();

    subgraph->inputs(kernelTask);
    subgraph->edges(kernelTask, collectorTask);
    subgraph->edges(collectorTask, phaseTransTask);
    subgraph->outputs(phaseTransTask);

    return subgraph;
}

/// Dual-output orchestrator for CorrFinal: receives BarrierData, runs global
/// work, then scatters MeshData for VelocityBCEdges AND passes BarrierData
/// through for RTE_SOURCE_CORRECTION (no scatter-collect overhead on RTE path).
class CorrFinalOrchDualTask
    : public hh::AbstractTask<1, BarrierData, MeshData<>, BarrierData> {
public:
    CorrFinalOrchDualTask(bool ccIBM)
        : hh::AbstractTask<1, BarrierData, MeshData<>, BarrierData>(
              "CorrFinalOrch", 1),
          ccIBM_(ccIBM) {}

    void execute(std::shared_ptr<BarrierData> data) override {
        auto t0 = std::chrono::steady_clock::now();
        fds_reset_wall_counter();
        if (ccIBM_) { fds_cc_end_step(data->meshes[0]->t, data->meshes[0]->dt, 0); }
        fds_mesh_exchange(6);
        auto t1 = std::chrono::steady_clock::now();
        totalTime_ += std::chrono::duration<double>(t1 - t0).count();
        ++invocations_;
        this->batchAddResult(data->meshes);  // scatter MeshData for VelBCEdges
        this->addResult(data);               // pass BarrierData for RTE chain
    }

    [[nodiscard]] std::string extraPrintingInformation() const override {
        std::ostringstream oss;
        oss << "RESET_WALL\\nCC_END_STEP\\nMESH_EXCHANGE(6)\\n"
            << std::fixed << std::setprecision(3) << totalTime_ << "s"
            << " / " << invocations_ << " calls";
        if (invocations_ > 0)
            oss << " / avg " << std::setprecision(3)
                << (totalTime_ * 1000.0 / invocations_) << "ms";
        return oss.str();
    }

private:
    bool ccIBM_;
    double totalTime_ = 0.0;
    int invocations_ = 0;
};

/// Build the CorrFinal sub-graph.
///
/// States converted to tasks:
///   - CorrFinalOrchestrator → CollectorTask + CorrFinalOrchDualTask
///   - CorrFinalCollector → CollectorTask + CorrFinalDumpTask (reduce+dump schedule)
///
/// Fork: RTE_SOURCE_CORRECTION runs in parallel with VelocityBCEdges.
/// RTE result is only consumed next timestep by DIVERGENCE_PART_1.
/// OrchDual emits MeshData (scatter) for VelBC and BarrierData (chain) for RTE,
/// avoiding unnecessary scatter-collect on the RTE path.
///
///   Collector → OrchDual → MeshData: VelBCEdges ─┐
///                        → BarrierData: RTE ──────┤→ ForkJoin → Collector → Dump
inline auto buildCorrFinalSubgraph(int nmeshes, size_t kernelThreads) {
    using SubGraphType = hh::Graph<1, MeshData<>, MeshData<>, BarrierData>;
    auto subgraph = std::make_shared<SubGraphType>("CorrFinal");

    bool ccIBM = fds_is_cc_ibm() != 0;

    // CorrFinalOrchestrator → CollectorTask + CorrFinalOrchDualTask
    auto orchCollectorTask = std::make_shared<CollectorTask>(nmeshes, "CorrFinalOrchCollector");
    auto orchDualTask = std::make_shared<CorrFinalOrchDualTask>(ccIBM);

    // Branch 1: VelocityBCEdges (per-mesh parallel)
    auto kernelTask = std::make_shared<VelocityBCEdgesTask>(
        kernelThreads, /*applyToEstimated=*/0, /*doIBEdges=*/1, /*isCorrFinal=*/true);

    // Branch 2: RTE_SOURCE_CORRECTION (global, single thread)
    // Receives BarrierData, runs RTE, scatters MeshData for join.
    // Runs in parallel with VelocityBCEdges. Result only needed next timestep.
    auto rteBarrierTask = makeBarrierTask("RTESourceCorr",
        "RTE_SOURCE_CORR",
        [](auto& meshes) {
            fds_rte_source_correction();
        });

    // ForkJoin: 2 branches per mesh → 1 token per mesh
    auto forkJoinTask = std::make_shared<ForkJoinTask>(nmeshes, 2, 1, "CorrFinalForkJoin");

    // CorrFinalCollector → CollectorTask + CorrFinalDumpTask
    auto dumpCollectorTask = std::make_shared<CollectorTask>(nmeshes, "CorrFinalDumpCollector");
    auto dumpTask = std::make_shared<CorrFinalDumpTask>();

    subgraph->inputs(orchCollectorTask);
    subgraph->edges(orchCollectorTask, orchDualTask);
    // Fork: VelocityBCEdges (MeshData) || RTE_SOURCE_CORRECTION (BarrierData)
    subgraph->edges(orchDualTask, kernelTask);
    subgraph->edges(orchDualTask, rteBarrierTask);
    // Join
    subgraph->edges(kernelTask, forkJoinTask);
    subgraph->edges(rteBarrierTask, forkJoinTask);
    subgraph->edges(forkJoinTask, dumpCollectorTask);
    subgraph->edges(dumpCollectorTask, dumpTask);
    subgraph->outputs(dumpTask);

    return subgraph;
}

#endif // VELOCITY_BC_SUBGRAPH_H
