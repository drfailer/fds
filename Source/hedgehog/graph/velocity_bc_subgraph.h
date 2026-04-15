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

/// Merged collector + orchestrator for CorrFinal: collects N MeshData<>,
/// runs global work (reset wall counter, CC_END_STEP, exchange(6)),
/// then dual-outputs: N MeshData for VelocityBCEdges + 1 BarrierData for RTE chain.
/// Replaces separate CollectorTask + CorrFinalOrchDualTask.
class CorrFinalOrchTask
    : public hh::AbstractTask<1, MeshData<>, MeshData<>, BarrierData> {
public:
    CorrFinalOrchTask(int nmeshes, bool ccIBM)
        : hh::AbstractTask<1, MeshData<>, MeshData<>, BarrierData>(
              "CorrFinalOrch", 1),
          nmeshes_(nmeshes), ccIBM_(ccIBM),
          nmOffset_(fds_get_lower_mesh_index()) {
        collected_.resize(nmeshes, nullptr);
    }

    void execute(std::shared_ptr<MeshData<>> data) override {
        collected_[data->nm - nmOffset_] = data;
        if (++count_ == nmeshes_) {
            auto t0 = std::chrono::steady_clock::now();
            fds_reset_wall_counter();
            if (ccIBM_) { fds_cc_end_step(collected_[0]->t, collected_[0]->dt, 0); }
            fds_mesh_exchange(6);
            auto t1 = std::chrono::steady_clock::now();
            totalTime_ += std::chrono::duration<double>(t1 - t0).count();
            ++invocations_;

            // Scatter MeshData for VelBCEdges
            this->batchAddResult(collected_);

            // Pass BarrierData for RTE chain
            auto bd = std::make_shared<BarrierData>();
            bd->meshes = collected_;
            this->addResult(bd);

            for (auto &md : collected_) { md = nullptr; }
            count_ = 0;
        }
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
    int nmeshes_, nmOffset_, count_ = 0;
    bool ccIBM_;
    double totalTime_ = 0.0;
    int invocations_ = 0;
    std::vector<std::shared_ptr<MeshData<>>> collected_;
};

/// Merged join + collector + dump for CorrFinal fork:
/// Collects N MeshData<> (from VelocityBCEdges) + 1 BarrierData (from RTESourceCorr),
/// then reduces HRR/MASS, checks dump schedule, and emits:
///   - MeshData<> (N tokens, only when dump needed)
///   - BarrierData (always, for DumpGlobalTask)
/// Replaces separate CorrFinalJoinCollectorTask + CorrFinalDumpTask.
class CorrFinalDumpTask
    : public hh::AbstractTask<2, MeshData<>, BarrierData, MeshData<>, BarrierData> {
public:
    explicit CorrFinalDumpTask(int nmeshes)
        : hh::AbstractTask<2, MeshData<>, BarrierData, MeshData<>, BarrierData>(
              "CorrFinalDump", 1),
          nmeshes_(nmeshes),
          nmOffset_(fds_get_lower_mesh_index()) {
        collected_.resize(nmeshes, nullptr);
    }

    void execute(std::shared_ptr<MeshData<>> data) override {
        collected_[data->nm - nmOffset_] = data;
        if (++count_ == nmeshes_ + 1) fire();
    }

    void execute(std::shared_ptr<BarrierData>) override {
        if (++count_ == nmeshes_ + 1) fire();
    }

    [[nodiscard]] std::string extraPrintingInformation() const override {
        std::ostringstream oss;
        oss << "REDUCE_HRR_MASS\\n"
            << "CHECK_DUMP_SCHEDULE\\n"
            << std::fixed << std::setprecision(3) << totalTime_ << "s"
            << " / " << invocations_ << " calls";
        if (invocations_ > 0)
            oss << " / avg " << std::setprecision(3)
                << (totalTime_ * 1000.0 / invocations_) << "ms";
        return oss.str();
    }

private:
    void fire() {
        auto t0 = std::chrono::steady_clock::now();

        // Reduce per-mesh accumulators into globals
        fds_reduce_hrr_mass(collected_[0]->dt);

        // Check if any mesh needs dump I/O this timestep
        bool anyDump = false;
        for (auto &md : collected_) {
            bool dump = false;
            fds_check_dump_schedule(md->t, md->nm, &dump);
            if (dump) { anyDump = true; break; }
        }

        // Build BarrierData (always emitted for DumpGlobalTask)
        auto bd = std::make_shared<BarrierData>();
        bd->meshes = collected_;  // copy shared_ptrs (TimestepState needs them)

        if (anyDump) {
            bd->skipMeshDump = false;
            this->batchAddResult(collected_);
        } else {
            bd->skipMeshDump = true;
        }

        this->addResult(bd);

        auto t1 = std::chrono::steady_clock::now();
        totalTime_ += std::chrono::duration<double>(t1 - t0).count();
        ++invocations_;

        for (auto &md : collected_) { md = nullptr; }
        count_ = 0;
    }

    int nmeshes_, nmOffset_, count_ = 0;
    double totalTime_ = 0.0;
    int invocations_ = 0;
    std::vector<std::shared_ptr<MeshData<>>> collected_;
};

#endif // VELOCITY_BC_SUBGRAPH_H
