#ifndef TIMESTEP_STATE_H
#define TIMESTEP_STATE_H

#include <hedgehog/hedgehog.h>
#include <chrono>
#include <iomanip>
#include <memory>
#include <sstream>
#include <vector>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../fds_fortran_interface.h"

/// Merged dump-join + timestep-loop + INSERT_PARTICLES state.
///
/// Three input paths:
///   1. MeshData<MeshState::Init> — graph initial injection (first iteration only).
///      Collects N tokens, runs INSERT_ALL_PARTICLES, emits MeshData<> to Predictor.
///   2. MeshData<> + BarrierData — dump fork-join (subsequent iterations).
///      Joins N MeshData<> from DumpMeshOutputsTask + 1 BarrierData from DumpGlobalTask.
///      After joining: STOP_CHECK, adjust DT, INSERT_ALL_PARTICLES, emit MeshData<>.
///
/// When skipMeshDump is set on BarrierData, no MeshData<> tokens are expected.
///
/// Termination: emits BarrierData → graph output for clean shutdown.
class TimestepState : public hh::AbstractState<3, MeshData<>, BarrierData, MeshData<MeshState::Init>,
                                                  MeshData<>, BarrierData> {
public:
    TimestepState(int nmeshes, double tEnd, std::shared_ptr<int> icyc)
        : nmeshes_(nmeshes), tEnd_(tEnd), icyc_(std::move(icyc)) {
        initCollected_.reserve(nmeshes);
    }

    /// Collect MeshData<> from DumpMeshOutputsTask (dump fork-join).
    void execute(std::shared_ptr<MeshData<>> /*data*/) override {
        ++meshCount_;
        tryFinalize();
    }

    /// Collect BarrierData from DumpGlobalTask (dump fork-join).
    void execute(std::shared_ptr<BarrierData> data) override {
        globalBarrier_ = data;
        if (data->skipMeshDump) meshCount_ = nmeshes_;
        tryFinalize();
    }

    /// Collect MeshData<Init> from graph input (first iteration).
    void execute(std::shared_ptr<MeshData<MeshState::Init>> data) override {
        initCollected_.push_back(data);
        if (static_cast<int>(initCollected_.size()) == nmeshes_) {
            // INSERT_ALL_PARTICLES moved to PredStep1KernelTask (parallel per-mesh)
            for (auto &md : initCollected_) {
                this->bufferResult(md->template retag<MeshState::Default>());
            }
            this->template flushResults<MeshData<>>();
            initCollected_.clear();
        }
    }

    [[nodiscard]] bool isDone() const { return done_; }

    [[nodiscard]] std::string info() const {
        std::ostringstream oss;
        oss << "INSERT_PART\\nSTOP_CHECK + cycle\\n"
            << std::fixed << std::setprecision(3) << totalTime_ << "s"
            << " / " << invocations_ << " calls";
        if (invocations_ > 0)
            oss << " / avg " << std::setprecision(3)
                << (totalTime_ * 1000.0 / invocations_) << "ms";
        if (skipped_ > 0)
            oss << "\\n" << skipped_ << " skip-dump";
        return oss.str();
    }

private:
    void tryFinalize() {
        if (meshCount_ < nmeshes_ || !globalBarrier_) return;

        auto t0 = std::chrono::steady_clock::now();

        double t = globalBarrier_->meshes[0]->t;
        double dt = globalBarrier_->meshes[0]->dt;

        if (globalBarrier_->skipMeshDump) ++skipped_;

        fds_stop_check(1, t, dt);

        int stopStatus = fds_get_stop_status();
        if (t >= tEnd_ || stopStatus != 0) {
            done_ = true;
            auto bd = std::make_shared<BarrierData>();
            bd->done = true;
            this->addResult(bd);
        } else {
            fds_set_predictor(1);
            fds_set_first_pass(1);
            double newDt = fds_adjust_dt(t, dt);
            ++(*icyc_);
            fds_set_icyc(*icyc_);

            for (auto &md : globalBarrier_->meshes) {
                md->phase = 0;
                md->dt = newDt;
                md->firstPass = true;
                md->dt_bc = 0.0;
                md->call_ht_1d = 0;
            }

            // INSERT_ALL_PARTICLES moved to PredStep1KernelTask (parallel per-mesh)

            for (auto &md : globalBarrier_->meshes) {
                this->bufferResult(md);
            }
            this->template flushResults<MeshData<>>();
        }

        meshCount_ = 0;
        globalBarrier_ = nullptr;

        auto t1 = std::chrono::steady_clock::now();
        totalTime_ += std::chrono::duration<double>(t1 - t0).count();
        ++invocations_;
    }

    int nmeshes_;
    double tEnd_;
    std::shared_ptr<int> icyc_;
    bool done_ = false;
    int meshCount_ = 0;
    std::shared_ptr<BarrierData> globalBarrier_ = nullptr;
    std::vector<std::shared_ptr<MeshData<MeshState::Init>>> initCollected_;
    double totalTime_ = 0.0;
    int invocations_ = 0;
    int skipped_ = 0;
};

/// StateManager for TimestepState.
/// Overrides canTerminate() to break the cycle when the simulation is done.
class TimestepStateManager
    : public hh::StateManager<3, MeshData<>, BarrierData, MeshData<MeshState::Init>,
                               MeshData<>, BarrierData> {
public:
    TimestepStateManager(std::shared_ptr<TimestepState> const &state,
                         std::string const &name)
        : hh::StateManager<3, MeshData<>, BarrierData, MeshData<MeshState::Init>,
                            MeshData<>, BarrierData>(
              state, name) {}

    [[nodiscard]] bool canTerminate() const override {
        this->state()->lock();
        auto ret = std::dynamic_pointer_cast<TimestepState>(
            this->state())->isDone();
        this->state()->unlock();
        return ret;
    }

    [[nodiscard]] std::string extraPrintingInformation() const override {
        this->state()->lock();
        auto ret = std::dynamic_pointer_cast<TimestepState>(
            this->state())->info();
        this->state()->unlock();
        return ret;
    }
};

#endif // TIMESTEP_STATE_H
