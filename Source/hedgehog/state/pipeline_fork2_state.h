#ifndef PIPELINE_FORK2_STATE_H
#define PIPELINE_FORK2_STATE_H

#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../data/pipeline_fork2_data.h"
#include "../fds_fortran_interface.h"

/// Fork state for Corrector Fork 2: RADIATION || DIV_P1.
/// Collects N MeshData tokens, runs InitDivIntegrals (zero DSUM/PSUM/USUM),
/// then dispatches N Fork2RadWork + N Fork2DivP1Work for concurrent execution.
class PipelineFork2State
    : public hh::AbstractState<1, MeshData, Fork2RadWork, Fork2DivP1Work> {
public:
    explicit PipelineFork2State(int nmeshes)
        : nmeshes_(nmeshes) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);

        if (static_cast<int>(collected_.size()) == nmeshes_) {
            // Zero divergence integrals (DSUM/PSUM/USUM) before DIV_P1
            fds_initialize_divergence_integrals();

            // Dispatch both branches for each mesh
            for (auto &md : collected_) {
                this->addResult(std::make_shared<Fork2RadWork>(md));
                this->addResult(std::make_shared<Fork2DivP1Work>(md));
            }

            collected_.clear();
            collected_.reserve(nmeshes_);
        }
    }

private:
    int nmeshes_;
    std::vector<std::shared_ptr<MeshData>> collected_;
};

/// Join state for Corrector Fork 2.
/// Collects Fork2RadBarrier and Fork2DivP1Barrier (one of each).
/// When both arrive, emits the BarrierData for MeshExchange(2).
class PipelineJoin2State
    : public hh::AbstractState<2, Fork2RadBarrier, Fork2DivP1Barrier, BarrierData> {
public:
    PipelineJoin2State() = default;

    void execute(std::shared_ptr<Fork2RadBarrier> radBarrier) override {
        radBarrier_ = radBarrier->barrierData;
        tryEmit();
    }

    void execute(std::shared_ptr<Fork2DivP1Barrier> divBarrier) override {
        divBarrier_ = divBarrier->barrierData;
        tryEmit();
    }

private:
    void tryEmit() {
        if (radBarrier_ && divBarrier_) {
            // Use the radiation barrier's data (contains mesh list)
            this->addResult(radBarrier_);
            radBarrier_.reset();
            divBarrier_.reset();
        }
    }

    std::shared_ptr<BarrierData> radBarrier_;
    std::shared_ptr<BarrierData> divBarrier_;
};

#endif // PIPELINE_FORK2_STATE_H
