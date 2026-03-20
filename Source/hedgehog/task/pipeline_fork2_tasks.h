#ifndef PIPELINE_FORK2_TASKS_H
#define PIPELINE_FORK2_TASKS_H

#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/pipeline_fork2_data.h"
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../fds_fortran_interface.h"

/// Branch D kernel task: COMBUSTION_BC + DIV_P1 (SKIP_QR, WORK_BRANCH=2).
/// Runs per-mesh in parallel. QR addition happens after the join.
class Fork2DivP1KernelTask
    : public hh::AbstractTask<1, Fork2DivP1Work, Fork2DivP1Work> {
public:
    explicit Fork2DivP1KernelTask(size_t numThreads)
        : hh::AbstractTask<1, Fork2DivP1Work, Fork2DivP1Work>(
              "Fork2DivP1Kernel", numThreads) {}

    void execute(std::shared_ptr<Fork2DivP1Work> work) override {
        auto &data = work->meshData;
        fds_combustion_bc_kernel(data->nm);
        fds_divergence_part_1_kernel_skip_qr_b(data->nm, data->t, data->dt);
        this->addResult(work);
    }

    std::shared_ptr<hh::AbstractTask<1, Fork2DivP1Work, Fork2DivP1Work>>
    copy() override {
        return std::make_shared<Fork2DivP1KernelTask>(this->numberThreads());
    }
};

/// Branch D collector: gathers N Fork2DivP1Work results, emits Fork2DivP1Barrier.
/// Uses indexed placement for deterministic mesh ordering.
class Fork2DivP1CollectorState
    : public hh::AbstractState<1, Fork2DivP1Work, Fork2DivP1Barrier> {
public:
    explicit Fork2DivP1CollectorState(int nmeshes)
        : nmeshes_(nmeshes), nmOffset_(fds_get_lower_mesh_index()) {
        collected_.resize(nmeshes, nullptr);
    }

    void execute(std::shared_ptr<Fork2DivP1Work> work) override {
        int idx = work->meshData->nm - nmOffset_;
        collected_[idx] = work->meshData;
        ++count_;
        if (count_ == nmeshes_) {
            auto bd = std::make_shared<BarrierData>();
            bd->meshes.reserve(nmeshes_);
            for (auto &md : collected_) {
                bd->meshes.push_back(md);
            }
            this->addResult(std::make_shared<Fork2DivP1Barrier>(bd));
            std::fill(collected_.begin(), collected_.end(), nullptr);
            count_ = 0;
        }
    }

private:
    int nmeshes_;
    int nmOffset_;
    int count_ = 0;
    std::vector<std::shared_ptr<MeshData>> collected_;
};

/// QR addition task: adds RTRM*QR to divergence after MeshExchange(2).
/// Uses WORK_BRANCH=2 (reads RTRM from WORK1_B where DIV_P1 stored it).
class DivP1QRAdditionTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit DivP1QRAdditionTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "DivP1QRAddition", numThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_divergence_part_1_add_qr_b(data->nm);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<DivP1QRAdditionTask>(this->numberThreads());
    }
};

#endif // PIPELINE_FORK2_TASKS_H
