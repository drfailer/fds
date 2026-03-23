#ifndef PIPELINE_FORK2_TASKS_H
#define PIPELINE_FORK2_TASKS_H

#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../fds_fortran_interface.h"

/// Branch D kernel task: COMBUSTION_BC + DIV_P1 (SKIP_QR, WORK_BRANCH=2).
/// Runs per-mesh in parallel. QR addition happens after the join.
class Fork2DivP1KernelTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit Fork2DivP1KernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "Fork2DivP1Kernel", numThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_combustion_bc_kernel(data->nm);
        fds_divergence_part_1_kernel_skip_qr_b(data->nm, data->t, data->dt);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<Fork2DivP1KernelTask>(this->numberThreads());
    }
};

/// Branch D collector: gathers N MeshData results, emits BarrierData.
class Fork2DivP1CollectorState
    : public hh::AbstractState<1, MeshData, BarrierData> {
public:
    explicit Fork2DivP1CollectorState(int nmeshes)
        : nmeshes_(nmeshes), nmOffset_(fds_get_lower_mesh_index()) {
        collected_.resize(nmeshes, nullptr);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        int idx = data->nm - nmOffset_;
        collected_[idx] = data;
        ++count_;
        if (count_ == nmeshes_) {
            auto bd = std::make_shared<BarrierData>();
            bd->meshes.reserve(nmeshes_);
            for (auto &md : collected_) {
                bd->meshes.push_back(md);
            }
            this->addResult(bd);
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
