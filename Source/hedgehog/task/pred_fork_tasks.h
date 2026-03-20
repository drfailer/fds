#ifndef PRED_FORK_TASKS_H
#define PRED_FORK_TASKS_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// DIV_P1 pre-fork task: PHASE=1 (zero DP, PREDICT_NORMAL_VELOCITY).
/// Runs per-mesh in parallel before the fork.
class DivP1PreforkTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit DivP1PreforkTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "DivP1Prefork", numThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_divergence_part_1_prefork(data->nm, data->t, data->dt);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<DivP1PreforkTask>(this->numberThreads());
    }
};

/// DIV_P1 early task: PHASE=2, WORK_BRANCH=2 (species diffusion, heat, thermal).
/// Runs per-mesh in parallel in Branch B after WALL_BC.
class DivP1EarlyTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit DivP1EarlyTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "DivP1Early", numThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_divergence_part_1_early_b(data->nm, data->t, data->dt);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<DivP1EarlyTask>(this->numberThreads());
    }
};

/// DIV_P1 late task: PHASE=3, WORK_BRANCH=2 (advection, RTRM, sources, zone sums).
/// Runs per-mesh in parallel after the join. The C wrapper copies
/// WORK1_B → WORK1 so DIV_P2 can find RTRM in the expected location.
class DivP1LateTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit DivP1LateTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "DivP1Late", numThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_divergence_part_1_late_b(data->nm, data->t, data->dt);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<DivP1LateTask>(this->numberThreads());
    }
};

/// Particle momentum kernel task (mesh-level, predictor fork Branch A).
class PredPartMomKernelTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit PredPartMomKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "PredPartMomKernel", numThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_particle_momentum_kernel(data->nm, data->dt);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<PredPartMomKernelTask>(this->numberThreads());
    }
};

#endif // PRED_FORK_TASKS_H
