#ifndef PRED_WALL_DIV_KERNEL_TASK_H
#define PRED_WALL_DIV_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel kernel task for predictor wall+div.
/// Calls PARTICLE_MOMENTUM_KERNEL + DIVERGENCE_PART_1_KERNEL per mesh.
class PredWallDivKernelTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit PredWallDivKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "PredWallDivKernel", numThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_particle_momentum_kernel(data->nm, data->dt);
        fds_divergence_part_1_kernel(data->nm, data->t, data->dt);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<PredWallDivKernelTask>(
            this->numberThreads());
    }
};

#endif // PRED_WALL_DIV_KERNEL_TASK_H
