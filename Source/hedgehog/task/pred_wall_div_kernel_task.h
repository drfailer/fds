#ifndef PRED_WALL_DIV_KERNEL_TASK_H
#define PRED_WALL_DIV_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/pred_wall_div_data.h"
#include "../fds_fortran_interface.h"

/// Parallel kernel task for predictor wall+div.
/// Calls PARTICLE_MOMENTUM_KERNEL + DIVERGENCE_PART_1_KERNEL per mesh.
class PredWallDivKernelTask : public hh::AbstractTask<1, PredWallDivWork, PredWallDivWork> {
public:
    explicit PredWallDivKernelTask(size_t numThreads)
        : hh::AbstractTask<1, PredWallDivWork, PredWallDivWork>("PredWallDivKernel", numThreads) {}

    void execute(std::shared_ptr<PredWallDivWork> work) override {
        fds_particle_momentum_kernel(work->nm, work->dt);
        fds_divergence_part_1_kernel(work->nm, work->t, work->dt);
        this->addResult(work);
    }

    std::shared_ptr<hh::AbstractTask<1, PredWallDivWork, PredWallDivWork>> copy() override {
        return std::make_shared<PredWallDivKernelTask>(this->numberThreads());
    }
};

#endif // PRED_WALL_DIV_KERNEL_TASK_H
