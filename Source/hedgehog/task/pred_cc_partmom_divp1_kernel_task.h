#ifndef PRED_CC_PARTMOM_DIVP1_KERNEL_TASK_H
#define PRED_CC_PARTMOM_DIVP1_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// CC_IBM predictor parallel kernel: particle momentum + divergence part 1.
/// Extracted from the CC_IBM "WallDiv+DivExch" barrier's Loop 1.
/// Both operations are per-mesh and thread-safe.
class PredCCPartMomDivP1KernelTask
    : public hh::AbstractTask<1, MeshData<>, MeshData<>> {
public:
    explicit PredCCPartMomDivP1KernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData<>, MeshData<>>(
              "PredCCPartMomDivP1Kernel", numThreads) {}

    void execute(std::shared_ptr<MeshData<>> data) override {
        fds_particle_momentum_kernel(data->nm, data->dt);
        fds_divergence_part_1_kernel(data->nm, data->t, data->dt);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<>, MeshData<>>>
    copy() override {
        return std::make_shared<PredCCPartMomDivP1KernelTask>(
            this->numberThreads());
    }
};

#endif // PRED_CC_PARTMOM_DIVP1_KERNEL_TASK_H
