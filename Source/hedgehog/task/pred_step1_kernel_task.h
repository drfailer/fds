#ifndef PRED_STEP1_KERNEL_TASK_H
#define PRED_STEP1_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel kernel task for predictor step 1.
/// Calls COMPUTE_VISCOSITY_KERNEL + MASS_FINITE_DIFFERENCES_NEW_KERNEL +
/// DENSITY_KERNEL per mesh (matches corrector's CorrStep1KernelTask pattern).
class PredStep1KernelTask
    : public hh::AbstractTask<1, MeshData<>, MeshData<>> {
public:
    explicit PredStep1KernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData<>, MeshData<>>(
              "PredStep1Kernel", numThreads) {}

    void execute(std::shared_ptr<MeshData<>> data) override {
        fds_insert_particles(data->t, data->nm);
        fds_compute_viscosity_kernel(data->nm, 0);  // estimated=0 for predictor
        fds_mass_finite_differences_kernel(data->nm);
        fds_density_kernel(data->nm, data->t, data->dt);
        fds_cc_density_ts(data->nm, data->t, data->dt);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<>, MeshData<>>>
    copy() override {
        return std::make_shared<PredStep1KernelTask>(
            this->numberThreads());
    }
};

#endif // PRED_STEP1_KERNEL_TASK_H
