#ifndef DIV_SETUP_KERNEL_TASK_H
#define DIV_SETUP_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel task that calls the thread-safe velocity flux kernel.
/// Each thread processes one mesh independently.
/// Uses data->phase to select predictor (0) or corrector (1) arrays.
class DivSetupKernelTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit DivSetupKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "DivSetupKernel", numThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_set_baroclinic_false(data->nm);
        fds_viscosity_bc_kernel(data->nm, data->phase);
        fds_velocity_flux_kernel(data->nm, data->t, data->dt,
                                 data->phase);
        if (data->phase)
            fds_agglomeration(data->dt, data->nm);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<DivSetupKernelTask>(
            this->numberThreads());
    }
};

#endif // DIV_SETUP_KERNEL_TASK_H
