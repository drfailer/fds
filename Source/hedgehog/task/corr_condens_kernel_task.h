#ifndef CORR_CONDENS_KERNEL_TASK_H
#define CORR_CONDENS_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel kernel task for corrector condensation + particle mass/energy.
/// Calls CONDENSATION_EVAPORATION_KERNEL + PARTICLE_MASS_ENERGY_KERNEL per mesh.
/// These are consecutive mesh-independent operations merged to eliminate queue overhead.
class CorrCondensKernelTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit CorrCondensKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "CorrCondensPartMEKernel", numThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_condensation_kernel(data->nm, data->dt);
        fds_particle_mass_energy_kernel(data->nm, data->t, data->dt);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<CorrCondensKernelTask>(
            this->numberThreads());
    }
};

#endif // CORR_CONDENS_KERNEL_TASK_H
