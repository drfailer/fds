#ifndef PIPELINE_FORK1_TASKS_H
#define PIPELINE_FORK1_TASKS_H

#include <hedgehog/hedgehog.h>
#include "../data/pipeline_fork1_data.h"
#include "../fds_fortran_interface.h"

/// Branch A task: VELOCITY_FLUX (mesh-level fallback when block decomposition unavailable).
/// Accepts Fork1VFluxWork, runs VFLUX kernel, emits Fork1VFluxResult.
class Fork1VFluxKernelTask
    : public hh::AbstractTask<1, Fork1VFluxWork, Fork1VFluxResult> {
public:
    explicit Fork1VFluxKernelTask(size_t numThreads)
        : hh::AbstractTask<1, Fork1VFluxWork, Fork1VFluxResult>(
              "Fork1VFluxKernel", numThreads) {}

    void execute(std::shared_ptr<Fork1VFluxWork> work) override {
        auto data = work->meshData;
        fds_set_baroclinic_false(data->nm);
        fds_viscosity_bc_kernel(data->nm, data->phase);
        fds_velocity_flux_kernel(data->nm, data->t, data->dt, data->phase);
        if (data->phase)
            fds_agglomeration(data->dt, data->nm);
        this->addResult(std::make_shared<Fork1VFluxResult>(data));
    }

    std::shared_ptr<hh::AbstractTask<1, Fork1VFluxWork, Fork1VFluxResult>>
    copy() override {
        return std::make_shared<Fork1VFluxKernelTask>(this->numberThreads());
    }
};

/// Branch B task: COMBUSTION.
/// Accepts Fork1CombWork, runs combustion kernel, emits Fork1CombResult.
class Fork1CombKernelTask
    : public hh::AbstractTask<1, Fork1CombWork, Fork1CombResult> {
public:
    explicit Fork1CombKernelTask(size_t numThreads)
        : hh::AbstractTask<1, Fork1CombWork, Fork1CombResult>(
              "Fork1CombKernel", numThreads) {}

    void execute(std::shared_ptr<Fork1CombWork> work) override {
        auto data = work->meshData;
        fds_combustion_kernel(data->nm, data->t, data->dt);
        this->addResult(std::make_shared<Fork1CombResult>(data));
    }

    std::shared_ptr<hh::AbstractTask<1, Fork1CombWork, Fork1CombResult>>
    copy() override {
        return std::make_shared<Fork1CombKernelTask>(this->numberThreads());
    }
};

#endif // PIPELINE_FORK1_TASKS_H
