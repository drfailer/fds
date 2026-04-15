#ifndef PRED_FORK_TASKS_H
#define PRED_FORK_TASKS_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// DIV_P1 late task: PHASE=3, WORK_BRANCH=2 (advection, RTRM, sources, zone sums).
/// Runs per-mesh in parallel after the join.
class DivP1LateTask
    : public hh::AbstractTask<1, MeshData<>, MeshData<>> {
public:
    explicit DivP1LateTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData<>, MeshData<>>(
              "DivP1LateKernel", numThreads) {}

    void execute(std::shared_ptr<MeshData<>> data) override {
        fds_divergence_part_1_late_b(data->nm, data->t, data->dt);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<>, MeshData<>>>
    copy() override {
        return std::make_shared<DivP1LateTask>(this->numberThreads());
    }
};

/// Merged predictor fork Branch A task: DivSetup + ParticleMomentum in one kernel.
/// Eliminates the PredFork-BranchA sub-graph (was DivSetupKernel -> PredPartMomKernel).
///
/// Calls (per mesh):
///   1. SET_BAROCLINIC_FALSE
///   2. VISCOSITY_BC_KERNEL
///   3. CC_VELOCITY_BC_TS (DO_IBEDGES=FALSE)
///   4. VELOCITY_FLUX_KERNEL
///   5. PARTICLE_MOMENTUM_KERNEL
class PredDivSetupPartMomTask
    : public hh::AbstractTask<1, MeshData<>, MeshData<>> {
public:
    explicit PredDivSetupPartMomTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData<>, MeshData<>>(
              "PredDivSetupPartMomKernel", numThreads) {}

    void execute(std::shared_ptr<MeshData<>> data) override {
        // DivSetup: baroclinic + viscosity BC + CC velocity BC + velocity flux
        fds_set_baroclinic_false(data->nm);
        fds_viscosity_bc_kernel(data->nm, data->phase);
        fds_cc_velocity_bc_ts(data->t, data->nm, data->phase, 0);
        fds_velocity_flux_kernel(data->nm, data->t, data->dt, data->phase);
        // Particle momentum
        fds_particle_momentum_kernel(data->nm, data->dt);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<>, MeshData<>>>
    copy() override {
        return std::make_shared<PredDivSetupPartMomTask>(this->numberThreads());
    }
};

#endif // PRED_FORK_TASKS_H
