#ifndef PREDICTOR_TASKS_H
#define PREDICTOR_TASKS_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Task 1: Insert particles, compute viscosity, mass finite differences
class PredStep1Task : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit PredStep1Task(size_t nThreads = 1)
        : hh::AbstractTask<1, MeshData, MeshData>("PredStep1", nThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_insert_particles(data->t, data->nm);
        fds_compute_viscosity(data->nm, 0);
        fds_mass_finite_differences(data->nm);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>> copy() override {
        return std::make_shared<PredStep1Task>(this->numberThreads());
    }
};

/// Task 2: Density prediction
class DensityPredTask : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit DensityPredTask(size_t nThreads = 1)
        : hh::AbstractTask<1, MeshData, MeshData>("DensityPred", nThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_density(data->t, data->dt, data->nm);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>> copy() override {
        return std::make_shared<DensityPredTask>(this->numberThreads());
    }
};

/// Task 3: Viscosity BC + velocity flux setup for divergence
class PredDivSetupTask : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit PredDivSetupTask(size_t nThreads = 1)
        : hh::AbstractTask<1, MeshData, MeshData>("PredDivSetup", nThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_set_baroclinic_false(data->nm);
        fds_viscosity_bc(data->nm, 0);
        fds_velocity_flux(data->t, data->dt, data->nm, 0);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>> copy() override {
        return std::make_shared<PredDivSetupTask>(this->numberThreads());
    }
};

/// Task 4: Wall BC + particle momentum + divergence part 1
class PredWallDivTask : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit PredWallDivTask(size_t nThreads = 1)
        : hh::AbstractTask<1, MeshData, MeshData>("PredWallDiv", nThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_wall_bc(data->t, data->dt, data->nm);
        fds_particle_momentum(data->dt, data->nm);
        fds_divergence_part_1(data->t, data->dt, data->nm);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>> copy() override {
        return std::make_shared<PredWallDivTask>(this->numberThreads());
    }
};

/// Task 5: Divergence part 2
class DivPart2PredTask : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit DivPart2PredTask(size_t nThreads = 1)
        : hh::AbstractTask<1, MeshData, MeshData>("DivPart2Pred", nThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_divergence_part_2(data->dt, data->nm);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>> copy() override {
        return std::make_shared<DivPart2PredTask>(this->numberThreads());
    }
};

/// Task 6: Velocity predictor
class VelPredictorTask : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit VelPredictorTask(size_t nThreads = 1)
        : hh::AbstractTask<1, MeshData, MeshData>("VelPredictor", nThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        // Initialize CHANGE_TIME_STEP_INDEX and DT_NEW before velocity prediction
        fds_init_change_time_step(data->dt);
        fds_velocity_predictor(data->t + data->dt, data->dt, data->nm);
        // Check for numerical instability
        fds_stop_check_zero();
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>> copy() override {
        return std::make_shared<VelPredictorTask>(this->numberThreads());
    }
};

/// Task 7: Match velocity + velocity BC (end of predictor)
class PredFinalTask : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit PredFinalTask(size_t nThreads = 1)
        : hh::AbstractTask<1, MeshData, MeshData>("PredFinal", nThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_match_velocity(data->nm);
        fds_synthetic_turbulence(data->dt, data->t, data->nm);
        fds_velocity_bc(data->t, data->nm, 1); // estimated=true for predictor end
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>> copy() override {
        return std::make_shared<PredFinalTask>(this->numberThreads());
    }
};

#endif // PREDICTOR_TASKS_H
