#ifndef CHANGE_TIMESTEP_TASKS_H
#define CHANGE_TIMESTEP_TASKS_H

#include <hedgehog/hedgehog.h>
#include "../data/change_timestep_data.h"
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../fds_fortran_interface.h"

/// Initial check: determines if retry is needed.
/// Always emits RetrySequenceData into the pipeline. If no retry is needed,
/// sets done=true so all downstream tasks pass it through without processing.
class CheckRetryTask : public hh::AbstractTask<1, BarrierData, RetrySequenceData> {
public:
    CheckRetryTask()
        : hh::AbstractTask<1, BarrierData, RetrySequenceData>("CheckRetry", 1) {}

    void execute(std::shared_ptr<BarrierData> barrier) override {
        fds_stop_check_zero();

        int needRetry = 0;
        double newDt = 0.0;
        fds_check_change_time_step(&needRetry, &newDt);

        if (!needRetry) {
            // No retry needed: mark as done, pass through pipeline to exit
            auto retryData = std::make_shared<RetrySequenceData>(
                barrier->meshes, barrier->t(), barrier->dt(), -1, true);
            this->addResult(retryData);
        } else {
            // Retry needed: start retry sequence
            auto retryData = std::make_shared<RetrySequenceData>(
                barrier->meshes, barrier->t(), newDt, 0, false);
            this->addResult(retryData);
        }
    }
};

/// Restore UVW and compute density for all meshes (per-mesh operations)
class RetryDensityTask : public hh::AbstractTask<1, RetrySequenceData, RetrySequenceData> {
public:
    RetryDensityTask()
        : hh::AbstractTask<1, RetrySequenceData, RetrySequenceData>("RetryDensity", 1) {}

    void execute(std::shared_ptr<RetrySequenceData> data) override {
        if (data->done) { this->addResult(data); return; }

        fds_set_first_pass(0);

        for (auto &md : data->meshes) {
            md->dt = data->dt;
            md->firstPass = false;
        }

        for (auto &md : data->meshes) {
            fds_cc_restore_uvw_unlinked(md->nm);
            fds_density(data->t, data->dt, md->nm);
        }

        this->addResult(data);
    }
};

/// CC_DENSITY global operation
class RetryCCDensityTask : public hh::AbstractTask<1, RetrySequenceData, RetrySequenceData> {
public:
    RetryCCDensityTask()
        : hh::AbstractTask<1, RetrySequenceData, RetrySequenceData>("RetryCCDensity", 1) {}

    void execute(std::shared_ptr<RetrySequenceData> data) override {
        if (data->done) { this->addResult(data); return; }
        fds_cc_density(data->t, data->dt);
        fds_mesh_exchange(1);
        this->addResult(data);
    }
};

/// Velocity flux computation (per-mesh)
class RetryVelocityFluxTask : public hh::AbstractTask<1, RetrySequenceData, RetrySequenceData> {
public:
    RetryVelocityFluxTask()
        : hh::AbstractTask<1, RetrySequenceData, RetrySequenceData>("RetryVelocityFlux", 1) {}

    void execute(std::shared_ptr<RetrySequenceData> data) override {
        if (data->done) { this->addResult(data); return; }
        for (auto &md : data->meshes) {
            fds_set_baroclinic_false(md->nm);
            fds_viscosity_bc(md->nm, 0);
            fds_velocity_flux(data->t, data->dt, md->nm, 0);
        }
        this->addResult(data);
    }
};

/// HVAC calculation (global barrier)
class RetryHvacTask : public hh::AbstractTask<1, RetrySequenceData, RetrySequenceData> {
public:
    RetryHvacTask()
        : hh::AbstractTask<1, RetrySequenceData, RetrySequenceData>("RetryHvac", 1) {}

    void execute(std::shared_ptr<RetrySequenceData> data) override {
        if (data->done) { this->addResult(data); return; }
        fds_hvac_calc(data->t, data->dt, 0);
        this->addResult(data);
    }
};

/// Initialize divergence integrals (global)
class RetryInitDivTask : public hh::AbstractTask<1, RetrySequenceData, RetrySequenceData> {
public:
    RetryInitDivTask()
        : hh::AbstractTask<1, RetrySequenceData, RetrySequenceData>("RetryInitDiv", 1) {}

    void execute(std::shared_ptr<RetrySequenceData> data) override {
        if (data->done) { this->addResult(data); return; }
        fds_initialize_divergence_integrals();
        this->addResult(data);
    }
};

/// Wall BC, particle momentum, and divergence part 1 (per-mesh)
class RetryDivergencePart1Task : public hh::AbstractTask<1, RetrySequenceData, RetrySequenceData> {
public:
    RetryDivergencePart1Task()
        : hh::AbstractTask<1, RetrySequenceData, RetrySequenceData>("RetryDivPart1", 1) {}

    void execute(std::shared_ptr<RetrySequenceData> data) override {
        if (data->done) { this->addResult(data); return; }
        for (auto &md : data->meshes) {
            fds_wall_bc(data->t, data->dt, md->nm);
            fds_particle_momentum(data->dt, md->nm);
            fds_divergence_part_1(data->t, data->dt, md->nm);
        }
        this->addResult(data);
    }
};

/// Exchange divergence info (global)
class RetryDivExchangeTask : public hh::AbstractTask<1, RetrySequenceData, RetrySequenceData> {
public:
    RetryDivExchangeTask()
        : hh::AbstractTask<1, RetrySequenceData, RetrySequenceData>("RetryDivExchange", 1) {}

    void execute(std::shared_ptr<RetrySequenceData> data) override {
        if (data->done) { this->addResult(data); return; }
        fds_exchange_divergence_info();
        this->addResult(data);
    }
};

/// Divergence part 2 (per-mesh)
class RetryDivergencePart2Task : public hh::AbstractTask<1, RetrySequenceData, RetrySequenceData> {
public:
    RetryDivergencePart2Task()
        : hh::AbstractTask<1, RetrySequenceData, RetrySequenceData>("RetryDivPart2", 1) {}

    void execute(std::shared_ptr<RetrySequenceData> data) override {
        if (data->done) { this->addResult(data); return; }
        for (auto &md : data->meshes) {
            fds_divergence_part_2(data->dt, md->nm);
        }
        this->addResult(data);
    }
};

/// Pressure iteration (global)
class RetryPressureTask : public hh::AbstractTask<1, RetrySequenceData, RetrySequenceData> {
public:
    RetryPressureTask()
        : hh::AbstractTask<1, RetrySequenceData, RetrySequenceData>("RetryPressure", 1) {}

    void execute(std::shared_ptr<RetrySequenceData> data) override {
        if (data->done) { this->addResult(data); return; }
        fds_pressure_iteration(data->t, data->dt);
        fds_init_change_time_step(data->dt);
        this->addResult(data);
    }
};

/// Velocity predictor (per-mesh)
class RetryVelocityPredictorTask : public hh::AbstractTask<1, RetrySequenceData, RetrySequenceData> {
public:
    RetryVelocityPredictorTask()
        : hh::AbstractTask<1, RetrySequenceData, RetrySequenceData>("RetryVelocityPredictor", 1) {}

    void execute(std::shared_ptr<RetrySequenceData> data) override {
        if (data->done) { this->addResult(data); return; }
        for (auto &md : data->meshes) {
            fds_velocity_predictor(data->t + data->dt, data->dt, md->nm);
        }

        fds_stop_check_zero();
        this->addResult(data);
    }
};

#endif // CHANGE_TIMESTEP_TASKS_H
