#ifndef PRESSURE_PARALLEL_TASK_H
#define PRESSURE_PARALLEL_TASK_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Packed parallel "thread pool" task for the pressure iteration pipeline.
///
/// Merges 3 per-mesh parallel kernels into one multi-threaded task:
///   1. Baroclinic:  PredPressure/CorrPressure/Pressure → correction  → PreSolveExch
///   2. Solve:       SolvePhase                         → solve+resid → PostSolveExch
///   3. VelError:    VelErrorPhase                      → vel error   → Pressure
///
/// Between phases, exchange operations (barriers or exchange graph) collect
/// N tokens, perform global exchange, and re-emit with the next phase's tag.
/// No canTerminate() needed: the sequential cycle partners (convergence state,
/// exchange barriers/router) receive TerminationData and terminate first.
///
/// One thread pool serves all 3 phases, saving threads vs. separate tasks.
class PressureParallelTask : public hh::AbstractTask<5,
    MeshData<MeshState::PredictorPressure>,   // initial entry from predictor
    MeshData<MeshState::CorrectorPressure>,   // initial entry from corrector
    MeshData<MeshState::Pressure>,            // cycle-back from convergence
    MeshData<MeshState::SolvePhase>,          // from exchange → solve phase
    MeshData<MeshState::VelErrorPhase>,       // from exchange → vel error phase
    MeshData<MeshState::PreSolveExch>,        // → pre-solve exchange
    MeshData<MeshState::PostSolveExch>,       // → post-solve exchange
    MeshData<MeshState::Pressure>>            // → convergence barrier
{
    using TaskBase = hh::AbstractTask<5,
        MeshData<MeshState::PredictorPressure>,
        MeshData<MeshState::CorrectorPressure>,
        MeshData<MeshState::Pressure>,
        MeshData<MeshState::SolvePhase>,
        MeshData<MeshState::VelErrorPhase>,
        MeshData<MeshState::PreSolveExch>,
        MeshData<MeshState::PostSolveExch>,
        MeshData<MeshState::Pressure>>;

public:
    explicit PressureParallelTask(size_t numThreads, int presFlag)
        : TaskBase("PressureParallel", numThreads),
          presFlag_(presFlag) {}

    /// Phase 1a: Baroclinic from predictor entry
    void execute(std::shared_ptr<MeshData<MeshState::PredictorPressure>> md) override {
        doBaroclinic(retag<MeshState::Pressure>(md));
    }

    /// Phase 1b: Baroclinic from corrector entry
    void execute(std::shared_ptr<MeshData<MeshState::CorrectorPressure>> md) override {
        doBaroclinic(retag<MeshState::Pressure>(md));
    }

    /// Phase 1c: Baroclinic from convergence cycle-back
    void execute(std::shared_ptr<MeshData<MeshState::Pressure>> md) override {
        doBaroclinic(md);
    }

    /// Phase 2: Pressure solve (from exchange)
    void execute(std::shared_ptr<MeshData<MeshState::SolvePhase>> spd) override {
        auto md = retag<MeshState::Pressure>(spd);
        doSolve(md);
        this->addResult(retag<MeshState::PostSolveExch>(md));
    }

    /// Phase 3: Velocity error (from exchange)
    void execute(std::shared_ptr<MeshData<MeshState::VelErrorPhase>> vepd) override {
        auto md = retag<MeshState::Pressure>(vepd);
        fds_compute_velocity_error_kernel(md->nm, md->dt);
        if (fds_is_cc_ibm()) {
            fds_cc_compute_velocity_error(md->dt, md->nm);
        }
        this->addResult(md);
    }

    std::shared_ptr<TaskBase> copy() override {
        return std::make_shared<PressureParallelTask>(
            this->numberThreads(), presFlag_);
    }

private:
    void doBaroclinic(std::shared_ptr<MeshData<MeshState::Pressure>> md) {
        if (fds_pressure_iteration_needs_baroclinic()) {
            fds_baroclinic_correction(md->t, md->nm);
        }
        if (fds_is_cc_ibm()) {
            fds_cc_no_flux(md->dt, md->nm, 1); // FORCE_FLG=TRUE
            fds_cc_exchange_prepare_fn(md->nm);
        }
        md->exchangeRound = 0;  // pre-solve exchange
        this->addResult(retag<MeshState::PreSolveExch>(md));
    }

    void doSolve(std::shared_ptr<MeshData<MeshState::Pressure>> md) {
        if (fds_pressure_iteration_needs_baroclinic() ||
            fds_get_pressure_iterations() == 1) {
            if (fds_is_cc_ibm()) {
                fds_cc_match_velocity_flux(md->nm);
            } else {
                fds_match_velocity_flux_kernel(md->nm);
            }
        }
        fds_no_flux_kernel(md->nm, md->dt);
        if (fds_is_cc_ibm()) {
            fds_cc_no_flux(md->dt, md->nm, 0); // FORCE_FLG=FALSE
        }
        if (fds_get_pressure_iterations() == 1) {
            fds_pressure_iteration_zero_wall_work1(md->nm);
        }
        fds_pressure_solver_compute_rhs_kernel(md->nm, md->t, md->dt);
        if (presFlag_ == ULMAT_PRES_FLAG) {
            fds_ulmat_solver_kernel(md->nm, md->t, md->dt);
            fds_ulmat_check_residuals_kernel(md->nm);
        } else {
            fds_pressure_solver_fft_kernel(md->nm);
            fds_pressure_check_residuals_kernel(md->nm);
        }
        md->exchangeRound = 1;  // post-solve exchange
    }

    static constexpr int ULMAT_PRES_FLAG = 3;
    int presFlag_;
};

/// Pass-through retag task for single-process mode.
/// Converts PreSolveExch/PostSolveExch → Pressure for the MeshExchangeGraph.
/// No canTerminate() needed — predecessor (PressureParallelTask) has it,
/// so this task terminates via default behavior once predecessor terminates.
class ExchangeInputRetagTask : public hh::AbstractTask<2,
    MeshData<MeshState::PreSolveExch>,
    MeshData<MeshState::PostSolveExch>,
    MeshData<MeshState::Pressure>>
{
public:
    ExchangeInputRetagTask()
        : hh::AbstractTask<2,
              MeshData<MeshState::PreSolveExch>,
              MeshData<MeshState::PostSolveExch>,
              MeshData<MeshState::Pressure>>("ExchangeInputRetag", 1) {}

    void execute(std::shared_ptr<MeshData<MeshState::PreSolveExch>> md) override {
        this->addResult(retag<MeshState::Pressure>(md));
    }

    void execute(std::shared_ptr<MeshData<MeshState::PostSolveExch>> md) override {
        this->addResult(retag<MeshState::Pressure>(md));
    }
};

#endif // PRESSURE_PARALLEL_TASK_H
