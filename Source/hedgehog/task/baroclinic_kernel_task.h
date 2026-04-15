#ifndef BAROCLINIC_KERNEL_TASK_H
#define BAROCLINIC_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel baroclinic correction kernel task.
///
/// Entry point for the pressure iteration subgraph. Accepts three input types:
///   - MeshData<PredictorPressure>: initial entry from predictor pipeline
///   - MeshData<CorrectorPressure>: initial entry from corrector pipeline
///   - MeshData<Pressure>: cycle-back from convergence state
///
/// All inputs are converted to MeshData<Pressure> for the internal pipeline.
///
/// Multi-threaded: each clone processes one mesh independently.
class BaroclinicKernelTask
    : public hh::AbstractTask<3,
          MeshData<MeshState::PredictorPressure>,
          MeshData<MeshState::CorrectorPressure>,
          MeshData<MeshState::Pressure>,
          MeshData<MeshState::Pressure>> {
public:
    explicit BaroclinicKernelTask(size_t kernelThreads)
        : hh::AbstractTask<3,
              MeshData<MeshState::PredictorPressure>,
              MeshData<MeshState::CorrectorPressure>,
              MeshData<MeshState::Pressure>,
              MeshData<MeshState::Pressure>>(
              "BaroclinicKernel", kernelThreads) {}

    void execute(std::shared_ptr<MeshData<MeshState::PredictorPressure>> md) override {
        doWork(md->retag<MeshState::Pressure>());
    }

    void execute(std::shared_ptr<MeshData<MeshState::CorrectorPressure>> md) override {
        doWork(md->retag<MeshState::Pressure>());
    }

    void execute(std::shared_ptr<MeshData<MeshState::Pressure>> md) override {
        doWork(md);
    }

    std::shared_ptr<hh::AbstractTask<3,
        MeshData<MeshState::PredictorPressure>,
        MeshData<MeshState::CorrectorPressure>,
        MeshData<MeshState::Pressure>,
        MeshData<MeshState::Pressure>>> copy() override {
        return std::make_shared<BaroclinicKernelTask>(this->numberThreads());
    }

private:
    void doWork(std::shared_ptr<MeshData<MeshState::Pressure>> md) {
        if (fds_pressure_iteration_needs_baroclinic()) {
            fds_baroclinic_correction(md->t, md->nm);
        }
        if (fds_is_cc_ibm()) {
            fds_cc_no_flux(md->dt, md->nm, 1); // FORCE_FLG=TRUE
            fds_cc_exchange_prepare_fn(md->nm); // Set FN_OMESH before exchange
        }
        md->exchangeRound = 0;  // pre-solve exchange
        this->addResult(md);
    }
};

#endif // BAROCLINIC_KERNEL_TASK_H
