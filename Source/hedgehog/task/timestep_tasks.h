#ifndef TIMESTEP_TASKS_H
#define TIMESTEP_TASKS_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../fds_fortran_interface.h"

/// Merged dump task: global file I/O (BarrierData) + per-mesh file I/O (MeshData<>).
///
/// Two execute methods handle independent work in parallel:
///   BarrierData path (1 per timestep):
///     SET_DIAGNOSTICS, EXCHANGE_GLOBAL_OUTPUTS, UPDATE_CONTROLS,
///     DUMP_GLOBAL_OUTPUTS, WRITE_STRINGS, WRITE_DIAGNOSTICS
///   MeshData path (N per timestep, skipped on non-dump steps):
///     DUMP_MESH_OUTPUTS_TS per mesh
///
/// Global and per-mesh files are independent — safe to run concurrently.
class DumpTask : public hh::AbstractTask<2, BarrierData, MeshData<>, BarrierData, MeshData<>> {
public:
    DumpTask(size_t numThreads, std::shared_ptr<int> icyc)
        : hh::AbstractTask<2, BarrierData, MeshData<>, BarrierData, MeshData<>>(
              "Dump", numThreads),
          icyc_(std::move(icyc)) {}

    void execute(std::shared_ptr<BarrierData> data) override {
        double t = data->t();
        double dt = data->dt();
        fds_set_diagnostics(*icyc_, t, dt);
        fds_exchange_global_outputs(t, dt);
        fds_update_controls(t, dt);
        fds_dump_global_outputs(t, dt);
        fds_write_strings(t, dt);
        fds_write_diagnostics(t, dt);
        this->addResult(data);
    }

    void execute(std::shared_ptr<MeshData<>> data) override {
        fds_dump_mesh_outputs_ts(data->t, data->dt, data->nm);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<2, BarrierData, MeshData<>, BarrierData, MeshData<>>>
    copy() override {
        return std::make_shared<DumpTask>(this->numberThreads(), icyc_);
    }

private:
    std::shared_ptr<int> icyc_;
};

#endif // TIMESTEP_TASKS_H
