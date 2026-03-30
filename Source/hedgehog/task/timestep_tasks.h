#ifndef TIMESTEP_TASKS_H
#define TIMESTEP_TASKS_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../fds_fortran_interface.h"

/// Global pre-dump operations for the timestep.
///
/// Receives BarrierData (all meshes collected after corrector), performs
/// global operations that must run before per-mesh dump I/O, then emits
/// individual MeshData tokens for parallel dump processing.
///
/// Operations: SET_DIAGNOSTICS, EXCHANGE_GLOBAL_OUTPUTS, UPDATE_CONTROLS
class TimestepGlobalTask : public hh::AbstractTask<1, BarrierData, MeshData> {
public:
    explicit TimestepGlobalTask(std::shared_ptr<int> icyc)
        : hh::AbstractTask<1, BarrierData, MeshData>("TimestepGlobal", 1),
          icyc_(std::move(icyc)) {}

    void execute(std::shared_ptr<BarrierData> data) override {
        double t = data->t();
        double dt = data->dt();

        fds_set_diagnostics(*icyc_, t, dt);
        fds_exchange_global_outputs(t, dt);
        fds_update_controls(t, dt);

        for (auto &md : data->meshes) {
            this->addResult(md);
        }
    }

private:
    std::shared_ptr<int> icyc_;
};

/// Per-mesh dump I/O task.
///
/// Calls DUMP_MESH_OUTPUTS for a single mesh. Each mesh writes to its own
/// independent files (SLCF, BNDF, PRT5, etc.).
///
/// Uses 1 thread because DUMP_MESH_OUTPUTS internally calls POINT_TO_MESH
/// which sets global module-level pointers (not thread-safe). Future work:
/// create a kernel version that uses M => MESHES(NM) directly to enable
/// multi-threaded dump I/O.
class DumpMeshOutputsTask : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    DumpMeshOutputsTask()
        : hh::AbstractTask<1, MeshData, MeshData>("DumpMeshOutputs", 1) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_dump_mesh_outputs_ts(data->t, data->dt, data->nm);
        this->addResult(data);
    }
};

#endif // TIMESTEP_TASKS_H
