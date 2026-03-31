#ifndef TIMESTEP_TASKS_H
#define TIMESTEP_TASKS_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../fds_fortran_interface.h"

/// Global computation and global file I/O task.
///
/// Runs in parallel with per-mesh dumps (DumpMeshOutputsTask). Performs:
///   1. SET_DIAGNOSTICS — compute diagnostic quantities
///   2. EXCHANGE_GLOBAL_OUTPUTS — accumulate HRR/mass/device values (MPI)
///   3. UPDATE_CONTROLS — evaluate control logic
///   4. DUMP_GLOBAL_OUTPUTS — write HRR.csv, mass.csv, devc.csv, ctrl.csv
///   5. WRITE_STRINGS — write Smokeview .smv metadata
///   6. WRITE_DIAGNOSTICS — write diagnostic output
///
/// All global files are independent of per-mesh files (.sf, .bf, .prt5, .iso).
class DumpGlobalTask : public hh::AbstractTask<1, BarrierData, BarrierData> {
public:
    explicit DumpGlobalTask(std::shared_ptr<int> icyc)
        : hh::AbstractTask<1, BarrierData, BarrierData>("DumpGlobal", 1),
          icyc_(std::move(icyc)) {}

    void execute(std::shared_ptr<BarrierData> data) override {
        double t = data->t();
        double dt = data->dt();

        // Global computation
        fds_set_diagnostics(*icyc_, t, dt);
        fds_exchange_global_outputs(t, dt);
        fds_update_controls(t, dt);

        // Global file I/O
        fds_dump_global_outputs(t, dt);
        fds_write_strings(t, dt);
        fds_write_diagnostics(t, dt);

        this->addResult(data);
    }

private:
    std::shared_ptr<int> icyc_;
};

/// Per-mesh dump I/O task (thread-safe, multi-threaded).
///
/// Calls fds_dump_mesh_outputs_ts for a single mesh. Each mesh writes to its
/// own independent files (SLCF, BNDF, PRT5, etc.). Uses M => MESHES(NM)
/// internally — no POINT_TO_MESH, fully thread-safe.
///
/// Runs in parallel with DumpGlobalTask (different files).
/// Receives no data when skipMeshDump=true (idle on non-dump timesteps).
class DumpMeshOutputsTask : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit DumpMeshOutputsTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData, MeshData>("DumpMeshOutputs", numThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_dump_mesh_outputs_ts(data->t, data->dt, data->nm);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<DumpMeshOutputsTask>(this->numberThreads());
    }
};

#endif // TIMESTEP_TASKS_H
