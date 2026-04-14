#ifndef WALLBC_KERNEL_TASK_H
#define WALLBC_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel task that calls WALL_BC_PREPROCESSING_KERNEL + WALL_BC_PROCESS_CELLS_KERNEL +
/// WALL_BC_FINALIZE. All three phases are thread-safe per-mesh operations:
/// - Preprocessing: ghost value assignment, near-surface gas variables, heat transfer coeff
/// - Cell processing: ~90% of wall cells (no cross-mesh dependencies)
/// - Finalize: remaining ~10% (INTERPOLATED_BC, HAS_BACK_MESH, particle off-gassing)
///   All finalize callees write only to local mesh (OMESH reads are copies, back-mesh reads are read-only).
///
/// Reads dt_bc and call_ht_1d from MeshData<> (set by upstream barrier).
class WallBCKernelTask
    : public hh::AbstractTask<1, MeshData<>, MeshData<>> {
public:
    explicit WallBCKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData<>, MeshData<>>(
              "WallBCKernel", numThreads) {}

    void execute(std::shared_ptr<MeshData<>> data) override {
        // Thread-safe preprocessing: ASSIGN_GHOST_VALUE_KERNEL + NEAR_SURFACE_GAS_VARIABLES + HTC
        fds_wall_bc_preprocessing_kernel(
            data->nm, data->t, data->dt_bc, data->call_ht_1d);
        // Process ~90% of wall cells (no cross-mesh dependencies)
        fds_wall_bc_process_cells_kernel(
            data->nm, data->t, data->dt, data->dt_bc, data->call_ht_1d);
        // Finalize: INTERPOLATED_BC, HAS_BACK_MESH, thin walls, particle off-gassing
        fds_wall_bc_finalize(data->nm, data->t, data->dt_bc, data->call_ht_1d);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<>, MeshData<>>>
    copy() override {
        return std::make_shared<WallBCKernelTask>(
            this->numberThreads());
    }
};

#endif // WALLBC_KERNEL_TASK_H
