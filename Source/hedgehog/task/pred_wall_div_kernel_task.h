#ifndef PRED_WALL_DIV_KERNEL_TASK_H
#define PRED_WALL_DIV_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Merged predictor fork Branch B task: WALL_BC + DIV_P1_early in one kernel.
/// Eliminates the PredFork-BranchB sub-graph overhead (was WallBCKernel -> DivP1EarlyKernel).
///
/// Calls (per mesh):
///   1. WALL_BC_PREPROCESSING_KERNEL  — ghost values, near-surface gas vars, HTC
///   2. WALL_BC_PROCESS_CELLS_KERNEL  — ~90% of wall cells (no cross-mesh deps)
///   3. WALL_BC_FINALIZE              — interpolated BC, back mesh, particle off-gassing
///   4. DIVERGENCE_PART_1_EARLY_B     — species diffusion, heat, thermal (PHASE=2, WORK_BRANCH=2)
///
/// Reads dt_bc and call_ht_1d from MeshData<> (predictor defaults: 0).
class PredWallBCDivEarlyTask
    : public hh::AbstractTask<1, MeshData<>, MeshData<>> {
public:
    explicit PredWallBCDivEarlyTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData<>, MeshData<>>(
              "PredWallBCDivEarlyKernel", numThreads) {}

    void execute(std::shared_ptr<MeshData<>> data) override {
        // WallBC: preprocessing + cell processing + finalize
        fds_wall_bc_preprocessing_kernel(
            data->nm, data->t, data->dt_bc, data->call_ht_1d);
        fds_wall_bc_process_cells_kernel(
            data->nm, data->t, data->dt, data->dt_bc, data->call_ht_1d);
        fds_wall_bc_finalize(data->nm, data->t, data->dt_bc, data->call_ht_1d);
        // DivP1 early: species diffusion, heat, thermal
        fds_divergence_part_1_early_b(data->nm, data->t, data->dt);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<>, MeshData<>>>
    copy() override {
        return std::make_shared<PredWallBCDivEarlyTask>(
            this->numberThreads());
    }
};

#endif // PRED_WALL_DIV_KERNEL_TASK_H
