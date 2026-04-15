#ifndef DIV_PART2_PREPROCESSING_KERNEL_TASK_H
#define DIV_PART2_PREPROCESSING_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel kernel task for DIVERGENCE_PART_2_PREPROCESSING.
/// Extracted from "DivExch+ZoneOps" barriers (non-CC_IBM only).
///
/// Zone ops (USUM adjustment) are idempotent across meshes: USUM_ADD is computed
/// identically by all meshes (global inputs + uniform-per-zone PBAR), so after the
/// first mesh adjusts USUM, subsequent meshes compute USUM_ADD=0. Concurrent
/// identical writes to USUM are benign (64-bit atomic). D_PBAR_DT is per-mesh.
///
/// NOT used for CC_IBM — GET_LINKED_VELOCITIES has cross-mesh writes.
class DivPart2PreprocessingKernelTask
    : public hh::AbstractTask<1, MeshData<>, MeshData<>> {
public:
    explicit DivPart2PreprocessingKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData<>, MeshData<>>(
              "DivP2PreprocessingKernel", numThreads) {}

    void execute(std::shared_ptr<MeshData<>> data) override {
        fds_divergence_part_2_preprocessing(data->nm, data->dt);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<>, MeshData<>>>
    copy() override {
        return std::make_shared<DivPart2PreprocessingKernelTask>(
            this->numberThreads());
    }
};

#endif // DIV_PART2_PREPROCESSING_KERNEL_TASK_H
