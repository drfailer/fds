#ifndef QR_ADD_COPY_KERNEL_TASK_H
#define QR_ADD_COPY_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel kernel task for QR addition + WORK1 copy.
/// Extracted from groupC barrier (corrector fork2 join) to run per-mesh in parallel.
///
/// Calls divergence_part_1_add_qr_b (adds QR contribution to divergence)
/// and copy_work1_b_to_work1 (moves branch-2 RTRM to WORK1 for DivP2).
class QRAddCopyKernelTask
    : public hh::AbstractTask<1, MeshData<>, MeshData<>> {
public:
    explicit QRAddCopyKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData<>, MeshData<>>(
              "QRAddCopyKernel", numThreads) {}

    void execute(std::shared_ptr<MeshData<>> data) override {
        fds_divergence_part_1_add_qr_b(data->nm);
        fds_copy_work1_b_to_work1(data->nm);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<>, MeshData<>>>
    copy() override {
        return std::make_shared<QRAddCopyKernelTask>(this->numberThreads());
    }
};

#endif // QR_ADD_COPY_KERNEL_TASK_H
