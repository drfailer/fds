#ifndef DIVERGENCE_PART2_KERNEL_TASK_H
#define DIVERGENCE_PART2_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel task that calls the thread-safe divergence part 2 block kernel.
/// Each thread processes one mesh independently.
///
/// Template parameter OutS controls the output MeshData state tag.
/// When OutS != Default, the task retags its output for type-based routing
/// (e.g. PredictorPressure to route to the shared pressure subgraph).
///
/// IMPORTANT: fds_divergence_part_2_preprocessing must be called sequentially
/// for all meshes in the preceding barrier BEFORE this task runs.
/// The preprocessing handles global zone ops (USUM modification, D_PBAR_DT
/// computation) which are not thread-safe.
template<MeshState OutS = MeshState::Default>
class DivergencePart2KernelTask
    : public hh::AbstractTask<1, MeshData<>, MeshData<OutS>> {
public:
    explicit DivergencePart2KernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData<>, MeshData<OutS>>(
              "DivPart2Kernel", numThreads) {}

    void execute(std::shared_ptr<MeshData<>> data) override {
        int kbar = fds_get_kbar(data->nm);
        fds_divergence_part_2_block_kernel(data->nm, data->dt, 1, kbar);
        if constexpr (OutS == MeshState::Default) {
            this->addResult(data);
        } else {
            this->addResult(data->template retag<OutS>());
        }
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<>, MeshData<OutS>>>
    copy() override {
        return std::make_shared<DivergencePart2KernelTask<OutS>>(
            this->numberThreads());
    }
};

#endif // DIVERGENCE_PART2_KERNEL_TASK_H
