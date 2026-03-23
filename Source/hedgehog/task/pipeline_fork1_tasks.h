#ifndef PIPELINE_FORK1_TASKS_H
#define PIPELINE_FORK1_TASKS_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Branch B task: COMBUSTION.
/// Runs combustion kernel per mesh in parallel.
class Fork1CombKernelTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit Fork1CombKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "Fork1CombKernel", numThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_combustion_kernel(data->nm, data->t, data->dt);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<Fork1CombKernelTask>(this->numberThreads());
    }
};

#endif // PIPELINE_FORK1_TASKS_H
