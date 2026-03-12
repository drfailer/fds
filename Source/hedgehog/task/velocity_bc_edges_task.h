#ifndef VELOCITY_BC_EDGES_TASK_H
#define VELOCITY_BC_EDGES_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/velocity_bc_data.h"
#include "../fds_fortran_interface.h"

/// Parallel task that calls MATCH_VELOCITY_KERNEL + VELOCITY_BC_PREPROCESSING +
/// VELOCITY_BC_PROCESS_EDGES_KERNEL for one mesh.
/// All three routines are thread-safe: explicit M% access, no POINT_TO_MESH.
class VelocityBCEdgesTask
    : public hh::AbstractTask<1, VelocityBCWork, VelocityBCWork> {
public:
    explicit VelocityBCEdgesTask(size_t numThreads)
        : hh::AbstractTask<1, VelocityBCWork, VelocityBCWork>(
              "VelocityBCEdges", numThreads) {}

    void execute(std::shared_ptr<VelocityBCWork> work) override {
        fds_match_velocity_kernel(work->nm, work->applyToEstimated);
        fds_velocity_bc_preprocessing(
            work->nm, work->t, work->applyToEstimated);
        fds_velocity_bc_process_edges_kernel(
            work->nm, work->t, work->applyToEstimated);
        this->addResult(work);
    }

    std::shared_ptr<hh::AbstractTask<1, VelocityBCWork, VelocityBCWork>>
    copy() override {
        return std::make_shared<VelocityBCEdgesTask>(this->numberThreads());
    }
};

#endif // VELOCITY_BC_EDGES_TASK_H
