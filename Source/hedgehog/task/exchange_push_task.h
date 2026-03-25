#ifndef EXCHANGE_PUSH_TASK_H
#define EXCHANGE_PUSH_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../tool/mesh_dependency_graph.h"
#include "../fds_fortran_interface.h"

/// Parallel task that pushes a mesh's flux data to all same-rank targets.
///
/// When mesh NM completes its baroclinic kernel, this task copies NM's
/// FVX/FVY/FVZ/H into each target's OMESH(NM) via fds_flux_copy_neighbor_ts.
/// The copy reads NM's pre-solve data, which is safe because the pressure
/// solve hasn't started yet (the downstream gate state holds NM until all
/// dependencies are satisfied).
///
/// After all copies complete, emits MeshData as a "push done" signal.
///
/// Thread safety: Multiple threads can push different source meshes
/// concurrently.  Each push writes to MESHES(target)%OMESH(source) — since
/// different source meshes write to different OMESH entries, there are no
/// write conflicts.  Reads from source meshes are concurrent-safe.
class ExchangePushTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    ExchangePushTask(size_t numThreads,
                     std::shared_ptr<MeshDependencyGraph> depGraph)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "ExchangePush", numThreads),
          depGraph_(std::move(depGraph)),
          myRank_(fds_mesh_process(depGraph_->lowerMesh())) {}

    void execute(std::shared_ptr<MeshData> md) override {
        for (int target : depGraph_->sendTargets(md->nm)) {
            if (fds_mesh_process(target) == myRank_) {
                fds_flux_copy_neighbor_ts(md->nm, target);
            }
        }
        this->addResult(md);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>> copy() override {
        return std::make_shared<ExchangePushTask>(
            this->numberThreads(), depGraph_);
    }

private:
    std::shared_ptr<MeshDependencyGraph> depGraph_;
    int myRank_;
};

#endif // EXCHANGE_PUSH_TASK_H
