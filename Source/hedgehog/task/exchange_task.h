#ifndef EXCHANGE_TASK_H
#define EXCHANGE_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/exchange_work.h"
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel task that performs mesh exchange copies.
///
/// Receives ExchangePullWork — all same-rank senders have arrived at the
/// orchestrator state, so their data is ready to be copied.  For each sender
/// in the work token, calls fds_flux_copy_neighbor_ts(sender, destNM) to
/// copy the sender's FVX/FVY/FVZ/H into destNM's OMESH(sender).
/// Then emits destNM's MeshData downstream.
///
/// Thread safety: Multiple threads can execute simultaneously because each
/// PullWork targets a different destination mesh.  The copies write to
/// disjoint OMESH entries: MESHES(dest1)%OMESH(*) vs MESHES(dest2)%OMESH(*).
/// Reads from source meshes are concurrent-safe (read-only).
class ExchangeTask
    : public hh::AbstractTask<1, ExchangePullWork, MeshData> {
public:
    explicit ExchangeTask(size_t numThreads)
        : hh::AbstractTask<1, ExchangePullWork, MeshData>(
              "ExchangeTask", numThreads) {}

    void execute(std::shared_ptr<ExchangePullWork> work) override {
        // Copy data from each same-rank sender into destination's OMESH
        for (int sender : work->sameRankSenders) {
            fds_flux_copy_neighbor_ts(sender, work->mesh->nm);
        }
        // (Future: MPI_Wait/unpack for cross-rank senders would go here)

        // Destination mesh is ready to proceed
        this->addResult(work->mesh);
    }

    std::shared_ptr<hh::AbstractTask<1, ExchangePullWork, MeshData>>
    copy() override {
        return std::make_shared<ExchangeTask>(this->numberThreads());
    }
};

#endif // EXCHANGE_TASK_H
