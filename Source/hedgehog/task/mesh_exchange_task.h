#ifndef FLUX_EXCHANGE_TASK_H
#define FLUX_EXCHANGE_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_exchange_data.h"
#include "../fds_fortran_interface.h"

/// Parallel task that performs pull-only flux copies for a mesh.
///
/// Part of a cycle with MeshDependenciesManagerState:
///   State emits MeshExchangeData -> this task -> MeshExchangeData back to state.
///
/// Pull model: each mesh pulls data from ALL its same-rank neighbors.
/// fds_flux_copy_neighbor_ts(NOM, NM) copies NOM's source data into NM's
/// OMESH buffer.  Each mesh writes only to its own OMESH entries, so
/// different meshes can be exchanged in parallel without races.
class FluxExchangeTask
    : public hh::AbstractTask<1, MeshExchangeData, MeshExchangeData> {
public:
    FluxExchangeTask(size_t numThreads)
        : hh::AbstractTask<1, MeshExchangeData, MeshExchangeData>(
              "FluxExchange", numThreads) {}

    void execute(std::shared_ptr<MeshExchangeData> data) override {
        int nm = data->mesh->nm;
        for (int nom : data->neighbors) {
            fds_flux_copy_neighbor_ts(nom, nm);
        }
        this->addResult(data);
    }
};

#endif // FLUX_EXCHANGE_TASK_H
