#ifndef FLUX_EXCHANGE_TASK_H
#define FLUX_EXCHANGE_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_exchange_data.h"
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Single-threaded task that performs bidirectional flux copies for a mesh.
///
/// Receives MeshExchangeData from MeshDependenciesManagerState.  For each
/// neighbor in the list, performs both pull (NOM → NM) and push (NM → NOM)
/// copies.  Neighbors already exchanged by a previous mesh are not in the
/// list, so no redundant copies are made.
///
/// Must be single-threaded (numThreads=1) to prevent races: if mesh A and
/// mesh B are emitted concurrently, A's push to B's OMESH could race with
/// B's solve reading that OMESH.  Single-threaded ensures exchanges are
/// serialized.
///
/// This task will be replaced by a CommunicatorTask for multi-process support.
class FluxExchangeTask
    : public hh::AbstractTask<1, MeshExchangeData, MeshData> {
public:
    FluxExchangeTask()
        : hh::AbstractTask<1, MeshExchangeData, MeshData>(
              "FluxExchange", 1) {}

    void execute(std::shared_ptr<MeshExchangeData> data) override {
        int nm = data->mesh->nm;
        for (int nom : data->neighbors) {
            // Pull: copy NOM's data into NM's ghost buffer
            fds_flux_copy_neighbor_ts(nom, nm);
            // Push: copy NM's data into NOM's ghost buffer
            fds_flux_copy_neighbor_ts(nm, nom);
        }
        this->addResult(data->mesh);
    }
};

#endif // FLUX_EXCHANGE_TASK_H
