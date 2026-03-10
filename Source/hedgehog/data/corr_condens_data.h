#ifndef CORR_CONDENS_DATA_H
#define CORR_CONDENS_DATA_H

#include <memory>
#include "mesh_data.h"

/// Work token for the corrector condensation sub-graph.
struct CorrCondensWork {
    int nm;
    double t;
    double dt;
    std::shared_ptr<MeshData> originalMeshData;

    CorrCondensWork(int nm_, double t_, double dt_, std::shared_ptr<MeshData> md)
        : nm(nm_), t(t_), dt(dt_), originalMeshData(md) {}
};

#endif // CORR_CONDENS_DATA_H
