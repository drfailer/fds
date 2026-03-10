#ifndef CORR_PARTICLE_DATA_H
#define CORR_PARTICLE_DATA_H

#include <memory>
#include "mesh_data.h"

/// Work token for the corrector particle sub-graph.
struct CorrParticleWork {
    int nm;
    double t;
    double dt;
    std::shared_ptr<MeshData> originalMeshData;

    CorrParticleWork(int nm_, double t_, double dt_, std::shared_ptr<MeshData> md)
        : nm(nm_), t(t_), dt(dt_), originalMeshData(md) {}
};

#endif // CORR_PARTICLE_DATA_H
