#ifndef DIV_SETUP_DATA_H
#define DIV_SETUP_DATA_H

#include "mesh_data.h"
#include <memory>

/// Work token for parallel velocity flux kernel execution.
struct DivSetupWork {
    int nm;         ///< Mesh index
    double t;       ///< Simulation time
    double dt;      ///< Time step
    int estimated;  ///< 0=predictor (use U,V,W), 1=corrector (use US,VS,WS)

    /// Preserve original MeshData for downstream routing
    std::shared_ptr<MeshData> originalMeshData;

    DivSetupWork(int nm_, double t_, double dt_, int estimated_,
                 std::shared_ptr<MeshData> md)
        : nm(nm_), t(t_), dt(dt_), estimated(estimated_),
          originalMeshData(md) {}
};

#endif // DIV_SETUP_DATA_H
