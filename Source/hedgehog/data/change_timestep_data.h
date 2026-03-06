#ifndef CHANGE_TIMESTEP_DATA_H
#define CHANGE_TIMESTEP_DATA_H

#include <memory>
#include <vector>
#include "mesh_data.h"

/// Data token flowing through the retry sequence
struct RetrySequenceData {
    std::vector<std::shared_ptr<MeshData>> meshes;
    double t;
    double dt;
    int iteration;     ///< Retry iteration count (for debugging)
    bool done;         ///< True when retry loop should exit

    RetrySequenceData(const std::vector<std::shared_ptr<MeshData>>& m, double t_, double dt_,
                      int iter = 0, bool d = false)
        : meshes(m), t(t_), dt(dt_), iteration(iter), done(d) {}

    int nm_count() const { return static_cast<int>(meshes.size()); }
};

#endif // CHANGE_TIMESTEP_DATA_H
