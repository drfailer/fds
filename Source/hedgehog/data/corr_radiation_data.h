#ifndef CORR_RADIATION_DATA_H
#define CORR_RADIATION_DATA_H

#include "mesh_data.h"
#include <memory>

struct CorrRadiationWork {
    int nm;
    double t;
    int radIter;
    double radQSumPartial;
    double kfst4SumPartial;
    std::shared_ptr<MeshData> originalMeshData;

    CorrRadiationWork(int nm_, double t_, int radIter_,
                      std::shared_ptr<MeshData> md)
        : nm(nm_), t(t_), radIter(radIter_),
          radQSumPartial(0.0), kfst4SumPartial(0.0),
          originalMeshData(md) {}
};

#endif // CORR_RADIATION_DATA_H
