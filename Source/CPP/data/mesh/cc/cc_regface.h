#ifndef CC_REGFACE_H
#define CC_REGFACE_H
#include "../../constants.h"

namespace mesh::cc {

struct CCRegface {
   int PRES_ZONE=-1; ///< Pressure Zone this regular face is located in.
   int IJK[MAX_DIM]; ///< Indexes [I J K] of X, Y or Z face.
   int JD[dim(1, 2) * dim(1, 2)];  ///< Indexes in H matrix where this face contributes coefficients.
};

} // end namespace mesh::cc

#endif
