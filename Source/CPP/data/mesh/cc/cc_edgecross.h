#ifndef CC_EDGECROSS_H
#define CC_EDGECROSS_H
#include <cstddef>
#include "../../constants.h"

namespace mesh::cc {

struct CCEdgeCross {
   size_t NCROSS;  ///< Number of BODINT_PLANE segments - Cartesian edge crossings.
   double SVAR[dim(1, CC_MAXCROSS_EDGE)]; ///< Locations along x2 axis of NCROSS intersections.
   size_t ISVAR[dim(1, CC_MAXCROSS_EDGE)]; ///< Type of intersection (i.e. SG, GS or GG).
   size_t IJK[MAX_DIM + 1]; ///< [I J K X2AXIS]
};

} // end namespace mesh::cc

#endif
