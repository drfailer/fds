#ifndef CC_CUTEDGE_H
#define CC_CUTEDGE_H
#include "../../constants.h"
#include <cstddef>

namespace mesh::cc {

struct CCCutedge {
   size_t   IE=0; ///< Index to entry in structured mesh EDGE array.
   size_t  NVERT; ///< Number of vertices.
   size_t  NEDGE; ///< Number of cut-edges in this cartesian entry.
   size_t STATUS; ///< Location of cut edge, either CC_GASPHASE or CC_INBOUNDARY.
   size_t NEDGE1; ///< Auxiliary number of cut-edges in this cartesian entry.
   size_t IJK[MAX_DIM+2]; ///< Position in mesh of this entry [ I J K X2AXIS cetype]
   double *XYZVERT; ///< Locations of vertices. (IAXIS:KAXIS,1:NVERT)
   int *   CEELEM;  ///< Cut-Edge connectivities.
   int *   INDSEG;  ///< Boundary segment triangles and geometry [ntr tr1 tr2 ibod]
   int *VERT_LIST;  ///< Vertex list for XYZVERT [VERT_TYPE I J K]
   int *FACE_LIST;  ///< [1:3, -2:2, JEC] Cut-face connected to edge.
   double *      DXX; ///< [DXX(1,JEC) DXX(2,JEC)]
   double *   DUIDXJ; ///< Unstructured velocity Grad components.
   double *MU_DUIDXJ; ///< Unstructured viscosity * velocity Grad components.
   int *NOD_PERM; ///< Permutation array for INSERT_FACE_VERT.
};

} // end namespace mesh::cc

#endif
