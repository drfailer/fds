#ifndef MESH_BOUNDARIES_BOUNDARY_COORDINATES_H
#define MESH_BOUNDARIES_BOUNDARY_COORDINATES_H
#include <cstddef>

namespace mesh::boundaries {

struct BoundaryCoordinates {
   size_t II;    ///< Ghost cell x index
   size_t JJ;    ///< Ghost cell y index
   size_t KK;    ///< Ghost cell z index
   size_t IIG;   ///< Gas cell x index
   size_t JJG;   ///< Gas cell y index
   size_t KKG;   ///< Gas cell z index
   size_t IOR=0; ///< Index of orientation of the WALL cell

   double *X;  ///< x coordinate of boundary cell center
   double *Y;  ///< y coordinate of boundary cell center
   double *Z;  ///< z coordinate of boundary cell center
   double *X1; ///< Lower x extent of boundary cell (m)
   double *X2; ///< Upper x extent of boundary cell (m)
   double *Y1; ///< Lower y extent of boundary cell (m)
   double *Y2; ///< Upper y extent of boundary cell (m)
   double *Z1; ///< Lower z extent of boundary cell (m)
   double *Z2; ///< Upper z extent of boundary cell (m)

   double NVEC[3] = {0}; ///< Normal vector
};

} // end namespace mesh::boundaries

#endif
