#ifndef MESH_BOUNDARIES_BOUNDARY_3D_H
#define MESH_BOUNDARIES_BOUNDARY_3D_H
#include <cstddef>

namespace mesh::boundaries {

struct Boundary3D {
  struct INTERNAL_NODE_TYPE {
     size_t *ALTERNATE_WALL_INDEX;  ///< Index of WALL cell in one of the two alternate directions
     size_t *ALTERNATE_WALL_NODE;   ///< Interior node of alternate WALL cell
     size_t *ALTERNATE_WALL_MESH;   ///< MESH number of alternate WALL cell
     size_t *ALTERNATE_WALL_TYPE;   ///< Type of alternate WALL cell (thin or not thin)
     size_t *ALTERNATE_WALL_IOR;    ///< Orientation index of alternate WALL cell
     double *ALTERNATE_WALL_WEIGHT; ///< Weight factor of alternate WALL cell
     int ALTERNATE_WALL_COUNT=0;    ///< Number of WALL cells that overlap the primary
     int I=-1;                      ///< I index of the node
     int J=-1;                      ///< J index of the node
     int K=-1;                      ///< K index of the node
     int MESH_NUMBER=-1;            ///< MESH number of the node
     bool HT3D=true;                ///< Indicates that this cell is meant to be updated in 3D
  } node; ///< Index of the interior solid cell
};

} // end namespace mesh::boundaries

#endif
