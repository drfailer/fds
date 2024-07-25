#ifndef EDGE_H
#define EDGE_H
#include <cstddef>

namespace edge {

struct Edge {
  size_t I = 0,J = 0,K = 0;     ///< Indices of the edge
  size_t AXIS = 0;          ///< Edge axis, 1 = x, 2 = y, 3 = z
  size_t CELL_INDEX_MM = 0; ///< Index of adjacent cell, (I,J,K)
  size_t CELL_INDEX_PM = 0; ///< Index of adjacent cell, (I+1,J,K) or (I,J+1,K)
  size_t CELL_INDEX_MP = 0; ///< Index of adjacent cell, (I,J+1,K) or (I,J,K+1)
  size_t CELL_INDEX_PP = 0; ///< Index of adjacent cell, (I,J+1,K+1) or (I+1,J,K+1) or (I+1,J+1,K)
  size_t NOM_1 = 0;         ///< Number of Other Mesh, first direction
  size_t IIO_1 = 0;         ///< I index in Other mesh, first direction
  size_t JJO_1 = 0;         ///< J index in Other mesh, first direction
  size_t KKO_1 = 0;         ///< K index in Other mesh, first direction
  size_t NOM_2 = 0;         ///< Number of Other Mesh, second direction
  size_t IIO_2 = 0;         ///< I index in Other mesh, second direction
  size_t JJO_2 = 0;         ///< J index in Other mesh, second direction
  size_t KKO_2 = 0;         ///< K index in Other mesh, second direction
  double U_AVG = -1.E6;  ///< Stored velocity component at edge (m/s)
  double V_AVG = -1.E6;  ///< Stored velocity component at edge (m/s)
  double W_AVG = -1.E6;  ///< Stored velocity component at edge (m/s)
  bool EXTERNAL = false; ///< Edge is at the edge of the mesh
  // DIMENSION(-2:2) :: OMEGA = -1.E6;
  // DIMENSION(-2:2) :: TAU = -1.E6;
  // DIMENSION(2) :: EDGE_INTERPOLATION_FACTOR = 1.;
  double OMEGA[5] = { -1.E6, -1.E6, -1.E6, -1.E6, -1.E6 }; ///< Vorticity at boundary with one of four orientations
  double TAU[5] = { -1.E6, -1.E6, -1.E6, -1.E6, -1.E6 };   ///< Strain at boundary with one of four orientations
  double EDGE_INTERPOLATION_FACTOR[2] = { 1., 1. };        ///< Strain at boundary with one of four orientations
};

} // end namespace edge

#endif
