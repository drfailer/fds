#ifndef CELL_H
#define CELL_H
#include <cstddef>

namespace cell {

struct Cell {
   bool SOLID=false;         ///< Indicates if grid cell is solid or not
   bool EXTERIOR=false;      ///< Indicates if the grid cell is outside the mesh
   bool EXTERIOR_EDGE=false; ///< Indicates if the grid cell is at an outside edge of the mesh
   size_t OBST_INDEX=0;      ///< Index of obstruction that fills cell (0 if none)
   size_t I,J,K;             ///< Indices of cell
   // DIMENSION(12) :: EDGE_INDEX=0
   // DIMENSION(-3:3) :: WALL_INDEX=0
   // DIMENSION(-3:3) :: SURF_INDEX=0
   // DIMENSION(-3:3,1:3) :: THIN_WALL_INDEX=0
   // DIMENSION(-3:3,1:3) :: THIN_SURF_INDEX=0
   // DIMENSION(-3:3,1:3) :: THIN_OBST_INDEX=0
   size_t EDGE_INDEX[12]={0}; ///< Indices of 12 cell edges
   size_t WALL_INDEX[7]={0};  ///< Indices of 6 adjacent wall cells
   size_t SURF_INDEX[7]={0};  ///< SURF indices of 6 adjacent wall cells
   size_t THIN_WALL_INDEX[21]={0}; ///< Indices of 6 adjacent thin wall cells
   size_t THIN_SURF_INDEX[21]={0}; ///< SURF indices of 6 adjacent thin wall cells
   size_t THIN_OBST_INDEX[21]={0}; ///< OBST indices of 6 adjacent thin wall cells
};

} // end namespace cell

#endif
