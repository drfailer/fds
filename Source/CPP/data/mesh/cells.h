#ifndef MESH_CELLS_H
#define MESH_CELLS_H
#include <cstddef>
#include "../cell/cell.h"

namespace mesh {

struct Cells {
  size_t *CELL_INDEX; ///< Unique integer identifier for grid cell (I,J,K)
  cell::Cell *CELL;   ///< Variables associated with nearby solids
  double *CELL_ILW;
  int *CELL_INTEGERS;  ///< 1-D array to hold CELL variables
  bool *CELL_LOGICALS; ///< 1-D array to hold CELL variables
};

} // end namespace mesh

#endif
