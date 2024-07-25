#ifndef MESH_WALL_H
#define MESH_WALL_H
#include "wall/external_wall.h"
#include "wall/thin_wall.h"
#include "wall/wall.h"
#include <cstddef>

namespace mesh {

struct Wall {
  size_t N_WALL_CELLS = 0, N_WALL_CELLS_DIM = 0, N_INTERNAL_WALL_CELLS = 0,
         N_EXTERNAL_WALL_CELLS = 0;
  size_t N_THIN_WALL_CELLS = 0, N_THIN_WALL_CELLS_DIM = 0;
  wall::Wall *WALL;
  wall::ThinWall *THIN_WALL;
  wall::ExternalWall *EXTERNAL_WALL;
};

} // end namespace mesh

#endif
