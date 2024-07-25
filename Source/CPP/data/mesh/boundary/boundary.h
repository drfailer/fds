#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <cstddef>
#include "boundary_1d.h"
#include "boundary_3d.h"
#include "boundary_coordinates.h"
#include "boundary_prop1.h"
#include "boundary_prop2.h"
#include "boundary_radia.h"

namespace mesh {

struct Boundary {
  boundaries::Boundary1D *BOUNDARY_ONE_D;
  boundaries::Boundary3D *BOUNDARY_THR_D;
  boundaries::BoundaryCoordinates *BOUNDARY_COORD;
  boundaries::BoundaryProp1 *BOUNDARY_PROP1;
  boundaries::BoundaryProp2 *BOUNDARY_PROP2;
  boundaries::BoundaryRadia *BOUNDARY_RADIA;
  size_t N_BOUNDARY_COORD_DIM=0,
         N_BOUNDARY_ONE_D_DIM=0,
         N_BOUNDARY_THR_D_DIM=0,
         N_BOUNDARY_PROP1_DIM=0,
         N_BOUNDARY_PROP2_DIM=0,
         N_BOUNDARY_RADIA_DIM=0;
  size_t NEXT_AVAILABLE_BOUNDARY_COORD_SLOT=1,
         NEXT_AVAILABLE_BOUNDARY_ONE_D_SLOT=1,
         NEXT_AVAILABLE_BOUNDARY_THR_D_SLOT=1,
         NEXT_AVAILABLE_BOUNDARY_PROP1_SLOT=1,
         NEXT_AVAILABLE_BOUNDARY_PROP2_SLOT=1,
         NEXT_AVAILABLE_BOUNDARY_RADIA_SLOT=1;
  size_t *BOUNDARY_COORD_OCCUPANCY,
         *BOUNDARY_ONE_D_OCCUPANCY,
         *BOUNDARY_THR_D_OCCUPANCY,
         *BOUNDARY_PROP1_OCCUPANCY,
         *BOUNDARY_PROP2_OCCUPANCY,
         *BOUNDARY_RADIA_OCCUPANCY,
         *PARTICLE_LAST;
};

} // end namespace mesh

#endif
