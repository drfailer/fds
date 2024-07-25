#ifndef PARTICLE_DRAG_H
#define PARTICLE_DRAG_H

namespace mesh {

struct ParticuleDrag {
  double *FVX_D; ///< Particle drag
  double *FVY_D; ///< Particle drag
  double *FVZ_D; ///< Particle drag
};

} // end namespace mesh

#endif
