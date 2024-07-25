#ifndef LUMPED_SPECIES_H
#define LUMPED_SPECIES_H

namespace mesh {

struct LumpedSpecies {
  double *ZZ;  ///< Lumped species, current time step, \f$ Z_{\alpha,ijk}^n \f$
  double *ZZS; ///< Lumped species, next time step, \f$ Z_{\alpha,ijk}^* \f$
};

} // end namespace mesh

#endif
