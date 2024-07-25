#ifndef VISCOSITY_H
#define VISCOSITY_H

namespace mesh {

struct Viscosity {
  double *MU;     ///< Turbulent viscosity (kg/m/s), \f$ \mu_{{\rm t},ijk} \f$
  double *MU_DNS; ///< Laminar viscosity (kg/m/s)
};

} // end namespace mesh

#endif
