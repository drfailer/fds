#ifndef DENSITY_H
#define DENSITY_H

namespace mesh {

struct Density {
  double *RHO;  ///< Density (kg/m3) at current time step, \f$ \rho_{ijk}^n \f$
  double *RHOS; ///< Density (kg/m3) at next time step, \f$ \rho_{ijk}^* \f$
  double *UII;  ///< Integrated intensity, \f$ U_{ijk}=\sum_{l=1}^N I_{ijk}^l\delta\Omega^l\f$
};

} // end namespace mesh

#endif
