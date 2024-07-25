#ifndef AMBIENT_H
#define AMBIENT_H

namespace mesh {

struct Ambient {
  double *P_0;   ///< Ambient pressure profile, \f$ \overline{p}_0(z) \f$ (Pa)
  double *RHO_0; ///< Ambient density profile, \f$ \overline{\rho}_0(z) \f$ (kg/m\f$^3\f$)
  double *TMP_0; ///< Ambient temperature profile, \f$ \overline{T}_0(z) \f$ (K)
};

} // end nemespace mesh

#endif
