#ifndef MESH_RADIATION_H
#define MESH_RADIATION_H

namespace mesh {

struct Radiation {
  double *KAPPA_GAS; ///< Radiation absorption coefficient by gas, \f$ \kappa_{ijk} \f$
  double *CHI_R;     ///< Radiative fraction, \f$ \chi_{{\rm r},ijk} \f$
  double *QR;        ///< Radiation source term, \f$ -\nabla \cdot \dot{\mathbf{q}}_{\rm r}'' \f$
  double *QR_W;      ///< Radiation source term, particles and droplets
};

}

#endif
