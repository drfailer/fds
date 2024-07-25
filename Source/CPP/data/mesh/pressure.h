#ifndef MESH_PRESSURE_H
#define MESH_PRESSURE_H
#include <cstddef>

namespace mesh {

struct Pressure {
  double *PBAR;   ///< Background pressure, current, \f$ \overline{p}_m^n(z,t) \f$ (Pa)
  double *PBAR_S; ///< Background pressure, estimated, \f$ \overline{p}_m^*(z,t) \f$ (Pa)
  size_t *PRESSURE_ZONE; ///< Index of the pressure zone for cell (I,J,K)
  size_t *MUNKH; ///< Stores pressure unknown indexes for the mesh.
  double *PP_RESIDUAL; ///< Pressure Poisson residual (debug)
  double *PRHS; ///< Right hand side of Poisson (pressure) equation
};

}

#endif
