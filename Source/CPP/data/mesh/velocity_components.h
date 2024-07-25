#ifndef VELOCITY_COMPONENTS_H
#define VELOCITY_COMPONENTS_H

namespace mesh {

struct VelocityComponents {
  double *U; ///< Velocity component at current time step, \f$u_{ijk}^n\f$
  double *V; ///< Velocity component at current time step, \f$v_{ijk}^n\f$
  double *W; ///< Velocity component at current time step, \f$w_{ijk}^n\f$
  double *US; ///< Velocity component estimated at next time step, \f$u_{ijk}^*\f$
  double *VS; ///< Velocity component estimated at next time step, \f$v_{ijk}^*\f$
  double *WS; ///< Velocity component estimated at next time step, \f$w_{ijk}^*\f$
};

} // end namespace mesh

#endif
