#ifndef MESH_DIVERGENCE_H
#define MESH_DIVERGENCE_H

namespace mesh {

struct Divergence {
  double *DDDT; ///< \f$(\partial D/\partial t)_{ijk}\f$ where \f$D_{ijk}\f$ is the divergence
  double *D;    ///< Divergence at current time step, \f$D_{ijk}^n\f$
  double *DS;   ///< Divergence estimate next time step, \f$D_{ijk}^*\f$
  double *D_SOURCE; ///< Source terms in the expression for the divergence
};

} // end namespace mesh

#endif
