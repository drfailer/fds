#ifndef MOMENTUM_FLUX_TERMS_H
#define MOMENTUM_FLUX_TERMS_H

namespace mesh {

struct MomentumFluxTerms {
  double *FVX;   ///< Momentum equation flux terms, \f$ F_{{\rm A},x,ijk}+F_{{\rm B},x,ijk} \f$
  double *FVY;   ///< Momentum equation flux terms, \f$ F_{{\rm A},y,ijk}+F_{{\rm B},y,ijk} \f$
  double *FVZ;   ///< Momentum equation flux terms, \f$ F_{{\rm A},z,ijk}+F_{{\rm B},z,ijk} \f$
  double *FVX_B; ///< Momentum equation flux terms, \f$ F_{{\rm B},x,ijk} \f$
  double *FVY_B; ///< Momentum equation flux terms, \f$ F_{{\rm B},y,ijk} \f$
  double *FVZ_B; ///< Momentum equation flux terms, \f$ F_{{\rm B},z,ijk} \f$
};

} // end namespace mesh

#endif
