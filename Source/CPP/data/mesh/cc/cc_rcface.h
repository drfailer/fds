#ifndef CC_RCFACE_H
#define CC_RCFACE_H
#include "../../constants.h"

namespace mesh::cc {

struct CCRCface {
  bool SHAREDZ=false; ///< Notes if RC face connects cells with same UNKZ.
  int IWC=0; ///< WALL CELL index (if present) in location of RC Face.
  int PRES_ZONE=-1; ///< Pressure zone where RC face is.
  int IJK[MAX_DIM+1]; ///< Location indexes and axis of RC face. [ I J K X1AXIS]
  int UNKF=0; ///< Momentum unknown number if face is linked.
  int UNKZ[dim(LOW_IND, HIGH_IND)]; ///< Scalar transport unknown numbers in connected cells.
  int UNKH[dim(LOW_IND, HIGH_IND)]; ///< Pressure unknown numbers in connected cells.
  double XCEN[MAX_DIM * dim(LOW_IND, HIGH_IND)]; ///< Centroid location of connected cells.
  int JDH[dim(1, 2) * dim(1, 2)]; ///< Index matrix for H Poisson matrix coefficients.
  int CELL_LIST[(MAX_DIM+1) * dim(LOW_IND, HIGH_IND)]; ///< Cell type/location of connected cells. [RC_TYPE I J K ]
  double TMP_FACE=0.; ///< Temperature in RC face.
  double *ZZ_FACE; ///< Species mass fractions in RC face.
  double *RHO_D_DZDN; ///< Species diffusive mass fluxes in RC face.
  double *H_RHO_D_DZDN; ///< Species heat fluxes due to RHO_D_DZDN in RC face.
  int *NOMICF; ///< OMESH face info for mesh boundary RC face.
};

} // end namespace mesh::cc

#endif
