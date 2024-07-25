#ifndef CC_REGFACEZ_H
#define CC_REGFACEZ_H
#include "../../constants.h"

namespace mesh::cc {

struct CCRegfaceZ {
   bool DO_LO_IND=false; ///< Defines if fluxes are to be assigned in low side cell.
   bool DO_HI_IND=false; ///< Defines if fluxes are to be assigned in high side cell.
   int IJK[MAX_DIM]; ///< Location indexes [ I J K ] of X,Y,Z REG face.
   int IWC=0; ///< WALL CELL index (if present) in location of REG Face.
   int *NOMICF; ///< OMESH face info for mesh boundary REG face.
   double FN_H_S=0.; ///< Stores components FX_H_S, FY_H_S, FZ_H_S as needed.
   double *RHOZZ_U; ///< Stores computed FX(I,J,K,N)*UU(I,J,K), etc. as needed.
   double *FN_ZZ; ///< Stores computed FX_ZZ(I,J,K), etc. as needed.
   double *RHO_D_DZDN; ///< Species diffusive mass fluxes in REG face.
   double *H_RHO_D_DZDN; ///< OMESH face info for mesh boundary REG face.
};

} // end namespace mesh::cc


#endif
