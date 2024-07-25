#ifndef CC_CUTCELL_H
#define CC_CUTCELL_H
#include "../../constants.h"

namespace mesh::cc {

struct CCCutcell {
   int IJK[MAX_DIM]; ///< Cartesian cell [ I J K ] location for this entry.
   int         NCELL=0; ///< Num cut-cells in this cartesian cell. Now fixed at 1 by blocking.
   int    NFACE_CELL=0; ///< Number of faces on cut-cells in this entry.
   int NFACE_DROPPED=0; ///< Parameter used for blocking cut-cells.
   double ALPHA_CC=0.; ///< Volume factor. ( SUM(VOLUME(1:NCELL))/(DX*DY*DZ) )
   int *CCELEM; ///< Cut-cells faces in FACE_LIST. (1:NFACE_CELL+1,1:NCELL)
   int *FACE_LIST; ///< List of cut-reg faces boundary of cut-cells for this entry.
   int *FACE_LIST_DROPPED; ///< Face list used for blocking cut-cells.
   int *IJK_LINK; ///< Cell/cut-cell each cut-cell is linked to.
   int *LINK_LEV; ///< Level in local Linking Hierarchy tree. (1:NCELL)
   double *VOLUME; ///< Cut-cell volumes. (1:NCELL)
   double *XYZCEN; ///< Cut-cell centroid locations. (IAXIS:KAXIS,1:NCELL)
   int *UNKZ; ///< Cut-cells unknown number for scalars.
   int *NOADVANCE; ///< Array to define if cut-cell should be blocked. (1:NCELL)
   int N_NOMICC=0; ///< Number of entries in NOMICC
   int *NOMICC; ///< OMESH cut-cells array. (1:2,1:N_NOMICC)

   double *RHO; ///< Corrector cut-cell densities. (1:NCELL)
   double *RHOS; ///< Predictor cut-cell densities.
   double *RSUM; ///< Cut-cells RSUM container. (1:NCELL)
   double *TMP; ///< Cut-cells temperatures. (1:NCELL)
   double *D; ///< Corrector cut-cells thermodynamic divg. (1:NCELL)
   double *DS; ///< Predictor cut-cells thermodynamic divg.
   double *DVOL; ///< Cut-cells thermodynamic divg*vol.
   double *DVOL_PR; ///< Predictor cut-cell thermodynamic divg*vol.
   double *Q; ///< Cut-cells volumetric heat source. (1:NCELL)
   double *QR; ///< Cut-cells volumetric radiative heat source.
   double *D_SOURCE; ///< Cut-cells divergence source terms.
   double *CHI_R; ///< Cut-cells radiative fraction.
   double *MIX_TIME; ///< Cut-cells species mixing time.
   double *Q_REAC; ///< Cut-cells volumetric heat source due to reaction.
   double *REAC_SOURCE_TERM; ///< Cut-cells species source term. (1:N_TOTAL_SCALARS,1:NCELL)
   double *ZZ; ///< Corrector cut-cells mass fractions.
   double *ZZS; ///< Predictor cut-cells mass fractions.
   double *M_DOT_PPP; ///< Cut-cells mass source term.

   int *UNKH; ///< Cut-cells unknown number for pressure H. (1:NCELL)
   double *KRES; ///< Cut-cells turbulent kinetic energy.
   double *H; ///< Cut-cells predictor pressure values.
   double *HS; ///< Cut-cells corrector pressure values.
   double *RTRM; ///< Cut-cells 1/(rho*c_p*T).
   double *R_H_G; ///< Cut-cells 1/(c_p*T)
   double *RHO_0; ///< Cut-cells background density.
   double *WVEL; ///< Cut-cells centroid vertical velocity.
   double *DDDTVOL; ///< Cut-cells dD/dT * vol.
   double *DELTA_RHO; ///< Cut-cells density change used in check mass density.
   double *DELTA_RHO_ZZ; ///< Cut-cells species density change used in check mass density.

   // Here: VIND=0, EP=1:INT_N_EXT_PTS, interpolation stencils for WVEL, RHO_0 into cut-cell centroids.
   int *INT_IJK;          ///< Interpolation (IAXIS:KAXIS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   double *INT_COEF;         ///< Interpolation (INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   double *INT_XYZBF,*INT_NOUT; ///< Interpolation (IAXIS:KAXIS,0:NCELL)
   int *INT_INBFC;        ///< Interpolation (1:3,0:NCELL)
   int *INT_NPE;          ///< Interpolation (LOW:HIGH,0,EP,0:NCELL)
   double *INT_XN,*INT_CN;    ///< Interpolation (0:INT_N_EXT_PTS,0:NCELL)  // 0 is interp pt.
   double *INT_CCVARS;       ///< Interpolation (1:N_INT_CCVARS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   int *INT_NOMIND;       ///< Interp. (LOW_IND:HIGH_IND,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)

   double *DEL_RHO_D_DEL_Z_VOL; ///< Cut-cells DEL_RHO_D_DEL_Z * VOL
   double *U_DOT_DEL_RHO_Z_VOL; ///< Cut-cells U_DOT_DEL_RHO_Z * VOL
};

} // end namespace mesh::cc


#endif
