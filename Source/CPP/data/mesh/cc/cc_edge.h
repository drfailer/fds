#ifndef CC_EDGE_H
#define CC_EDGE_H
#include "../../constants.h"

namespace mesh::cc {

struct CCEdge {
   int IJK[MAX_DIM+1];            ///< Edge I,J,K,AXIS indexesin mesh. [ I J K X1AXIS]
   int IE=0;           ///< Index to entry in structured mesh EDGE array.
   int *FACE_LIST;      ///< [1:3, -2:2] Reg/Cut faces connected to edge.
   // Here: VIND=IAXIS:KAXIS, EP=1:INT_N_EXT_PTS,
   // INT_VEL_IND = 1; INT_VELS_IND = 2; INT_FV_IND = 3; INT_DHDX_IND = 4; N_INT_FVARS = 4;
   // INT_NPE_LO = INT_NPE(LOW,VIND,EP,IFACE); INT_NPE_LO = INT_NPE(HIGH,VIND,EP,IEDGE).
   // IEDGE = 0, Cartesian GASPHASE EDGE.
   int *INT_IJK;        ///< Interpolation (IAXIS:KAXIS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   double *INT_COEF;       ///< Interpolation (INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   double *INT_DCOEF;      ///< Interpolation (IAXIS:KAXIS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   double *INT_XYZBF;      ///< Interpolation (IAXIS:KAXIS,IEDGE)
   double *INT_NOUT;       ///< Interpolation (IAXIS:KAXIS,IEDGE)
   int *INT_INBFC;      ///< Interpolation (1:3,IEDGE)
   int *INT_NPE;        ///< Interpolation (LOW:HIGH,VIND,EP,IEDGE)
   double *INT_XN,*INT_CN;  ///< Interpolation (0:INT_N_EXT_PTS,IEDGE)
   double *INT_FVARS;      ///< Interpolation (1:N_INT_FVARS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   double *INT_CVARS;      ///< Interpolation (1:N_INT_CVARS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   int *INT_NOMIND;     ///< Interp. (LOW_IND:HIGH_IND,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   double *XB_IB;          ///< For edges inside geometry, distance on each face to boundary.
   double *DUIDXJ;         ///< Wall modeled velocity Grad components.
   double *MU_DUIDXJ;      ///< Wall modeled viscosity * velocity Grad components.
   int *SURF_INDEX;     ///< Surface type for boundary on each face connected to edge in geom.
   bool *PROCESS_EDGE_ORIENTATION; ///< Computation logical.
   bool *EDGE_IN_MESH;   ///< Logical stating edge inside mesh (not boundary).
   bool *SIDE_IN_GEOM;   ///< Logical stating edge side inside geometry or actual boundary.
};

} // end namespace mesh::cc

#endif
