#ifndef CC_CUTFACE_H
#define CC_CUTFACE_H
#include "../../constants.h"

namespace mesh::cc {

struct CCCutface {
  int IWC = 0;        ///< Wall cell index if wall cell present.
  int PRES_ZONE = -1; ///< Pressure zone cut-face is immersed in.
  int NVERT = 0;      ///< Number of gas phase vertices in cartesian face/cell.
  int NSVERT = 0; ///< Number of solid phase Vertices in cartesian face/cell. For Plotting
  int NFACE = 0;  ///< Number of gas phase cut-faces in cartesian face/cell.
  int NSFACE = 0; ///< Number of solid phase cut-faces in cartesian face/cell. For Plotting
  int STATUS;     ///< Status of these cut-faces, CC_GASPHASE (cart faces), CC_INBOUNDARY (cart cells).

  int IJK[MAX_DIM + 1]; ///< Cartesian face/cell location for these cut-faces. [I J K X1AXIS]
  double *XYZVERT;      ///< Locations of vertices. (IAXIS:KAXIS,1:NVERT)
  int *CFELEM;          ///< Cut-faces connectivities.
  int *CEDGES;          ///< Cut-Edges per cut-face. Points to EDGE_LIST.
  double *AREA;         ///< Cut-faces areas. (1:NFACE+NSFACE)
  double *AREA_ADJUST;  ///< Cut-faces area adjustment factors.
  double *XYZCEN; ///< Cut-faces centroid locations. (IAXIS:KAXIS,1:NFACE+NSFACE)
  bool *SHARED;
  bool *BLK_TAG; ///< Array used in cut-cell blocking.
  int *LINK_LEV; ///< Level in gas cut-face linking hierarchy, used for momentum stability.
  double *INXAREA;   ///< Integral used in cut-cell volume computations.
  double *INXSQAREA; ///< Integral used in cut-cell x centroid computations.
  double *JNYSQAREA; ///< Integral used in cut-cell y centroid computations.
  double *KNZSQAREA; ///< Integral used in cut-cell z centroid computations.
  int *BODTRI;       ///< GEOMETRY and triangle associated with CC_INBOUNDARY cut-faces (CFACEs).
  int *EDGE_LIST; ///< Edges list. [CE_TYPE IEC JEC] or [RG_TYPE SIDE_LOHI AXIS]
  int *CELL_LIST; ///< Connected cut-cells list. [RC_TYPE I J K  ]
  int *CFACE_INDEX;  ///< In boundary cut-faces only, index in CFACE(:).
  int *SURF_INDEX;   ///< In boundary cut-faces only.
  int *CFACE_ORIGIN; ///< In boundary face origin (built, blocked small cell, blocked split cell).
  int *NOMICF; ///< For external boundary boundary cutfaces, NOM and cut-face list.

  int * UNKH; ///< Low and high side cut-cell H unknown number. (LOW:HIGH,1:NFACE)
  int * UNKZ; ///< Low and high side cut-cell Z unknown number. (LOW:HIGH,1:NFACE)
  int *UNKF; ///< Momentum unknown number, used for face linking. (1:NFACE)
  double *XCENLOW;      ///< Centroid position for cut-cells in low side. (IAXIS:KAXIS,1:NFACE)
  double *XCENHIGH;     ///< Centroid position for cut-cells in high side. (IAXIS:KAXIS,1:NFACE)
  double *ZZ_FACE;      ///< Scalar values interpolated to cut-faces.
  double *TMP_FACE;     ///< Gas phase cut-face temperature array. (1:NFACE)
  double *RHO_D_DZDN;   ///< Diffusive mass flux for species and cut-faces.
  double *H_RHO_D_DZDN; ///< Heat flux due to diffusive mass flux for species and cut-faces.
  double *VEL;          ///< Corrector velocity normal to cut-faces. (1:NFACE)
  double *VELS;         ///< Predictor velocity normal to cut-faces. (1:NFACE)
  double *FN;   ///< Momentum RHS in cut-faces (Advective+diffusive+body+...).
  double *FN_B; ///< Baroclinic Force RHS in cut-faces.
  double *VEL_SAVE;  ///< Saved unlinked velocities container for cut-faces.
  double *VEL_LNK;   ///< Linked velocity container for cut-faces.
  double *VEL_OMESH; ///< OMESH Corrector velocity normal to cut-faces of MESHES(NOM).
  double *VELS_OMESH; ///< OMESH Predictor velocity normal to cut-faces. (1:NFACE)
  double *VEL_LNK_OMESH; ///< OMESH Linked velocities.
  double *FN_OMESH;      ///< OMESH Momentum RHS.
  int *JDH;              ///< Index matrix per cutface in H Poisson matrix.
  double FV = 0., FV_B = 0.; ///< Momentum RHS and baroclinic torque in Cartesian face.
  double ALPHA_CF = 1.; ///< Area fraction for all gas cut-faces in a given cartesian face.
  double VEL_CF = 0., VEL_CRT = 0.; ///< Average cut-face velocity, cartesian velocity containers.

  // Here: VIND=IAXIS:KAXIS, EP=1:INT_N_EXT_PTS, INT_VEL_IND = 1; INT_VELS_IND = 2; INT_FV_IND = 3; INT_DHDX_IND = 4; N_INT_FVARS = 4; INT_NPE_LO = INT_NPE(LOW,VIND,EP,IFACE); INT_NPE_HI = INT_NPE(HIGH,VIND,EP,IFACE). Interpolation of variables for boundary cut-faces (CFACEs).
  int *INT_IJK;      ///< Interpolation (IAXIS:KAXIS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
  double *INT_COEF;  ///< Interpolation (INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
  double *INT_DCOEF; ///< Interpolation (IAXIS:KAXIS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
  double *INT_XYZBF, *INT_NOUT; ///< Interpolation (IAXIS:KAXIS,0:NFACE)
  int *INT_INBFC;               ///< Interpolation (1:3,0:NFACE)
  int *INT_NPE;                 ///< Interpolation (LOW:HIGH,VIND,EP,0:NFACE)
  double *INT_XN, *INT_CN; ///< Interpolation (0:INT_N_EXT_PTS,0:NFACE) //0 is interp point.
  double *INT_FVARS; ///< Interp (1:N_INT_FVARS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
  int *INT_NOMIND; ///< Interp (LOW_IND:HIGH_IND,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
  // Fields used in INBOUNDARY faces: Here: VIND=0, EP=1:INT_N_EXT_PTS, INT_H_IND=1, etc. N_INT_CVARS=INT_P_IND+N_TRACKED_SPECIES
  double *INT_CVARS; ///< Interpolation (1:N_INT_CVARS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
};

} // end namespace mesh::cc

#endif
