#ifndef CFACE_H
#define CFACE_H

namespace mesh {

struct CFace {
  int CFACE_INDEX=0;     ///< Index of itself -- used to determine if the CFACE cell has been assigned
  int BC_INDEX=0;        ///< Derived type carrying coordinate variables
  int OD_INDEX=0;        ///< Derived type carrying 1-D solid info
  int B1_INDEX=0;        ///< Derived type carrying most of the surface boundary conditions
  int B2_INDEX=0;        ///< Derived type carrying most of the surface boundary conditions
  int BR_INDEX=0;        ///< Derived type carrying angular-specific radiation intensities
  int SURF_INDEX=0;
  int VENT_INDEX=0;
  int BOUNDARY_TYPE=0;
  int CUT_FACE_IND1=-11; ///< First index pointing to CUT_FACE array for this CFACE.
  int CUT_FACE_IND2=-11; ///< Second index pointing to CUT_FACE array for this CFACE.
  int N_REALS=0;         ///< Number of reals to pack into restart or send/recv buffer
  int N_INTEGERS=0;      ///< Number of integers to pack into restart or send/recv buffer
  int N_LOGICALS=0;      ///< Number of logicals to pack into restart or send/recv buffer
  double AREA=0.;        ///< CFACE area. From CUT_FACE(CUT_FACE_IND1)%AREA(CUT_FACE_IND2)
  double DUNDT=0.;
  double RSUM_G=0.;      ///< \f$ R_0 \sum_\alpha Z_\alpha/W_\alpha \f$ in first gas phase cell
  double MU_G=0.1;       ///< Viscosity, \f$ \mu \f$, in adjacent gas phase cell
  double PRES_BXN;       ///< Pressure H boundary condition for this CFACE (used in unstructured solvers).
};

} // end namespace mesh

#endif
