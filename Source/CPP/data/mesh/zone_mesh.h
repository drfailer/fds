#ifndef ZONE_MESH_H
#define ZONE_MESH_H

#ifdef WITH_MKL
#include "../mkl/mkl_pardiso_handle.h"
#endif

namespace mesh {

struct ZoneMesh {
#ifdef WITH_MKL
  mkl::MKLPardisoHandle *T_H; ///< Internal solver memory pointer
#else
  int *PT_H;
#endif /* WITH_MKL */
  int NUNKH = 0;                 ///< Number of unknowns in pressure solution for a given ZONE_MESH
  int NCVLH = 0;                 ///< Number of pressure control volumes for a given ZONE_MESH
  int ICVL = 0;                  ///< Control volume counter for parent ZONE
  int IROW = 0;                  ///< Parent ZONE matrix row index
  int NUNKH_CART = 0;            ///< Number of unknowns in Cartesian cells of ZONE_MESH
  int NCVLH_CART = 0;            ///< Number of pressure CVs in Cartesian cells of ZONE_MESH
  int MTYPE = 0;                 ///< Matrix type (symmetric indefinite, or symm positive definite)
  int CONNECTED_ZONE_PARENT = 0; ///< Index of first zone in a connected zone list
  bool USE_FFT = true;           ///< Flag for use of FFT solver
  int *MESH_IJK;     ///< I,J,K positions of cell with unknown row IROW (1:3,1:NUNKH)
  double *A_H;       ///< Matrix coefficients for ZONE_MESH, up triang part, CSR format.
  int *IA_H, *JA_H;  ///< Matrix indexes for ZONE_MESH, up triang part, CSR format.
  double *F_H, *X_H; ///< RHS and Solution containers for the ZONE_MESH.
};

} // end namespace mesh

#endif
