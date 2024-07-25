#ifndef EXTERNAL_WALL_H
#define EXTERNAL_WALL_H

namespace mesh::wall {

struct ExternalWall {
   int NOM;               ///< Number of the adjacent (Other) Mesh
   int NIC_MIN;           ///< Start of indices for the cell in the other mesh
   int NIC_MAX;           ///< End of indices for the cell in the other mesh
   int NIC;               ///< NIC_MAX-NIC_MIN
   int IIO_MIN;           ///< Minimum I index of adjacent cell in other mesh
   int IIO_MAX;           ///< Maximum I index of adjacent cell in other mesh
   int JJO_MIN;           ///< Minimum J index of adjacent cell in other mesh
   int JJO_MAX;           ///< Maximum J index of adjacent cell in other mesh
   int KKO_MIN;           ///< Minimum K index of adjacent cell in other mesh
   int KKO_MAX;           ///< Maximum K index of adjacent cell in other mesh
   int PRESSURE_BC_TYPE;  ///< Poisson boundary condition, NEUMANN or DIRICHLET
   int SURF_INDEX_ORIG=0; ///< Original SURFace index for this cell
   double AREA_RATIO;     ///< Ratio of face areas of adjoining cells
   double DUNDT=0;        ///< \f$ \partial u_n / \partial t \f$
   double *FVN;           ///< Flux-limited \f$ \int \rho Y_\alpha u_n \f$
   double *FVNS;          ///< Estimated value of FVN at next time step
   double *RHO_D_DZDN;    ///< Species diffusive flux as computed in divg.f90
   double *RHO_D_DZDNS;   ///< RHO_D_DZDN estimated at next time step
};

} // end namespace mesh::wall

#endif
