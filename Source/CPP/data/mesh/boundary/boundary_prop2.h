#ifndef MESH_BOUNDARIES_BOUNDARY_PROP2_H
#define MESH_BOUNDARIES_BOUNDARY_PROP2_H

namespace mesh::boundaries {

struct BoundaryProp2 {
   double *A_LP_MPUA; ///< Accumulated liquid droplet mass per unit area (kg/m2)
   double *LP_CPUA;   ///< Liquid droplet cooling rate unit area (W/m2)
   double *LP_EMPUA;  ///< Ember masss generated per unit area (kg/m2)
   double *LP_MPUA;   ///< Liquid droplet mass per unit area (kg/m2)
   double *LP_TEMP;   ///< Liquid droplet mean temperature (K)

   double U_TAU=0.;         ///< Friction velocity (m/s)
   double Y_PLUS=1.;        ///< Dimensionless boundary layer thickness unit
   double Z_STAR=1.;        ///< Dimensionless boundary layer unit
   double PHI_LS=-1.;       ///< Level Set value for output only
   double WORK1=0.;         ///< Work array
   double WORK2=0.;         ///< Work array
   double WORK3=0.;         ///< Work array
   double K_SUPPRESSION=0.; ///< Suppression coefficent (m2/kg/s)
   double V_DEP=0.;         ///< Deposition velocity (m/s)

   int SURF_INDEX=-1;          ///< Surface index
   int HEAT_TRANSFER_REGIME=0; ///< 1=Forced convection, 2=Natural convection, 3=Impact convection, 4=Resolved
};

} // end namespace mesh::boundaries

#endif

