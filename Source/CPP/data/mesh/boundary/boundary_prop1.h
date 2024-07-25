#ifndef MESH_BOUNDARIES_BOUNDARY_PROP1_H
#define MESH_BOUNDARIES_BOUNDARY_PROP1_H

namespace mesh::boundaries {

struct BoundaryProp1 {
  double *M_DOT_G_PP_ACTUAL; ///< (1:N_TRACKED_SPECIES) Actual mass production rate per unit area
  double *M_DOT_G_PP_ADJUST; ///< (1:N_TRACKED_SPECIES) Adjusted mass production rate per unit area
  double *ZZ_F;              ///< (1:N_TRACKED_SPECIES) Species mixture mass fraction at surface
  double *ZZ_G;              ///< (1:N_TRACKED_SPECIES) Species mixture mass fraction in gas cell
  double *RHO_D_F;           ///< (1:N_TRACKED_SPECIES) Diffusion at surface, \f$ \rho D_\alpha \f$
  double *RHO_D_DZDN_F;      ///< \f$ \rho D_\alpha \partial Z_\alpha / \partial n \f$
  double *AWM_AEROSOL;       ///< Accumulated aerosol mass per unit area (kg/m2)
  double *QDOTPP_INT;        ///< Integrated HRRPUA for the S_Pyro method (kJ/m2)
  double *Q_IN_SMOOTH_INT;   ///< Integrated smoothed incident flux for the S_Pyro method (kJ)


   int SURF_INDEX=-1;   ///< SURFACE index
   int PRESSURE_ZONE=0; ///< Pressure ZONE of the adjacent gas phase cell
   int NODE_INDEX=0;    ///< HVAC node index associated with surface
   int N_SUBSTEPS=1;    ///< Number of substeps in the 1-D conduction/reaction update

   double AREA=0.;            ///< Face area (m2)
   double HEAT_TRANS_COEF=0.; ///< Heat transfer coefficient (W/m2/K)
   double Q_CON_F=0.;         ///< Convective heat flux at surface (W/m2)
   double Q_RAD_IN=0.;        ///< Incoming radiative flux (W/m2)
   double Q_RAD_OUT=0.;       ///< Outgoing radiative flux (W/m2)
   double EMISSIVITY=1.;      ///< Surface emissivity
   double AREA_ADJUST=1.;     ///< Ratio of actual surface area to grid cell face area
   double T_IGN=0.;           ///< Ignition time (s)
   double TMP_F;              ///< Surface temperature (K)
   double TMP_G;              ///< Near-surface gas temperature (K)
   double TMP_F_OLD;          ///< Holding value for surface temperature (K)
   double TMP_B;              ///< Back surface temperature (K)
   double U_NORMAL=0.;        ///< Normal component of velocity (m/s) at surface, start of time step
   double U_NORMAL_S=0.;      ///< Estimated normal component of velocity (m/s) at next time step
   double U_NORMAL_0=0.;      ///< Initial or specified normal component of velocity (m/s) at surface
   double U_TANG=0.;          ///< Tangential velocity (m/s) near surface
   double U_IMPACT=0.;        ///< Impact velocity from stagnation pressure
   double RHO_F;              ///< Gas density at the wall (kg/m3)
   double RHO_G;              ///< Gas density in near wall cell (kg/m3)
   double RDN=1.;             ///< \f$ 1/ \delta n \f$ at the surface (1/m)
   double K_G=0.025;          ///< Thermal conductivity of gas in adjacent gas phase cell near wall
   double Q_DOT_G_PP=0.;      ///< Heat release rate per unit area (W/m2)
   double Q_DOT_O2_PP=0.;     ///< Heat release rate per unit area (W/m2) due to oxygen consumption
   double Q_CONDENSE=0.;      ///< Heat release rate per unit area (W/m2) due to gas condensation
   double BURN_DURATION=0.;   ///< Duration of a specified fire (s)
   double Q_IN_SMOOTH=0.;     ///< Smoothed incident flux for the S_Pyro method (kJ)
   double T_MATL_PART=0.;     ///< Time interval for current value in PART_MASS and PART_ENTHALPY arrays (s)
   double B_NUMBER=0.;        ///< B number for droplet or wall
   double M_DOT_PART_ACTUAL;  ///< Mass flux of all particles (kg/m2/s)
   double Q_LEAK=0.;          ///< Heat production of leaking gas (W/m3)
   double VEL_ERR_NEW=0.;     ///< Velocity mismatch at mesh or solid boundary (m/s)

   bool BURNAWAY=false;      ///< Indicater if cell can burn away when fuel is exhausted
   bool LAYER_REMOVED=false; ///< Indicator that at least one layer has been removed during the time step
};

} // end namespace mesh::boundaries

#endif
