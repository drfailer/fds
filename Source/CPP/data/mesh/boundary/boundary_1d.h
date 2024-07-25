#ifndef MESH_BOUDARIES_BOUNDARY_1D_H
#define MESH_BOUDARIES_BOUNDARY_1D_H
#include <cstddef>

namespace mesh::boundaries {

struct Boundary1D {
  double *M_DOT_S_PP;         ///< (1:SF\%N_MATL) Mass production rate of solid species
  double *X;                  ///< (0:NWP) Depth (m), \f$ x_{{\rm s},i} \f$
  double *TMP;                ///< Temperature in center of each solid cell, \f$ T_{{\rm s},i} \f$
  double *DELTA_TMP;          ///< Temperature change (K)
  double *LAYER_THICKNESS;    ///< (1:ONE_D\%N_LAYERS) Thickness of layer (m)
  double *MINIMUM_LAYER_THICKNESS; ///< (1:ONE_D\%N_LAYERS) Minimum thickness of layer (m)
  double *MIN_DIFFUSIVITY;    ///< (1:ONE_D\%N_LAYERS) Min diffusivity of all matls in layer (m2/s)
  double *RHO_C_S;            ///< Solid density times specific heat (J/m3/K)
  double *K_S;                ///< Solid conductivity (W/m/K)
  double *DDSUM;              ///< Scaling factor to get minimum cell size
  double *SMALLEST_CELL_SIZE; ///< Minimum cell size (m)
  double *PART_MASS;          ///< Accumulated mass of particles waiting to be injected (kg/m2)
  double *PART_ENTHALPY;      ///< Accumulated enthalpy of particles waiting to be injected (kJ/m2)
  double *CELL_SIZE_FACTOR;   ///< Amount to resize smallest grid cell near the surface
  double *CELL_SIZE;          ///< Specified constant cell size (m)
  double *STRETCH_FACTOR;     ///< Amount to stretch cells away from the surface
  double *HEAT_SOURCE;        ///< Heat addition within solid (W/m3)

  struct MatlCompType {
    double *MASS_FRACTION; ///< (1:N_LAYERS) Mass Fraction
    double *RHO;           ///< (1:NWP) Solid density (kg/m3)
  } *MATL_COMP; //< (1:SF\%N_MATL) Material component

   size_t *N_LAYER_CELLS;       ///< (1:SF\%N_LAYERS) Number of cells in the layer
   size_t *MATL_INDEX;          ///< (1:ONE_D\%N_MATL) Number of materials
   size_t *N_LAYER_CELLS_MAX;   ///< (1:SF\%N_LAYERS) Maximum possible number of cells in the layer
   size_t *RAMP_IHS_INDEX;      ///< (1:SF\%N_LAYERS) RAMP index for HEAT_SOURCE

   int SURF_INDEX=-1;    ///< SURFACE index
   size_t N_CELLS_MAX=0; ///< Maximum number of interior cells
   size_t N_CELLS_INI=0; ///< Initial number of interior cells
   size_t N_LAYERS=0;    ///< Number of material layers
   size_t N_MATL=0;      ///< Number of materials
   size_t N_LPC=0;       ///< Number of Lagrangian Particle Classes
   size_t BACK_INDEX=0;  ///< WALL index of back side of obstruction or exterior wall cell
   size_t BACK_MESH=0;   ///< Mesh number on back side of obstruction or exterior wall cell
   size_t BACK_SURF=0;   ///< SURF_INDEX on back side of obstruction or exterior wall cell

   bool *HT3D_LAYER; ///< (1:ONE_D\%N_LAYERS) Indicator that layer in 3D
};

} // end namespace mesh::boundaries

#endif
