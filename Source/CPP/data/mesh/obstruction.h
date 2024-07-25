#ifndef OBSTRUCTION_H
#define OBSTRUCTION_H
#include "../constants.h"

namespace mesh {

struct Obstruction {
  char DEVC_ID[LABEL_LENGTH] = "null"; ///< Name of controlling device
  char CTRL_ID[LABEL_LENGTH] = "null"; ///< Name of controller
  char ID[LABEL_LENGTH] = "null";      ///< Name of obstruction

  size_t SURF_INDEX[dim(-3, 3)] = {0}; ///< SURFace properties for each face
  size_t SURF_INDEX_INTERIOR = 0;      ///< SURFace properties for a newly exposed interior cell
  unsigned char RGB[3] = {0,0,0};      ///< Color indices for Smokeview

  double TRANSPARENCY = 1.;             ///< Transparency index for Smokeview, 0 = invisible, 1 = solid
  double VOLUME_ADJUST = 1.;            ///< Effective volume divided by user specified volume
  double BULK_DENSITY = -1.;            ///< Mass per unit volume (kg/m3) of specified OBST
  double X1 = 0.;                       ///< Lower specified \f$ x \f$ boundary (m)
  double X2 = 1.;                       ///< Upper specified \f$ x \f$ boundary (m)
  double Y1 = 0.;                       ///< Lower specified \f$ y \f$ boundary (m)
  double Y2 = 1.;                       ///< Upper specified \f$ y \f$ boundary (m)
  double Z1 = 0.;                       ///< Lower specified \f$ z \f$ boundary (m)
  double Z2 = 1.;                       ///< Upper specified \f$ z \f$ boundary (m)
  double MASS = 1.E6;                   ///< Actual mass of the obstruction (kg)
  double HEAT_SOURCE = 0.;              ///< Internal heat source (W/m3) used in HT3D applications
  double STRETCH_FACTOR = -1.;          ///< Node stretching factor used in HT3D applications
  double CELL_SIZE = -1.;               ///< Cell size (m) used in HT3D applications
  double CELL_SIZE_FACTOR = -1.;        ///< Cell size factor used in HT3D applications

  double INPUT_AREA[3]  =  { -1, -1, -1 };              ///< Specified area of x, y, and z faces (m2)
  double UNDIVIDED_INPUT_AREA[3]  =  { -1, -1, -1 };    ///< Area of x, y, z faces (m2) unbroken by mesh boundaries
  double UNDIVIDED_INPUT_LENGTH[3]  =  { -1, -1, -1 };  ///< Length in x, y, z direction (m) unbroken by mesh boundaries
  double SHAPE_AREA[3]  =  {0};               ///< Area of idealized top, sides, bottom (m2)
  double TEXTURE[3]  =  {0};                  ///< Origin of texture map (m)
  double FDS_AREA[3]  =  { -1., -1, -1 };                ///< Effective areas of x, y, and z faces (m2)
  double MATL_MASS_FRACTION[MAX_MATERIALS]  =  {0};      ///< Mass fraction of material component

  int I1 = -1;                ///< Lower I node
  int I2 = -1;                ///< Upper I node
  int J1 = -1;                ///< Lower J node
  int J2 = -1;                ///< Upper J node
  int K1 = -1;                ///< Lower K node
  int K2 = -1;                ///< Upper K node
  int COLOR_INDICATOR = -1;   ///< Coloring code: -3 = use specified color, -2 = invisible, -1 = no color specified
  int TYPE_INDICATOR = -1;    ///< Smokeview code: 2 = outline, -1 = solid
  int ORDINAL = 0;            ///< Order of OBST in input file
  int SHAPE_TYPE = -1;        ///< Indicator of shape carved out of larger obstruction
  int DEVC_INDEX = -1;        ///< Index of controlling device
  int CTRL_INDEX = -1;        ///< Index of controlling controller
  int DEVC_INDEX_O = -1;      ///< Original DEVC_INDEX
  int CTRL_INDEX_O = -1;      ///< Original CTRL_INDEX
  int MULT_INDEX = -1;        ///< Index of multiplier function
  int N_LAYER_CELLS_MAX = -1; ///< Maximum number of cells allowed in the layer, used in HT3D applications
  int RAMP_IHS_INDEX = 0;     ///< Index for internal heat source RAMP
  // WARN: this array is not initialized properly
  int MATL_INDEX[MAX_MATERIALS] = {-1};      ///< Index of material

  bool SHOW_BNDF[dim(-3, 3)] = {true,true,true,true,true,true,true}; ///< Show boundary quantities in Smokeview
  bool THIN = false;                      ///< The obstruction is zero cells thick
  bool HIDDEN = false;                    ///< Hide obstruction in Smokeview and ignore in simulation
  bool PERMIT_HOLE = true;                ///< Allow the obstruction to have a hole cutout
  bool ALLOW_VENT = true;                 ///< Allow a VENT to sit on the OBST
  bool CONSUMABLE = false;                ///< The obstruction can burn away
  bool REMOVABLE = false;                 ///< The obstruction can be removed from the simulation
  bool HOLE_FILLER = false;               ///< The obstruction fills a HOLE
  bool OVERLAY = true;                    ///< The obstruction can have another obstruction overlap a surface
  bool SCHEDULED_FOR_REMOVAL = false;     ///< The obstruction is scheduled for removal during the current time step
  bool SCHEDULED_FOR_CREATION = false;    ///< The obstruction is scheduled for creation during the current time step

};

} // end namespace mesh

#endif
