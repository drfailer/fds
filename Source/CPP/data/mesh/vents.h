#ifndef VENTS_H
#define VENTS_H
#include "../constants.h"

namespace mesh {

struct Vents {
  int I1 = -1, I2 = -1, J1 = -1, J2 = -1, K1 = -1, K2 = -1, BOUNDARY_TYPE = 0,
      IOR = 0, SURF_INDEX = 0, DEVC_INDEX = -1, CTRL_INDEX = -1,
      COLOR_INDICATOR = 99, TYPE_INDICATOR = 0, ORDINAL = 0,
      PRESSURE_RAMP_INDEX = 0, TMP_EXTERIOR_RAMP_INDEX = 0, NODE_INDEX = -1,
      OBST_INDEX = 0;
  int RGB[3] = {-1, -1, -1};
  double TRANSPARENCY = 1.;
  double TEXTURE[3] = {0};
  double X1 = 0., X2 = 0., Y1 = 0., Y2 = 0., Z1 = 0., Z2 = 0., FDS_AREA = 0.,
         X1_ORIG = 0., X2_ORIG = 0., Y1_ORIG = 0., Y2_ORIG = 0., Z1_ORIG = 0.,
         Z2_ORIG = 0., X0 = -9.E6, Y0 = -9.E6, Z0 = -9.E6, FIRE_SPREAD_RATE,
         UNDIVIDED_INPUT_AREA = 0., INPUT_AREA = 0., TMP_EXTERIOR = -1000.,
         DYNAMIC_PRESSURE = 0., RADIUS = -1.;
  double UVW[3] = {-1.E12, -1.E12, -1.E12};
  bool ACTIVATED = true, GEOM = false;
  char DEVC_ID[LABEL_LENGTH] = "null";
  char CTRL_ID[LABEL_LENGTH] = "null";
  char ID[LABEL_LENGTH] = "null";
  // turbulent inflow (experimental)
  int N_EDDY = 0;
  double R_IJ[3 * 3] = {0}, A_IJ[3 * 3] = {0}, SIGMA_IJ[3 * 3] = {0};
  double EDDY_BOX_VOLUME = 0.;
  double X_EDDY_MIN = 0., X_EDDY_MAX = 0., Y_EDDY_MIN = 0., Y_EDDY_MAX = 0.,
         Z_EDDY_MIN = 0., Z_EDDY_MAX = 0.;
  double *U_EDDY, *V_EDDY, *W_EDDY;
  double *X_EDDY, *Y_EDDY, *Z_EDDY, *CU_EDDY, *CV_EDDY, *CW_EDDY;
};

} // end namespace mesh

#endif
