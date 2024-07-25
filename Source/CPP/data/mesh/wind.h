#ifndef WIND_H
#define WIND_H

namespace mesh {

struct Wind {
  double *U_WIND; ///< Component of wind in x direction
  double *V_WIND; ///< Component of wind in y direction
  double *W_WIND; ///< Component of wind in z direction
};

}

#endif
