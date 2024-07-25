#ifndef RAD_FILE_H
#define RAD_FILE_H

namespace mesh {

struct RadFile {
  int I1,I2,J1,J2,K1,K2,I_STEP=1,J_STEP=1,K_STEP=1,N_POINTS;
  double *IL_SAVE;
};

} // end namespace mesh

#endif
