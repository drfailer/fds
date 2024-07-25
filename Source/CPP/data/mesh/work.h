#ifndef MESH_WORK_H
#define MESH_WORK_H
#include "storage.h"
#include <cstddef>
#include <complex.h>

namespace mesh {

struct Work {
  double *SCALAR_WORK1;
  double *SCALAR_WORK2;
  double *SCALAR_WORK3;
  double *SCALAR_WORK4;
  double *WORK1, WORK2, WORK3, WORK4, WORK5, WORK6, WORK7, WORK8, WORK9;
  size_t *IWORK1;
  double *PWORK1, PWORK2, PWORK3, PWORK4;
  _Complex double *PWORK5, PWORK6, PWORK7, PWORK8;
  double *TURB_WORK1, TURB_WORK2, TURB_WORK3, TURB_WORK4;
  double *TURB_WORK5, TURB_WORK6, TURB_WORK7, TURB_WORK8;
  double *TURB_WORK9, TURB_WORK10;
  double *CCVELDIV, CARTVELDIV;
  double *WALL_WORK1, WALL_WORK2, FACE_WORK1, FACE_WORK2, FACE_WORK3;
  double *QQ, QQ2;
  double *PP, PPN, BNDF_TIME_INTEGRAL;
  size_t *IBK;
  size_t *IBLK;
  Storage WALL_STORAGE, CFACE_STORAGE;
};

}

#endif
