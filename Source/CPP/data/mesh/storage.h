#ifndef MESH_STORAGE_H
#define MESH_STORAGE_H

namespace mesh {

struct Storage {
  int N_ITEMS=0;        ///< Number of WALL cells, CFACEs, or PARTICLEs listed in ITEM_INDEX
  int N_ITEMS_DIM=0;    ///< Dimension of 1-D arrays ITEM_INDEX and SURF_INDEX
  int *ITEM_INDEX;      ///< Array of indices of the WALL cells, CFACEs, or PARTICLEs
  int *SURF_INDEX;      ///< Array of SURF indices of the WALL cells, CFACEs, or PARTICLEs
  int N_REALS_DIM=0;    ///< Dimension of the array REALS
  int N_INTEGERS_DIM=0; ///< Dimension of the array INTEGERS
  int N_LOGICALS_DIM=0; ///< Dimension of the array LOGICALS
  int N_REALS=0;        ///< Number of reals stored in REALS
  int N_INTEGERS=0;     ///< Number of integers stored in INTEGERS
  int N_LOGICALS=0;     ///< Number of logicals stored in LOGICALS
  double *REALS;        ///< Array of reals
  int *INTEGERS;        ///< Array of integers
  bool *LOGICALS;       ///< Array of logicals
};

};

#endif
