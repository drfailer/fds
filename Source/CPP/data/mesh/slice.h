#ifndef MESH_SLICE_H
#define MESH_SLICE_H
#include "../constants.h"

namespace mesh {

struct Slice {
  int I1,I2,J1,J2,K1,K2,GEOM_INDEX=-1,INDEX,INDEX2=0,Z_INDEX=-999,Y_INDEX=-999,MATL_INDEX=-999,
             PART_INDEX=0,VELO_INDEX=0,PROP_INDEX=0,REAC_INDEX=0,SLCF_INDEX;
  double MINMAX[2];
  float RLE_MIN, RLE_MAX;
  double AGL_SLICE;
  bool TERRAIN_SLICE=false,CELL_CENTERED=false,RLE=false,DEBUG=false,THREE_D=false;
  char SLICETYPE[LABEL_LENGTH]="STRUCTURED",SMOKEVIEW_LABEL[LABEL_LENGTH];
  char SMOKEVIEW_BAR_LABEL,ID[LABEL_LENGTH]="null",MATL_ID[LABEL_LENGTH]="null";
};

} // end namespace mesh

#endif
