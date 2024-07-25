#ifndef MESH_PATCH_H
#define MESH_PATCH_H

namespace mesh {

struct Patch {
  int I1, I2, J1, J2, K1, K2, IG1, IG2, JG1, JG2, KG1, KG2, IOR, OBST_INDEX;
};

} // end namespace mesh

#endif
