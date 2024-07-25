#ifndef CC_IBM_ARRAYS_H
#define CC_IBM_ARRAYS_H
#include "cc_cutcell.h"
#include "cc_cutedge.h"
#include "cc_cutface.h"
#include "cc_edgecross.h"
#include <cstddef>

namespace mesh::cc {

struct CCIBMArrays {
   size_t N_EDGE_CROSS=0,  N_CUTEDGE_MESH=0, N_CUTFACE_MESH=0, N_CUTCELL_MESH=0;
   size_t N_BBCUTFACE_MESH=0, N_GCCUTFACE_MESH=0, N_GCCUTCELL_MESH=0;
   int FINEST_LINK_LEV=0, NUNK_F=0;
   int *VERTVAR, CCVAR;
   int *ECVAR, FCVAR;
   CCEdgeCross *EDGE_CROSS;
   CCCutedge *CUT_EDGE;
   CCCutface *CUT_FACE;
   CCCutcell *CUT_CELL;
};

} // end namespace mesh::cc

#endif
