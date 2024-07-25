#ifndef CC_INBCF_AREA_H
#define CC_INBCF_AREA_H
#include <cstddef>

namespace mesh::cc {

struct CCINBCFArea {
   size_t NCELL=0;
   double *AINB, *NEW_AINB;
   int *SURF_INDEX;
   struct CCINBCFIJCF {
     size_t NWFACE=0;
     int *NEW_ICFINB, *NEW_JCFINB;
   } *IJCF;
};

} // end namespace mesh::cc

#endif
