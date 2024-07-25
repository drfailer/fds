#ifndef BOUNDARY_RADIA_H
#define BOUNDARY_RADIA_H

namespace mesh::boundaries {

using BandType = double *;

struct BoundaryRadia {
   double *IL;    ///< (1:NSB) Radiance (W/m2/sr); output only
   BandType BAND; ///< (1:NSB) Radiation wavelength band
};

} // end namespace mesh::boundaries

#endif
