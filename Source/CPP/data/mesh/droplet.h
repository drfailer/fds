#ifndef DROPLET_H
#define DROPLET_H

namespace mesh {

struct Droplet {
  double *AVG_DROP_DEN;  ///< Droplet mass per unit volume for a certain droplet type
  double *AVG_DROP_TMP;  ///< Average temperature for a certain droplet type
  double *AVG_DROP_RAD;  ///< Average radius for a certain droplet type
  double *AVG_DROP_AREA; ///< Average area for a certain droplet type
};

} // end namespace mesh

#endif
