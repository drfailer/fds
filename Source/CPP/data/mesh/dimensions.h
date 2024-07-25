#ifndef DIMENSIONS_H
#define DIMENSIONS_H
#include <cstddef>

namespace mesh {

struct Dimensions {
  size_t IBAR; ///< Number of cells in the \f$ x \f$ direction, \f$ I \f$
  size_t JBAR; ///< Number of cells in the \f$ y \f$ direction, \f$ J \f$
  size_t KBAR; ///< Number of cells in the \f$ z \f$ direction, \f$ K \f$
  size_t IBM1; ///< IBAR minus 1
  size_t JBM1; ///< JBAR minus 1
  size_t KBM1; ///< KBAR minus 1
  size_t IBP1; ///< IBAR plus 1
  size_t JBP1; ///< JBAR plus 1
  size_t KBP1; ///< KBAR plus 1
  double *LES_FILTER_WIDTH; ///< Characteristic cell dimension (m)
};

}

#endif
