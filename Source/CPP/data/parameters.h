#ifndef LOOP_FDS_PARAMETERS_H
#define LOOP_FDS_PARAMETERS_H

#include <cstddef>
#include <memory>

struct Parameters {
  Parameters(size_t ibp1, size_t jpb1, size_t kbp1, size_t iBlockSize, size_t jBlockSize,
             size_t kBlockSize) :
          ibp1(ibp1), jbp1(jpb1), kbp1(kbp1),
          iBlockSize(iBlockSize),
          jBlockSize(jBlockSize),
          kBlockSize(kBlockSize),
          // matrix: (0:IBPX) where IBPX = XBAR + 1
          nbBlocksI((ibp1 + 1) / iBlockSize + (((ibp1 + 1) % iBlockSize == 0)? 0 : 1)),
          nbBlocksJ((jpb1 + 1) / jBlockSize + (((jbp1 + 1) % jBlockSize == 0)? 0 : 1)),
          nbBlocksK((kbp1 + 1) / kBlockSize + (((kbp1 + 1) % kBlockSize == 0)? 0 : 1)) {}

  Parameters(Parameters const &other):
    Parameters(other.ibp1, other.jbp1, other.kbp1, other.iBlockSize, other.jBlockSize, other.kBlockSize) {}

  explicit Parameters(std::shared_ptr<Parameters> const &other):
    Parameters(other->ibp1, other->jbp1, other->kbp1, other->iBlockSize, other->jBlockSize, other->kBlockSize) {}

  explicit Parameters(std::shared_ptr<Parameters> &&other) noexcept :
    Parameters(other->ibp1, other->jbp1, other->kbp1, other->iBlockSize, other->jBlockSize, other->kBlockSize) {}

  size_t ibp1 = 0;
  size_t jbp1 = 0;
  size_t kbp1 = 0;
  size_t nedge = 0;
  size_t iBlockSize = 0;
  size_t jBlockSize = 0;
  size_t kBlockSize = 0;
  size_t nbBlocksI = 0;
  size_t nbBlocksJ = 0;
  size_t nbBlocksK = 0;
};

#endif //LOOP_FDS_PARAMETERS_H
