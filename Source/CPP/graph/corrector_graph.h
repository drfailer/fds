#ifndef CORRECTOR_GRAPH_H
#define CORRECTOR_GRAPH_H
#include "../data/parameters.h"
#include <hedgehog/hedgehog.h>

#define CorrectorGraphInNb 1
#define CorrectorGraphIn Parameters
#define CorrectorGraphOut Parameters

class CorrectorGraph : public hh::Graph<CorrectorGraphInNb, CorrectorGraphIn,
                                        CorrectorGraphOut> {
public:
  CorrectorGraph()
      : hh::Graph<CorrectorGraphInNb, CorrectorGraphIn, CorrectorGraphOut>(
            "Corrector Graph") {}
};

#endif

