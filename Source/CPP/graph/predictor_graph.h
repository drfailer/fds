#ifndef PREDICTOR_GRAPH_H
#define PREDICTOR_GRAPH_H
#include "../data/parameters.h"
#include <hedgehog/hedgehog.h>

#define PredictorGraphInNb 1
#define PredictorGraphIn Parameters
#define PredictorGraphOut Parameters

class PredictorGraph : public hh::Graph<PredictorGraphInNb, PredictorGraphIn,
                                        PredictorGraphOut> {
public:
  PredictorGraph()
      : hh::Graph<PredictorGraphInNb, PredictorGraphIn, PredictorGraphOut>(
            "Predictor Graph") {}
};

#endif
