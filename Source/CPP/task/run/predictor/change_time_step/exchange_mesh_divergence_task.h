#ifndef EXCHANGE_MESH_DIVERGENCE_TASK_H
#define EXCHANGE_MESH_DIVERGENCE_TASK_H
#include "../../../../../Fortran/include/fds.h"
#include "../../../../data/parameters.h"
#include <hedgehog/hedgehog.h>

#define ExchangeMeshDivergenceTaskInNb 1
#define ExchangeMeshDivergenceTaskIn Parameters<ParameterIds::None>
#define ExchangeMeshDivergenceTaskOut Parameters<ParameterIds::None>

class ExchangeMeshDivergenceTask
    : public hh::AbstractTask<ExchangeMeshDivergenceTaskInNb,
                              ExchangeMeshDivergenceTaskIn,
                              ExchangeMeshDivergenceTaskOut> {
public:
  ExchangeMeshDivergenceTask(size_t nbThreads)
      : hh::AbstractTask<ExchangeMeshDivergenceTaskInNb,
                         ExchangeMeshDivergenceTaskIn,
                         ExchangeMeshDivergenceTaskOut>(
            "Exchange Mesh Divergence Task", nbThreads) {}

  void
  execute(std::shared_ptr<Parameters<ParameterIds::None>> parameters) override {
    INFO_GRP("exchangeMeshDivergence", TIME_STEP_LOOP_INFO_GRP);
    exchangeMeshDivergence();
    this->addResult(parameters);
  }

  std::shared_ptr<hh::AbstractTask<ExchangeMeshDivergenceTaskInNb,
                                   ExchangeMeshDivergenceTaskIn,
                                   ExchangeMeshDivergenceTaskOut>>
  copy() override {
    return std::make_shared<ExchangeMeshDivergenceTask>(this->numberThreads());
  }
};

#endif
