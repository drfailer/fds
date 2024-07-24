#ifndef PROCESS_EXTERN_VARIABLE_H
#define PROCESS_EXTERN_VARIABLE_H
#include "../../../../Fortran/include/fds.h"
#include "../../../data/parameters.h"
#include <hedgehog/hedgehog.h>

#define PEVTaskInNb 1
#define PEVTaskIn Parameters<ParameterIds::None>
#define PEVTaskOut Parameters<ParameterIds::None>

class ProcessExternVariablesTask
    : public hh::AbstractTask<PEVTaskInNb, PEVTaskIn, PEVTaskOut> {
public:
  ProcessExternVariablesTask(size_t nbThreads)
      : hh::AbstractTask<PEVTaskInNb, PEVTaskIn, PEVTaskOut>(
            "Process External Variables", nbThreads) {}

  void
  execute(std::shared_ptr<Parameters<ParameterIds::None>> parameters) override {
    processExternalVariables();
    this->addResult(parameters);
  }

  std::shared_ptr<hh::AbstractTask<PEVTaskInNb, PEVTaskIn, PEVTaskOut>>
  copy() override {
    return std::make_shared<ProcessExternVariablesTask>(this->numberThreads());
  }
};

#endif
