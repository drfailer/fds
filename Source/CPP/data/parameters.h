#ifndef LOOP_FDS_PARAMETERS_H
#define LOOP_FDS_PARAMETERS_H
#include <memory>
#include <string>
#include "parameter_ids.h"
#include "../log.h"

template <ParameterIds Id>
struct Parameters {
  Parameters() : fnInput("null", 250) {}

  template <ParameterIds OtherId>
  Parameters(Parameters<OtherId> &&parameters)
      : fnInput(std::move(parameters.fnInput)) {}

  template <ParameterIds OtherId>
  Parameters(std::shared_ptr<Parameters<OtherId>> const &parameters)
      : fnInput(parameters->fnInput) {}

  template <ParameterIds OtherId>
  Parameters(std::shared_ptr<Parameters<OtherId>> &&parameters)
      : fnInput(std::move(parameters->fnInput)) {}

  std::string fnInput = "";
};

#endif // LOOP_FDS_PARAMETERS_H
