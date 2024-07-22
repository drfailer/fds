#ifndef LOOP_FDS_PARAMETERS_H
#define LOOP_FDS_PARAMETERS_H
#include <memory>
#include <string>

struct Parameters {
  Parameters() : fnInput("null", 250) {}

  Parameters(Parameters &&parameters)
      : fnInput(std::move(parameters.fnInput)) {}

  Parameters(std::shared_ptr<Parameters> parameters)
      : fnInput(parameters->fnInput) {}

  Parameters(std::shared_ptr<Parameters> &&parameters)
      : fnInput(std::move(parameters->fnInput)) {}

  std::string fnInput = "";
};

#endif // LOOP_FDS_PARAMETERS_H
