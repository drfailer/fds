#ifndef EXCHANGE_VALUES_TASK_H
#define EXCHANGE_VALUES_TASK_H
#include "../../../../../Fortran/include/fds.h"
#include "../../../../data/parameters.h"
#include <hedgehog/hedgehog.h>

#define ExchangeValuesTaskInNb 1
#define ExchangeValuesTaskIn Parameters<ParameterIds::None>
#define ExchangeValuesTaskOut Parameters<ParameterIds::None>

class ExchangeValuesTask
    : public hh::AbstractTask<ExchangeValuesTaskInNb, ExchangeValuesTaskIn,
                              ExchangeValuesTaskOut> {
public:
  ExchangeValuesTask(size_t nbThreads)
      : hh::AbstractTask<ExchangeValuesTaskInNb, ExchangeValuesTaskIn,
                         ExchangeValuesTaskOut>("Exchange Values Task",
                                                 nbThreads) {}

  void
  execute(std::shared_ptr<Parameters<ParameterIds::None>> parameters) override {
    INFO("exchangeValues");
    exchangeValues();
    this->addResult(parameters);
  }

  std::shared_ptr<hh::AbstractTask<
      ExchangeValuesTaskInNb, ExchangeValuesTaskIn, ExchangeValuesTaskOut>>
  copy() override {
    return std::make_shared<ExchangeValuesTask>(this->numberThreads());
  }
};

#endif
