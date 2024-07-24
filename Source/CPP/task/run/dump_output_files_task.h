#ifndef DUMP_OUTPUT_FILES_TASK_H
#define DUMP_OUTPUT_FILES_TASK_H
#include "../../../Fortran/include/fds.h"
#include "../../data/parameters.h"
#include <hedgehog/hedgehog.h>

#define DumpOutputFilesTaskInNb 1
#define DumpOutputFilesTaskIn Parameters<ParameterIds::None>
#define DumpOutputFilesTaskOut Parameters<ParameterIds::None>

class DumpOutputFilesTask
    : public hh::AbstractTask<DumpOutputFilesTaskInNb, DumpOutputFilesTaskIn,
                              DumpOutputFilesTaskOut> {
public:
  DumpOutputFilesTask(size_t nbThreads)
      : hh::AbstractTask<DumpOutputFilesTaskInNb, DumpOutputFilesTaskIn,
                         DumpOutputFilesTaskOut>("Dump Output Files Task",
                                                 nbThreads) {}

  void
  execute(std::shared_ptr<Parameters<ParameterIds::None>> parameters) override {
    INFO_GRP("dump output files", MAIN_LOOP_GRP);
    dumpOutputFiles();
    this->addResult(parameters);
  }

  std::shared_ptr<hh::AbstractTask<
      DumpOutputFilesTaskInNb, DumpOutputFilesTaskIn, DumpOutputFilesTaskOut>>
  copy() override {
    return std::make_shared<DumpOutputFilesTask>(this->numberThreads());
  }
};

#endif
