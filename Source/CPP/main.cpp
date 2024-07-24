#include "data/parameters.h"
#include "graph/fds_graph.h"
#include <cstring>
#include "log.h"

int main(int argc, char **argv) {
  auto parameters = std::make_shared<Parameters<ParameterIds::None>>();
  if (argc == 2) {
    parameters->fnInput = argv[1];
  }

  INFO("info log active");
  WARN("warn log active");
  ERROR("error log active");
  TODO("todo log active");

  FDSGraph fdsGraph;

  fdsGraph.executeGraph(true);

  fdsGraph.pushData(parameters);
  fdsGraph.finishPushingData();
  fdsGraph.waitForTermination();

  fdsGraph.createDotFile("./graphFDS.dot", hh::ColorScheme::EXECUTION,
                         hh::StructureOptions::QUEUE);

  return 0;
}
