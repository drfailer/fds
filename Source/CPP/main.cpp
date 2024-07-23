#include "data/parameters.h"
#include "graph/fds_graph.h"
#include <cstring>

int main(int argc, char **argv) {
  auto parameters = std::make_shared<Parameters>();
  if (argc == 2) {
    parameters->fnInput = argv[1];
  }

  FDSGraph fdsGraph;

  fdsGraph.executeGraph(true);

  fdsGraph.pushData(parameters);
  fdsGraph.finishPushingData();
  fdsGraph.waitForTermination();

  fdsGraph.createDotFile("./graphFDS.dot", hh::ColorScheme::EXECUTION,
                         hh::StructureOptions::QUEUE);

  return 0;
}
