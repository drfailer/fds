/// @file main_hh.cpp
/// @brief C++ main entry point for FDS with Hedgehog dataflow graph.
///
/// This replaces the Fortran PROGRAM FDS main time-stepping loop with a
/// Hedgehog dataflow graph. All Fortran physics subroutines remain unchanged;
/// they are called through ISO_C_BINDING wrappers.
///
/// Phase 1: Sequential execution (numThreads=1) for correctness verification.
/// Phase 2: Parallel mesh processing (numThreads=nmeshes).

#include <iostream>
#include <memory>
#include <string>
#include <thread>
#include <algorithm>
#include <cstdlib>

#include "fds_fortran_interface.h"
#include "data/mesh_data.h"
#include "data/mesh_dim.h"
#include "data/termination_data.h"
#include "exchange/pack_data.h"
#include "service/fds_comm_service.h"
#include "tool/thread_budget.h"
#include "graph/fds_graph.h"

int main(int argc, char *argv[]) {
    double t_init = 0.0, dt_init = 0.0;
    int nmeshes = 0;

    // Parse command-line arguments.
    // Usage: fds_hh [options] <input_file>
    std::string inputFile;
    MeshDim meshDim;
    int pressureSubgraphOverride = -1;  // -1=auto, 0=off, 1=on

    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "--mesh-dim" && i + 3 < argc) {
            meshDim.i = std::atoi(argv[++i]);
            meshDim.j = std::atoi(argv[++i]);
            meshDim.k = std::atoi(argv[++i]);
        } else if (arg == "--pressure-subgraph" && i + 1 < argc) {
            std::string val(argv[++i]);
            if (val == "on")        pressureSubgraphOverride = 1;
            else if (val == "off")  pressureSubgraphOverride = 0;
            else if (val == "auto") pressureSubgraphOverride = -1;
            else {
                std::cerr << "[FDS-HH] Invalid --pressure-subgraph value: " << val
                          << " (expected on|off|auto)" << std::endl;
                return 1;
            }
        } else if (arg[0] != '-') {
            inputFile = arg;
        }
    }

    // Step 0: Pass the input file name to Fortran before initialization.
    // The Fortran runtime's GET_COMMAND_ARGUMENT may not work when main is C++.
    if (!inputFile.empty()) {
        fds_set_input_file(inputFile.c_str(), static_cast<int>(inputFile.size()));
    }

    // Step 0b: Set target mesh dimensions for automatic re-decomposition.
    // Each user-configured mesh will be split into sub-meshes of approximately
    // meshDim cells per dimension.
    if (meshDim.enabled()) {
        fds_set_target_mesh_dims(meshDim.i, meshDim.j, meshDim.k);
        std::cout << "[FDS-HH] Mesh re-decomposition: target dims="
                  << meshDim.i << "x" << meshDim.j << "x" << meshDim.k << std::endl;
    }

    // Step 0c: Set pressure subgraph override (before initialization so it's
    // available when the graph is built after fds_initialize_all).
    if (pressureSubgraphOverride != -1) {
        fds_set_pressure_subgraph(pressureSubgraphOverride);
        std::cout << "[FDS-HH] Pressure subgraph: "
                  << (pressureSubgraphOverride ? "on (forced)" : "off (forced)") << std::endl;
    }

    // Step 1: Run the entire Fortran initialization sequence.
    // This initializes MPI, reads the input file, sets up meshes, etc.
    fds_initialize_all(&t_init, &dt_init, &nmeshes);

    // Get the end time and mesh indices for this MPI process
    double tEnd = fds_get_t_end();
    int lower_mesh_index = fds_get_lower_mesh_index();
    int upper_mesh_index = fds_get_upper_mesh_index();
    int local_nmeshes = upper_mesh_index - lower_mesh_index + 1;

    // Clamp initial DT to T_END (matches main.f90 top-of-MAIN_LOOP logic).
    // Without this, simulations where initial DT > T_END overshoot on the first step.
    double t = t_init;
    double dt = fds_adjust_dt(t_init, dt_init);

    std::cout << "[FDS-HH] Initialization complete." << std::endl;
    std::cout << "[FDS-HH] Total meshes=" << nmeshes
              << " Local meshes=" << local_nmeshes
              << " (range: " << lower_mesh_index << "-" << upper_mesh_index << ")"
              << " t=" << t << " dt=" << dt << " tEnd=" << tEnd << std::endl;
    std::cout << "[FDS-HH] Hardware threads=" << std::thread::hardware_concurrency()
              << std::endl;

    // Step 2: Build the Hedgehog dataflow graph.
    // Compute thread budget from hardware capabilities
    int hwThreads = static_cast<int>(std::thread::hardware_concurrency());
    if (hwThreads <= 0) hwThreads = local_nmeshes;  // fallback
    auto budget = ThreadBudget::compute(hwThreads, local_nmeshes);
    budget.print(std::cout);

    // Initialize communicator service (reuses FDS's already-initialized MPI)
    FDSMPIService commService(true);
    std::cout << "[FDS-HH] CommService: rank=" << commService.rank()
              << " nbProcesses=" << commService.nbProcesses() << std::endl;

    // Build mesh dependency graph for exchange overlap
    auto depGraph = std::make_shared<MeshDependencyGraph>(
        lower_mesh_index, upper_mesh_index);
    std::cout << "[FDS-HH] " << depGraph->dump();

    // Initialize max slab sizes for MPI receive buffers
    int maxSlab1 = fds_exchange_max_slab_size(1);
    PackData<MeshState::MeshExch1>::maxSlabSize_ = maxSlab1;

    int maxSlab3 = fds_exchange_max_slab_size(3);
    PackData<MeshState::MeshExch3>::maxSlabSize_ = maxSlab3;

    int maxSlab4 = fds_exchange_max_slab_size(4);
    PackData<MeshState::MeshExch4>::maxSlabSize_ = maxSlab4;

    int maxSlab5 = fds_exchange_max_slab_size(5);
    PackData<MeshState::PreSolveExch>::maxSlabSize_ = maxSlab5;
    PackData<MeshState::PostSolveExch>::maxSlabSize_ = maxSlab5;
    PackData<MeshState::MeshExch5>::maxSlabSize_ = maxSlab5;

    int maxSlab7 = fds_exchange_max_slab_size(7);
    PackData<MeshState::MeshExch7>::maxSlabSize_ = maxSlab7;

    // Only pass commService for multi-process runs
    hh::comm::CommService *commPtr =
        (commService.nbProcesses() > 1) ? &commService : nullptr;

    auto graph = buildFDSGraph(local_nmeshes, t, dt, tEnd, budget,
                               depGraph, commPtr);

    // Step 3: Execute the graph (spawns threads).
    graph->executeGraph();

    // Step 4: Push initial MeshData<Init> tokens (one per LOCAL mesh) into the graph.
    // Set PREDICTOR=TRUE and FIRST_PASS=TRUE for the first time step.
    // TimestepState collects these, runs INSERT_ALL_PARTICLES, then emits MeshData<>
    // to the Predictor subgraph.
    // IMPORTANT: Only push tokens for meshes owned by this MPI process
    fds_set_predictor(1);
    fds_set_first_pass(1);
    fds_set_icyc(1);

    for (int nm = lower_mesh_index; nm <= upper_mesh_index; ++nm) {
        auto md = std::make_shared<MeshData<MeshState::Init>>(nm, t, dt, 0);
        graph->pushData(md);
    }

    // Step 5: Wait for graph output (simulation complete).
    // getBlockingResult() blocks until the graph produces BarrierData output,
    // which happens when TimestepState emits done=true BarrierData.
    std::cout << "[FDS-HH] Waiting for graph termination..." << std::endl;

    graph->getBlockingResult();

    graph->pushData(std::make_shared<TerminationData>());
    graph->finishPushingData();

    commService.terminate();

    graph->waitForTermination();

    std::cout << "[FDS-HH] Graph terminated." << std::endl;

    // Close all per-mesh output files and flush remaining I/O buffers
    fds_close_all_mesh_output_files();
    fds_flush_output_files();

    // Step 7: Generate dot file for visualization
    // Single process: fds_hh_graph.dot; MPI: fds_hh_graph_{rank}.dot
    std::string dotFile = (commService.nbProcesses() > 1)
        ? "fds_hh_graph_" + std::to_string(commService.rank()) + ".dot"
        : "fds_hh_graph.dot";
    graph->createDotFile(
        dotFile,
        hh::ColorScheme::EXECUTION,
        hh::StructureOptions::QUEUE);

    std::cout << "[FDS-HH] Graph dot file written to " << dotFile << std::endl;

    // Step 8: Finalize FDS (deallocate solvers, MPI_Finalize, etc.)
    fds_finalize_all(t, dt);

    return 0;
}
