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

#include "fds_fortran_interface.h"
#include "data/mesh_data.h"
#include "graph/fds_graph.h"

int main(int argc, char *argv[]) {
    double t_init = 0.0, dt_init = 0.0;
    int nmeshes = 0;

    // Step 0: Pass the input file name to Fortran before initialization.
    // The Fortran runtime's GET_COMMAND_ARGUMENT may not work when main is C++.
    if (argc > 1) {
        std::string inputFile(argv[1]);
        fds_set_input_file(inputFile.c_str(), static_cast<int>(inputFile.size()));
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
    // kernelThreads controls inter-mesh parallelism. Block decomposition of velocity/momentum
    // kernels is available but provides negligible benefit (<0.2% of runtime) — the heavy kernels
    // (WallBC, DivPart1, DivSetup, VelocityBCEdges) are mesh-level due to wall/zone loops.
    // Using hw_threads would enable intra-mesh block decomposition but adds overhead without
    // measurable speedup. Keep kernelThreads = nmeshes for now.
    size_t kernelThreads = local_nmeshes;
    std::cout << "[FDS-HH] Kernel threads=" << kernelThreads << std::endl;
    auto graph = buildFDSGraph(local_nmeshes, t, dt, tEnd, kernelThreads);

    // Step 3: Execute the graph (spawns threads).
    graph->executeGraph();

    // Step 4: Push initial MeshData tokens (one per LOCAL mesh) into the graph.
    // Set PREDICTOR=TRUE and FIRST_PASS=TRUE for the first time step
    // IMPORTANT: Only push tokens for meshes owned by this MPI process
    fds_set_predictor(1);
    fds_set_first_pass(1);
    fds_set_icyc(1);

    for (int nm = lower_mesh_index; nm <= upper_mesh_index; ++nm) {
        auto md = std::make_shared<MeshData>(nm, t, dt, 0);  // phase=0 (predictor)
        graph->pushData(md);
    }

    // Step 5: Signal that no more data will be pushed from outside.
    graph->finishPushingData();

    // Step 6: Wait for the graph to complete.
    // TimestepLoopStateManager::canTerminate() breaks the main cycle when done.
    std::cout << "[FDS-HH] Waiting for graph termination..." << std::endl;

    graph->waitForTermination();

    std::cout << "[FDS-HH] Graph terminated." << std::endl;

    // Flush Fortran I/O buffers to ensure all outputs are written to disk
    fds_flush_output_files();

    // Step 7: Generate dot file for visualization
    graph->createDotFile(
        "fds_hh_graph.dot",
        hh::ColorScheme::EXECUTION,
        hh::StructureOptions::QUEUE);

    std::cout << "[FDS-HH] Graph dot file written to fds_hh_graph.dot" << std::endl;

    // Step 8: Finalize FDS (deallocate solvers, MPI_Finalize, etc.)
    fds_finalize_all(t, dt);

    return 0;
}
