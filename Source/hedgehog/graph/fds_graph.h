#ifndef FDS_GRAPH_H
#define FDS_GRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../data/velocity_corrector_data.h"
#include "../data/velocity_predictor_data.h"
#include "../data/divergence_part2_data.h"
#include "../data/corr_step1_data.h"
#include "../data/density_pred_data.h"
#include "../data/corr_div_part1_data.h"
#include "../data/div_setup_data.h"
#include "../data/pred_step1_data.h"
#include "../data/corr_condens_data.h"
#include "../data/pred_wall_div_data.h"
#include "../data/corr_particle_data.h"
#include "../task/predictor_tasks.h"
#include "../task/corrector_tasks.h"
#include "../task/barrier_tasks.h"
#include "../task/velocity_predictor_kernel_task.h"
#include "../task/velocity_corrector_kernel_task.h"
#include "../task/divergence_part2_kernel_task.h"
#include "../task/corr_step1_kernel_task.h"
#include "../task/density_pred_kernel_task.h"
#include "../task/corr_div_part1_kernel_task.h"
#include "../task/div_setup_kernel_task.h"
#include "../task/pred_step1_kernel_task.h"
#include "../task/corr_condens_kernel_task.h"
#include "../task/pred_wall_div_kernel_task.h"
#include "../task/corr_particle_kernel_task.h"
#include "../state/collector_state.h"
#include "../state/mesh_barrier_state.h"
#include "../state/timestep_state.h"
#include "../state/velocity_predictor_state.h"
#include "../state/velocity_corrector_state.h"
#include "../state/divergence_part2_state.h"
#include "../state/corr_step1_state.h"
#include "../state/density_pred_state.h"
#include "../state/corr_div_part1_state.h"
#include "../state/div_setup_state.h"
#include "../state/pred_step1_state.h"
#include "../state/corr_condens_state.h"
#include "../state/pred_wall_div_state.h"
#include "../state/corr_particle_state.h"
#include "change_timestep_subgraph.h"

/// Build the FDS Hedgehog dataflow graph.
///
/// The graph implements the FDS time-stepping loop as a dataflow pipeline:
///   Predictor tasks -> barriers -> Corrector tasks -> barriers -> cycle back
///
/// Barrier pattern: each synchronization point is a CollectorState (pure
/// data-flow: collects N MeshData -> emits 1 BarrierData) followed by a
/// barrier task (computation: receives BarrierData -> emits N MeshData).
///
/// @param nmeshes Number of meshes
/// @param t Initial simulation time
/// @param dt Initial time step
/// @param tEnd End time
/// @param kernelThreads Number of threads for parallel kernel tasks (1 for sequential)
/// @return Shared pointer to the constructed graph
inline auto buildFDSGraph(int nmeshes, double t, double dt, double tEnd, size_t kernelThreads) {

    using GraphType = hh::Graph<1, MeshData, BarrierData>;
    auto graph = std::make_shared<GraphType>("FDS Hedgehog Graph");

    // --- Create predictor tasks (all sequential) ---
    // NOTE: predStep1 replaced by pred step 1 sub-graph (see below)
    // NOTE: densityPred replaced by density predictor sub-graph (see below)
    // NOTE: predDivSetup replaced by div setup sub-graph (see below)
    // NOTE: predWallDiv replaced by pred wall div sub-graph (see below)
    // NOTE: divPart2Pred replaced by divergence part 2 sub-graph (see below)
    // NOTE: velPredictor replaced by velocity predictor sub-graph (see below)
    auto predFinal       = std::make_shared<PredFinalTask>(1);

    // --- Create corrector tasks (all sequential) ---
    // NOTE: corrStep1 replaced by corrector step 1 sub-graph (see below)
    // NOTE: corrDivSetup replaced by div setup sub-graph (see below)
    // NOTE: corrCondens replaced by corr condens sub-graph (see below)
    // NOTE: corrParticle replaced by corr particle sub-graph (see below)
    auto corrWallBC      = std::make_shared<CorrWallBCTask>(1);
    auto corrRadiation   = std::make_shared<CorrRadiationTask>(1);
    // NOTE: corrDivPart1 replaced by corr div part 1 sub-graph (see below)
    // NOTE: corrDivPart2 replaced by divergence part 2 sub-graph (see below)
    // NOTE: corrVelocity replaced by velocity corrector sub-graph (see below)
    auto corrFinal       = std::make_shared<CorrFinalTask>(1);

    // --- Create velocity predictor sub-graph components ---
    // Only the kernel task is parallelized; orchestrator and collector are always sequential
    auto velPredOrchSM = std::make_shared<hh::StateManager<1, MeshData, VelocityPredictorWork>>(
        std::make_shared<VelocityPredictorOrchestrator>(nmeshes), "VelPredOrch");
    auto velPredKernelTask = std::make_shared<VelocityPredictorKernelTask>(kernelThreads);
    auto velPredCollectorSM = std::make_shared<hh::StateManager<1, VelocityPredictorWork, MeshData>>(
        std::make_shared<VelocityPredictorCollector>(nmeshes), "VelPredCollector");

    // --- Create velocity corrector sub-graph components ---
    // Only the kernel task is parallelized; orchestrator and collector are always sequential
    auto velCorrOrchSM = std::make_shared<hh::StateManager<1, MeshData, VelocityCorrectorWork>>(
        std::make_shared<VelocityCorrectorOrchestrator>(nmeshes), "VelCorrOrch");
    auto velCorrKernelTask = std::make_shared<VelocityCorrectorKernelTask>(kernelThreads);
    auto velCorrCollectorSM = std::make_shared<hh::StateManager<1, VelocityCorrectorWork, MeshData>>(
        std::make_shared<VelocityCorrectorCollector>(nmeshes), "VelCorrCollector");

    // --- Create divergence part 2 sub-graph components (predictor instance) ---
    auto predDivP2OrchSM = std::make_shared<hh::StateManager<1, MeshData, DivergencePart2Work>>(
        std::make_shared<DivergencePart2Orchestrator>(nmeshes), "PredDivP2Orch");
    auto predDivP2KernelTask = std::make_shared<DivergencePart2KernelTask>(kernelThreads);
    auto predDivP2CollectorSM = std::make_shared<hh::StateManager<1, DivergencePart2Work, MeshData>>(
        std::make_shared<DivergencePart2Collector>(nmeshes), "PredDivP2Collector");

    // --- Create divergence part 2 sub-graph components (corrector instance) ---
    auto corrDivP2OrchSM = std::make_shared<hh::StateManager<1, MeshData, DivergencePart2Work>>(
        std::make_shared<DivergencePart2Orchestrator>(nmeshes), "CorrDivP2Orch");
    auto corrDivP2KernelTask = std::make_shared<DivergencePart2KernelTask>(kernelThreads);
    auto corrDivP2CollectorSM = std::make_shared<hh::StateManager<1, DivergencePart2Work, MeshData>>(
        std::make_shared<DivergencePart2Collector>(nmeshes), "CorrDivP2Collector");

    // --- Create corrector step 1 sub-graph components ---
    // Bundles COMPUTE_VISCOSITY + MASS_FINITE_DIFFERENCES + DENSITY kernels
    auto corrStep1OrchSM = std::make_shared<hh::StateManager<1, MeshData, CorrStep1Work>>(
        std::make_shared<CorrStep1Orchestrator>(nmeshes), "CorrStep1Orch");
    auto corrStep1KernelTask = std::make_shared<CorrStep1KernelTask>(kernelThreads);
    auto corrStep1CollectorSM = std::make_shared<hh::StateManager<1, CorrStep1Work, MeshData>>(
        std::make_shared<CorrStep1Collector>(nmeshes), "CorrStep1Collector");

    // --- Create density predictor sub-graph components ---
    auto densPredOrchSM = std::make_shared<hh::StateManager<1, MeshData, DensityPredWork>>(
        std::make_shared<DensityPredOrchestrator>(nmeshes), "DensPredOrch");
    auto densPredKernelTask = std::make_shared<DensityPredKernelTask>(kernelThreads);
    auto densPredCollectorSM = std::make_shared<hh::StateManager<1, DensityPredWork, MeshData>>(
        std::make_shared<DensityPredCollector>(nmeshes), "DensPredCollector");

    // --- Create corrector divergence part 1 sub-graph components ---
    // Sequential COMBUSTION_BC (cross-mesh) in orchestrator, parallel DIVERGENCE_PART_1_KERNEL
    auto corrDivP1OrchSM = std::make_shared<hh::StateManager<1, MeshData, CorrDivPart1Work>>(
        std::make_shared<CorrDivPart1Orchestrator>(nmeshes), "CorrDivP1Orch");
    auto corrDivP1KernelTask = std::make_shared<CorrDivPart1KernelTask>(kernelThreads);
    auto corrDivP1CollectorSM = std::make_shared<hh::StateManager<1, CorrDivPart1Work, MeshData>>(
        std::make_shared<CorrDivPart1Collector>(nmeshes), "CorrDivP1Collector");

    // --- Create predictor div setup sub-graph components ---
    // Sequential VISCOSITY_BC (cross-mesh) in orchestrator, parallel VELOCITY_FLUX_KERNEL
    auto predDivSetupOrchSM = std::make_shared<hh::StateManager<1, MeshData, DivSetupWork>>(
        std::make_shared<PredDivSetupOrchestrator>(nmeshes), "PredDivSetupOrch");
    auto predDivSetupKernelTask = std::make_shared<DivSetupKernelTask>(kernelThreads);
    auto predDivSetupCollectorSM = std::make_shared<hh::StateManager<1, DivSetupWork, MeshData>>(
        std::make_shared<DivSetupCollector>(nmeshes), "PredDivSetupCollector");

    // --- Create corrector div setup sub-graph components ---
    // Sequential VISCOSITY_BC + AGGLOMERATION in orchestrator, parallel VELOCITY_FLUX_KERNEL
    auto corrDivSetupOrchSM = std::make_shared<hh::StateManager<1, MeshData, DivSetupWork>>(
        std::make_shared<CorrDivSetupOrchestrator>(nmeshes), "CorrDivSetupOrch");
    auto corrDivSetupKernelTask = std::make_shared<DivSetupKernelTask>(kernelThreads);
    auto corrDivSetupCollectorSM = std::make_shared<hh::StateManager<1, DivSetupWork, MeshData>>(
        std::make_shared<DivSetupCollector>(nmeshes), "CorrDivSetupCollector");

    // --- Create predictor step 1 sub-graph components ---
    // Sequential INSERT_ALL_PARTICLES in orchestrator, parallel COMPUTE_VISCOSITY + MASS_FD kernels
    auto predStep1OrchSM = std::make_shared<hh::StateManager<1, MeshData, PredStep1Work>>(
        std::make_shared<PredStep1Orchestrator>(nmeshes), "PredStep1Orch");
    auto predStep1KernelTask = std::make_shared<PredStep1KernelTask>(kernelThreads);
    auto predStep1CollectorSM = std::make_shared<hh::StateManager<1, PredStep1Work, MeshData>>(
        std::make_shared<PredStep1Collector>(nmeshes), "PredStep1Collector");

    // --- Create corrector condensation sub-graph components ---
    // Pure kernel (Pattern A): CONDENSATION_EVAPORATION_KERNEL
    auto corrCondensOrchSM = std::make_shared<hh::StateManager<1, MeshData, CorrCondensWork>>(
        std::make_shared<CorrCondensOrchestrator>(nmeshes), "CorrCondensOrch");
    auto corrCondensKernelTask = std::make_shared<CorrCondensKernelTask>(kernelThreads);
    auto corrCondensCollectorSM = std::make_shared<hh::StateManager<1, CorrCondensWork, MeshData>>(
        std::make_shared<CorrCondensCollector>(nmeshes), "CorrCondensCollector");

    // --- Create predictor wall+div sub-graph components ---
    // Sequential WALL_BC in orchestrator, parallel PARTICLE_MOMENTUM + DIVERGENCE_PART_1 kernels
    auto predWallDivOrchSM = std::make_shared<hh::StateManager<1, MeshData, PredWallDivWork>>(
        std::make_shared<PredWallDivOrchestrator>(nmeshes), "PredWallDivOrch");
    auto predWallDivKernelTask = std::make_shared<PredWallDivKernelTask>(kernelThreads);
    auto predWallDivCollectorSM = std::make_shared<hh::StateManager<1, PredWallDivWork, MeshData>>(
        std::make_shared<PredWallDivCollector>(nmeshes), "PredWallDivCollector");

    // --- Create corrector particle sub-graph components ---
    // Sequential PARTICLE_MASS_ENERGY + MOVE_PARTICLES in orchestrator, parallel PARTICLE_MOMENTUM_KERNEL
    auto corrParticleOrchSM = std::make_shared<hh::StateManager<1, MeshData, CorrParticleWork>>(
        std::make_shared<CorrParticleOrchestrator>(nmeshes), "CorrParticleOrch");
    auto corrParticleKernelTask = std::make_shared<CorrParticleKernelTask>(kernelThreads);
    auto corrParticleCollectorSM = std::make_shared<hh::StateManager<1, CorrParticleWork, MeshData>>(
        std::make_shared<CorrParticleCollector>(nmeshes), "CorrParticleCollector");

    // --- Create barrier collector state managers + barrier tasks ---

    // Predictor: MESH_EXCHANGE(1) after density
    auto collector1SM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "Collector(1)");
    auto meshExchange1 = std::make_shared<MeshExchangeTask>(1);

    // Predictor: HVAC barrier
    auto predHvacCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "PredHvacCollector");
    auto predHvacTask = std::make_shared<HvacTask>(1);  // first=1

    // Predictor: INITIALIZE_DIVERGENCE_INTEGRALS
    auto predInitDivCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "PredInitDivCollector");
    auto predInitDivTask = std::make_shared<InitDivIntegralsTask>();

    // Predictor: EXCHANGE_DIVERGENCE_INFO
    auto predDivCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "PredDivCollector");
    auto predDivExchangeTask = std::make_shared<DivergenceExchangeTask>(/*corrector=*/false);

    // Predictor: PRESSURE_ITERATION
    auto predPressureCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "PredPressureCollector");
    auto predPressureTask = std::make_shared<PressureIterationTask>(/*predictor=*/true);

    // Predictor: CHANGE_TIME_STEP_LOOP (subgraph with cycle for CFL retry)
    auto changeTimeStepCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "ChangeTimeStepCollector");
    auto changeTimeStepSubgraph = buildChangeTimeStepSubgraph(tEnd);

    // Predictor: MESH_EXCHANGE(3) after CFL check
    auto collector3SM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "Collector(3)");
    auto meshExchange3 = std::make_shared<MeshExchangeTask>(3);

    // Predictor->Corrector phase transition
    auto phaseTransCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "PhaseTransCollector");
    auto phaseTransTask = std::make_shared<PhaseTransitionTask>();

    // Corrector: MESH_EXCHANGE(4)
    auto collector4SM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "Collector(4)");
    auto meshExchange4 = std::make_shared<MeshExchangeTask>(4);

    // Corrector: COMBUSTION barrier
    auto combustionCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "CombustionCollector");
    auto combustionTask = std::make_shared<CombustionTask>();

    // Corrector: HVAC barrier
    auto corrHvacCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "CorrHvacCollector");
    auto corrHvacTask = std::make_shared<HvacTask>(1);  // first=1

    // Corrector: MESH_EXCHANGE(7) particles
    auto collector7SM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "Collector(7)");
    auto meshExchange7 = std::make_shared<MeshExchangeTask>(7);

    // Corrector: MESH_EXCHANGE(6) after wall BC
    auto collector6aSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "Collector(6a)");
    auto meshExchange6a = std::make_shared<MeshExchangeTask>(6);

    // Corrector: MESH_EXCHANGE(2) after radiation
    auto collector2SM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "Collector(2)");
    auto meshExchange2 = std::make_shared<MeshExchangeTask>(2);

    // Corrector: INITIALIZE_DIVERGENCE_INTEGRALS
    auto corrInitDivCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "CorrInitDivCollector");
    auto corrInitDivTask = std::make_shared<InitDivIntegralsTask>();

    // Corrector: EXCHANGE_DIVERGENCE_INFO + RTE
    auto corrDivCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "CorrDivCollector");
    auto corrDivExchangeTask = std::make_shared<DivergenceExchangeTask>(/*corrector=*/true);

    // Corrector: PRESSURE_ITERATION
    auto corrPressureCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "CorrPressureCollector");
    auto corrPressureTask = std::make_shared<PressureIterationTask>(/*predictor=*/false);

    // Corrector: MESH_EXCHANGE(6) after velocity
    auto collector6bSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "Collector(6b)");
    auto meshExchange6b = std::make_shared<MeshExchangeTask>(6);

    // Passthrough barriers (pure data flow, no computation — kept as direct states)
    auto predStep1BarrierSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
        std::make_shared<PassthroughBarrierState>(nmeshes), "PredStep1Barrier");
    auto corrCondensBarrierSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
        std::make_shared<PassthroughBarrierState>(nmeshes), "CorrCondensBarrier");

    // Timestep: collector -> computation task -> loop state (cycle management)
    auto timestepCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "TimestepCollector");
    auto timestepTask = std::make_shared<TimestepTask>(tEnd);
    auto timestepLoopSM = std::make_shared<TimestepLoopStateManager>(
        std::make_shared<TimestepLoopState>(), "TimestepLoop");

    // Termination sink: receives final BarrierData and outputs to graph
    auto terminationSinkSM = std::make_shared<hh::StateManager<1, BarrierData, BarrierData>>(
        std::make_shared<TerminationSinkState>(), "TerminationSink");

    // --- Wire the graph ---

    // Graph input goes to predStep1 orchestrator
    graph->inputs(predStep1OrchSM);

    // Predictor step 1 sub-graph (sequential INSERT_ALL_PARTICLES + parallel kernels)
    graph->edges(predStep1OrchSM, predStep1KernelTask);       // Orchestrator -> Kernel (parallel)
    graph->edges(predStep1KernelTask, predStep1CollectorSM);  // Kernel -> Collector
    graph->edges(predStep1CollectorSM, predStep1BarrierSM);   // Collector -> Passthrough barrier
    // Density predictor sub-graph (parallel multi-mesh execution)
    graph->edges(predStep1BarrierSM, densPredOrchSM);         // Barrier -> DensPred Orchestrator
    graph->edges(densPredOrchSM, densPredKernelTask);          // Orchestrator -> Kernel (parallel)
    graph->edges(densPredKernelTask, densPredCollectorSM);     // Kernel -> Collector
    graph->edges(densPredCollectorSM, collector1SM);           // Collector -> MESH_EXCHANGE(1)
    graph->edges(collector1SM, meshExchange1);                // Do MESH_EXCHANGE(1)
    // Predictor div setup sub-graph (sequential VISCOSITY_BC + parallel VELOCITY_FLUX_KERNEL)
    graph->edges(meshExchange1, predDivSetupOrchSM);           // MeshExch -> DivSetup Orchestrator
    graph->edges(predDivSetupOrchSM, predDivSetupKernelTask);  // Orchestrator -> Kernel (parallel)
    graph->edges(predDivSetupKernelTask, predDivSetupCollectorSM); // Kernel -> Collector
    graph->edges(predDivSetupCollectorSM, predHvacCollectorSM); // Collector -> HVAC
    graph->edges(predHvacCollectorSM, predHvacTask);          // Do HVAC_CALC
    graph->edges(predHvacTask, predInitDivCollectorSM);       // Collect for INIT_DIV
    graph->edges(predInitDivCollectorSM, predInitDivTask);    // Do INIT_DIV_INTEGRALS
    // Predictor wall+div sub-graph (sequential WALL_BC + parallel PARTICLE_MOMENTUM + DIV_PART_1 kernels)
    graph->edges(predInitDivTask, predWallDivOrchSM);          // InitDiv -> PredWallDiv Orchestrator
    graph->edges(predWallDivOrchSM, predWallDivKernelTask);    // Orchestrator -> Kernel (parallel)
    graph->edges(predWallDivKernelTask, predWallDivCollectorSM); // Kernel -> Collector
    graph->edges(predWallDivCollectorSM, predDivCollectorSM);  // Collector -> DIV_EXCHANGE
    graph->edges(predDivCollectorSM, predDivExchangeTask);    // Do EXCHANGE_DIV_INFO
    // Divergence part 2 sub-graph (predictor, parallel multi-mesh execution)
    graph->edges(predDivExchangeTask, predDivP2OrchSM);       // DivExchange -> DivP2 Orchestrator
    graph->edges(predDivP2OrchSM, predDivP2KernelTask);       // Orchestrator -> Kernel (parallel)
    graph->edges(predDivP2KernelTask, predDivP2CollectorSM);  // Kernel -> Collector
    graph->edges(predDivP2CollectorSM, predPressureCollectorSM); // Collector -> PRESSURE
    graph->edges(predPressureCollectorSM, predPressureTask);  // Do PRESSURE_ITERATION
    // Velocity predictor sub-graph (parallel multi-mesh execution)
    graph->edges(predPressureTask, velPredOrchSM);            // Pressure → VelPredOrchestrator
    graph->edges(velPredOrchSM, velPredKernelTask);           // Orchestrator → Kernel (parallel)
    graph->edges(velPredKernelTask, velPredCollectorSM);      // Kernel → Collector
    graph->edges(velPredCollectorSM, changeTimeStepCollectorSM); // Collector → CFL check
    graph->edges(changeTimeStepCollectorSM, changeTimeStepSubgraph);  // Do CHANGE_TIME_STEP_LOOP
    graph->edges(changeTimeStepSubgraph, collector3SM);               // Collect for MESH_EXCHANGE(3)
    graph->edges(collector3SM, meshExchange3);                 // Do MESH_EXCHANGE(3)
    graph->edges(meshExchange3, predFinal);
    graph->edges(predFinal, phaseTransCollectorSM);           // Collect for phase transition
    graph->edges(phaseTransCollectorSM, phaseTransTask);      // Do phase transition
    // Corrector step 1 sub-graph (parallel multi-mesh execution)
    graph->edges(phaseTransTask, corrStep1OrchSM);            // PhaseTrans -> CorrStep1 Orchestrator
    graph->edges(corrStep1OrchSM, corrStep1KernelTask);       // Orchestrator -> Kernel (parallel)
    graph->edges(corrStep1KernelTask, corrStep1CollectorSM);  // Kernel -> Collector

    // Corrector pipeline
    graph->edges(corrStep1CollectorSM, collector4SM);         // Collect for MESH_EXCHANGE(4)
    graph->edges(collector4SM, meshExchange4);                // Do MESH_EXCHANGE(4)
    // Corrector div setup sub-graph (sequential VISCOSITY_BC + AGGLOMERATION + parallel VELOCITY_FLUX_KERNEL)
    graph->edges(meshExchange4, corrDivSetupOrchSM);           // MeshExch -> DivSetup Orchestrator
    graph->edges(corrDivSetupOrchSM, corrDivSetupKernelTask);  // Orchestrator -> Kernel (parallel)
    graph->edges(corrDivSetupKernelTask, corrDivSetupCollectorSM); // Kernel -> Collector
    graph->edges(corrDivSetupCollectorSM, combustionCollectorSM); // Collector -> COMBUSTION
    graph->edges(combustionCollectorSM, combustionTask);      // Do COMBUSTION
    graph->edges(combustionTask, corrHvacCollectorSM);        // Collect for HVAC
    graph->edges(corrHvacCollectorSM, corrHvacTask);          // Do HVAC_CALC
    // Corrector condensation sub-graph (parallel CONDENSATION_EVAPORATION_KERNEL)
    graph->edges(corrHvacTask, corrCondensOrchSM);             // HVAC -> CorrCondens Orchestrator
    graph->edges(corrCondensOrchSM, corrCondensKernelTask);    // Orchestrator -> Kernel (parallel)
    graph->edges(corrCondensKernelTask, corrCondensCollectorSM); // Kernel -> Collector
    graph->edges(corrCondensCollectorSM, corrCondensBarrierSM); // Collector -> Passthrough barrier
    // Corrector particle sub-graph (sequential MASS_ENERGY + MOVE + parallel MOMENTUM kernel)
    graph->edges(corrCondensBarrierSM, corrParticleOrchSM);    // Barrier -> CorrParticle Orchestrator
    graph->edges(corrParticleOrchSM, corrParticleKernelTask);  // Orchestrator -> Kernel (parallel)
    graph->edges(corrParticleKernelTask, corrParticleCollectorSM); // Kernel -> Collector
    graph->edges(corrParticleCollectorSM, collector7SM);       // Collector -> MESH_EXCHANGE(7)
    graph->edges(collector7SM, meshExchange7);                // Do MESH_EXCHANGE(7)
    graph->edges(meshExchange7, corrWallBC);
    graph->edges(corrWallBC, collector6aSM);                  // Collect for MESH_EXCHANGE(6)
    graph->edges(collector6aSM, meshExchange6a);              // Do MESH_EXCHANGE(6)
    graph->edges(meshExchange6a, corrRadiation);
    graph->edges(corrRadiation, collector2SM);                // Collect for MESH_EXCHANGE(2)
    graph->edges(collector2SM, meshExchange2);                // Do MESH_EXCHANGE(2)
    graph->edges(meshExchange2, corrInitDivCollectorSM);      // Collect for INIT_DIV
    graph->edges(corrInitDivCollectorSM, corrInitDivTask);    // Do INIT_DIV_INTEGRALS
    // Corrector divergence part 1 sub-graph (sequential COMBUSTION_BC + parallel kernel)
    graph->edges(corrInitDivTask, corrDivP1OrchSM);            // InitDiv -> CorrDivP1 Orchestrator
    graph->edges(corrDivP1OrchSM, corrDivP1KernelTask);        // Orchestrator -> Kernel (parallel)
    graph->edges(corrDivP1KernelTask, corrDivP1CollectorSM);   // Kernel -> Collector
    graph->edges(corrDivP1CollectorSM, corrDivCollectorSM);    // Collector -> DIV_EXCHANGE + RTE
    graph->edges(corrDivCollectorSM, corrDivExchangeTask);    // Do EXCHANGE_DIV_INFO + RTE
    // Divergence part 2 sub-graph (corrector, parallel multi-mesh execution)
    graph->edges(corrDivExchangeTask, corrDivP2OrchSM);       // DivExchange -> DivP2 Orchestrator
    graph->edges(corrDivP2OrchSM, corrDivP2KernelTask);       // Orchestrator -> Kernel (parallel)
    graph->edges(corrDivP2KernelTask, corrDivP2CollectorSM);  // Kernel -> Collector
    graph->edges(corrDivP2CollectorSM, corrPressureCollectorSM); // Collector -> PRESSURE
    graph->edges(corrPressureCollectorSM, corrPressureTask);  // Do PRESSURE_ITERATION
    // Velocity corrector sub-graph (parallel multi-mesh execution)
    graph->edges(corrPressureTask, velCorrOrchSM);            // Pressure → VelCorrOrchestrator
    graph->edges(velCorrOrchSM, velCorrKernelTask);           // Orchestrator → Kernel (parallel)
    graph->edges(velCorrKernelTask, velCorrCollectorSM);      // Kernel → Collector
    graph->edges(velCorrCollectorSM, collector6bSM);          // Collector → MESH_EXCHANGE(6)
    graph->edges(collector6bSM, meshExchange6b);              // Do MESH_EXCHANGE(6)
    graph->edges(meshExchange6b, corrFinal);

    // End of time step: corrFinal -> collector -> timestep task -> loop state -> cycle
    graph->edges(corrFinal, timestepCollectorSM);
    graph->edges(timestepCollectorSM, timestepTask);
    graph->edges(timestepTask, timestepLoopSM);

    // Cycle: timestep loop state -> back to predictor (MeshData output)
    graph->edges(timestepLoopSM, predStep1OrchSM);

    // Termination path: timestep loop state -> termination sink (BarrierData output when done)
    graph->edges(timestepLoopSM, terminationSinkSM);

    // Graph output: termination sink emits final BarrierData for clean shutdown
    graph->outputs(terminationSinkSM);

    return graph;
}

#endif // FDS_GRAPH_H
