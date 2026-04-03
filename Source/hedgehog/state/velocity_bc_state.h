#ifndef VELOCITY_BC_STATE_H
#define VELOCITY_BC_STATE_H

// This file previously contained PredFinalCollector, CorrFinalCollector,
// and CorrFinalOrchestrator states. These have been converted to tasks:
//
//   PredFinalCollector  → CollectorState + PhaseTransitionTask (barrier_tasks.h)
//   CorrFinalOrchestrator → CollectorState + BarrierTask (barrier_state.h)
//   CorrFinalCollector  → CollectorState + CorrFinalDumpTask (barrier_tasks.h)
//
// See velocity_bc_subgraph.h for the new wiring.

#endif // VELOCITY_BC_STATE_H
