#ifndef TERMINATION_DATA_H
#define TERMINATION_DATA_H

#include <hedgehog/hedgehog.h>

/// Empty marker type pushed through the graph to signal termination.
///
/// Pushed by the main thread after getBlockingResult() returns (simulation
/// complete). States that manage cycles receive this via execute(), set
/// done_=true, and their canTerminate() returns true to break the cycle.
struct TerminationData {};

/// Sink state for TerminationData when no cycle states need it.
/// Used in predictor/corrector subgraphs when parallel pressure is disabled.
class TerminationDataSink : public hh::AbstractState<1, TerminationData, TerminationData> {
public:
    void execute(std::shared_ptr<TerminationData>) override {}
};

#endif // TERMINATION_DATA_H
