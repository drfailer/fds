#ifndef TERMINATION_SIGNAL_H
#define TERMINATION_SIGNAL_H

#include <atomic>

/// Shared termination signal for sub-graph cycle termination.
///
/// The main timestep loop owns this signal and calls terminate() when the
/// simulation is done. Sub-graph states with internal cycles check
/// isTerminated() in their canTerminate() override. This decouples cycle
/// termination from data-driven conditions (like time), centralizing the
/// decision in the timestep loop.
struct TerminationSignal {
    std::atomic<bool> terminated{false};

    void terminate() { terminated.store(true, std::memory_order_release); }
    bool isTerminated() const { return terminated.load(std::memory_order_acquire); }
};

#endif // TERMINATION_SIGNAL_H
