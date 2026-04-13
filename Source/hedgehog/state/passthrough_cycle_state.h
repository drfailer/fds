#ifndef PASSTHROUGH_CYCLE_STATE_H
#define PASSTHROUGH_CYCLE_STATE_H

#include <hedgehog/hedgehog.h>
#include "../data/termination_data.h"

/// Generic passthrough state that breaks cycles via TerminationData.
///
/// Passes T data through unchanged. Receives TerminationData to set
/// done_=true, enabling canTerminate() to return true and break the cycle.
///
/// Used at the fds_graph level between the shared pressure subgraph output
/// and the predictor/corrector input to break the bidirectional cycle:
///   predictor → pressure → predReturnSM → predictor
template<typename T>
class PassthroughCycleState
    : public hh::AbstractState<2, T, TerminationData, T> {
public:
    void execute(std::shared_ptr<T> data) override {
        this->addResult(data);
    }

    void execute(std::shared_ptr<TerminationData>) override {
        done_ = true;
    }

    [[nodiscard]] bool isDone() const { return done_; }

private:
    bool done_ = false;
};

/// State manager wrapper with canTerminate() for cycle breaking.
template<typename T>
class PassthroughCycleManager
    : public hh::StateManager<2, T, TerminationData, T> {
public:
    PassthroughCycleManager(
        std::shared_ptr<PassthroughCycleState<T>> const &state,
        std::string const &name)
        : hh::StateManager<2, T, TerminationData, T>(state, name) {}

    [[nodiscard]] bool canTerminate() const override {
        this->state()->lock();
        auto s = std::dynamic_pointer_cast<PassthroughCycleState<T>>(
            this->state());
        bool ret = s->isDone();
        this->state()->unlock();
        return ret;
    }
};

#endif // PASSTHROUGH_CYCLE_STATE_H
