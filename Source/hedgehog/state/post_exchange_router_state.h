#ifndef POST_EXCHANGE_ROUTER_STATE_H
#define POST_EXCHANGE_ROUTER_STATE_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/pressure_iteration_data.h"
#include "../data/termination_data.h"

/// Pass-through router: dispatches each MeshData immediately (no collection)
/// based on exchangeRound to SolvePhaseData or VelErrorPhaseData.
///
/// exchangeRound % 2 == 0 -> SolvePhaseData (pre-solve exchange done)
/// exchangeRound % 2 == 1 -> VelErrorPhaseData (post-solve exchange done)
///
/// Also receives TerminationData so the StateManager's canTerminate()
/// can break the cycle: Exchange -> Router -> SolveKernel -> Exchange.
class PostExchangeRouterState
    : public hh::AbstractState<2, MeshData, TerminationData,
                               SolvePhaseData, VelErrorPhaseData> {
public:
    void execute(std::shared_ptr<MeshData> md) override {
        if (md->exchangeRound % 2 == 0) {
            this->addResult(std::make_shared<SolvePhaseData>(std::move(md)));
        } else {
            this->addResult(std::make_shared<VelErrorPhaseData>(std::move(md)));
        }
    }

    void execute(std::shared_ptr<TerminationData>) override {
        done_ = true;
    }

    [[nodiscard]] bool isDone() const { return done_; }

private:
    bool done_ = false;
};

/// StateManager for PostExchangeRouterState with canTerminate to break cycle.
class PostExchangeRouterManager
    : public hh::StateManager<2, MeshData, TerminationData,
                              SolvePhaseData, VelErrorPhaseData> {
public:
    PostExchangeRouterManager(
        std::shared_ptr<PostExchangeRouterState> const &state,
        std::string const &name)
        : hh::StateManager<2, MeshData, TerminationData,
                           SolvePhaseData, VelErrorPhaseData>(
              state, name) {}

    [[nodiscard]] bool canTerminate() const override {
        this->state()->lock();
        auto s = std::dynamic_pointer_cast<PostExchangeRouterState>(
            this->state());
        bool ret = s->isDone();
        this->state()->unlock();
        return ret;
    }
};

#endif // POST_EXCHANGE_ROUTER_STATE_H
