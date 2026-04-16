#ifndef POST_EXCHANGE_ROUTER_STATE_H
#define POST_EXCHANGE_ROUTER_STATE_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/termination_data.h"

/// Pass-through router: dispatches each MeshData<Pressure> immediately
/// (no collection) based on exchangeRound to SolvePhase or VelErrorPhase.
///
/// exchangeRound % 2 == 0 -> MeshData<SolvePhase> (pre-solve exchange done)
/// exchangeRound % 2 == 1 -> MeshData<VelErrorPhase> (post-solve exchange done)
///
/// Also receives TerminationData so the StateManager's canTerminate()
/// can break the cycle: Exchange -> Router -> SolveKernel -> Exchange.
class PostExchangeRouterState
    : public hh::AbstractState<2,
          MeshData<MeshState::Pressure>, TerminationData,
          MeshData<MeshState::SolvePhase>, MeshData<MeshState::VelErrorPhase>> {
public:
    void execute(std::shared_ptr<MeshData<MeshState::Pressure>> md) override {
        if (md->exchangeRound % 2 == 0) {
            this->addResult(retag<MeshState::SolvePhase>(md));
        } else {
            this->addResult(retag<MeshState::VelErrorPhase>(md));
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
    : public hh::StateManager<2,
          MeshData<MeshState::Pressure>, TerminationData,
          MeshData<MeshState::SolvePhase>, MeshData<MeshState::VelErrorPhase>> {
public:
    PostExchangeRouterManager(
        std::shared_ptr<PostExchangeRouterState> const &state,
        std::string const &name)
        : hh::StateManager<2,
              MeshData<MeshState::Pressure>, TerminationData,
              MeshData<MeshState::SolvePhase>, MeshData<MeshState::VelErrorPhase>>(
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
