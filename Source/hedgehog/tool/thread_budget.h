#ifndef THREAD_BUDGET_H
#define THREAD_BUDGET_H

#include <algorithm>
#include <cstdlib>
#include <cstddef>
#include <iostream>
#include <vector>

/// Thread budget for the FDS Hedgehog graph.
///
/// The graph is divided into sections separated by barriers. Within each
/// section, concurrent tasks share the available threads proportionally
/// to their compute weight (per-element cost):
///   HEAVY (4): > 500us/element
///   MEDIUM (2): 100-500us/element
///   LIGHT (1): < 100us/element
///
/// Standalone sections: single task gets a fraction of cap based on weight.
/// Fork sections: concurrent tasks share cap by weight.
/// Pressure iteration: pipeline tasks share cap by weight.
struct ThreadBudget {
    int cap_;       ///< available hardware threads
    int nmeshes_;

    // --- Predictor standalone sections ---
    size_t predStep1;           // PredStep1KernelTask          (HEAVY) — CC_IBM path only
    size_t predPreforkDiv;      // PredPreforkDivTask            (HEAVY) — merged PredStep1+Prefork+Fork(DivSetup+PartMom || WallBC+DivEarly)
                                //   Real OS threads = 2 × predPreforkDiv (each HH thread owns an AsyncWorker)
    size_t predDivPart2;        // DivergencePart2KernelTask    (LIGHT) — CC_IBM path only
    size_t velPredictor;        // VelocityPredictorKernelTask  (LIGHT) — CC_IBM path only
    size_t predSynTurbVelBC;    // PredSynTurbVelBCTask          (MEDIUM) — merged VelPred+SynTurb+VelBC
    size_t retryMomDiv;         // RetryMomentumDivKernelTask   (LIGHT)

    // --- Corrector standalone sections ---
    size_t corrDivSetupCombPart;// CorrDivSetupCombPartTask      (HEAVY) — merged CorrStep1+Fork1(DivSetup||Comb)+ParticleOps
                                //   Real OS threads = 2 × corrDivSetupCombPart (each HH thread owns AsyncWorker)
    size_t corrDivPart2;        // DivergencePart2KernelTask    (MEDIUM)
    size_t corrFinal;           // CorrFinalKernelTask          (MEDIUM) — merged VelCorr+VelBCEdges+QRAdd+DivPart2

    // --- Shared: DivExchangeTask (pred & corr) ---
    size_t divExchange;         // DivExchangeTask pool threads  (LIGHT) — DivP2Pre is ~283us/elem

    // --- Two-lane parallel compute (pred+corr unified) ---
    size_t forkableLane;        // ForkableParallelComputeLane   (HEAVY) — P1+P2+C1+C2+C3, has AsyncWorker
                                //   Real OS threads = 2 × forkableLane during fork phases (P2, C2)
    size_t parallelLane;        // ParallelComputeLane           (MEDIUM) — P3+P4+P5+C4+C5+C6+C7+C8
    bool useAsyncWorker = true; // ForkableLane uses AsyncWorker for P2/C2 fork phases

    // --- Pressure iteration (shared by pred & corr instances) ---
    size_t pressureParallel;    // PressureParallelTask          (HEAVY) — merged Baroclinic+Solve+VelError

    /// Compute standalone thread count for a given weight (1-4).
    /// Useful for tasks not in the named fields (e.g., CC_IBM path).
    [[nodiscard]] size_t standalone(int weight) const {
        return static_cast<size_t>(
            std::max(1, std::min(nmeshes_, cap_ * weight / 4)));
    }

    /// Compute thread budget from hardware capabilities and mesh count.
    static ThreadBudget compute(int hwThreads, int nmeshes) {
        ThreadBudget b{};
        b.nmeshes_ = nmeshes;
        b.cap_ = std::max(1, hwThreads);
        int cap = b.cap_;
        int meshCap = std::min(nmeshes, hwThreads);

        // --- Standalone: fraction of meshCap based on weight/4 ---
        // Clamped by nmeshes (only nmeshes items in flight per phase)
        auto solo = [meshCap, nmeshes](int weight) -> size_t {
            return static_cast<size_t>(
                std::max(1, std::min(nmeshes, meshCap * weight / 4)));
        };

        // --- Fork: distribute meshCap among concurrent tasks by weight ---
        auto distribute = [meshCap, nmeshes](std::initializer_list<int> weights) {
            int total = 0;
            for (int w : weights) total += w;
            std::vector<size_t> r;
            for (int w : weights) {
                r.push_back(static_cast<size_t>(
                    std::max(1, std::min(nmeshes, meshCap * w / total))));
            }
            return r;
        };

        // --- Predictor standalone ---
        b.predStep1      = solo(4);  // 1.4ms/elem
        // Merged prefork+fork: each HH thread owns an AsyncWorker.
        // Full cap so sequential parts (prefork, post-join) keep full parallelism.
        b.predPreforkDiv = solo(4);
        b.predDivPart2   = solo(1);  // 283us/elem
        b.velPredictor   = solo(1);  // 90us/elem
        b.predSynTurbVelBC = solo(2);  // merged SynTurb(light)+VelBC(medium)
        b.retryMomDiv    = solo(1);  // rarely used

        // --- Corrector standalone ---
        // Merged CorrStep1+Fork1+ParticleOps: each HH thread owns an AsyncWorker.
        // Full cap so ParticleOps (sequential after join) keeps full parallelism.
        // Brief 2× oversubscription during fork phase is acceptable.
        b.corrDivSetupCombPart = solo(4);
        b.corrDivPart2      = solo(2);  // 604us/elem
        b.corrFinal         = solo(2);  // merged VelCorr+VelBCEdges+QRAdd+DivPart2
        b.divExchange       = solo(1);  // DivExchangeTask: DivP2Pre pool (LIGHT, ~283us/elem)

        // --- Two-lane parallel compute ---
        // Both lanes can be active concurrently during Fork2 pipeline.
        // Oversubscription is acceptable — most phases are sequential (only 1 lane active).
        // No mesh count clamp: pipelining can use more threads than meshes.
        {
            int total = 3 + 2;
            b.forkableLane = static_cast<size_t>(std::max(1, cap * 3 / total));
            b.parallelLane = static_cast<size_t>(std::max(1, cap * 2 / total));
        }

        // Override with env vars if set
        if (auto *v = std::getenv("HH_FORKABLE_THREADS"))
            b.forkableLane = static_cast<size_t>(std::max(1, std::atoi(v)));
        if (auto *v = std::getenv("HH_PARALLEL_THREADS"))
            b.parallelLane = static_cast<size_t>(std::max(1, std::atoi(v)));
        if (auto *v = std::getenv("HH_NO_ASYNC_WORKER"))
            b.useAsyncWorker = (std::atoi(v) == 0);

        // --- Pressure iteration ---
        b.pressureParallel = solo(4);

        return b;
    }

    void print(std::ostream &os) const {
        os << "[FDS-HH] Thread budget (cap=" << cap_ << "):\n"
           << "  Predictor:  step1=" << predStep1
           << " preforkDiv=" << predPreforkDiv << "(+aw)"
           << " divP2=" << predDivPart2
           << " velPred=" << velPredictor
           << " synTurbVelBC=" << predSynTurbVelBC
           << " retry=" << retryMomDiv << "\n"
           << "  Corrector:  divSetupCombPart=" << corrDivSetupCombPart << "(+aw)"
           << " divP2=" << corrDivPart2
           << " corrFinal=" << corrFinal << "\n"
           << "  DivExchange: pool=" << divExchange << "\n"
           << "  TwoLane:    forkable=" << forkableLane
           << (useAsyncWorker ? "(+aw)" : "")
           << " parallel=" << parallelLane << "\n"
           << "  Pressure:   parallel=" << pressureParallel << std::endl;
    }
};

#endif // THREAD_BUDGET_H
