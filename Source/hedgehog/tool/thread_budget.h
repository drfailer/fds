#ifndef THREAD_BUDGET_H
#define THREAD_BUDGET_H

#include <algorithm>
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
    int cap_;       ///< min(nmeshes, hwThreads) — max useful parallelism
    int nmeshes_;

    // --- Predictor standalone sections ---
    size_t predStep1;           // PredStep1KernelTask          (HEAVY)
    size_t predPreforkDiv;      // PredPreforkDivTask            (HEAVY) — merged Prefork+Fork(DivSetup+PartMom || WallBC+DivEarly)
                                //   Real OS threads = 2 × predPreforkDiv (each HH thread owns an AsyncWorker)
    size_t predDivPart2;        // DivergencePart2KernelTask    (LIGHT) — CC_IBM path only
    size_t velPredictor;        // VelocityPredictorKernelTask  (LIGHT)
    size_t predSynTurbVelBC;    // PredSynTurbVelBCTask          (MEDIUM) — merged SynTurb+VelBC
    size_t retryMomDiv;         // RetryMomentumDivKernelTask   (LIGHT)

    // --- Corrector standalone sections ---
    size_t corrStep1;           // CorrStep1KernelTask          (HEAVY)
    size_t corrDivSetupCombPart;// CorrDivSetupCombPartTask      (HEAVY) — merged Fork1(DivSetup||Comb)+ParticleOps
                                //   Real OS threads = 2 × corrDivSetupCombPart (each HH thread owns AsyncWorker)
    size_t corrDivPart2;        // DivergencePart2KernelTask    (MEDIUM)
    size_t corrFinal;           // CorrFinalKernelTask          (MEDIUM) — merged VelCorr+VelBCEdges+RTE
    size_t corrWallBC;          // WallBCKernelTask             (MEDIUM) — includes finalize
    size_t corrDivParallel;     // CorrDivParallelTask           (MEDIUM) — merged QRAddCopy+DivP2Pre+DivPart2


    // --- Corrector fork2: {Radiation} || {Fork2DivP1} ---
    size_t corrFork2Radiation;  // CorrRadiationKernelTask (LIGHT)
    size_t corrFork2DivP1;      // Fork2DivP1KernelTask    (MEDIUM)

    // --- Pressure iteration (shared by pred & corr instances) ---
    size_t pressureParallel;    // PressureParallelTask          (HEAVY) — merged Baroclinic+Solve+VelError
    size_t exchangePush;        // ExchangePushBufferTask       (LIGHT)
    size_t exchangePull;        // ExchangePullBufferTask       (LIGHT)

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
        b.cap_ = std::min(nmeshes, hwThreads);
        int cap = b.cap_;

        // --- Standalone: fraction of cap based on weight/4 ---
        auto solo = [cap, nmeshes](int weight) -> size_t {
            return static_cast<size_t>(
                std::max(1, std::min(nmeshes, cap * weight / 4)));
        };

        // --- Fork: distribute cap among concurrent tasks by weight ---
        auto distribute = [cap, nmeshes](std::initializer_list<int> weights) {
            int total = 0;
            for (int w : weights) total += w;
            std::vector<size_t> r;
            for (int w : weights) {
                r.push_back(static_cast<size_t>(
                    std::max(1, std::min(nmeshes, cap * w / total))));
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
        b.corrStep1         = solo(4);  // 1.5ms/elem
        // Merged Fork1+ParticleOps: each HH thread owns an AsyncWorker.
        // Full cap so ParticleOps (sequential after join) keeps full parallelism.
        // Brief 2× oversubscription during fork phase is acceptable.
        b.corrDivSetupCombPart = solo(4);
        b.corrDivPart2      = solo(2);  // 604us/elem
        b.corrFinal         = solo(2);  // merged VelCorr(light)+VelBCEdges(medium)+RTE
        b.corrWallBC        = solo(2);  // 513us/elem (includes finalize)
        b.corrDivParallel   = solo(2);  // merged QRAddCopy+DivP2Pre+DivPart2 (heaviest is MEDIUM)

        // --- Corrector fork2: A{Radiation(1)} || B{Fork2DivP1(2)} ---
        {
            auto t = distribute({1, 2});
            b.corrFork2Radiation = t[0];
            b.corrFork2DivP1     = t[1];
        }

        // --- Pressure iteration pipeline ---
        // PressureParallel ↔ Exchange(Push+Pull) ↔ PressureParallel
        {
            auto t = distribute({4, 2, 2});
            b.pressureParallel = t[0];
            b.exchangePush     = t[1];
            b.exchangePull     = t[2];
        }

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
           << "  Corrector:  step1=" << corrStep1
           << " divSetupCombPart=" << corrDivSetupCombPart << "(+aw)"
           << " divP2=" << corrDivPart2
           << " corrFinal=" << corrFinal
           << " wallBC=" << corrWallBC
           << " divParallel=" << corrDivParallel << "\n"
           << "  Corr fork2: radiation=" << corrFork2Radiation
           << " divP1=" << corrFork2DivP1 << "\n"
           << "  Pressure:   parallel=" << pressureParallel
           << " push=" << exchangePush
           << " pull=" << exchangePull << std::endl;
    }
};

#endif // THREAD_BUDGET_H
