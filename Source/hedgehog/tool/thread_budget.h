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
    size_t predDivPrefork;      // DivP1PreforkKernelTask       (LIGHT) — split from barrier
    size_t predDivPart2;        // DivergencePart2KernelTask    (LIGHT) — CC_IBM path only
    size_t velPredictor;        // VelocityPredictorKernelTask  (LIGHT)
    size_t predSynTurbVelBC;    // PredSynTurbVelBCTask          (MEDIUM) — merged SynTurb+VelBC
    size_t retryMomDiv;         // RetryMomentumDivKernelTask   (LIGHT)
    size_t predDivParallel;     // PredDivParallelTask           (LIGHT) — merged DivP1Late+DivP2Pre+DivPart2

    // --- Predictor fork: {DivSetup+PartMom} || {WallBC+DivP1Early} ---
    size_t predForkDivSetupPartMom; // PredDivSetupPartMomTask   (HEAVY) — merged DivSetup+PartMom
    size_t predForkWallBCDivEarly;  // PredWallBCDivEarlyTask    (HEAVY) — merged WallBC+DivP1Early

    // --- Corrector standalone sections ---
    size_t corrStep1;           // CorrStep1KernelTask          (HEAVY)
    size_t corrParticleOps;     // ParticleOpsKernelTask         (MEDIUM) — split from barrier
    size_t corrDivPart2;        // DivergencePart2KernelTask    (MEDIUM)
    size_t velCorrector;        // VelocityCorrectorKernelTask  (LIGHT)
    size_t corrFinalVelBC;      // VelocityBCEdgesTask          (MEDIUM)
    size_t corrWallBC;          // WallBCKernelTask             (MEDIUM) — includes finalize
    size_t corrDivParallel;     // CorrDivParallelTask           (MEDIUM) — merged QRAddCopy+DivP2Pre+DivPart2

    // --- Corrector fork1: {DivSetup} || {Fork1Comb} ---
    size_t corrFork1DivSetup;   // DivSetupKernelTask     (MEDIUM)
    size_t corrFork1Comb;       // Fork1CombKernelTask    (MEDIUM)

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
        b.predDivPrefork = solo(1);  // split from barrier (light)
        b.predDivPart2   = solo(1);  // 283us/elem
        b.velPredictor   = solo(1);  // 90us/elem
        b.predSynTurbVelBC = solo(2);  // merged SynTurb(light)+VelBC(medium)
        b.retryMomDiv    = solo(1);  // rarely used
        b.predDivParallel = solo(1); // merged DivP1Late+DivP2Pre+DivPart2 thread pool

        // --- Predictor fork: A{DivSetupPartMom(4)} || B{WallBCDivEarly(4)} ---
        {
            auto t = distribute({4, 4});
            b.predForkDivSetupPartMom = t[0];
            b.predForkWallBCDivEarly  = t[1];
        }

        // --- Corrector standalone ---
        b.corrStep1         = solo(4);  // 1.5ms/elem
        b.corrParticleOps   = solo(2);  // split from barrier (medium)
        b.corrDivPart2      = solo(2);  // 604us/elem
        b.velCorrector      = solo(1);  // 89us/elem
        b.corrFinalVelBC    = solo(2);  // 472us/elem
        b.corrWallBC        = solo(2);  // 513us/elem (includes finalize)
        b.corrDivParallel   = solo(2);  // merged QRAddCopy+DivP2Pre+DivPart2 (heaviest is MEDIUM)

        // --- Corrector fork1: A{DivSetup(2)} || B{Fork1Comb(2)} ---
        {
            auto t = distribute({2, 2});
            b.corrFork1DivSetup = t[0];
            b.corrFork1Comb     = t[1];
        }

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
           << " divPrefork=" << predDivPrefork
           << " divP2=" << predDivPart2
           << " velPred=" << velPredictor
           << " synTurbVelBC=" << predSynTurbVelBC
           << " retry=" << retryMomDiv
           << " divParallel=" << predDivParallel << "\n"
           << "  Pred fork:  divSetupPartMom=" << predForkDivSetupPartMom
           << " wallBCDivEarly=" << predForkWallBCDivEarly << "\n"
           << "  Corrector:  step1=" << corrStep1
           << " particleOps=" << corrParticleOps
           << " divP2=" << corrDivPart2
           << " velCorr=" << velCorrector
           << " finalVBC=" << corrFinalVelBC
           << " wallBC=" << corrWallBC
           << " divParallel=" << corrDivParallel << "\n"
           << "  Corr fork1: divSetup=" << corrFork1DivSetup
           << " comb=" << corrFork1Comb << "\n"
           << "  Corr fork2: radiation=" << corrFork2Radiation
           << " divP1=" << corrFork2DivP1 << "\n"
           << "  Pressure:   parallel=" << pressureParallel
           << " push=" << exchangePush
           << " pull=" << exchangePull << std::endl;
    }
};

#endif // THREAD_BUDGET_H
