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
    size_t predDivPart2;        // DivergencePart2KernelTask    (LIGHT)
    size_t velPredictor;        // VelocityPredictorKernelTask  (LIGHT)
    size_t predSynTurb;         // SyntheticTurbulenceKernelTask(LIGHT) — split from barrier
    size_t predFinalVelBC;      // VelocityBCEdgesTask          (MEDIUM)
    size_t retryMomDiv;         // RetryMomentumDivKernelTask   (LIGHT)
    size_t predDivP1Late;       // DivP1LateKernelTask          (LIGHT) — split from barrier

    // --- Predictor fork: {DivSetup+PartMom} || {WallBC+DivP1Early} ---
    size_t predForkDivSetup;    // DivSetupKernelTask     (HEAVY)
    size_t predForkPartMom;     // PredPartMomKernelTask  (LIGHT)
    size_t predForkWallBC;      // WallBCKernelTask       (MEDIUM)
    size_t predForkDivP1Early;  // DivP1EarlyTask         (MEDIUM)

    // --- Corrector standalone sections ---
    size_t corrStep1;           // CorrStep1KernelTask          (HEAVY)
    size_t corrParticleOps;     // ParticleOpsKernelTask         (MEDIUM) — split from barrier
    size_t corrDivPart2;        // DivergencePart2KernelTask    (MEDIUM)
    size_t velCorrector;        // VelocityCorrectorKernelTask  (LIGHT)
    size_t corrFinalVelBC;      // VelocityBCEdgesTask          (MEDIUM)
    size_t corrWallBC;          // WallBCKernelTask             (MEDIUM)
    size_t corrWallBCFinalize;  // WallBCFinalizeKernelTask     (LIGHT) — split from barrier
    size_t corrQRAddCopy;       // QRAddCopyKernelTask          (LIGHT) — split from barrier

    // --- Corrector fork1: {DivSetup} || {Fork1Comb} ---
    size_t corrFork1DivSetup;   // DivSetupKernelTask     (MEDIUM)
    size_t corrFork1Comb;       // Fork1CombKernelTask    (MEDIUM)

    // --- Corrector fork2: {Radiation} || {Fork2DivP1} ---
    size_t corrFork2Radiation;  // CorrRadiationKernelTask (LIGHT)
    size_t corrFork2DivP1;      // Fork2DivP1KernelTask    (MEDIUM)

    // --- Pressure iteration (shared by pred & corr instances) ---
    size_t baroclinic;          // BaroclinicKernelTask    (LIGHT)
    size_t fluxExchange;        // FluxExchangeTask        (LIGHT)
    size_t pressureSolve;       // PressureSolveKernelTask (MEDIUM)
    size_t velError;            // VelocityErrorTask       (LIGHT)

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
        b.predSynTurb    = solo(1);  // split from barrier (light)
        b.predFinalVelBC = solo(2);  // 615us/elem
        b.retryMomDiv    = solo(1);  // rarely used
        b.predDivP1Late  = solo(1);  // split from barrier (light)

        // --- Predictor fork: A{DivSetup(4)+PartMom(1)} || B{WallBC(2)+DivP1Early(2)} ---
        {
            auto t = distribute({4, 1, 2, 2});
            b.predForkDivSetup  = t[0];
            b.predForkPartMom   = t[1];
            b.predForkWallBC    = t[2];
            b.predForkDivP1Early = t[3];
        }

        // --- Corrector standalone ---
        b.corrStep1         = solo(4);  // 1.5ms/elem
        b.corrParticleOps   = solo(2);  // split from barrier (medium)
        b.corrDivPart2      = solo(2);  // 604us/elem
        b.velCorrector      = solo(1);  // 89us/elem
        b.corrFinalVelBC    = solo(2);  // 472us/elem
        b.corrWallBC        = solo(2);  // 513us/elem
        b.corrWallBCFinalize = solo(1); // split from barrier (light)
        b.corrQRAddCopy     = solo(1);  // split from barrier (light)

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
        // Baroclinic → ExchPre → Solve → ExchPost → VelError
        {
            auto t = distribute({1, 1, 4, 2});
            b.baroclinic    = t[0];
            b.fluxExchange  = t[1];
            b.pressureSolve = t[2];
            b.velError      = t[3];
        }

        return b;
    }

    void print(std::ostream &os) const {
        os << "[FDS-HH] Thread budget (cap=" << cap_ << "):\n"
           << "  Predictor:  step1=" << predStep1
           << " divPrefork=" << predDivPrefork
           << " divP2=" << predDivPart2
           << " velPred=" << velPredictor
           << " synTurb=" << predSynTurb
           << " finalVBC=" << predFinalVelBC
           << " retry=" << retryMomDiv
           << " divP1Late=" << predDivP1Late << "\n"
           << "  Pred fork:  divSetup=" << predForkDivSetup
           << " partMom=" << predForkPartMom
           << " wallBC=" << predForkWallBC
           << " divP1Early=" << predForkDivP1Early << "\n"
           << "  Corrector:  step1=" << corrStep1
           << " particleOps=" << corrParticleOps
           << " divP2=" << corrDivPart2
           << " velCorr=" << velCorrector
           << " finalVBC=" << corrFinalVelBC
           << " wallBC=" << corrWallBC
           << " wallBCFin=" << corrWallBCFinalize
           << " qrAddCopy=" << corrQRAddCopy << "\n"
           << "  Corr fork1: divSetup=" << corrFork1DivSetup
           << " comb=" << corrFork1Comb << "\n"
           << "  Corr fork2: radiation=" << corrFork2Radiation
           << " divP1=" << corrFork2DivP1 << "\n"
           << "  Pressure:   baroclinic=" << baroclinic
           << " exchange=" << fluxExchange
           << " solve=" << pressureSolve
           << " velError=" << velError << std::endl;
    }
};

#endif // THREAD_BUDGET_H
