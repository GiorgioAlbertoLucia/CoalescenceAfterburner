#ifndef COALESCENCEENGINEHE3_H
#define COALESCENCEENGINEHE3_H

#include "Particle.h"
#include "SourceSize.h"
#include "WignerDensityA.h"
#include "Event.h"

#include "TRandom3.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include <memory>
#include <vector>

struct QAHistograms {
    TH1D* fHPtNucleon;
    TH1D* fHYNucleon;
    TH1D* fHPhiNucleon;
    TH1D* fHRadiusNucleon;
    TH1D* fHPhiCoordinateNucleon;
    TH1D* fHCosThetaCoordinateNucleon;
    TH1D* fHPtNucleus;
    TH1D* fHYNucleus;
    TH1D* fHPhiNucleus;
};

// ─────────────────────────────────────────────────────────────────────────────
// CoalescenceEngineHe3
//
// Event-by-event He3 (2p+2n) coalescence via the A-body Wigner function.
//
// For each unique (p,p,n,n) quartet in an event:
//   1. Boost quartet to CM frame
//   2. Transform to Jacobi coordinates
//   3. Look up D_A = prod_j D_1(k_j, r_j) from a precomputed TH2D map
//      (falls back to analytic evaluation if no map is set)
//   4. Accumulate SA * D_A as the event weight
//
// The lookup map is a TH2D of D_1(|k|, |r|) — the single Jacobi-pair
// Wigner density — shared as a non-owning pointer from IntegratorA3.
// ─────────────────────────────────────────────────────────────────────────────

class CoalescenceEngineHe3 {
public:

    CoalescenceEngineHe3(const SourceSize&               srcNucleon,
                         std::shared_ptr<WignerDensityA> wigner,
                         unsigned int                    seed = 42,
                         unsigned int                    threadIndex = 0);

    // Process one event. Returns the summed SA*D_A weight over all quartets.
    // Positions must already be assigned to ev.protons and ev.neutrons
    // by the caller (via SourceSampler) before calling this.
    float processEvent(const Event& ev) const;

    // Inject a precomputed D_1(|k| [GeV/c], |r| [fm]) lookup map.
    // Non-owning: the TH2D must remain alive for the engine's lifetime.
    // If not set, analytic GaussianWigner evaluation is used (slower).
    void setWignerMap(const TH2D* map) { fWignerMap = map; }

    void setSeed(unsigned int seed);

    void setWignerDensity(std::shared_ptr<WignerDensityA> wigner) {
        fWigner = std::move(wigner);
    }
    const WignerDensityA& wignerDensity() const { return *fWigner; }

    QAHistograms getQAHistograms() const {
        return QAHistograms{
            fHPtNucleon,
            fHYNucleon,
            fHPhiNucleon,
            fHRadiusNucleon,
            fHPhiCoordinateNucleon,
            fHCosThetaCoordinateNucleon,
            fHPtNucleus,
            fHYNucleus,
            fHPhiNucleus
        };
    }

private:

    // Look up or compute D_1 for a single Jacobi pair (k, r)
    double wignerD1(const TVector3& k_GeV, const TVector3& r_fm) const;

    SourceSize                              fSrcNucleon;
    std::shared_ptr<WignerDensityA>         fWigner;
    std::shared_ptr<SingleGaussianWigner>   fWigner1;
    const TH2D*                             fWignerMap{nullptr}; // non-owning lookup table
    mutable TRandom3                        fRng;
    int                                    fThreadIndex{0}; // for QA histogram naming

    // QA histograms
    TH1D* fHPtNucleon{nullptr};
    TH1D* fHYNucleon{nullptr};
    TH1D* fHPhiNucleon{nullptr};
    TH1D* fHRadiusNucleon{nullptr};
    TH1D* fHPhiCoordinateNucleon{nullptr};
    TH1D* fHCosThetaCoordinateNucleon{nullptr};
    TH1D* fHPtNucleus{nullptr};
    TH1D* fHYNucleus{nullptr};
    TH1D* fHPhiNucleus{nullptr};
};

#endif // COALESCENCEENGINEHE3_H