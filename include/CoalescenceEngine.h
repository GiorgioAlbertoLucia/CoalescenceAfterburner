#ifndef COALESCENCEENGINE_H
#define COALESCENCEENGINE_H

#include "Particle.h"
#include "SourceSize.h"
#include "SourceSampler.h"
#include "WignerDensity.h"
#include "EventReader.h"

#include "TRandom3.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include <memory>
#include <vector>

// ─────────────────────────────────────────────────────────────────────────────
// CoalescenceResult
//
// Output of processing one event: the list of accepted Li4 candidates,
// each carrying the 4-momentum of the He3+p system, plus per-pair
// diagnostics (q, r) for ALL pairs evaluated (accepted or not).
// ─────────────────────────────────────────────────────────────────────────────

struct PairDiagnostic {
    double q_GeV; // relative momentum magnitude in PRF [GeV/c]
    double r_fm;  // relative distance magnitude in PRF [fm]
};

struct CoalescenceResult {
    std::vector<Particle>       li4s;  // accepted Li4 candidates
    std::vector<PairDiagnostic> pairs; // all evaluated pairs (for diagnostics)

    void clear() { li4s.clear(); pairs.clear(); }
    int nLi4() const { return static_cast<int>(li4s.size()); }
};


// ─────────────────────────────────────────────────────────────────────────────
// CoalescenceEngine
//
// Implements event-by-event Li4 (He3 + p) coalescence via the Wigner function
// formalism described in Mahlein et al. (arXiv:2302.12696).
//
// For each He3-proton pair in an event the engine:
//   1. Assigns emission positions via Gaussian sampling (SourceSampler)
//   2. Boosts the pair into its rest frame (PRF)
//   3. Computes the relative momentum q and relative position r in the PRF
//   4. Evaluates the Wigner density D(q, r) via the WignerDensity instance
//   5. Accepts the pair as Li4 with probability P = D(q,r) / D_max
//      using statistical rejection sampling
//
// The Li4 4-momentum is set to the sum of the He3 and proton 4-momenta.
//
// Ownership:
//   - WignerDensity is held via shared_ptr (polymorphic, swappable)
//   - SourceSize objects are owned by the engine
//   - SourceSampler objects hold const refs to the owned SourceSize objects
// ─────────────────────────────────────────────────────────────────────────────

class CoalescenceEngine {
public:

    // srcP    : SourceSize for protons
    // srcHe3  : SourceSize for He3
    // wigner  : Wigner density implementation (e.g. GaussianWigner)
    // seed    : RNG seed for rejection sampling and position sampling
    CoalescenceEngine(const SourceSize&                  srcP,
                      const SourceSize&                  srcHe3,
                      std::shared_ptr<WignerDensity>     wigner,
                      unsigned int                       seed = 42);

    // Process one event: sample positions, loop over pairs, apply coalescence.
    // Returns a CoalescenceResult containing accepted Li4 candidates.
    CoalescenceResult processEvent(const Event& ev);

    // Precompute the coalescence probability P(r, q) and store in a TH2D.
    // Must be called once before processEvent().
    // rMax_fm  : upper edge of the r axis [fm]
    // qMax_GeV : upper edge of the q axis [GeV/c]
    // nBinsR, nBinsQ : number of bins on each axis
    void precompute(double rMax_fm  = 20., double qMax_GeV = 2.,
                    int nBinsR = 200,      int nBinsQ      = 200);

    // Access the precomputed histogram (e.g. to write it to a ROOT file)
    const TH2D* coalescenceProbability() const { return fProbMap.get(); }

    // Reseed all internal RNGs
    void setSeed(unsigned int seed);

    // Access the Wigner density (e.g. to swap wavefunction hypothesis)
    void setWignerDensity(std::shared_ptr<WignerDensity> wigner) {
        fWigner = std::move(wigner);
    }
    const WignerDensity& wignerDensity() const { return *fWigner; }

    void setProbMap(const TH2D* map) { fProbMap.reset(); fSharedProbMap = map; }

private:

    // Sample positions for all particles in a list
    void samplePositions(std::vector<Particle>& particles,
                         SourceSampler&         sampler) const;

    // Boost pair to PRF, return relative momentum [GeV/c] and
    // relative position [fm] as out-parameters.
    void boostToPRF(const Particle& a,
                    const Particle& b,
                    TVector3&       q_out,
                    TVector3&       r_out) const;
    void boostToPRF(const Particle& a,
                    const Particle& b,
                    TVector3&       q_out) const;

    // Build the Li4 candidate from a He3+p pair
    static Particle makeLi4(const Particle& he3, const Particle& p);

    // ── Owned objects ────────────────────────────────────────────────────────
    SourceSize   fSrcProton;   // source size parametrisation for protons
    SourceSize   fSrcHe3;      // source size parametrisation for He3
    PairSpaceDistribution fPairDist; // pair space distribution for the two sources

    std::shared_ptr<WignerDensity> fWigner; // Wigner density (polymorphic)

    mutable TRandom3          fRng;     // for rejection sampling
    std::unique_ptr<TH2D>     fProbMap; // precomputed P(r [fm], q [GeV/c])
    const TH2D*               fSharedProbMap{nullptr}; // non-owning, from main thread
};

#endif // COALESCENCEENGINE_H
