#include "CMFrameBooster.h"
#include "CoalescenceEngineLi4.h"
#include "JacobiTransform.h"
 
#include <cmath>
#include <stdexcept>
 
// ─────────────────────────────────────────────────────────────────────────────
// Constructor
// ─────────────────────────────────────────────────────────────────────────────
 
CoalescenceEngineLi4::CoalescenceEngineLi4(const SourceSize&                srcNucleon,
                                            std::shared_ptr<WignerDensityA> wigner,
                                            unsigned int                    seed)
    : fSrcNucleon(srcNucleon)
    , fWigner(std::move(wigner))
    , fRng(seed + 2)
{
    if (!fWigner)
        throw std::invalid_argument(
            "CoalescenceEngineLi4: WignerDensityA must not be null");
    fWigner1 = std::make_shared<SingleGaussianWigner>(fWigner->getD());
}
 
// ─────────────────────────────────────────────────────────────────────────────
// wignerD1
//
// Returns D_1(k, r) for a single Jacobi pair.
// Uses the precomputed TH2D map if available (fast path),
// otherwise falls back to analytic evaluation (slow path).
// ─────────────────────────────────────────────────────────────────────────────

double CoalescenceEngineLi4::wignerD1(const TVector3& k_GeV,
                                       const TVector3& r_fm) const
{
    if (fWignerMap) {
        const double kMag = k_GeV.Mag();
        const double rMag = r_fm.Mag();
 
        // TH2::Interpolate is only valid within the bin-centre range
        const double kMax = fWignerMap->GetXaxis()->GetBinCenter(
                                fWignerMap->GetNbinsX());
        const double rMax = fWignerMap->GetYaxis()->GetBinCenter(
                                fWignerMap->GetNbinsY());
 
        if (kMag < kMax && rMag < rMax)
            return fWignerMap->Interpolate(kMag, rMag);
        return 0.; // exponentially suppressed outside the map range
    }
 
    // Analytic fallback: extract D_1 from the A-body Wigner by evaluating
    // with a single-pair input
    return fWigner->evaluate({k_GeV}, {r_fm});
}

// ─────────────────────────────────────────────────────────────────────────────
// processEvent
// ─────────────────────────────────────────────────────────────────────────────
 
float CoalescenceEngineLi4::processEvent(const Event& ev) const {
    
    float result = 0.f;
    const int A = 4; // number of nucleons in Li4
    const float SA = 3.f / 16.f; // spin-isospin factor for Li4 (3/16)
 
    if (ev.nProtons() < 3 || ev.nNeutrons() < 1) return result;
 
    // ── 1. Copy particles so we can assign positions ─────────────────────────
    std::vector<Particle> protons = ev.protons;
    std::vector<Particle> neutrons = ev.neutrons;
 
    // Loop over all 4 nucleon combinations (3p + 1n) and apply coalescence.
    // All particles can coalesce any time (no "used" tracking), so we just sum the weights.
    // Always use different protons (no self-pairs)

    for (std::size_t iP1 = 0; iP1 < protons.size(); ++iP1) {
        for (std::size_t iP2 = iP1 + 1; iP2 < protons.size(); ++iP2) {
            for (std::size_t iP3 = iP2 + 1; iP3 < protons.size(); ++iP3) {
                for (std::size_t iN = 0; iN < neutrons.size(); ++iN) {
                    
                    const Particle& p1 = protons[iP1];
                    const Particle& p2 = protons[iP2];
                    const Particle& p3 = protons[iP3];
                    const Particle& n  = neutrons[iN];

                    // Sample positions for the four nucleons
                    std::vector<Particle> nucleons{p1, p2, p3, n};
                    std::vector<Particle> boostedNucleons = CMFrameBooster::boostParticles(nucleons);

                    // Transform to Jacobi coordinates
                    std::vector<Particle> jacobiNucleons = JacobiTransform::transform(boostedNucleons);
                    std::vector<TVector3> r_jacobiForIntegration, k_jacobiForIntegration;
                    r_jacobiForIntegration.reserve(A-1);
                    k_jacobiForIntegration.reserve(A-1);
                    for (int j = 1; j < A; ++j) {
                        r_jacobiForIntegration.push_back(jacobiNucleons[j].pos);
                        k_jacobiForIntegration.push_back(jacobiNucleons[j].mom.Vect());
                    }

                    // Evaluate A-body Wigner density
                    double D = 1.;
                    for (std::size_t j = 0; j < r_jacobiForIntegration.size(); ++j)
                        D *= wignerD1(k_jacobiForIntegration[j], r_jacobiForIntegration[j]);
            
                    result += static_cast<float>(SA * D);
                }
            }
        }
    }
 
    return result;
}

void CoalescenceEngineLi4::setSeed(unsigned int seed) {
    fRng.SetSeed(seed + 2);
}
