#include "CMFrameBooster.h"
#include "CoalescenceEngineHe3.h"
#include "JacobiTransform.h"

#include "TString.h"
 
#include <cmath>
#include <stdexcept>
 
// ─────────────────────────────────────────────────────────────────────────────
// Constructor
// ─────────────────────────────────────────────────────────────────────────────
 
CoalescenceEngineHe3::CoalescenceEngineHe3(const SourceSize&                srcNucleon,
                                            std::shared_ptr<WignerDensityA> wigner,
                                            unsigned int                    seed,
                                            unsigned int                    threadIndex)
    : fSrcNucleon(srcNucleon)
    , fWigner(std::move(wigner))
    , fRng(seed + 2)
    , fThreadIndex(threadIndex)
{
    if (!fWigner)
        throw std::invalid_argument(
            "CoalescenceEngineHe3: WignerDensityA must not be null");
    fWigner1 = std::make_shared<SingleGaussianWigner>(fWigner->getD());

    fHPtNucleon = new TH1D(
        Form("hPtNucleon_CEHe3_thread%d", fThreadIndex),
        "Sampled nucleon p_{T};p_{T} (GeV/c);Counts",
        100, 0., 10.);
    fHYNucleon = new TH1D(
        Form("hYNucleon_CEHe3_thread%d", fThreadIndex),
        "Sampled nucleon y;y (GeV/c);Counts",
        50, -1., 1.);
    fHPhiNucleon = new TH1D(
        Form("hPhiNucleon_CEHe3_thread%d", fThreadIndex),
        "Sampled nucleon #phi emission direction;#phi (rad);Counts",
        50, 0., 2. * M_PI);
    fHRadiusNucleon = new TH1D(
        Form("hRadiusNucleon_CEHe3_thread%d", fThreadIndex),
        "Sampled nucleon radius;r (fm);Counts",
        50, 0., 10.);
    fHPhiCoordinateNucleon = new TH1D(
        Form("hPhiCoordinateNucleon_CEHe3_thread%d", fThreadIndex),
        "Sampled nucleon #phi coordinate;#phi (rad);Counts",
        50, -M_PI, M_PI);
    fHCosThetaCoordinateNucleon = new TH1D(
        Form("hCosThetaCoordinateNucleon_CEHe3_thread%d", fThreadIndex),
        "Sampled nucleon #theta coordinate;#theta (rad);Counts",
        50, -1., 1.);
    fHPtNucleus = new TH1D(
        Form("hPtNucleus_CEHe3_thread%d", fThreadIndex),
        "Sampled nucleus p_{T};p_{T} (GeV/c);Counts",
        100, 0., 10.);
    fHYNucleus = new TH1D(
        Form("hYNucleus_CEHe3_thread%d", fThreadIndex),
        "Sampled nucleus y;y (GeV/c);Counts",
        250, -2.5, 2.5);
    fHPhiNucleus = new TH1D(
        Form("hPhiNucleus_CEHe3_thread%d", fThreadIndex),
        "Sampled nucleus #phi;#phi (rad);Counts",
        50, 0., 2. * M_PI);
}
 
// ─────────────────────────────────────────────────────────────────────────────
// wignerD1
//
// Returns D_1(k, r) for a single Jacobi pair.
// Uses the precomputed TH2D map if available (fast path),
// otherwise falls back to analytic evaluation (slow path).
// ─────────────────────────────────────────────────────────────────────────────

double CoalescenceEngineHe3::wignerD1(const TVector3& k_GeV,
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
 
float CoalescenceEngineHe3::processEvent(const Event& ev) const {
    
    float result = 0.f;
    const int A = 3; // number of nucleons in He3
    const float SA = 1.f / 12.f; // spin-isospin factor for He3 (1/12)
 
    if (ev.nProtons() < 2 || ev.nNeutrons() < 1) return result;
 
    // ── 1. Copy particles so we can assign positions ─────────────────────────
    std::vector<Particle> protons = ev.protons;
    std::vector<Particle> neutrons = ev.neutrons;
 
    // Loop over all 3 nucleon combinations (2p + 1n) and apply coalescence.
    // All particles can coalesce any time (no "used" tracking), so we just sum the weights.
    // Always use different protons (no self-pairs)

    for (std::size_t iP1 = 0; iP1 < protons.size(); ++iP1) {
        for (std::size_t iP2 = iP1 + 1; iP2 < protons.size(); ++iP2) {
            for (std::size_t iN1 = 0; iN1 < neutrons.size(); ++iN1) {
                    
                const Particle& p1 = protons[iP1];
                const Particle& p2 = protons[iP2];
                const Particle& n1 = neutrons[iN1];

                if (true) {
                    fHPtNucleon->Fill(p1.mom.Pt());
                    fHPtNucleon->Fill(p2.mom.Pt());
                    fHPtNucleon->Fill(n1.mom.Pt());
                    fHYNucleon->Fill(p1.mom.Rapidity());
                    fHYNucleon->Fill(p2.mom.Rapidity());
                    fHYNucleon->Fill(n1.mom.Rapidity());
                    fHPhiNucleon->Fill(p1.mom.Phi());
                    fHPhiNucleon->Fill(p2.mom.Phi());
                    fHPhiNucleon->Fill(n1.mom.Phi());

                    fHRadiusNucleon->Fill(p1.pos.Mag());
                    fHRadiusNucleon->Fill(p2.pos.Mag());
                    fHRadiusNucleon->Fill(n1.pos.Mag());
                    fHPhiCoordinateNucleon->Fill(p1.pos.Phi());
                    fHPhiCoordinateNucleon->Fill(p2.pos.Phi());
                    fHPhiCoordinateNucleon->Fill(n1.pos.Phi());
                    fHCosThetaCoordinateNucleon->Fill(p1.pos.CosTheta());
                    fHCosThetaCoordinateNucleon->Fill(p2.pos.CosTheta());
                    fHCosThetaCoordinateNucleon->Fill(n1.pos.CosTheta());
                }

                std::vector<Particle> nucleons{p1, p2, n1};
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
                if (true) {
                    TLorentzVector pTot;
                    for (const auto& p : nucleons)
                        pTot += p.mom;
                    fHPtNucleus->Fill(pTot.Pt(), SA*D);
                    fHYNucleus->Fill(pTot.Rapidity(), SA*D);
                    fHPhiNucleus->Fill(pTot.Phi(), SA*D);
                }
            }
        }
    }
 
    return result;
}

void CoalescenceEngineHe3::setSeed(unsigned int seed) {
    fRng.SetSeed(seed + 2);
}
