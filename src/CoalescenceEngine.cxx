#include "CoalescenceEngine.h"

#include <cmath>
#include <stdexcept>

// ─────────────────────────────────────────────────────────────────────────────
// Constructor
// ─────────────────────────────────────────────────────────────────────────────

CoalescenceEngine::CoalescenceEngine(const SourceSize&              srcP,
                                     const SourceSize&              srcHe3,
                                     std::shared_ptr<WignerDensity> wigner,
                                     unsigned int                   seed)
    : fSrcProton(srcP)
    , fSrcHe3(srcHe3)
    , fWigner(std::move(wigner))
    , fRng(seed + 2)                        // different seed for rejection
    , fPairDist(srcP, srcHe3)               // pair space distribution for the two sources
{
    if (!fWigner)
        throw std::invalid_argument(
            "CoalescenceEngine: WignerDensity must not be null");
}

// ─────────────────────────────────────────────────────────────────────────────
// precompute
// ─────────────────────────────────────────────────────────────────────────────

void CoalescenceEngine::precompute(double rMax_fm, double qMax_GeV,
                                   int nBinsR,     int nBinsQ)
{
    fProbMap = std::make_unique<TH2D>(
        "hCoalescenceProb",
        "Coalescence probability;r (fm);q (GeV/c)",
        nBinsR, 0., rMax_fm,
        nBinsQ, 0., qMax_GeV);
    fProbMap->SetDirectory(nullptr);

    const double rMin = 0., rMax = rMax_fm;
    const double RpairMin = 0., RpairMax = rMax_fm; // range of pair separation to consider
    const double isospinFactor = 3./16.;
    const double spinFactor = 5./16. + 3./16. + 1./16. + 3./16.; // considering excited states
    
    for (int iq = 1; iq <= nBinsQ; ++iq) {
        const double q = fProbMap->GetYaxis()->GetBinCenter(iq);
        const TVector3 q_vec(q, 0., 0.);

        for (int ir0 = 1; ir0 <= nBinsR; ++ir0) {
            const double r0 = fProbMap->GetXaxis()->GetBinCenter(ir0);
            PairSpaceDistribution pairDist(r0);

            int nSamples = 0;
            double sumP = 0.;

            for (double rIter = rMin; rIter <= rMax; rIter += (rMax - rMin) / nBinsR) {

                const TVector3 r_vec(rIter, 0., 0.);
                double D = fWigner->evaluate(q_vec, r_vec);
                
                double H_pHe3 = pairDist.evaluateIntegrateGlobal(rIter);
                sumP +=  D * H_pHe3 * 4 * M_PI * rIter * rIter; // weight by spherical volume element
                nSamples++;
            }
            
            double P = spinFactor * isospinFactor * sumP / nSamples;
            fProbMap->SetBinContent(ir0, iq, P);
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// processEvent
// ─────────────────────────────────────────────────────────────────────────────
 
CoalescenceResult CoalescenceEngine::processEvent(const Event& ev) {
    CoalescenceResult result;
 
    if (ev.nProtons() == 0 || ev.nHe3() == 0) return result;
 
    // ── 1. Copy particles so we can assign positions ─────────────────────────
    std::vector<Particle> protons = ev.protons;
    std::vector<Particle> he3s    = ev.he3s;
 
    // ── 3. Loop over all He3-proton pairs ────────────────────────────────────
    // Each particle can only coalesce once (greedy first-match).
    // Track which particles have already been used.
    std::vector<bool> protonUsed(protons.size(), false);
    std::vector<bool> he3Used   (he3s.size(),    false);
 
    for (std::size_t iHe3 = 0; iHe3 < he3s.size(); ++iHe3) {
        if (he3Used[iHe3]) continue;
 
        for (std::size_t iP = 0; iP < protons.size(); ++iP) {
            if (protonUsed[iP]) continue;
 
            const Particle& he3    = he3s[iHe3];
            const Particle& proton = protons[iP];
 
            // ── 3a. Boost to PRF, get q and r ────────────────────────────────
            TVector3 q_GeV;
            boostToPRF(he3, proton, q_GeV);
 
            // ── 3b. Store pair diagnostics (all pairs, before rejection) ─────
            const double qMag = q_GeV.Mag();
            const double rMag = fPairDist.pairR0();
            result.pairs.push_back({qMag, rMag});
 
            // ── 3c. Look up precomputed coalescence probability ───────────────
            // Prefer the shared (non-owning) map injected by the main thread;
            // fall back to the owned map, then to direct Wigner evaluation.
            const TH2D* map = fSharedProbMap ? fSharedProbMap : fProbMap.get();
            double P = 0.;
            if (map) {
                const double rMax = map->GetXaxis()->GetBinCenter(map->GetNbinsX());
                const double qMax = map->GetYaxis()->GetBinCenter(map->GetNbinsY());
                if (rMag < rMax && qMag < qMax)
                    P = map->Interpolate(rMag, qMag);
            } else {
                P = fWigner->evaluate(q_GeV, TVector3(rMag,0,0)) / fWigner->maxValue();
            }
 
            // ── 3d. Rejection sampling ────────────────────────────────────────
            const double u = fRng.Uniform(0., 1.);
            if (u < P) {
                result.li4s.push_back(makeLi4(he3, proton));
                protonUsed[iP]   = true;
                he3Used[iHe3]    = true;
                break; // move to next He3
            }
        }
    }
 
    return result;
}

// ─────────────────────────────────────────────────────────────────────────────
// samplePositions
// ─────────────────────────────────────────────────────────────────────────────

void CoalescenceEngine::samplePositions(std::vector<Particle>& particles,
                                        SourceSampler&         sampler) const {
    for (Particle& p : particles)
        sampler.sample(p);
}

// ─────────────────────────────────────────────────────────────────────────────
// boostToPRF
//
// Boosts the He3-proton pair into its rest frame and computes:
//   q = relative 3-momentum  = (p_a - p_b) / 2   [GeV/c]  in PRF
//   r = relative 3-position  = pos_a - pos_b      [fm]     in PRF
//
// The position boost uses the non-relativistic approximation for r,
// consistent with the equal-time approximation used in the paper (Sec. 2).
// ─────────────────────────────────────────────────────────────────────────────

void CoalescenceEngine::boostToPRF(const Particle& a,
                                   const Particle& b,
                                   TVector3&       q_out,
                                   TVector3&       r_out) const
{
    // 4-momentum of the pair (= Li4 4-momentum before binding correction)
    TLorentzVector pPair = a.mom + b.mom;

    // Beta vector of the pair centre-of-mass
    TVector3 beta = pPair.BoostVector();

    // Copy 4-momenta and boost into PRF
    TLorentzVector momA = a.mom;
    TLorentzVector momB = b.mom;
    momA.Boost(-beta);
    momB.Boost(-beta);

    // Relative momentum in PRF: q = (p_a - p_b) / 2
    q_out = (momA.Vect() - momB.Vect()) * 0.5;

    // Relative position in the lab frame, then Lorentz-contract along boost
    // Under the equal-time approximation we work with the 3D distance directly.
    // We apply the same boost to position 4-vectors (x, y, z, 0).
    TLorentzVector posA(a.pos, 0.);
    TLorentzVector posB(b.pos, 0.);
    posA.Boost(-beta);
    posB.Boost(-beta);

    r_out = posA.Vect() - posB.Vect();
}

void CoalescenceEngine::boostToPRF(const Particle& a,
                                   const Particle& b,
                                   TVector3&       q_out) const
{
    // 4-momentum of the pair (= Li4 4-momentum before binding correction)
    TLorentzVector pPair = a.mom + b.mom;

    // Beta vector of the pair centre-of-mass
    TVector3 beta = pPair.BoostVector();

    // Copy 4-momenta and boost into PRF
    TLorentzVector momA = a.mom;
    TLorentzVector momB = b.mom;
    momA.Boost(-beta);
    momB.Boost(-beta);

    // Relative momentum in PRF: q = (p_a - p_b) / 2
    q_out = (momA.Vect() - momB.Vect()) * 0.5;

    // Relative position in the lab frame, then Lorentz-contract along boost
    // Under the equal-time approximation we work with the 3D distance directly.
    // We apply the same boost to position 4-vectors (x, y, z, 0).
    TLorentzVector posA(a.pos, 0.);
    TLorentzVector posB(b.pos, 0.);
    posA.Boost(-beta);
    posB.Boost(-beta);
}

// ─────────────────────────────────────────────────────────────────────────────
// makeLi4
// ─────────────────────────────────────────────────────────────────────────────

Particle CoalescenceEngine::makeLi4(const Particle& he3, const Particle& p) {
    Particle li4;
    li4.pdg = PDG::kLi4;
    li4.mom = he3.mom + p.mom; // 4-momentum sum (binding energy neglected,
                                // same semi-classical approximation as paper)
    return li4;
}

// ─────────────────────────────────────────────────────────────────────────────
// setSeed
// ─────────────────────────────────────────────────────────────────────────────

void CoalescenceEngine::setSeed(unsigned int seed) {
    fRng.SetSeed(seed + 2);
}
