#ifndef MOMENTUMSAMPLER_H
#define MOMENTUMSAMPLER_H

#include "Particle.h"

#include "TH1.h"
#include "TRandom3.h"
#include "TLorentzVector.h"

#include <string>
#include <stdexcept>
#include <cmath>

// ─────────────────────────────────────────────────────────────────────────────
// MomentumSampler
//
// Samples the full 3-momentum of a particle from an input pT histogram,
// following the ToMCCA approach:
//
//   1. pT   is sampled from the input TH1 (treated as a probability density)
//   2. y    is drawn uniformly from [-yMax, +yMax]  (flat rapidity)
//   3. phi  is drawn uniformly from [0, 2*pi]       (azimuthal symmetry)
//
// The input histogram is expected to be a d^2N/(dy dpT) spectrum (or any
// shape proportional to it) already corrected for bin widths. It is used
// as a probability density via ROOT's TH1::GetRandom().
//
// A particle is fully constructed with:
//   px = pT * cos(phi)
//   py = pT * sin(phi)
//   pz = mT * sinh(y)       where mT = sqrt(m^2 + pT^2)
//   E  = mT * cosh(y)
//
// The histogram is cloned internally so the caller retains ownership of
// the original.
// ─────────────────────────────────────────────────────────────────────────────

class MomentumSampler {
public:

    // pdg     : PDG code of the species to sample
    // mass    : rest mass [GeV/c^2]
    // hPt     : input pT histogram (will be cloned)
    // yMax    : half-width of the flat rapidity window (default 0.5)
    // seed    : TRandom3 seed
    MomentumSampler(int           pdg,
                    double        mass_GeV,
                    const TH1*    hPt,
                    double        yMax = 0.5,
                    unsigned int  seed = 42)
        : fPdg(pdg)
        , fMass(mass_GeV)
        , fYMax(yMax)
        , fRng(seed)
    {
        if (!hPt)
            throw std::invalid_argument("MomentumSampler: hPt is null");
        if (mass_GeV <= 0.)
            throw std::invalid_argument("MomentumSampler: mass must be positive");
        if (yMax <= 0.)
            throw std::invalid_argument("MomentumSampler: yMax must be positive");

        fHPt = hPt;
    }

    // Sample one particle. pos is left at (0,0,0) for SourceSampler to fill.
    Particle sample() const {
        const double pT  = fHPt->GetRandom();
        const double phi = fRng.Uniform(0., 2. * M_PI);
        const double y   = fRng.Uniform(-fYMax, fYMax);

        const double mT  = std::sqrt(fMass * fMass + pT * pT);
        const double px  = pT * std::cos(phi);
        const double py  = pT * std::sin(phi);
        const double pz  = mT * std::sinh(y);
        const double E   = mT * std::cosh(y);

        Particle p;
        p.pdg = fPdg;
        p.mom.SetPxPyPzE(px, py, pz, E);
        // pos left at (0,0,0) — filled later by SourceSampler
        return p;
    }

    // Sample n particles at once
    std::vector<Particle> sampleN(int n) const {
        std::vector<Particle> out;
        out.reserve(n);
        for (int i = 0; i < n; ++i)
            out.push_back(sample());
        return out;
    }

    void setSeed(unsigned int seed) { fRng.SetSeed(seed); }

    int    pdg()  const { return fPdg;  }
    double mass() const { return fMass; }
    double yMax() const { return fYMax; }

private:
    int                    fPdg;
    double                 fMass;
    double                 fYMax;
    mutable TRandom3       fRng;
    const TH1*   fHPt{nullptr};
};

#endif // MOMENTUMSAMPLER_H
