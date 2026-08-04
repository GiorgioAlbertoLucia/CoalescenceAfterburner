#ifndef MOMENTUMSAMPLER_H
#define MOMENTUMSAMPLER_H

#include "Particle.h"

#include "TH1.h"
#include "TRandom3.h"
#include "TLorentzVector.h"

#include <string>
#include <stdexcept>
#include <cmath>

/** 
 * Samples the full 3-momentum of a particle from an input pT histogram,
 * following the ToMCCA approach:
 *
 *   1. pT   is sampled from the input TH1 (treated as a probability density)
 *   2. y    is drawn uniformly from [-yMax, +yMax]  (flat rapidity)
 *   3. phi  is drawn uniformly from [0, 2*pi]       (azimuthal symmetry)
 *
 * The input histogram is expected to be a d^2N/(dy dpT) spectrum (or any
 * shape proportional to it) already corrected for bin widths. It is used
 * as a probability density via ROOT's TH1::GetRandom().
 *
 * A particle is fully constructed with:
 *   px = pT * cos(phi)
 *   py = pT * sin(phi)
 *   pz = mT * sinh(y)       where mT = sqrt(m^2 + pT^2)
 *   E  = mT * cosh(y)
 *
 * The histogram is cloned internally so the caller retains ownership of
 * the original.
 */
class MomentumSampler {
public:

    /**
     * @brief Construct a new Momentum Sampler object
     * @param pdg PDG code of the species to sample
     * @param mass rest mass [GeV/c^2]
     * @param hPt input pT histogram (will be cloned)
     * @param yMax half-width of the flat rapidity window
     * @param seed TRandom3 seed
     */
    MomentumSampler(const double        mass_GeV,
                    const TH1*          hPt,
                    const double        yMax = 1.5,
                    const unsigned int  seed = 42)
        : fRng(seed)
    {
        if (!hPt)
            throw std::invalid_argument("MomentumSampler: hPt is null");
        if (mass_GeV <= 0.)
            throw std::invalid_argument("MomentumSampler: mass must be positive");
        if (yMax <= 0.)
            throw std::invalid_argument("MomentumSampler: yMax must be positive");
            
        fMass = mass_GeV;
        fYMax = yMax;
        fHPt = hPt;
    }

    /**
     * @brief Sample a single particle with 3-momentum from the input pT spectrum
     * @return A Particle with pdg, mom, and pos (pos left at (0,0,0) to be filled later)
     */
    Particle sample() const {
        TVector3 mom;
        sampleMomentum(mom);

        Particle p;
        p.mom = TLorentzVector(mom, 0);
        p.mom.SetE(std::sqrt(fMass * fMass + mom.Mag2()));
        // pos left at (0,0,0) — filled later by SourceSampler
        return p;
    }

    /**
     * @brief Sample N particles with 3-momentum from the input pT spectrum
     * @param n number of particles to sample
     * @return A vector of Particles with pdg, mom, and pos (pos left at (0,0,0) to be filled later)
     */
    std::vector<Particle> sampleN(int n) const {
        std::vector<Particle> out;
        out.reserve(n);
        for (int i = 0; i < n; ++i)
            out.push_back(sample());
        return out;
    }

    /**
     * @brief Sample a 3-momentum vector from the input pT spectrum
     * @param mom output 3-momentum vector (px, py, pz)
     */
    void sampleMomentum(TVector3& mom) const {
        const double pT  = fHPt->GetRandom();
        const double phi = fRng.Uniform(0., 2. * M_PI);
        const double y   = fRng.Uniform(-fYMax, fYMax);

        const double mT  = std::sqrt(fMass * fMass + pT * pT);
        const double px  = pT * std::cos(phi);
        const double py  = pT * std::sin(phi);
        const double pz  = mT * std::sinh(y);

        mom.SetXYZ(px, py, pz);
    }

    void setSeed(unsigned int seed) { fRng.SetSeed(seed); }

    double mass() const { return fMass; }
    double yMax() const { return fYMax; }

private:
    double                 fMass;
    double                 fYMax;
    mutable TRandom3       fRng;
    const TH1*             fHPt{nullptr};
};

#endif // MOMENTUMSAMPLER_H