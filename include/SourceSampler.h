#ifndef SOURCESAMPLER_H
#define SOURCESAMPLER_H

#include "Particle.h"
#include "SourceSize.h"

#include "TRandom3.h"
#include "TVector3.h"

#include <stdexcept>

// ─────────────────────────────────────────────────────────────────────────────
// SourceSampler
//
// Assigns a 3D emission position to a Particle by sampling from a spherically
// symmetric Gaussian source:
//
//   P(r) ~ exp( -r^2 / (2 * r0^2) )
//
// where r0 = r0(mT) is provided by a SourceSize instance.
//
// Each spatial component is sampled independently as:
//   x_i ~ Gaussian(0, r0)
//
// so the 3D distance distribution has:
//   <r^2> = 3 * r0^2
//   sigma_1D = r0
//
// This is consistent with the Gaussian source model used in the paper
// (Eq. 16), where Hnp factorises into single-particle distributions h(r)
// each of width r0.
//
// The mean inter-particle distance r_mu relates to r0 through:
//   r_mu = (4/sqrt(pi)) * r0        (for the 3D distance distribution)
// which is the same relation used in the paper (Sec. 2).
//
// A single TRandom3 instance is owned by the sampler and seeded at
// construction. Pass seed=0 for a time-based random seed.
// ─────────────────────────────────────────────────────────────────────────────

class SourceSampler {
public:

    // srcSize : SourceSize instance for the species being sampled
    // seed    : seed for TRandom3 (0 = time-based)
    explicit SourceSampler(const SourceSize& srcSize, unsigned int seed = 42)
        : fSrcSize(srcSize)
        , fRng(seed)
    {}

    // Sample a 3D emission position for the given particle and assign it
    // to particle.pos [fm].
    // The mT of the particle is used to determine r0.
    void sample(Particle& particle) const {
        const double mT = particle.mT();
        if (mT <= 0.)
            throw std::runtime_error(
                "SourceSampler::sample : particle has non-positive mT");

        const double r0 = fSrcSize.r0();
        particle.pos = samplePosition(r0);
    }

    // Sample a position vector directly given an explicit r0 [fm].
    // Useful when the caller has already computed the source size.
    TVector3 samplePosition(double r0_fm) const {
        // Each component ~ Gaussian(0, r0)
        // TRandom3::Gaus(mean, sigma)
        const double x = fRng.Gaus(0., r0_fm);
        const double y = fRng.Gaus(0., r0_fm);
        const double z = fRng.Gaus(0., r0_fm);
        return TVector3(x, y, z);
    }

    // Pair source size derived from two single-particle source sizes.
    // r0_pair = sqrt(r0_a^2 + r0_b^2)
    // Convenience wrapper around SourceSize::pairR0.
    static double pairR0(double r0a_fm, double r0b_fm) {
        return SourceSize::pairR0(r0a_fm, r0b_fm);
    }

    // Reseed the internal RNG
    void setSeed(unsigned int seed) { fRng.SetSeed(seed); }

    // Access the underlying SourceSize
    const SourceSize& sourceSize() const { return fSrcSize; }

private:
    const SourceSize& fSrcSize; // reference: lifetime managed by caller
    mutable TRandom3  fRng;     // mutable: sampling is logically const
};

#endif // SOURCESAMPLER_H