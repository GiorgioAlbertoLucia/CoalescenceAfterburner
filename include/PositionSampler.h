#ifndef POSITIONSAMPLER_H
#define POSITIONSAMPLER_H

#include "Particle.h"

#include "TRandom3.h"
#include "TVector3.h"

#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

/**
 * Preset centrality classes for PbPb at sqrt(sNN) = 5.02 TeV
 *
 * Parameters are placeholders marked with TODO.
 * Replace R and alpha with the values from the ALICE Run 2 publication
 * (arXiv:2505.01061 or equivalent) for each centrality class.
 *
 * Usage:
 *   auto [src_p, src_He3] = SourceSizeParameters::PbPb502_0_10();
 */
namespace Source {

    enum class CentralityClass {
        PbPb502_0_10 = 0,
        PbPb502_10_30,
        PbPb502_30_50,
        AllCentralities
    };

    inline constexpr std::array<std::array<float, 3>, 3> parameters = {{
        {{3.59f, 2.78f, -1.91f}},
        {{3.38f, 2.21f, -3.16f}},
        {{2.50f, 1.28f, -2.13f}}
    }};

    inline float getR(const std::array<float, 3>& params, float mT) {
        return params[0] * std::pow(mT, params[1]) * std::exp(params[2] * mT);
    }
    inline float getR(CentralityClass centrality, float mT) {
        if (centrality == CentralityClass::AllCentralities) {
            throw std::invalid_argument("SourceSizeParameters: AllCentralities does not have a single R(mT)");
        }
        return getR(parameters[static_cast<int>(centrality)], mT);
    }

} // namespace Source

class PositionSampler {
public:

    /** 
     * @brief Construct a PositionSampler with a nominal radius and random seed
     * @param nominalRadius The nominal radius of the Gaussian source [fm]
     * @param seed The random seed for the internal RNG
     */
    explicit PositionSampler(const double nominalRadius, unsigned int seed = 42)
        : fNominalRadius(nominalRadius)
        , fRng(seed)
    {}

    /**
     * @brief Sample a Particle's position from a Gaussian source with width r0 [fm]
     * @param particle Particle whose position will be sampled
     * @param useMt If true, the source size will be scaled with the particle's transverse mass (mT)
     *               using the ALICE Run 2 parameterization. If false, the nominal radius will be used.
     *               Default is false.
     * @return void (the particle's position is modified in place)
     */
    void sample(Particle& particle, const bool useMt = false) const {
        double r0 = fNominalRadius;
        if (useMt) {
            const double mT = particle.mom.Mt();
            r0 = Source::getR(Source::CentralityClass::AllCentralities, mT);
        }
        particle.pos = samplePosition(r0);
    }

    Particle sample(const bool useMt = false) const {
        Particle p;
        sample(p, useMt);
        return p;
    }

    void sampleN(std::vector<Particle>& particles, const bool useMt = false) const {
        for (auto& p : particles) {
            sample(p, useMt);
        }
    }

    std::vector<Particle> sampleN(int n, const bool useMt = false) const {
        std::vector<Particle> particles;
        particles.reserve(n);
        for (int i = 0; i < n; ++i) {
            particles.push_back(sample(useMt));
        }
        return particles;
    }

    /** 
     * @brief Sample a 3D emission position from a Gaussian source with width r0 [fm]
     * @param r0_fm Gaussian source width [fm]
     * @return A TVector3 representing the sampled position [fm]
     */
    TVector3 samplePosition(double r0_fm = -1.) const {

        if (r0_fm <= 0.) {
            r0_fm = fNominalRadius;
        }

        const double x = fRng.Gaus(0., r0_fm);
        const double y = fRng.Gaus(0., r0_fm);
        const double z = fRng.Gaus(0., r0_fm);
        return TVector3(x, y, z);
    }

    void setSeed(unsigned int seed) { fRng.SetSeed(seed); }
    double getNominalRadius() const { return fNominalRadius; }

private:
    double fNominalRadius; 
    mutable TRandom3  fRng;
};

#endif // POSITIONSAMPLER_H