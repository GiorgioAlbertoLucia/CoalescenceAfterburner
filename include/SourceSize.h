#ifndef SOURCESIZE_H
#define SOURCESIZE_H

#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

// ─────────────────────────────────────────────────────────────────────────────
// Preset centrality classes for PbPb at sqrt(sNN) = 5.02 TeV
//
// Parameters are placeholders marked with TODO.
// Replace R and alpha with the values from the ALICE Run 2 publication
// (arXiv:2505.01061 or equivalent) for each centrality class.
//
// Usage:
//   auto [src_p, src_He3] = SourceSizeParameters::PbPb502_0_10();
// ─────────────────────────────────────────────────────────────────────────────

namespace SourceSizeParameters {

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

} // namespace SourceSizeParameters


// ─────────────────────────────────────────────────────────────────────────────
// SourceSize
//
// Encodes the ALICE Run 2 mT-scaling parametrisation for the single-particle
// femtoscopic source radius:
//
//   r0(mT) = R * (mT / mT_ref)^(-alpha)        [fm]
//
// where:
//   mT      = transverse mass of the particle [GeV/c^2]
//   mT_ref  = reference transverse mass, conventionally 1 GeV/c^2
//   R       = source size at mT = mT_ref [fm]
//   alpha   = power-law exponent (positive, typically ~0.3-0.5 for baryons)
//
// One instance is created per particle species (proton, He3).
// The pair source size is derived externally as:
//   r0_pair = sqrt(r0_p^2 + r0_He3^2)
//
// Parameters R and alpha are provided at construction time and should be
// taken from the ALICE Run 2 PbPb measurement at sqrt(sNN) = 5.02 TeV
// for the relevant centrality class.
// ─────────────────────────────────────────────────────────────────────────────

class SourceSize {
public:

    // species label is for bookkeeping only
    SourceSize(const SourceSizeParameters::CentralityClass& centrality,
               const double mT_ref_GeV = 1.0)
        : fCentrality(centrality),
          fMTRef(mT_ref_GeV)
    {
        computeR(mT_ref_GeV, centrality);
    }
    SourceSize() = default;

    // Returns r0(mT) in fm
    double r0() const {
        return fR;
    }

    double mtRef() const { return fMTRef; }
    SourceSizeParameters::CentralityClass centrality() const { return fCentrality; }

    void setMtRef(double mT_ref_GeV) {
        fMTRef = mT_ref_GeV;
        computeR(mT_ref_GeV, fCentrality);
    }

    void setCentrality(const SourceSizeParameters::CentralityClass& centrality) {
        fCentrality = centrality;
        computeR(fMTRef, centrality);
    }

    // Pair source size from two single-particle source sizes [fm]
    // r0_pair = sqrt(r0_a^2 + r0_b^2)
    static double pairR0(double r0a_fm, double r0b_fm) {
        return std::sqrt(r0a_fm * r0a_fm + r0b_fm * r0b_fm);
    }


private:

    void computeR(const double mT_ref_GeV, const SourceSizeParameters::CentralityClass& centrality) {
        fCentrality = centrality;
        fMTRef = mT_ref_GeV;
        const double a = SourceSizeParameters::parameters[static_cast<int>(fCentrality)][0];
        const double b = SourceSizeParameters::parameters[static_cast<int>(fCentrality)][1];
        const double c = SourceSizeParameters::parameters[static_cast<int>(fCentrality)][2];
        fR = a + b * std::pow(fMTRef, c);
    }

    double fR;      // [fm]
    SourceSizeParameters::CentralityClass fCentrality;
    double fMTRef; // [GeV/c^2]
};

class PairSpaceDistribution {
public:
    PairSpaceDistribution(const SourceSize& srcA, const SourceSize& srcB)
        : fSrcA(srcA), fSrcB(srcB) {
        fR0Pair = SourceSize::pairR0(fSrcA.r0(), fSrcB.r0());
    }
    PairSpaceDistribution(const double r0Pair_fm) : fR0Pair(r0Pair_fm) {}
    ~PairSpaceDistribution() = default;

    double pairR0() const {
        return fR0Pair;
    }

    /**
     * Evaluate the pair space distribution at a given relative distance r and pair source size R_pair.
     */
    double evaluate(const double r, const double R_pair) const {
        // Gaussian distribution in relative distance r with width R_pair
        return 1/(2* M_PI * fR0Pair * fR0Pair) * std::exp(-(r * r + R_pair * R_pair) / (4. * fR0Pair * fR0Pair));
    }

private:
    SourceSize fSrcA;
    SourceSize fSrcB;
    double fR0Pair; // [fm]
    
};



#endif // SOURCESIZE_H