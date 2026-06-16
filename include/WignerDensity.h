#ifndef WIGNERDENSITY_H
#define WIGNERDENSITY_H

#include "TVector3.h"
#include <cmath>
#include <stdexcept>

// ─────────────────────────────────────────────────────────────────────────────
// Abstract base class
// Subclasses implement the 3D Wigner density D(q, r) for a given wavefunction.
//
// Convention (following the paper, Eq. 13):
//   q  = relative momentum of the pair in the PRF [GeV/c]
//   r  = relative position of the pair in the PRF [fm]
//
// Note on units: q is in GeV/c and r is in fm.
// The conversion factor is: 1 GeV/c * 1 fm = 1 / (hbar*c) with
// hbar*c = 0.197327 GeV*fm, so q*r is dimensionless as required.
// ─────────────────────────────────────────────────────────────────────────────

class WignerDensity {
public:
    virtual ~WignerDensity() = default;

    // Returns D(q, r). Maximum value is 8 (see footnote 4 of the paper).
    virtual double evaluate(const TVector3& q_GeV,
                            const TVector3& r_fm) const = 0;

    // Maximum value of D, used for rejection sampling normalisation.
    // For all physical wavefunctions D_max = 8.
    virtual double maxValue() const { return 8.0; }

    virtual std::shared_ptr<WignerDensity> clone() const = 0;

};


// ─────────────────────────────────────────────────────────────────────────────
// Gaussian wavefunction implementation
//
// phi_d(r) = exp(-r^2 / 2d^2) / (pi * d^2)^(3/4)
//
// The corresponding Wigner density is analytic (paper Eq. 21):
//   D(q, r) = 8 * exp( -(d^4 * q^2 + r^2) / d^2 )
//
// where d is the Gaussian width in fm, related to the rms radius by:
//   r_rms = sqrt(3/2) * d
//   -> d = r_rms * sqrt(2/3)
//
// For Li4 with r_rms = 2 fm: d = 2 * sqrt(2/3) ≈ 1.6330 fm
// ─────────────────────────────────────────────────────────────────────────────

class GaussianWigner : public WignerDensity {
public:

    // d_fm: Gaussian width parameter in fm
    explicit GaussianWigner(double d_fm) : fD(d_fm) {
        if (d_fm <= 0.)
            throw std::invalid_argument("GaussianWigner: d must be positive");
        fD2 = d_fm * d_fm;       // d^2  [fm^2]
        fD4 = fD2 * fD2;         // d^4  [fm^4]
    }

    // D(q, r) = 8 * exp( -(d^4 * q_fm^2 + r^2) / d^2 )
    // q is given in GeV/c and converted to fm^-1 via hbar*c
    double evaluate(const TVector3& q_GeV,
                    const TVector3& r_fm) const override {

        // convert q from GeV/c to fm^-1: q[fm^-1] = q[GeV/c] / hbarc
        constexpr double hbarc = 0.197327; // GeV*fm
        const double q2_fm = q_GeV.Mag2() / (hbarc * hbarc); // [fm^-2]
        const double r2_fm = r_fm.Mag2();                     // [fm^2]

        return 8.0 * std::exp( -(fD4 * q2_fm + r2_fm) / fD2 );
    }

    double getD()  const { return fD;  }
    double getRms() const { return fD * std::sqrt(1.5); } // r_rms = sqrt(3/2)*d

    std::shared_ptr<WignerDensity> clone() const override {
        return std::make_shared<GaussianWigner>(fD);
    }

private:
    double fD;   // [fm]
    double fD2;  // [fm^2]
    double fD4;  // [fm^4]
};

#endif // WIGNERDENSITY_H