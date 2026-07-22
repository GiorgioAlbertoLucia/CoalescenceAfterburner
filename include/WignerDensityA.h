#ifndef WIGNERDENSITYA_H
#define WIGNERDENSITYA_H

#include "WignerDensity.h"
#include "TVector3.h"

#include <vector>
#include <cmath>
#include <stdexcept>

// ─────────────────────────────────────────────────────────────────────────────
// Gaussian wavefunction implementation
//
// phi(r) = exp(-r^2 / (2d^2)) / (pi * d^2)^(3/4)
//
// The corresponding Wigner density is analytic (paper Eq. 21):
//   D(q, r) = int d^3y phi(r + y) phi*(r - y) exp(-i q . y)
//           = exp( -(d^4 * q^2 + r^2) / d^2 )
//
// where d is the Gaussian width in fm, related to the rms radius by:
//   r_rms = sqrt(3/2) * d
//   -> d = r_rms * sqrt(2/3)
//
// For Li4 with r_rms = 2 fm: d = 2 * sqrt(2/3) ≈ 1.6330 fm
// ─────────────────────────────────────────────────────────────────────────────

class SingleGaussianWigner : public WignerDensity {
public:

    // d_fm: Gaussian width parameter in fm
    explicit SingleGaussianWigner(double d_fm) : fD(d_fm) {
        if (d_fm <= 0.)
            throw std::invalid_argument("SingleGaussianWigner: d must be positive");
        fD2 = d_fm * d_fm;       // d^2  [fm^2]
        fD4 = fD2 * fD2;         // d^4  [fm^4]
    }

    // D(q, r) = exp( -(d^4 * q_fm^2 / 4. + r^2) / d^2 )
    // q is given in GeV/c and converted to fm^-1 via hbar*c
    double evaluate(const TVector3& q_GeV,
                    const TVector3& r_fm) const override {

        // convert q from GeV/c to fm^-1: q[fm^-1] = q[GeV/c] / hbarc
        constexpr double hbarc = 0.197327; // GeV*fm
        const double q2_fm = q_GeV.Mag2() / (hbarc * hbarc); // [fm^-2]
        const double r2_fm = r_fm.Mag2();                     // [fm^2]

        return 8 * std::exp( -(fD4 * q2_fm + r2_fm) / fD2 );
    }

    double getD()  const { return fD;  }

    std::shared_ptr<WignerDensity> clone() const override {
        return std::make_shared<SingleGaussianWigner>(fD);
    }

private:
    double fD;   // [fm]
    double fD2;  // [fm^2]
    double fD4;  // [fm^4]
};


// ─────────────────────────────────────────────────────────────────────────────
// WignerDensityA
//
// A-body Wigner density for a Gaussian wavefunction, computed as a product
// of A-1 identical two-body Gaussian Wigner densities evaluated on the
// Jacobi coordinate pairs (k_j, r_j):
//
//   D_A({k_j}, {r_j}) = prod_{j=1}^{A-1} D_1(k_j, r_j)
//
// where D_1 is the two-body Gaussian Wigner density (Eq. 21 of the paper):
//
//   D_1(k, r) = exp( -(d^4*k^2 + r^2) / d^2 )
//
// The factorisation is exact for the Gaussian wavefunction.
//
// Inputs are A-1 Jacobi relative momenta [GeV/c] and A-1 Jacobi relative
// positions [fm], as returned by JacobiTransform::relative().
// ─────────────────────────────────────────────────────────────────────────────

class WignerDensityA {
public:

    // d_fm : Gaussian width parameter [fm], same as GaussianWigner
    explicit WignerDensityA(double d_fm)
        : fWigner1(d_fm)
    {}

    // Evaluate D_A given A-1 Jacobi momenta [GeV/c] and A-1 Jacobi positions [fm].
    // Both vectors must have the same length (= A-1).
    double evaluate(const std::vector<TVector3>& k_jacobi_GeV,
                    const std::vector<TVector3>& r_jacobi_fm) const {
        if (k_jacobi_GeV.size() != r_jacobi_fm.size())
            throw std::invalid_argument(
                "WignerDensityA::evaluate: k and r vectors must have equal size");

        double D = 1.0;
        for (std::size_t j = 0; j < k_jacobi_GeV.size(); ++j)
            D *= fWigner1.evaluate(k_jacobi_GeV[j], r_jacobi_fm[j]);
        return D;
    }

    // Maximum value: D_max = 8^(A-1)
    double maxValue(int A) const {
        return std::pow(fWigner1.maxValue(), A - 1);
    }

    double getD()   const { return fWigner1.getD(); }

    const SingleGaussianWigner& getWigner1() const { return fWigner1; }

    std::shared_ptr<WignerDensityA> clone() const {
        return std::make_shared<WignerDensityA>(fWigner1.getD());
    }

private:
    SingleGaussianWigner fWigner1; // two-body Gaussian Wigner density
};

#endif // WIGNERDENSITYA_H
