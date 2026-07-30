#ifndef WIGNERDENSITY_H
#define WIGNERDENSITY_H

#include "TVector3.h"

#include <vector>
#include <cmath>
#include <stdexcept>



/**
 * Abstract base class
 * Subclasses implement the 3D Wigner density D(q, r) for a given wavefunction.
 * Reference paper: https://doi.org/10.48550/arXiv.2302.12696
 *
 * Convention (following the paper, Eq. 13):
 *   q  = relative momentum of the pair in the PRF [GeV/c]
 *   r  = relative position of the pair in the PRF [fm]
 *
 * Note on units: q is in GeV/c and r is in fm.
 * The conversion factor is: 1 GeV/c * 1 fm = 1 / (hbar*c) with
 * hbar*c = 0.197327 GeV*fm, so q*r is dimensionless as required.
 * */
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


/** 
 * Gaussian wavefunction implementation
 *
 * phi_d(r) = exp(-r^2 / 2d^2) / (pi * d^2)^(3/4)
 *
 * The corresponding Wigner density is analytic (paper Eq. 21):
 * 
 *   D(q, r) = 8 * exp( -(d^4 * q^2 + r^2) / d^2 )
 *
 * where d is the Gaussian width in fm, related to the rms radius by:
 *   r_rms = sqrt(3/2) * d
 *   -> d = r_rms * sqrt(2/3)
 *
 * For Li4 with r_rms = 2 fm: d = 2 * sqrt(2/3) ≈ 1.6330 fm
 */
class GaussianWigner : public WignerDensity {
public:

    /**
     * @param d_fm Gaussian width parameter in fm. Must be positive.
     */
    explicit GaussianWigner(double d_fm) : fD(d_fm) {
        if (d_fm <= 0.)
            throw std::invalid_argument("GaussianWigner: d must be positive");
        fD2 = d_fm * d_fm;       // d^2  [fm^2]
        fD4 = fD2 * fD2;         // d^4  [fm^4]
    }

    /**
     * @brief Evaluate the Wigner density D(q, r) for given relative momentum and position
     */
    double evaluate(const TVector3& q_GeV,
                    const TVector3& r_fm) const override {

        constexpr double hbarc = 0.197327; // GeV*fm
        const double q2_fm = q_GeV.Mag2() / (hbarc * hbarc); // [fm^-2]
        const double r2_fm = r_fm.Mag2();                     // [fm^2]

        return 8.0 * std::exp( -(fD4 * q2_fm + r2_fm) / fD2 );
    }

    double getD()  const { return fD; }

    std::shared_ptr<WignerDensity> clone() const override {
        return std::make_shared<GaussianWigner>(fD);
    }

private:
    double fD;   // [fm]
    double fD2;  // [fm^2]
    double fD4;  // [fm^4]
};


/**
 * A-body Wigner density for a Gaussian wavefunction, computed as a product
 * of A-1 identical two-body Gaussian Wigner densities evaluated on the
 * Jacobi coordinate pairs (k_j, r_j):
 *
 *   D_A({k_j}, {r_j}) = prod_{j=1}^{A-1} D_1(k_j, r_j)
 *
 * where D_1 is the two-body Gaussian Wigner density (Eq. 21 of the paper):
 *
 *   D_1(k, r) = exp( -(d^4*k^2 + r^2) / d^2 )
 *
 * The factorisation is exact for the Gaussian wavefunction.
 *
 * Inputs are A-1 Jacobi relative momenta [GeV/c] and A-1 Jacobi relative
 * positions [fm], as returned by JacobiTransform::relative().
 */
class GaussianWignerDensityA {
public:

    /** \brief Constructor
     *  \param d_fm Gaussian width parameter [fm]
     */
    explicit GaussianWignerDensityA(double d_fm)
        : fWigner1(d_fm)
    {}

    /** \brief Evaluate the A-body Wigner density
     *  \param k_jacobi_GeV A-1 Jacobi relative momenta [GeV/c]
     *  \param r_jacobi_fm A-1 Jacobi relative positions [fm]
     *  \return The Wigner density D_A({k_j}, {r_j})
     */
    double evaluate(const std::vector<TVector3>& k_jacobi_GeV,
                    const std::vector<TVector3>& r_jacobi_fm) const {
        if (k_jacobi_GeV.size() != r_jacobi_fm.size())
            throw std::invalid_argument(
                "GaussianWignerDensityA::evaluate: k and r vectors must have equal size");

        double D = 1.0;
        for (std::size_t j = 0; j < k_jacobi_GeV.size(); ++j)
            D *= fWigner1.evaluate(k_jacobi_GeV[j], r_jacobi_fm[j]);
        return D;
    }

    /** \brief Get the maximum value of the A-body Wigner density
     *  \param A Number of particles
     *  \return The maximum value D_max = 8^(A-1)
     */
    double maxValue(int A) const {
        return std::pow(fWigner1.maxValue(), A - 1);
    }

    double getD()   const { return fWigner1.getD(); }

    const GaussianWigner& getWigner1() const { return fWigner1; }

    std::shared_ptr<GaussianWignerDensityA> clone() const {
        return std::make_shared<GaussianWignerDensityA>(fWigner1.getD());
    }

private:
    GaussianWigner fWigner1; // two-body Gaussian Wigner density
};

#endif // WIGNERDENSITY_H