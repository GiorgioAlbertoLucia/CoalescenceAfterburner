#ifndef WIGNERDENSITY_H
#define WIGNERDENSITY_H

#include "TH2.h"
#include "TFile.h"
#include "TString.h"
#include "TVector3.h"

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>



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

/**
 * Wigner density read from a tabulated TH2 (q, r) -> D.
 *
 * Use this when the internal wavefunction does not admit a closed-form
 * Wigner density (Hulthen, Argonne v18, chi-EFT ... see Appendix A.2/A.3 and
 * Fig. 10 of arXiv:2302.12696, where D(q,r) has a genuine q-r correlation
 * that GaussianWigner's factorized ansatz cannot capture). The table must be
 * produced externally (e.g. by numerically evaluating Eq. 13 on a grid of
 * |q|, |r|) and stored as a ROOT TH2.
 *
 * Convention (must match how the table was produced):
 *   x axis: |r|  relative position magnitude in the PRF [fm]
 *   y axis: |q|  relative momentum magnitude in the PRF [GeV/c]
 *   bin content: D(r, q), same normalisation as GaussianWigner
 *                (D(0,0) -> 8 for a properly normalised wavefunction)
 *
 * Only the magnitudes of q and r are used -- isotropic/s-wave assumption,
 * same as GaussianWigner and as Eq. 21 of the paper.
 *
 * Evaluation uses TH2::Interpolate (bilinear). A query outside the range
 * spanned by the table's bin centers returns 0 rather than extrapolating --
 * make sure the table's range comfortably covers what PositionSampler /
 * MomentumSampler actually populate, or the tails will be silently cut.
 */
class HistogramWigner : public WignerDensity {
public:

    /**
     * @param hDensity TH2 with x = r [fm], y = q [GeV/c], content = D(r, q).
     *                 A private clone is taken; the caller keeps ownership
     *                 of the original and may delete/close its file freely
     *                 afterwards.
     */
    explicit HistogramWigner(const TH2* hDensity) {
        if (!hDensity)
            throw std::invalid_argument("HistogramWigner: hDensity is null");
        init(hDensity);
    }

    /**
     * @brief Load the Wigner-density table from a ROOT file.
     * @param filePath path to the .root file
     * @param histName name of the TH2 inside the file
     */
    HistogramWigner(const std::string& filePath, const std::string& histName) {
        TFile* file = TFile::Open(filePath.c_str());
        if (!file || file->IsZombie())
            throw std::runtime_error("HistogramWigner: failed to open file: " + filePath);

        TH2* h = dynamic_cast<TH2*>(file->Get(histName.c_str()));
        if (!h) {
            file->Close();
            delete file;
            throw std::runtime_error("HistogramWigner: histogram not found or not a TH2: " + histName);
        }
        init(h);
        file->Close();
        delete file;
    }

    ~HistogramWigner() override {
        if (fHDensity) delete fHDensity;
    }

    // Owns a raw TH2* -- use clone() rather than copying.
    HistogramWigner(const HistogramWigner&) = delete;
    HistogramWigner& operator=(const HistogramWigner&) = delete;

    /**
     * @brief Evaluate the Wigner density at a given relative momentum and position.
     * @param q_GeV relative momentum in the PRF [GeV/c]
     * @param r_fm relative position in the PRF [fm]
     * @return the Wigner density D(r, q)
     */
    double evaluate(const TVector3& q_GeV, const TVector3& r_fm) const override {

        const double q = q_GeV.Mag();
        const double r = r_fm.Mag();

        if (q < fQMin || q > fQMax || r < fRMin || r > fRMax) {
            return 0.0; // outside the tabulated range: treat as non-coalescing
        }

        const double D = fHDensity->Interpolate(r, q);
        return D > 0. ? D : 0.0; // guard against small negative interpolation artefacts
    }

    double maxValue() const override { return fDMax; }

    std::shared_ptr<WignerDensity> clone() const override {
        return std::make_shared<HistogramWigner>(fHDensity);
    }

    const TH2* histogram() const { return fHDensity; }

private:

    void init(const TH2* hDensity) {
        fHDensity = static_cast<TH2*>(hDensity->Clone());
        fHDensity->SetDirectory(nullptr);

        // Interpolate() is only well-defined within [first, last] bin centers
        fRMin = fHDensity->GetXaxis()->GetBinCenter(1);
        fRMax = fHDensity->GetXaxis()->GetBinCenter(fHDensity->GetNbinsX());
        fQMin = fHDensity->GetYaxis()->GetBinCenter(1);
        fQMax = fHDensity->GetYaxis()->GetBinCenter(fHDensity->GetNbinsY());

        fDMax = fHDensity->GetMaximum();
        if (fDMax <= 0.)
            throw std::invalid_argument("HistogramWigner: histogram maximum is <= 0 -- is the table filled?");
    }

    TH2*   fHDensity{nullptr};
    double fQMin{0.}, fQMax{0.};
    double fRMin{0.}, fRMax{0.};
    double fDMax{0.};
};

#endif // WIGNERDENSITY_H