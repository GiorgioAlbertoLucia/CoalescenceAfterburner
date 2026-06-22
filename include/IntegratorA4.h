#ifndef INTEGRATORA4_H
#define INTEGRATORA4_H

#include "Particle.h"
#include "SourceSize.h"
#include "SourceSampler.h"
#include "MomentumSampler.h"
#include "JacobiTransform.h"
#include "WignerDensityA.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TRandom3.h"

#include <memory>
#include <vector>

// ─────────────────────────────────────────────────────────────────────────────
// IntegratorA4
//
// Computes the Li4 pT spectrum via Monte Carlo integration:
//
//   dN_A/dpT ~ (1/N_samples) * sum_samples D_A({k_j},{r_j}) * delta(pT - pT_Li4)
//
// For each sample:
//   1. Draw A=4 momenta from MomentumSampler (importance sampling over G_A)
//   2. Draw A=4 positions from SourceSampler (Gaussian source)
//   3. Transform both to Jacobi coordinates
//   4. Evaluate D_A = prod_{j=1}^{3} D_1(k_j, r_j)
//   5. Fill pT(Li4) = |sum p_{T,i}| with weight D_A
//
// G_A is already encoded in the sampling, so it must NOT be multiplied
// explicitly — standard importance sampling.
// ─────────────────────────────────────────────────────────────────────────────

class IntegratorA4 {
public:

    struct Config {
        int          A         = 4;     // number of nucleons
        long long    nSamples  = 1e7;   // MC samples
        unsigned int seed      = 42;
        double       yMax      = 0.5;   // rapidity window half-width
    };

    IntegratorA4(const SourceSize&    srcSize,
                 std::shared_ptr<WignerDensityA> wigner,
                 const TH1*           hPtNucleon,
                 const Config&        cfg);

    // Run the integration. Fills and returns internal histograms.
    void run();

    // ── Output histograms ─────────────────────────────────────────────────────
    const TH1D* hLi4Pt()           const { return fHLi4Pt.get();       }
    const TH1D* hNucleonPt()       const { return fHNucleonPt.get();   }
    const TH1D* hSourceDist()      const { return fHSourceDist.get();  }
    const TH2D* hJacobiKvsR()      const { return fHJacobiKvsR.get();  }

    // Write all histograms to the current ROOT directory
    void write() const;

private:

    // Draw one sample: fills positions and momenta for A nucleons
    void drawSample(std::vector<TVector3>& positions,
                    std::vector<TVector3>& momenta) const;

    // Compute Li4 pT from lab-frame momenta
    static double li4pT(const std::vector<TVector3>& momenta);

    // ── Configuration ─────────────────────────────────────────────────────────
    Config fCfg;

    // ── Physics components ────────────────────────────────────────────────────
    SourceSize                     fSrcSize;
    std::shared_ptr<WignerDensityA> fWigner;

    // ── Samplers ──────────────────────────────────────────────────────────────
    std::unique_ptr<MomentumSampler> fMomSampler;
    mutable TRandom3                 fRng;

    // ── Histograms ────────────────────────────────────────────────────────────
    std::unique_ptr<TH1D> fHLi4Pt;       // Li4 pT spectrum (weighted)
    std::unique_ptr<TH1D> fHNucleonPt;   // sampled nucleon pT (unweighted, diagnostic)
    std::unique_ptr<TH1D> fHSourceDist;  // sampled inter-nucleon distances [fm]
    std::unique_ptr<TH2D> fHJacobiKvsR;  // |k_j| vs |r_j| for each Jacobi pair
    std::unique_ptr<TH2D> fHWignerKvsR;  // |k_j| vs |r_j| for each Jacobi pair, Wigner evaluated
};

#endif // INTEGRATORA4_H