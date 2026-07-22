#ifndef INTEGRATORA3_H
#define INTEGRATORA3_H

#include "Particle.h"
#include "SourceSize.h"
#include "SourceSampler.h"
#include "MomentumSampler.h"
#include "MultiplicityModel.h"
#include "JacobiTransform.h"
#include "WignerDensityA.h"
#include "CoalescenceEngineHe3.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom3.h"

#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <stdexcept>
#include <thread>
#include <vector>

// ─────────────────────────────────────────────────────────────────────────────
// IntegratorA3
//
// Computes the Nucleus pT spectrum via Monte Carlo integration:
//
//   dN_A/dpT ~ (1/N_samples) * sum_samples D_A({k_j},{r_j}) * delta(pT - pT_Nucleus)
//
// For each sample:
//   1. Draw A=4 momenta from MomentumSampler (importance sampling over G_A)
//   2. Draw A=4 positions from SourceSampler (Gaussian source)
//   3. Transform both to Jacobi coordinates
//   4. Evaluate D_A = prod_{j=1}^{3} D_1(k_j, r_j)
//   5. Fill pT(Nucleus) = |sum p_{T,i}| with weight D_A
//
// G_A is already encoded in the sampling, so it must NOT be multiplied
// explicitly — standard importance sampling.
// ─────────────────────────────────────────────────────────────────────────────

class IntegratorA3 {
public:

    struct Config {
        int              A         = 3;     // number of nucleons
        long long        nEvents   = 1e7;   // MC samples
        unsigned int     seed      = 42;
        double           yMax      = 0.5;   // rapidity window half-width
        double           meanNProtons = -1.0;   // mean number of protons per event
        int              nThreads = 4;  // number of worker threads
        MultiplicityMode multMode = MultiplicityMode::Poissonian;
    };

    IntegratorA3(const SourceSize&    srcSize,
                 std::shared_ptr<WignerDensityA> wigner,
                 const TH1*           hPtNucleon,
                 const Config&        cfg);

    // Run the integration. Fills and returns internal histograms.
    void run(long long nEvents);
    void workerRun(long long      evStart,
                   long long      evEnd,
                   int            threadId,
                   const TH1*     hPtNucleonClone,
                   std::pair<float, float>& threadResult,
                   std::vector<double>& threadYields,
                   QAHistograms& qaClone) const;

    // ── Output histograms ─────────────────────────────────────────────────────
    const TH1D* hNucleusPt()       const { return fHNucleusPt.get();   }
    const TH1D* hNucleonPt()       const { return fHNucleonPt.get();   }
    const TH1D* hSourceDist()      const { return fHSourceDist.get();  }
    const TH2D* hJacobiKvsR()      const { return fHJacobiKvsR.get();  }
    const TH2D* hWignerKvsR()      const { return fHWignerKvsR.get();  }
    const TH1D* hNucleusYield()    const { return fHNucleusYield; }
    const TH1D* hNucleusYieldUncertainty() const { return fHNucleusYieldUncertainty; }
    const TH1D* hNucleusYieldDistribution() const { return fHNucleusYieldDistribution; }

    // Write all histograms to the current ROOT directory
    void write() const;

private:

    const WignerDensityA& wignerDensity() const { return *fWigner; }

    // Draw one sample: fills positions and momenta for A nucleons
    void drawSample(std::vector<TVector3>& positions,
                    std::vector<TVector3>& momenta) const;

    // Compute Nucleus pT from lab-frame momenta
    static double NucleuspT(const std::vector<TVector3>& momenta);

    void writeQAHistograms() const;

    // ── Configuration ─────────────────────────────────────────────────────────
    Config fCfg;

    // ── Physics components ────────────────────────────────────────────────────
    SourceSize                              fSrcSize;
    std::shared_ptr<WignerDensityA>         fWigner;
    std::shared_ptr<SingleGaussianWigner>   fWigner1;

    // ── Samplers ──────────────────────────────────────────────────────────────
    std::unique_ptr<MomentumSampler>        fMomSampler;
    std::unique_ptr<TH1>                    fHPtNucleon;
    mutable TRandom3                        fRng;

    // ── Histograms ────────────────────────────────────────────────────────────
    std::unique_ptr<TH1D> fHNucleusPt;       // Nucleus pT spectrum (weighted)
    std::unique_ptr<TH1D> fHNucleonPt;   // sampled nucleon pT (unweighted, diagnostic)
    std::unique_ptr<TH1D> fHSourceDist;  // sampled inter-nucleon distances [fm]
    std::unique_ptr<TH2D> fHJacobiKvsR;  // |k_j| vs |r_j| for each Jacobi pair
    std::unique_ptr<TH2D> fHWignerKvsR;  // |k_j| vs |r_j| for each Jacobi pair, Wigner evaluated

    TH1D* fHNucleusYield{nullptr}; // Nucleus yield histogram (unweighted)
    TH1D* fHNucleusYieldUncertainty{nullptr}; // Nucleus yield uncertainty
    TH1D* fHNucleusYieldDistribution{nullptr}; // Nucleus yield distribution histogram (unweighted)

    QAHistograms fQAHistograms; // QA histograms for the main thread
};

#endif // INTEGRATORA3_H