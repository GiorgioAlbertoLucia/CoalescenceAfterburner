#ifndef TOYMCENGINE_H
#define TOYMCENGINE_H

#include "Particle.h"
#include "MomentumSampler.h"
#include "MultiplicityModel.h"
#include "SourceSize.h"
#include "SourceSampler.h"
#include "WignerDensity.h"
#include "CoalescenceEngine.h" // reuse CoalescenceResult, PairDiagnostic

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"

#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <stdexcept>
#include <thread>
#include <vector>

// ─────────────────────────────────────────────────────────────────────────────
// ToyMCEngine
//
// Self-contained toy Monte Carlo coalescence afterburner inspired by ToMCCA.
// Unlike the FIST-based CoalescenceEngine which reads pre-generated events,
// ToyMCEngine generates proton and He3 momenta internally from input ROOT
// pT histograms, assigns emission positions via Gaussian source sampling,
// and runs the Wigner-function coalescence on every event.
//
// The input pT histograms are read from a ROOT file at construction time.
// The mean particle multiplicities per event are derived from the histogram
// integrals scaled by the rapidity window and a user-provided yield factor.
//
// Workflow per event:
//   1. Draw N_p  ~ Poisson(<N_p>)  protons
//   2. Draw N_He3 ~ Poisson(<N_He3>) He3s
//   3. Sample 3-momenta from hPtProton, hPtHe3
//   4. Sample emission positions via Gaussian source (SourceSampler)
//   5. Run coalescence: loop over He3-p pairs, evaluate Wigner density,
//      accept/reject as Li4
//
// The coalescence step is fully delegated to CoalescenceEngine, so all
// physics (PRF boost, Wigner density, rejection sampling) is shared
// between the two implementations.
// ─────────────────────────────────────────────────────────────────────────────

class ToyMCEngine {
public:

    struct Config {
        std::string rootFile;        // path to ROOT file with input histograms
        std::string hNameProton;     // name of proton pT histogram in file
        std::string hNameHe3;        // name of He3 pT histogram in file

        double meanNProtons = 1.0;   // mean number of protons per event
        double meanNHe3     = 0.01;  // mean number of He3 per event
                                     // (set from integrated yield * dY)

        double yMax         = 0.5;   // rapidity window half-width

        MultiplicityMode multMode = MultiplicityMode::Poissonian;

        unsigned int seed     = 42;
        int          nThreads = 4;  // number of worker threads
    };

    // Per-thread histogram bundle — filled independently, merged at the end.
    // Declared public so main can read the merged result directly.
    struct Histograms {
        TH1D* hLi4Pt          = nullptr;
        TH1D* hLi4Rapidity    = nullptr;
        TH1D* hProtonPt       = nullptr;
        TH1D* hHe3Pt          = nullptr;
        TH1D* hNLi4PerEvent   = nullptr;
        TH1D* hNProtonPerEvent= nullptr;
        TH1D* hNHe3PerEvent   = nullptr;
        TH2D* hQvsR           = nullptr;
        const TH2D* hProbMap  = nullptr; // from CoalescenceEngine precompute()

        void book(const std::string& suffix);
        void merge(const Histograms& other); // Add other into this
        void normalise(long long nEvents);
        void write() const;
    };

    // srcProton, srcHe3 : source size parametrisations (same as FIST impl.)
    // wigner            : Wigner density (shared with CoalescenceEngine)
    // cfg               : configuration (file paths, multiplicities, etc.)
    ToyMCEngine(const SourceSize&              srcProton,
                const SourceSize&              srcHe3,
                std::shared_ptr<WignerDensity> wigner,
                const Config&                  cfg);

    // Run nEvents toy events across nThreads worker threads.
    // Each thread fills its own Histograms; results are merged and returned.
    Histograms run(long long nEvents);

    // Access internal components for diagnostics
    const MultiplicityModel& multiplicityModel() const { return *fMultModel; }
    const CoalescenceEngine& coalescenceEngine() const { return *fCoalEngine; }

    void setSeed(unsigned int seed);

private:

    // Work executed by each thread: processes [evStart, evEnd) events.
    // Fills threadHists in place. Uses its own cloned engine and samplers.
    void workerRun(long long          evStart,
                   long long          evEnd,
                   int                threadId,
                   const TH2D*        sharedProbMap,
                   const TH1*         hPtProton,
                   const TH1*         hPtHe3,
                   Histograms&        threadHists) const;

    // Load histogram from ROOT file, clone it, return owned pointer
    static std::unique_ptr<TH1> loadHistogram(const std::string& filename,
                                               const std::string& hname);

    // ── Owned components ─────────────────────────────────────────────────────
    std::unique_ptr<MomentumSampler>   fSamplerProtonMom;
    std::unique_ptr<MomentumSampler>   fSamplerHe3Mom;
    std::unique_ptr<MultiplicityModel> fMultModel;
    std::unique_ptr<CoalescenceEngine> fCoalEngine;

    // Keep a copy of config and source sizes for cloning in worker threads
    Config     fCfg;
    SourceSize fSrcProton;
    SourceSize fSrcHe3;

    // Cloned input histograms for worker-thread momentum sampling
    std::unique_ptr<TH1> fHPtProton;
    std::unique_ptr<TH1> fHPtHe3;
};

#endif // TOYMCENGINE_H