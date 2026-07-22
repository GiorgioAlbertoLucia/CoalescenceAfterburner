#include "CMFrameBooster.h"
#include "CoalescenceEngineHe3.h"
#include "Event.h"
#include "IntegratorA3.h"
#include "MultiplicityModel.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

#include "TH1F.h"

// ─────────────────────────────────────────────────────────────────────────────
// Constructor
// ─────────────────────────────────────────────────────────────────────────────

IntegratorA3::IntegratorA3(const SourceSize&               srcSize,
                            std::shared_ptr<WignerDensityA> wigner,
                            const TH1*                      hPtNucleon,
                            const Config&                   cfg)
    : fCfg(cfg)
    , fSrcSize(srcSize)
    , fWigner(std::move(wigner))
    , fRng(cfg.seed)
    , fQAHistograms()
{
    if (!fWigner)
        throw std::invalid_argument("IntegratorA3: WignerDensityA must not be null");
    if (!hPtNucleon)
        throw std::invalid_argument("IntegratorA3: hPtNucleon must not be null");
    fWigner1 = std::make_shared<SingleGaussianWigner>(fWigner->getD());

    constexpr double kMassNucleon = 0.938272; // GeV/c^2
    fHPtNucleon = std::unique_ptr<TH1>(static_cast<TH1*>(hPtNucleon->Clone()));
    fHPtNucleon->SetDirectory(nullptr);

    fMomSampler = std::make_unique<MomentumSampler>(
        PDG::kProton, kMassNucleon, hPtNucleon, cfg.yMax, cfg.seed + 1);
    fCfg.meanNProtons = hPtNucleon->Integral() * 2. * cfg.yMax;

    // ── Book histograms ───────────────────────────────────────────────────────
    fHNucleusPt = std::make_unique<TH1D>(
        "hNucleusPt_A3",
        "^{3}He p_{T};p_{T} (GeV/c);d^{2}N/(dydp_{T})",
        50, 0., 12.);
    fHNucleusPt->SetDirectory(nullptr);

    fHNucleonPt = std::make_unique<TH1D>(
        "hNucleonPt_A3",
        "Sampled nucleon p_{T};p_{T} (GeV/c);Counts",
        50, 0., 5.);
    fHNucleonPt->SetDirectory(nullptr);

    fHSourceDist = std::make_unique<TH1D>(
        "hSourceDist_A3",
        "Inter-nucleon distance;r (fm);Counts",
        100, 0., 20.);
    fHSourceDist->SetDirectory(nullptr);

    fHJacobiKvsR = std::make_unique<TH2D>(
        "hJacobiKvsR_A3",
        "Jacobi |k_{j}| vs |r_{j}|;|k| (GeV/c);|r| (fm)",
        100, 0., 2., 100, 0., 20.);
    fHJacobiKvsR->SetDirectory(nullptr);

    fHWignerKvsR = std::make_unique<TH2D>(
        "hWignerKvsR",
        "Wigner |k_{j}| vs |r_{j}| for a single nucleon pair;|k| (GeV/c);|r| (fm)",
        200, 0., 2., 200, 0., 20.);
    fHWignerKvsR->SetDirectory(nullptr);
    for (int xbin = 1; xbin <= fHWignerKvsR->GetNbinsX(); ++xbin) {
        for (int ybin = 1; ybin <= fHWignerKvsR->GetNbinsY(); ++ybin) {
            const double k = fHWignerKvsR->GetXaxis()->GetBinCenter(xbin);
            const double r = fHWignerKvsR->GetYaxis()->GetBinCenter(ybin);
            const double w = fWigner1->evaluate({TVector3(k, 0, 0)}, {TVector3(r, 0, 0)});
            fHWignerKvsR->SetBinContent(xbin, ybin, w);
        }
    }

    fHNucleusYield = new TH1D(
        "hNucleusYield_A3",
        "^{3}He yield;Yield;Counts",
        1, 0., 1.);
    fHNucleusYield->SetDirectory(nullptr);
    fHNucleusYieldUncertainty = new TH1D(
        "hNucleusYieldUncertainty_A3",
        "^{3}He yield uncertainty;Yield;Counts",
        1, 0., 1.);
    fHNucleusYieldUncertainty->SetDirectory(nullptr);
}

void IntegratorA3::writeQAHistograms() const {
    fQAHistograms.fHPtNucleon->Write();
    fQAHistograms.fHYNucleon->Write();
    fQAHistograms.fHPhiNucleon->Write();
    fQAHistograms.fHRadiusNucleon->Write();
    fQAHistograms.fHPhiCoordinateNucleon->Write();
    fQAHistograms.fHCosThetaCoordinateNucleon->Write();
    fQAHistograms.fHPtNucleus->Write();
    fQAHistograms.fHYNucleus->Write();
    fQAHistograms.fHPhiNucleus->Write();
}

// ─────────────────────────────────────────────────────────────────────────────
// run
// ─────────────────────────────────────────────────────────────────────────────

void IntegratorA3::run(long long nEvents) {
 
    const int nThreads = fCfg.nThreads;
 
    // ── Divide events across threads ─────────────────────────────────────────
    const long long evPerThread = nEvents / nThreads;
    const long long remainder   = nEvents % nThreads;
 
    std::vector<std::thread>    workers;
    std::vector<std::pair<float, float>>     threadResults(nThreads);
 
    long long evStart = 0;
 
    // Pre-clone histograms on the main thread — one clone per worker.
    // This avoids concurrent TH1::Clone() calls which are not thread-safe
    // due to gROOT registration.
    std::vector<std::unique_ptr<TH1>> protonClones(nThreads);
    std::vector<std::vector<double>> threadYields(nThreads);
    std::vector<QAHistograms> qaClones(nThreads);
    for (int t = 0; t < nThreads; ++t) {
        protonClones[t].reset(static_cast<TH1*>(
            fHPtNucleon->Clone(("hPtNucleon_t" + std::to_string(t)).c_str())));
        protonClones[t]->SetDirectory(nullptr);
        qaClones[t] = QAHistograms{};
    }
 
    for (int t = 0; t < nThreads; ++t) {
        const long long evEnd = evStart + evPerThread + (t < remainder ? 1 : 0);
        workers.emplace_back(&IntegratorA3::workerRun, this,
                             evStart, evEnd, t,
                             protonClones[t].get(),
                             std::ref(threadResults[t]),
                             std::ref(threadYields[t]),
                             std::ref(qaClones[t]));
        evStart = evEnd;
    }
 
    for (auto& w : workers) w.join();
 
    // ── Merge per-thread histograms into thread 0's set ──────────────────────
    float totalYield = 0.f, totalEvents = 0.f;
    for (int t = 0; t < nThreads; ++t) {
        totalYield += threadResults[t].first;
        totalEvents += threadResults[t].second;
    }
    std::cout << "Total yield: " << totalYield << " from " << totalEvents << " events\n";
    std::cout << "Average yield per event: " << (totalYield / totalEvents) << "\n";

    const double yield = totalYield / totalEvents;
    double yieldUncertainty = 0.0;
    fHNucleusYieldDistribution = new TH1D(
        "hNucleusYieldDistribution_A3",
        "^{3}He yield distribution;Yield;Counts",
        1000, yield - 10 * yield, yield + 10 * yield);
    for (int t = 0; t < nThreads; ++t) {
        for (const auto& w : threadYields[t]) {
            const double diff = w - yield;
            yieldUncertainty += diff * diff;
            fHNucleusYieldDistribution->Fill(w);
        }
    }
    yieldUncertainty = std::sqrt(yieldUncertainty / (totalEvents - 1));

    fHNucleusYield->SetBinContent(1, yield);
    fHNucleusYieldUncertainty->SetBinContent(1, yieldUncertainty / std::sqrt(totalEvents));

    for (int t = 1; t < nThreads; ++t) {
        qaClones[0].fHPtNucleon->Add(qaClones[t].fHPtNucleon);
        qaClones[0].fHYNucleon->Add(qaClones[t].fHYNucleon);
        qaClones[0].fHPhiNucleon->Add(qaClones[t].fHPhiNucleon);
        qaClones[0].fHRadiusNucleon->Add(qaClones[t].fHRadiusNucleon);
        qaClones[0].fHPhiCoordinateNucleon->Add(qaClones[t].fHPhiCoordinateNucleon);
        qaClones[0].fHCosThetaCoordinateNucleon->Add(qaClones[t].fHCosThetaCoordinateNucleon);
        qaClones[0].fHPtNucleus->Add(qaClones[t].fHPtNucleus);
        qaClones[0].fHYNucleus->Add(qaClones[t].fHYNucleus);
        qaClones[0].fHPhiNucleus->Add(qaClones[t].fHPhiNucleus);
    }
    fQAHistograms = qaClones[0];

}

void IntegratorA3::workerRun(long long      evStart,
                             long long      evEnd,
                             int            threadId,
                             const TH1*     hPtProton,
                             std::pair<float, float>& threadResult,
                             std::vector<double>& threadYields,
                             QAHistograms& qaClone
                            ) const
{

    const long long nEvents   = fCfg.nEvents;
    const int       A          = fCfg.A;

    // SA for Nucleus (all excited states)
    constexpr double SA = 1./12.;
    constexpr double nucleonMass_GeV = 0.938272; // GeV/c^2

    MultiplicityModel multModel(fCfg.meanNProtons, 0.,
                                 fCfg.multMode,    fCfg.seed + 2);
    MomentumSampler samplerP  (PDG::kProton, nucleonMass_GeV,
                                hPtProton, fCfg.yMax, fCfg.seed + 3);
    //CoalescenceEngineNucleus engine(fSrcSize, wignerDensity().clone(), fCfg.seed + 4);
    CoalescenceEngineHe3 engine(fSrcSize, wignerDensity().clone(), fCfg.seed + 4, threadId);
    engine.setWignerMap(fHWignerKvsR.get()); // <-- add this line

    float yield = 0.f;

    // ── Event loop ────────────────────────────────────────────────────────────
    const long long nLocal   = evEnd - evStart;
    const long long printEvery = std::max(1LL, nLocal / 5);
 
    for (long long iEv = evStart; iEv < evEnd; ++iEv) {
        if ((iEv - evStart) % printEvery == 0)
            std::cout << " [Thread " << threadId << "] Event "
                      << (iEv - evStart) << " / " << nLocal
                      << "  (" << (5 * (iEv - evStart) / nLocal) << "%)\n";

        const int nProtons = multModel.drawNProtons();
        const int nNeutrons = multModel.drawNProtons();

        std::vector<Particle> protons(nProtons), neutrons(nNeutrons);
        for (int i = 0; i < nProtons; ++i) {
            protons[i] = samplerP.sample();
            if (i == 0) {
                protons[i].pos.SetXYZ(0., 0., 0.); // First proton at origin
            } else {
                fSrcSize.samplePosition(protons[i].pos);
            }
        }
        for (int i = 0; i < nNeutrons; ++i) {
            neutrons[i] = samplerP.sample();
            fSrcSize.samplePosition(neutrons[i].pos);
        }

        Event ev;
        ev.protons = std::move(protons);
        ev.neutrons = std::move(neutrons);
        const float weight = engine.processEvent(ev);
        yield += weight;
        threadYields.push_back(weight);
    }

    threadResult = std::make_pair(yield, static_cast<float>(nLocal));
    std::cout << " [Thread " << threadId << "] Finished " << nLocal
              << " events, yield = " << yield << "\n";
    qaClone = engine.getQAHistograms();
}

// ─────────────────────────────────────────────────────────────────────────────
// NucleuspT
// ─────────────────────────────────────────────────────────────────────────────

double IntegratorA3::NucleuspT(const std::vector<TVector3>& momenta) {
    double px = 0., py = 0.;
    for (const auto& p : momenta) {
        px += p.X();
        py += p.Y();
    }
    return std::sqrt(px*px + py*py);
}

// ─────────────────────────────────────────────────────────────────────────────
// write
// ─────────────────────────────────────────────────────────────────────────────

void IntegratorA3::write() const {
    fHNucleusPt ->Write();
    fHNucleonPt ->Write();
    fHSourceDist->Write();
    fHJacobiKvsR->Write();
    fHWignerKvsR->Write();
    fHNucleusYield->Write();
    fHNucleusYieldUncertainty->Write();
    fHNucleusYieldDistribution->Write();

    writeQAHistograms();
}