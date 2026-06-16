#include "ToyMCEngine.h"

#include <iostream>
#include <stdexcept>
#include <sstream>

// ─────────────────────────────────────────────────────────────────────────────
// Histograms helpers
// ─────────────────────────────────────────────────────────────────────────────

void ToyMCEngine::Histograms::book(const std::string& suffix) {
    auto n = [&](const char* base) { return std::string(base) + suffix; };

    hLi4Pt           = new TH1D(n("hLi4Pt").c_str(),
                                 "Li4 p_{T};p_{T} (GeV/c);d^{2}N/(dydp_{T})",
                                 50, 0., 10.);
    hLi4Rapidity     = new TH1D(n("hLi4Rapidity").c_str(),
                                 "Li4 rapidity;y;dN/dy",
                                 60, -3., 3.);
    hProtonPt        = new TH1D(n("hProtonPt").c_str(),
                                 "Proton p_{T};p_{T} (GeV/c);d^{2}N/(dydp_{T})",
                                 50, 0., 5.);
    hHe3Pt           = new TH1D(n("hHe3Pt").c_str(),
                                 "^{3}He p_{T};p_{T} (GeV/c);d^{2}N/(dydp_{T})",
                                 50, 0., 8.);
    hNLi4PerEvent    = new TH1D(n("hNLi4PerEvent").c_str(),
                                 "Li4 per event;N;Counts", 10, 0., 10.);
    hNProtonPerEvent = new TH1D(n("hNProtonPerEvent").c_str(),
                                 "Protons per event;N;Counts", 50, 0., 50.);
    hNHe3PerEvent    = new TH1D(n("hNHe3PerEvent").c_str(),
                                 "He3 per event;N;Counts", 20, 0., 20.);
    hQvsR            = new TH2D(n("hQvsR").c_str(),
                                 "PRF q vs r;q (GeV/c);r (fm)",
                                 100, 0., 2., 100, 0., 20.);

    // Detach from any global ROOT directory — ownership is ours
    for (auto* h : {hLi4Pt, hLi4Rapidity, hProtonPt, hHe3Pt,
                    hNLi4PerEvent, hNProtonPerEvent, hNHe3PerEvent})
        h->SetDirectory(nullptr);
    hQvsR->SetDirectory(nullptr);
}

void ToyMCEngine::Histograms::merge(const Histograms& other) {
    hLi4Pt          ->Add(other.hLi4Pt);
    hLi4Rapidity    ->Add(other.hLi4Rapidity);
    hProtonPt       ->Add(other.hProtonPt);
    hHe3Pt          ->Add(other.hHe3Pt);
    hNLi4PerEvent   ->Add(other.hNLi4PerEvent);
    hNProtonPerEvent->Add(other.hNProtonPerEvent);
    hNHe3PerEvent   ->Add(other.hNHe3PerEvent);
    hQvsR           ->Add(other.hQvsR);
}

void ToyMCEngine::Histograms::normalise(long long nEvents) {
    if (nEvents == 0) return;
    auto normPt = [&](TH1D* hist) {
        hist->Sumw2();
        hist->Scale(1. / nEvents);
        for (int ib = 1; ib <= hist->GetNbinsX(); ++ib) {
            const double w = hist->GetBinWidth(ib);
            hist->SetBinContent(ib, hist->GetBinContent(ib) / w);
            hist->SetBinError  (ib, hist->GetBinError(ib)   / w);
        }
    };
    normPt(hLi4Pt);
    normPt(hProtonPt);
    normPt(hHe3Pt);
}

void ToyMCEngine::Histograms::write() const {
    hLi4Pt          ->Write();
    hLi4Rapidity    ->Write();
    hProtonPt       ->Write();
    hHe3Pt          ->Write();
    hNLi4PerEvent   ->Write();
    hNProtonPerEvent->Write();
    hNHe3PerEvent   ->Write();
    hQvsR           ->Write();
    hProbMap        ->Write();
}

// ─────────────────────────────────────────────────────────────────────────────
// Constructor
// ─────────────────────────────────────────────────────────────────────────────

ToyMCEngine::ToyMCEngine(const SourceSize&              srcProton,
                         const SourceSize&              srcHe3,
                         std::shared_ptr<WignerDensity> wigner,
                         const Config&                  cfg)
    : fCfg(cfg)
    , fSrcProton(srcProton)
    , fSrcHe3(srcHe3)
{
    // ── Load and store input pT histograms ───────────────────────────────────
    fHPtProton = loadHistogram(cfg.rootFile, cfg.hNameProton);
    fHPtHe3    = loadHistogram(cfg.rootFile, cfg.hNameHe3);

    // ── Build the primary engine (used for precompute(); not run directly) ───
    constexpr double kMassProton = 0.938272;
    constexpr double kMassHe3   = 2.808391;
    fCfg.meanNProtons = fHPtProton->Integral() * 2. * cfg.yMax;
    fCfg.meanNHe3     = fHPtHe3->Integral()    * 2. * cfg.yMax;

    fSamplerProtonMom = std::make_unique<MomentumSampler>(
        PDG::kProton, kMassProton, fHPtProton.get(), cfg.yMax, cfg.seed);
    fSamplerHe3Mom = std::make_unique<MomentumSampler>(
        PDG::kHe3, kMassHe3, fHPtHe3.get(), cfg.yMax, cfg.seed + 1);
    fMultModel = std::make_unique<MultiplicityModel>(
        fCfg.meanNProtons, fCfg.meanNHe3, cfg.multMode, cfg.seed + 2);
    fCoalEngine = std::make_unique<CoalescenceEngine>(
        srcProton, srcHe3, std::move(wigner), cfg.seed + 3);
    
}

// ─────────────────────────────────────────────────────────────────────────────
// run  (multithreaded)
// ─────────────────────────────────────────────────────────────────────────────
 
ToyMCEngine::Histograms ToyMCEngine::run(long long nEvents) {
 
    const int nThreads = fCfg.nThreads;
 
    // ── Ensure probability map is precomputed before spawning threads ─────────
    // precompute() is idempotent — safe to call even if already done.
    fCoalEngine->precompute();
    const TH2D* sharedProbMap = fCoalEngine->coalescenceProbability();
 
    // ── Divide events across threads ─────────────────────────────────────────
    const long long evPerThread = nEvents / nThreads;
    const long long remainder   = nEvents % nThreads;
 
    std::vector<std::thread>    workers;
    std::vector<Histograms>     threadHists(nThreads);
 
    for (int t = 0; t < nThreads; ++t)
        threadHists[t].book("_t" + std::to_string(t));
 
    long long evStart = 0;
 
    // Pre-clone histograms on the main thread — one clone per worker.
    // This avoids concurrent TH1::Clone() calls which are not thread-safe
    // due to gROOT registration.
    std::vector<std::unique_ptr<TH1>> protonClones(nThreads);
    std::vector<std::unique_ptr<TH1>> he3Clones(nThreads);
    for (int t = 0; t < nThreads; ++t) {
        protonClones[t].reset(static_cast<TH1*>(
            fHPtProton->Clone(("hPtProton_t" + std::to_string(t)).c_str())));
        he3Clones[t].reset(static_cast<TH1*>(
            fHPtHe3->Clone(("hPtHe3_t" + std::to_string(t)).c_str())));
        protonClones[t]->SetDirectory(nullptr);
        he3Clones[t]->SetDirectory(nullptr);
    }
 
    for (int t = 0; t < nThreads; ++t) {
        const long long evEnd = evStart + evPerThread + (t < remainder ? 1 : 0);
        workers.emplace_back(&ToyMCEngine::workerRun, this,
                             evStart, evEnd, t,
                             sharedProbMap,
                             protonClones[t].get(),
                             he3Clones[t].get(),
                             std::ref(threadHists[t]));
        evStart = evEnd;
    }
 
    for (auto& w : workers) w.join();
 
    // ── Merge per-thread histograms into thread 0's set ──────────────────────
    for (int t = 1; t < nThreads; ++t)
        threadHists[0].merge(threadHists[t]);
    threadHists[0].hProbMap = sharedProbMap; // non-owning pointer for output
 
    threadHists[0].normalise(nEvents);
    return threadHists[0];
}

// ─────────────────────────────────────────────────────────────────────────────
// workerRun
// ─────────────────────────────────────────────────────────────────────────────
 
void ToyMCEngine::workerRun(long long      evStart,
                             long long      evEnd,
                             int            threadId,
                             const TH2D*    sharedProbMap,
                             const TH1*     hPtProton,
                             const TH1*     hPtHe3,
                             Histograms&    threadHists) const
{
    // ── Per-thread objects (no sharing with other threads) ────────────────────
    constexpr double kMassProton = 0.938272;
    constexpr double kMassHe3   = 2.808391;
 
    // Offset seed by threadId to ensure independent RNG sequences
    const unsigned int tSeed = fCfg.seed + threadId * 1000;
 
    // Use pre-cloned histograms — no TH1::Clone() call inside the thread
    MomentumSampler samplerP  (PDG::kProton, kMassProton,
                                hPtProton, fCfg.yMax, tSeed);
    MomentumSampler samplerHe3(PDG::kHe3,    kMassHe3,
                                hPtHe3,    fCfg.yMax, tSeed + 1);
    MultiplicityModel multModel(fCfg.meanNProtons, fCfg.meanNHe3,
                                 fCfg.multMode,    tSeed + 2);
 
    // Clone the engine with offset seed; inject the shared (read-only) probMap
    CoalescenceEngine engine(fSrcProton, fSrcHe3,
                              fCoalEngine->wignerDensity().clone(),
                              tSeed + 3);
    engine.setProbMap(sharedProbMap);
 
    // ── Event loop ────────────────────────────────────────────────────────────
    const long long nLocal   = evEnd - evStart;
    const long long printEvery = std::max(1LL, nLocal / 100);
 
    for (long long iEv = evStart; iEv < evEnd; ++iEv) {
        if ((iEv - evStart) % printEvery == 0)
            std::cout << " [Thread " << threadId << "] Event "
                      << (iEv - evStart) << " / " << nLocal
                      << "  (" << (100 * (iEv - evStart) / nLocal) << "%)\n";
 
        const int nP   = multModel.drawNProtons();
        const int nHe3 = multModel.drawNHe3();
 
        auto protons = samplerP.sampleN(nP);
        auto he3s    = samplerHe3.sampleN(nHe3);
 
        Event ev;
        ev.id      = static_cast<int>(iEv);
        ev.protons = std::move(protons);
        ev.he3s    = std::move(he3s);
 
        CoalescenceResult result = engine.processEvent(ev);
 
        // Fill thread-local histograms — no mutex needed
        for (const auto& p   : ev.protons) threadHists.hProtonPt->Fill(p.pT());
        for (const auto& he3 : ev.he3s)    threadHists.hHe3Pt->Fill(he3.pT());
 
        threadHists.hNProtonPerEvent->Fill(ev.nProtons());
        threadHists.hNHe3PerEvent   ->Fill(ev.nHe3());
        threadHists.hNLi4PerEvent   ->Fill(result.nLi4());
 
        for (const auto& li4 : result.li4s) {
            threadHists.hLi4Pt      ->Fill(li4.pT());
            threadHists.hLi4Rapidity->Fill(li4.y());
        }
        for (const auto& pd : result.pairs)
            threadHists.hQvsR->Fill(pd.q_GeV, pd.r_fm);
    }
}


// ─────────────────────────────────────────────────────────────────────────────
// loadHistogram
// ─────────────────────────────────────────────────────────────────────────────

std::unique_ptr<TH1> ToyMCEngine::loadHistogram(const std::string& filename,
                                                 const std::string& hname)
{
    TFile f(filename.c_str(), "READ");
    if (f.IsZombie())
        throw std::runtime_error(
            "ToyMCEngine: cannot open ROOT file: " + filename);
    TH1* h = dynamic_cast<TH1*>(f.Get(hname.c_str()));
    if (!h)
        throw std::runtime_error(
            "ToyMCEngine: histogram not found: " + hname + " in " + filename);
    auto owned = std::unique_ptr<TH1>(
        static_cast<TH1*>(h->Clone((hname + "_toymc").c_str())));
    owned->SetDirectory(nullptr);
    return owned;
}

// ─────────────────────────────────────────────────────────────────────────────
// setSeed
// ─────────────────────────────────────────────────────────────────────────────

void ToyMCEngine::setSeed(unsigned int seed) {
    fCfg.seed = seed;
    fSamplerProtonMom->setSeed(seed);
    fSamplerHe3Mom->setSeed(seed + 1);
    fMultModel->setSeed(seed + 2);
    fCoalEngine->setSeed(seed + 3);
}