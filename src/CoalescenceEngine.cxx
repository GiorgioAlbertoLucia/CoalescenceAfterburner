#include "CoalescenceEngine.h"

#include "TRandom3.h"
#include "TFile.h"
#include "TString.h"
#include "TLorentzVector.h"

#include <thread>

CoalescenceEngine::CoalescenceEngine(Config& config)
    : fConfig(config)
{

    TFile* file = TFile::Open(config.inputPtHistogramFile.c_str());
    if (!file) {
        throw std::runtime_error("Failed to open input pT histogram file");
    }
    fHPt = static_cast<TH1D*>(file->Get(config.inputPtHistogramName.c_str()));
    if (!fHPt) {
        std::string availableHists;
        TIter next(file->GetListOfKeys());
        while (TObject* obj = next()) {
            availableHists += obj->GetName();
            availableHists += " ";
        }
        throw std::runtime_error("Failed to get input pT histogram: " + config.inputPtHistogramName +"\n"
                                 "Available histograms in file: " + availableHists);
    }
    fHPt = static_cast<TH1D*>(fHPt->Clone());
    fHPt->SetDirectory(nullptr);
    file->Close();

    fAverageNucleons = fHPt->Integral() * (fYMax * 2.0);
}

bool CoalescenceEngine::checkNucleusWithinRapidity(const std::vector<Particle>& nucleons, double yMin, double yMax) const
{
    TLorentzVector pTot;
    for (const auto& p : nucleons) {
        pTot += p.mom;
    }
    pTot.SetE(std::sqrt(fMass * fMass + pTot.Vect().Mag2()));
    double yNucleus = pTot.Rapidity();
    return (yNucleus >= yMin && yNucleus <= yMax);
}

void CoalescenceEngine::run(BookKeeping& bookKeeping) {
    
    const int nThreads = fConfig.nThreads;
    const int nEvents = fConfig.nEvents;
    const long long evPerThread = nEvents / nThreads;
    const long long remainder   = nEvents % nThreads;
 
    std::vector<std::thread>                workers;
    std::vector<std::pair<float, float>>    threadResults(nThreads);
 
    long long evStart = 0;
 
    // Pre-clone histograms on the main thread — one clone per worker.
    // This avoids concurrent TH1::Clone() calls which are not thread-safe
    // due to gROOT registration.
    std::vector<std::vector<double>> threadYields(nThreads);
    std::vector<BookKeeping> bookKeepingClones(nThreads);
    std::vector<std::unique_ptr<TH1D>> hPtNucleonClones(nThreads);
    for (int t = 0; t < nThreads; ++t) {
        bookKeepingClones[t] = BookKeeping{};
        hPtNucleonClones[t] = std::unique_ptr<TH1D>(static_cast<TH1D*>(fHPt->Clone(Form("hInputPtNucleon_%d", t))));
    }
 
    for (int t = 0; t < nThreads; ++t) {
        const long long evEnd = evStart + evPerThread + (t < remainder ? 1 : 0);
        workers.emplace_back(&CoalescenceEngine::workerRun, this, evStart, evEnd, t, hPtNucleonClones[t].get(), 
                             std::ref(threadResults[t]), std::ref(threadYields[t]), std::ref(bookKeepingClones[t]));
        evStart = evEnd;
    }
 
    for (auto& w : workers) w.join();
    for (int t = 1; t < nThreads; ++t) {
        bookKeepingClones[0].fHPtNucleon->Add(bookKeepingClones[t].fHPtNucleon);
        bookKeepingClones[0].fHYNucleon->Add(bookKeepingClones[t].fHYNucleon);
        bookKeepingClones[0].fHPhiNucleon->Add(bookKeepingClones[t].fHPhiNucleon);
        bookKeepingClones[0].fHPtNucleus->Add(bookKeepingClones[t].fHPtNucleus);
        bookKeepingClones[0].fHYNucleus->Add(bookKeepingClones[t].fHYNucleus);
        bookKeepingClones[0].fHPhiNucleus->Add(bookKeepingClones[t].fHPhiNucleus);
    }
    fBookKeeping = bookKeepingClones[0];
 
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
    for (int t = 0; t < nThreads; ++t) {
        for (const auto& w : threadYields[t]) {
            const double diff = w - yield;
            yieldUncertainty += diff * diff;
            fBookKeeping.fHYieldDistributionNucleus->Fill(w);
        }
    }
    yieldUncertainty = std::sqrt(yieldUncertainty / (totalEvents - 1));

    fBookKeeping.fHYieldNucleus->SetBinContent(1, yield);
    fBookKeeping.fHYieldUncertaintyNucleus->SetBinContent(1, yieldUncertainty / std::sqrt(totalEvents));
    bookKeeping = fBookKeeping;
}


void CoalescenceEngine::generateEvent(Event& event, PositionSampler& positionSampler, MomentumSampler& momentumSampler, TRandom3& random) const {
    const int nProtons = static_cast<int>(random.Poisson(fAverageNucleons));
    const int nNeutrons = static_cast<int>(random.Poisson(fAverageNucleons));

    event.protons = momentumSampler.sampleN(nProtons);
    event.neutrons = momentumSampler.sampleN(nNeutrons);

    positionSampler.sampleN(event.protons);
    positionSampler.sampleN(event.neutrons); 
}
