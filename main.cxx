#include "EventReader.h"
#include "CoalescenceEngine.h"
#include "SourceSize.h"
#include "WignerDensity.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

#include <iostream>
#include <string>
#include <memory>
#include <stdexcept>
#include <cmath>

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

void printUsage(const char* progName) {
    std::cerr << "Usage: " << progName
              << " <input.dat> <output.root> <centrality> [seed]\n"
              << "\n"
              << "  input.dat    : FIST-sampler ASCII output\n"
              << "  output.root  : ROOT file for histograms\n"
              << "  centrality   : 0 = 0-10%, 1 = 10-30%, 2 = 30-50%\n"
              << "  seed         : RNG seed (optional, default 42)\n";
}

std::pair<SourceSize, SourceSize> makeSourceSizes(int centrality) {
    // Proton and He3 share the same mT-scaling within current uncertainties
    // Individual values for the mT of the species is used

    const double mTproton_GeV = 1.74; // GeV/c^2
    const double mTHe3_GeV    = 4.97; // GeV/c^2
    
    return {
        SourceSize(static_cast<SourceSizeParameters::CentralityClass>(centrality), mTproton_GeV),
        SourceSize(static_cast<SourceSizeParameters::CentralityClass>(centrality), mTHe3_GeV)
    };
}

// ─────────────────────────────────────────────────────────────────────────────
// Histogram booking
// ─────────────────────────────────────────────────────────────────────────────

struct Histograms {
    // Li4 spectra
    TH1D* hPt      = nullptr; // d^2N / (dy dpT)
    TH1D* hRapidity= nullptr;

    // Coalescence parameter BA (generalised to He3+p -> Li4)
    // BA = (1/2pi) * (d^2N_Li4/dy dpT) / [(1/2pi) * (d^2N_He3/dy dpT)]
    //                                   / [(1/2pi) * (d^2N_p  /dy dpT)]
    // stored differentially vs pT/A
    TH1D* hBA      = nullptr;

    // Input spectra (for BA computation)
    TH1D* hPtProton= nullptr;
    TH1D* hPtHe3   = nullptr;

    // Diagnostics
    TH1D* hNLi4PerEvent   = nullptr;
    TH1D* hNProtonPerEvent= nullptr;
    TH1D* hNHe3PerEvent   = nullptr;
    TH2D* hQvsR           = nullptr; // relative momentum vs distance in PRF
};

Histograms bookHistograms() {
    Histograms h;

    // pT binning for Li4: 0 to 10 GeV/c in 50 bins
    h.hPt       = new TH1D("hLi4Pt",
                            "Li4 p_{T} spectrum;p_{T} (GeV/c);d^{2}N/(dydp_{T})",
                            50, 0., 10.);
    h.hRapidity = new TH1D("hLi4Rapidity",
                            "Li4 rapidity;y;dN/dy",
                            60, -3., 3.);

    // Input spectra — same binning, for BA
    h.hPtProton = new TH1D("hProtonPt",
                            "Proton p_{T};p_{T} (GeV/c);d^{2}N/(dydp_{T})",
                            50, 0., 5.);
    h.hPtHe3    = new TH1D("hHe3Pt",
                            "^{3}He p_{T};p_{T} (GeV/c);d^{2}N/(dydp_{T})",
                            50, 0., 8.);

    // BA vs pT/A  (A=4 for Li4)
    h.hBA       = new TH1D("hBA",
                            "Coalescence parameter B_{A};p_{T}/A (GeV/c);B_{A}",
                            25, 0., 2.5);

    // Diagnostics
    h.hNLi4PerEvent    = new TH1D("hNLi4PerEvent",
                                   "Li4 per event;N_{Li4};Counts",
                                   10, 0., 10.);
    h.hNProtonPerEvent = new TH1D("hNProtonPerEvent",
                                   "Protons per event;N_{p};Counts",
                                   50, 0., 50.);
    h.hNHe3PerEvent    = new TH1D("hNHe3PerEvent",
                                   "He3 per event;N_{He3};Counts",
                                   20, 0., 20.);
    h.hQvsR            = new TH2D("hQvsR",
                                   "PRF relative momentum vs distance;"
                                   "q (GeV/c);r (fm)",
                                   100, 0., 2., 100, 0., 20.);

    return h;
}

// ─────────────────────────────────────────────────────────────────────────────
// main
// ─────────────────────────────────────────────────────────────────────────────

int main(int argc, char** argv) {

    constexpr double kLi4Radius_fm = 2.0; // fm, from ALICE measurement of Li4 rms radius

    // ── Parse arguments ──────────────────────────────────────────────────────
    if (argc < 4) { printUsage(argv[0]); return 1; }

    const std::string inputFile  = argv[1];
    const std::string outputFile = argv[2];
    const int         centrality = std::stoi(argv[3]);
    const unsigned int seed      = (argc >= 5) ? std::stoul(argv[4]) : 42u;

    std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
              << " Li4 coalescence afterburner\n"
              << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
              << " Input  : " << inputFile  << "\n"
              << " Output : " << outputFile << "\n"
              << " Centrality class : " << centrality << "\n"
              << " RNG seed         : " << seed       << "\n"
              << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";

    // ── Set up source sizes ──────────────────────────────────────────────────
    auto [srcProton, srcHe3] = makeSourceSizes(centrality);

    std::cout << " Source size (proton): R = " << srcProton.r0() << " fm\n"
              << " Source size (He3)   : R = " << srcHe3.r0()    << " fm\n"
              << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";

    // ── Set up Wigner density ────────────────────────────────────────────────
    // d = r_rms * sqrt(2/3), r_rms = 2 fm for Li4
    constexpr double rRms_fm = kLi4Radius_fm;
    const double     d_fm    = rRms_fm * std::sqrt(2./3.);

    std::cout << " Li4 Gaussian width: d = " << d_fm << " fm"
              << "  (r_rms = " << rRms_fm << " fm)\n"
              << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";

    auto wigner = std::make_shared<GaussianWigner>(d_fm);

    // ── Build engine ─────────────────────────────────────────────────────────
    CoalescenceEngine engine(srcProton, srcHe3, wigner, seed);

    // ── Open input file ──────────────────────────────────────────────────────
    EventReader reader(inputFile);
    if (!reader.isOpen()) {
        std::cerr << "ERROR: cannot open input file: " << inputFile << "\n";
        return 1;
    }

    // ── Book histograms ──────────────────────────────────────────────────────
    Histograms h = bookHistograms();

    // ── Event loop ───────────────────────────────────────────────────────────
    Event ev;
    long long nEvents = 0;
    long long nLi4Total = 0;

    constexpr int kPrintEvery = 10000;

    while (reader.nextEvent(ev)) {
        ++nEvents;
        if (nEvents % kPrintEvery == 0)
            std::cout << " Processed " << nEvents << " events ...\n";

        // Fill input spectra
        for (const auto& p : ev.protons)
            h.hPtProton->Fill(p.pT());
        for (const auto& he3 : ev.he3s)
            h.hPtHe3->Fill(he3.pT());

        // Diagnostics
        h.hNProtonPerEvent->Fill(ev.nProtons());
        h.hNHe3PerEvent   ->Fill(ev.nHe3());

        // Run coalescence
        CoalescenceResult result = engine.processEvent(ev);

        // Fill Li4 histograms
        for (const auto& li4 : result.li4s) {
            h.hPt      ->Fill(li4.pT());
            h.hRapidity->Fill(li4.y());
            ++nLi4Total;
        }

        // Fill q vs r diagnostic for all evaluated pairs
        for (const auto& pd : result.pairs)
            h.hQvsR->Fill(pd.q_GeV, pd.r_fm);

        h.hNLi4PerEvent->Fill(result.nLi4());
    }

    std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
              << " Done.\n"
              << " Total events processed : " << nEvents   << "\n"
              << " Total Li4 produced     : " << nLi4Total << "\n"
              << " Mean Li4 per event     : "
              << (nEvents > 0 ? double(nLi4Total)/nEvents : 0.) << "\n"
              << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";

    // ── Normalise spectra by number of events and bin width ──────────────────
    // d^2N / (dy dpT): divide by nEvents and pT bin width.
    // Rapidity window is filled as-is (not normalised to delta_y here).
    auto normalise = [&](TH1D* hist) {
        if (nEvents == 0) return;
        hist->Sumw2();
        hist->Scale(1. / nEvents);
        // Divide each bin by its width
        for (int ib = 1; ib <= hist->GetNbinsX(); ++ib) {
            const double width = hist->GetBinWidth(ib);
            hist->SetBinContent(ib, hist->GetBinContent(ib) / width);
            hist->SetBinError  (ib, hist->GetBinError(ib)   / width);
        }
    };

    normalise(h.hPt);
    normalise(h.hPtProton);
    normalise(h.hPtHe3);

    // ── Compute BA ───────────────────────────────────────────────────────────
    // BA(pT/A) = [d^2N_Li4 / (dy dpT)] evaluated at pT
    //          / [d^2N_He3 / (dy dpT)] evaluated at pT/A * (3/4)
    //          / [d^2N_p   / (dy dpT)] evaluated at pT/A * (1/4)
    // where pT/A is the per-nucleon transverse momentum.
    // Here we fill hBA bin-by-bin using the per-nucleon pT axis.
    // For a rough first estimate we use the global yield ratio:
    //   BA ~ (dN_Li4/dpT) / [(dN_He3/dpT)(dN_p/dpT)] * (2*pi*pT)^(A-1)
    // A proper differential BA requires matched pT bins — left as a post-
    // processing step in ROOT macros where finer control is easier.
    // The raw histograms are sufficient for that computation.

    // ── Write output ─────────────────────────────────────────────────────────
    TFile outFile(outputFile.c_str(), "RECREATE");
    if (outFile.IsZombie()) {
        std::cerr << "ERROR: cannot create output file: " << outputFile << "\n";
        return 1;
    }

    h.hPt            ->Write();
    h.hRapidity      ->Write();
    h.hPtProton      ->Write();
    h.hPtHe3         ->Write();
    h.hBA            ->Write();
    h.hNLi4PerEvent  ->Write();
    h.hNProtonPerEvent->Write();
    h.hNHe3PerEvent  ->Write();
    h.hQvsR          ->Write();

    outFile.Close();

    std::cout << " Results written to: " << outputFile << "\n";
    return 0;
}
