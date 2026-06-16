#include "ToyMCEngine.h"
#include "SourceSize.h"
#include "WignerDensity.h"
#include "CoalescenceEngine.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStopwatch.h"

#include <iostream>
#include <memory>
#include <cmath>
#include <string>

// ─────────────────────────────────────────────────────────────────────────────
// Usage
// ─────────────────────────────────────────────────────────────────────────────

void printUsage(const char* prog) {
    std::cerr << "Usage: " << prog
              << " <spectra.root> <output.root> <nEvents> <centrality>"
                 " [hNameProton] [hNameHe3] [seed]\n"
              << "\n"
              << "  spectra.root  : ROOT file with input pT histograms\n"
              << "  output.root   : output ROOT file\n"
              << "  nEvents       : number of toy events to generate\n"
              << "  centrality    : 0 = 0-10%, 1 = 10-30%, 2 = 30-50%\n"
              << "  hNameProton   : histogram name for proton pT (default: hProtonPt)\n"
              << "  hNameHe3      : histogram name for He3 pT    (default: hHe3Pt)\n"
              << "  seed          : RNG seed (default: 42)\n";
}

// ─────────────────────────────────────────────────────────────────────────────
// Source sizes  (same helper as main.cxx — TODO: move to a shared Config.h)
// ─────────────────────────────────────────────────────────────────────────────

std::pair<SourceSize, SourceSize> makeSourceSizes(int centrality) {
    const double mTproton_GeV = 1.74; // GeV/c^2
    const double mTHe3_GeV    = 4.97; // GeV/c^2
    
    return {
        SourceSize(static_cast<SourceSizeParameters::CentralityClass>(centrality), mTproton_GeV),
        SourceSize(static_cast<SourceSizeParameters::CentralityClass>(centrality), mTHe3_GeV)
    };
}

// ─────────────────────────────────────────────────────────────────────────────
// main
// ─────────────────────────────────────────────────────────────────────────────

int main(int argc, char** argv) {

    TStopwatch timer;
    constexpr double kLi4Radius_fm = 2.0; // fm, from ALICE measurement of Li4 rms radius

    if (argc < 5) { printUsage(argv[0]); return 1; }

    const std::string spectraFile  = argv[1];
    const std::string outputFile   = argv[2];
    const long long   nEvents      = std::stoll(argv[3]);
    const int         centrality   = std::stoi(argv[4]);
    const std::string hNameProton  = (argc >= 6) ? argv[5] : "hProtonPt";
    const std::string hNameHe3     = (argc >= 7) ? argv[6] : "hHe3Pt";
    const unsigned int seed        = (argc >= 8) ? std::stoul(argv[7]) : 42u;

    std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
              << " Li4 coalescence — ToyMC mode\n"
              << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
              << " Spectra file : " << spectraFile << "\n"
              << " Output       : " << outputFile  << "\n"
              << " nEvents      : " << nEvents     << "\n"
              << " Centrality   : " << centrality  << "\n"
              << " h(proton)    : " << hNameProton << "\n"
              << " h(He3)       : " << hNameHe3    << "\n"
              << " Seed         : " << seed        << "\n"
              << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";

    // ── Source sizes and Wigner density ──────────────────────────────────────
    auto [srcProton, srcHe3] = makeSourceSizes(centrality);

    constexpr double rRms_fm = kLi4Radius_fm;
    const double     d_fm    = rRms_fm * std::sqrt(2./3.);
    auto wigner = std::make_shared<GaussianWigner>(d_fm);

    std::cout << " Li4 Gaussian width: d = " << d_fm << " fm\n"
              << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";

    // ── Build ToyMC engine ────────────────────────────────────────────────────
    ToyMCEngine::Config cfg;
    cfg.rootFile     = spectraFile;
    cfg.hNameProton  = hNameProton;
    cfg.hNameHe3     = hNameHe3;
    cfg.meanNProtons = -1.; // Will be set by ToyMCEngine constructor
    cfg.meanNHe3     = -1.; // Will be set by ToyMCEngine constructor
    cfg.yMax         = 0.5;
    cfg.multMode     = MultiplicityMode::Poissonian;
    cfg.seed         = seed;
    cfg.nThreads     = 10;

    ToyMCEngine engine(srcProton, srcHe3, wigner, cfg);

    // ── Run ───────────────────────────────────────────────────────────────────
    timer.Start();    
    ToyMCEngine::Histograms h = engine.run(nEvents);
    timer.Stop();
    std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
              << " Finished in " << timer.RealTime() << " seconds\n"
              << " CPU time: " << timer.CpuTime() << " seconds\n"
              << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";

    // ── Write output ──────────────────────────────────────────────────────────
    TFile outFile(outputFile.c_str(), "RECREATE");
    if (outFile.IsZombie()) {
        std::cerr << "ERROR: cannot create " << outputFile << "\n";
        return 1;
    }
    h.write();
    outFile.Close();

    std::cout << " Results written to: " << outputFile << "\n";
    return 0;
}