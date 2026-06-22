#ifndef TOYMCENGINEA4_H
#define TOYMCENGINEA4_H

#include "IntegratorA4.h"
#include "SourceSize.h"
#include "WignerDensityA.h"

#include "TDirectory.h"
#include "TFile.h"
#include "TH1.h"

#include <memory>
#include <string>
#include <stdexcept>

// ─────────────────────────────────────────────────────────────────────────────
// ToyMCEngineA4
//
// Thin driver that wires together SourceSize, WignerDensityA, and
// IntegratorA4. Responsible for:
//   - Loading the input nucleon pT histogram from a ROOT file
//   - Constructing and running IntegratorA4
//   - Writing output to a ROOT file
// ─────────────────────────────────────────────────────────────────────────────

class ToyMCEngineA4 {
public:

    struct Config {
        std::string rootFile;       // input ROOT file with nucleon pT histogram
        std::string hNameNucleon;   // histogram name

        double rRms_fm   = 2.0;     // Li4 rms radius [fm]
        double yMax      = 0.5;     // rapidity window half-width

        IntegratorA4::Config intCfg; // passed through to IntegratorA4
    };

    ToyMCEngineA4(const SourceSize& srcSize,
                  const Config&     cfg);

    // Run integration and write results to outputFile
    void run(TDirectory* out);

private:

    static std::unique_ptr<TH1> loadHistogram(const std::string& filename,
                                               const std::string& hname);

    SourceSize                     fSrcSize;
    std::unique_ptr<TH1>           fHPtNucleon;
    std::shared_ptr<WignerDensityA> fWigner;
    Config                         fCfg;
};

#endif // TOYMCENGINEA4_H