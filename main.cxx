#include "Config.h"
#include "PositionSampler.h"
#include "ToyMcEngine.h"

#include "TFile.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include <iostream>
#include <string>

void printUsage(const char* prog) {
    std::cerr << "Usage: " << prog;
}

int main(int argc, char** argv) {
    
    //if (argc < 4) { printUsage(argv[0]); return 1; }

    ROOT::EnableThreadSafety();

    Config config = {
        .nucleusName = "He3",
        .nucleusRadius = 1.96 * std::sqrt(2.), // fm
        .sourceRadius = 4.6,   // fm
        .inputPtHistogramFile = "../input/spectra_0_10.root",
        .inputPtHistogramName = "hProton_0_10",
        .outputFile = "../output/output.root",
        .nEvents = 100,
        .nThreads = 20,
        .randomSeed = 42
    };

    std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
              << " Li4 A=4 integrator\n"
              << "━━━━━━━━━━━━━━━━━━━━━━════════════════════════════\n"
              << " Input      : " << config.inputPtHistogramFile << "\n"
              << " Output     : " << config.outputFile << "\n"
              << " Centrality : " << config.inputPtHistogramFile << "\n"
              << " Histogram  : " << config.inputPtHistogramName << "\n"
              << " N samples  : " << config.nEvents << "\n"
              << " Li4 radius : " << config.nucleusRadius << " fm\n"
              << " Seed       : " << config.randomSeed << "\n"
              << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";

    auto output = TFile::Open(config.outputFile.c_str(), "RECREATE");
    if (!output || output->IsZombie()) {
        std::cerr << "Error: cannot create output file " << config.outputFile << "\n";
        return 1;
    }

    TStopwatch stopwatch;

    ToyMcEngine engine(output, config);
    
    stopwatch.Start();
    engine.run(output);
    stopwatch.Stop();
    
    std::cout << "Elapsed time: " << stopwatch.RealTime() << " s\n";
    std::cout << "CPU time: " << stopwatch.CpuTime() << " s\n";

    output->Close();

    return 0;
}