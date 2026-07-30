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
    
    /* 
    // He3 configuration
    Config config = {
        .nucleusName = "He3",
        .sourceRadius = 4.6,   // fm
        .inputPtHistogramFile = "../input/spectra_0_10.root",
        .inputPtHistogramName = "hProton_0_10",
        .outputFile = "../output/output_he3.root",
        .nEvents = 100,
        .nThreads = 20,
        .randomSeed = 42
    };
    */

    // He4 configuration
    Config config = {
        .nucleusName = "He4",
        .sourceRadius = 4.6,   // fm
        .inputPtHistogramFile = "../input/spectra_0_10.root",
        .inputPtHistogramName = "hProton_0_10",
        .outputFile = "../output/output_he4.root",
        .nEvents = 100,
        .nThreads = 20,
        .randomSeed = 42
    };

    std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
              << " A integrator\n"
              << "━━━━━━━━━━━━━━━━━━━━━━════════════════════════════\n"
              << " Input spectrum : " << config.inputPtHistogramFile << "\n"
              << " Histogram      : " << config.inputPtHistogramName << "\n"
              << " Output         : " << config.outputFile << "\n"
              << " N samples      : " << config.nEvents << "\n"
              << " Seed           : " << config.randomSeed << "\n"
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