#include "ToyMCEngineA4.h"
#include "SourceSize.h"

#include <iostream>
#include <string>

void printUsage(const char* prog) {
    std::cerr << "Usage: " << prog
              << " <spectra.root> <output.root> <centrality>"
                 " [hName] [nSamples] [seed]\n"
              << "  centrality: 0=0-10%, 1=10-30%, 2=30-50%\n";
}

std::pair<SourceSize, SourceSize> makeSourceSizes(int centrality);  // defined in main_toymc.cxx — TODO: move to shared header

int main(int argc, char** argv) {
    if (argc < 4) { printUsage(argv[0]); return 1; }

    const std::string spectraFile = argv[1];
    const std::string outputFile  = argv[2];
    const int         centrality  = std::stoi(argv[3]);
    const std::string hName       = (argc >= 5) ? argv[4] : "hNucleonPt";
    const long long   nSamples    = (argc >= 6) ? std::stoll(argv[5]) : 10000000LL;
    const double      radius        = (argc >= 7) ? std::stod(argv[6]) : 3.;
    const unsigned int seed       = (argc >= 8) ? std::stoul(argv[7]) : 42u;

    std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
              << " Li4 A=4 integrator\n"
              << "━━━━━━━━━━━━━━━━━━━━━━════════════════════════════\n"
              << " Input      : " << spectraFile << "\n"
              << " Output     : " << outputFile  << "\n"
              << " Centrality : " << centrality  << "\n"
              << " Histogram  : " << hName       << "\n"
              << " N samples  : " << nSamples    << "\n"
              << " Li4 radius : " << radius      << " fm\n"
              << " Seed       : " << seed        << "\n"
              << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";

    const double mTproton_GeV = 1.74; // GeV/c^2
    const SourceSizeParameters::CentralityClass centClass = SourceSizeParameters::CentralityClass::PbPb502_0_10; 
    const int seedSourceSize = seed - 1; // different seed for source size sampling
    
    auto output = TFile::Open(outputFile.c_str(), "RECREATE");
    if (!output || output->IsZombie()) {
        std::cerr << "Error: cannot create output file " << outputFile << "\n";
        return 1;
    }
    const double radius_PbPb_Run3_0_10_fm = 4.6; // fm, for PbPb 0-10% centrality, Run 3
            
    SourceSize srcSize(centClass, mTproton_GeV, seedSourceSize);
    srcSize.setR0(radius_PbPb_Run3_0_10_fm);

    auto outputDir = output->mkdir(Form("d%.1f", radius));
    std::cout << " Running with rRms = " << radius << " fm\n";
    ToyMCEngineA4::Config cfg;
    cfg.rootFile       = spectraFile;
    cfg.hNameNucleon   = hName;
    cfg.rRms_fm        = radius;
    cfg.yMax           = 0.5;
    cfg.intCfg.nSamples = nSamples;
    cfg.intCfg.seed     = seed;
    cfg.intCfg.yMax     = cfg.yMax;

    ToyMCEngineA4 engine(srcSize, cfg);
    engine.run(outputDir);

    output->Close();

    return 0;
}