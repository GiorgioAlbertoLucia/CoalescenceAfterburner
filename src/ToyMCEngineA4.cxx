#include "ToyMCEngineA4.h"

#include <iostream>
#include <cmath>

ToyMCEngineA4::ToyMCEngineA4(const SourceSize& srcSize,
                               const Config&     cfg)
    : fSrcSize(srcSize)
    , fCfg(cfg)
{
    fHPtNucleon = loadHistogram(cfg.rootFile, cfg.hNameNucleon);

    //const double d_fm = cfg.rRms_fm * std::sqrt(2./3.);
    const double d_fm = cfg.rRms_fm;
    fWigner = std::make_shared<WignerDensityA>(d_fm);

    std::cout << " ToyMCEngineA4: d = " << d_fm
              << " fm  (r_rms = " << cfg.rRms_fm << " fm)\n";
}

void ToyMCEngineA4::run(TDirectory* out) {
    IntegratorA4 integrator(fSrcSize, fWigner, fHPtNucleon.get(), fCfg.intCfg);
    integrator.run();

    out->cd();
    integrator.write();
}

std::unique_ptr<TH1> ToyMCEngineA4::loadHistogram(const std::string& filename,
                                                    const std::string& hname)
{
    TFile f(filename.c_str(), "READ");
    if (f.IsZombie())
        throw std::runtime_error("ToyMCEngineA4: cannot open " + filename);

    TH1* h = dynamic_cast<TH1*>(f.Get(hname.c_str()));
    if (!h)
        throw std::runtime_error("ToyMCEngineA4: histogram " + hname
                                 + " not found in " + filename);

    auto owned = std::unique_ptr<TH1>(static_cast<TH1*>(h->Clone()));
    owned->SetDirectory(nullptr);
    return owned;
}