#include "ToyMCEngineA3.h"

#include <iostream>
#include <cmath>

ToyMCEngineA3::ToyMCEngineA3(const SourceSize& srcSize,
                               const Config&     cfg)
    : fSrcSize(srcSize)
    , fCfg(cfg)
{
    fHPtNucleon = loadHistogram(cfg.rootFile, cfg.hNameNucleon);

    //const double d_fm = cfg.rRms_fm * std::sqrt(2./3.);
    const double d_fm = cfg.rRms_fm;
    fWigner = std::make_shared<WignerDensityA>(d_fm);

    std::cout << " ToyMCEngineA3: d = " << d_fm
              << " fm  (r_rms = " << cfg.rRms_fm << " fm)\n";
}

void ToyMCEngineA3::run(TDirectory* out) {
    IntegratorA3 integrator(fSrcSize, fWigner, fHPtNucleon.get(), fCfg.intCfg);
    integrator.run(fCfg.intCfg.nEvents);

    out->cd();
    integrator.write();
}

std::unique_ptr<TH1> ToyMCEngineA3::loadHistogram(const std::string& filename,
                                                    const std::string& hname)
{
    TFile f(filename.c_str(), "READ");
    if (f.IsZombie())
        throw std::runtime_error("ToyMCEngineA3: cannot open " + filename);

    TH1* h = dynamic_cast<TH1*>(f.Get(hname.c_str()));
    if (!h)
        throw std::runtime_error("ToyMCEngineA3: histogram " + hname
                                 + " not found in " + filename);

    auto owned = std::unique_ptr<TH1>(static_cast<TH1*>(h->Clone()));
    owned->SetDirectory(nullptr);
    return owned;
}