#ifndef BOOKKEEPING_H
#define BOOKKEEPING_H

#include "TH1D.h"
#include "TString.h"

struct BookKeeping {
    TH1D* fHPtNucleon = nullptr;
    TH1D* fHYNucleon = nullptr;
    TH1D* fHPhiNucleon = nullptr;

    TH1D* fHPtNucleus = nullptr;
    TH1D* fHYNucleus = nullptr;
    TH1D* fHPhiNucleus = nullptr;

    TH1D* fHYieldNucleus = nullptr;
    TH1D* fHYieldUncertaintyNucleus = nullptr;
    TH1D* fHYieldDistributionNucleus = nullptr;

    void reset() {
        for (auto hist : {fHPtNucleon, fHYNucleon, fHPhiNucleon,
                      fHPtNucleus, fHYNucleus, fHPhiNucleus,
                      fHYieldNucleus, fHYieldUncertaintyNucleus, fHYieldDistributionNucleus}) {
            if (hist) {
                delete hist;
                hist = nullptr;
            }
        }
    }

    void init(const int index) {
        fHPtNucleon = new TH1D(Form("hPtNucleon_%d", index), "Sampled nucleon;p_{T} (GeV/c);Counts", 100, 0., 10.);
        fHYNucleon = new TH1D(Form("hYNucleon_%d", index), "Sampled nucleon;y (GeV/c);Counts", 50, -2., 2.);
        fHPhiNucleon = new TH1D(Form("hPhiNucleon_%d", index),"Sampled nucleon;#phi (rad);Counts", 50, 0., 2. * M_PI);
        
        fHPtNucleus = new TH1D(Form("hPtNucleus_%d", index), "Sampled nucleus;p_{T} (GeV/c);Counts", 100, 0., 10.);
        fHYNucleus = new TH1D(Form("hYNucleus_%d", index), "Sampled nucleus;y (GeV/c);Counts", 150, -1.5, 1.5);
        fHPhiNucleus = new TH1D(Form("hPhiNucleus_%d", index), "Sampled nucleus;#phi (rad);Counts", 50, 0., 2. * M_PI);

        fHYieldNucleus = new TH1D(Form("hYieldNucleus_%d", index), "Sampled nucleus yield;Yield;Counts", 1, 0., 1.);
        fHYieldUncertaintyNucleus = new TH1D(Form("hYieldUncertaintyNucleus_%d", index), "Sampled nucleus ;Yield uncertainty;Counts", 1, 0., 1.);
        fHYieldDistributionNucleus = new TH1D(Form("hYieldDistributionNucleus_%d", index), "Sampled nucleus ;Yield distribution;Counts", 1000, 0., 0.001);
    }

    void write(TDirectory* out) {
        if (!out)
            throw std::invalid_argument("BookKeeping::write: output directory is null");

        out->cd();
        for (auto hist : {fHPtNucleon, fHYNucleon, fHPhiNucleon,
                      fHPtNucleus, fHYNucleus, fHPhiNucleus,
                      fHYieldNucleus, fHYieldUncertaintyNucleus, fHYieldDistributionNucleus}) {
            if (hist) {
                hist->Write();
            }
        }
    }   

};

#endif // BOOKKEEPING_H