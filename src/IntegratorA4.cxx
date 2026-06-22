#include "IntegratorA4.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

#include "TH1F.h"

// ─────────────────────────────────────────────────────────────────────────────
// Constructor
// ─────────────────────────────────────────────────────────────────────────────

IntegratorA4::IntegratorA4(const SourceSize&               srcSize,
                            std::shared_ptr<WignerDensityA> wigner,
                            const TH1*                      hPtNucleon,
                            const Config&                   cfg)
    : fCfg(cfg)
    , fSrcSize(srcSize)
    , fWigner(std::move(wigner))
    , fRng(cfg.seed)
{
    if (!fWigner)
        throw std::invalid_argument("IntegratorA4: WignerDensityA must not be null");
    if (!hPtNucleon)
        throw std::invalid_argument("IntegratorA4: hPtNucleon must not be null");

    constexpr double kMassNucleon = 0.938272; // GeV/c^2

    fMomSampler = std::make_unique<MomentumSampler>(
        PDG::kProton, kMassNucleon, hPtNucleon, cfg.yMax, cfg.seed + 1);

    // ── Book histograms ───────────────────────────────────────────────────────
    fHLi4Pt = std::make_unique<TH1D>(
        "hLi4Pt_A4",
        "^{4}Li p_{T};p_{T} (GeV/c);d^{2}N/(dydp_{T})",
        50, 0., 12.);
    fHLi4Pt->SetDirectory(nullptr);

    fHNucleonPt = std::make_unique<TH1D>(
        "hNucleonPt_A4",
        "Sampled nucleon p_{T};p_{T} (GeV/c);Counts",
        50, 0., 5.);
    fHNucleonPt->SetDirectory(nullptr);

    fHSourceDist = std::make_unique<TH1D>(
        "hSourceDist_A4",
        "Inter-nucleon distance;r (fm);Counts",
        100, 0., 20.);
    fHSourceDist->SetDirectory(nullptr);

    fHJacobiKvsR = std::make_unique<TH2D>(
        "hJacobiKvsR_A4",
        "Jacobi |k_{j}| vs |r_{j}|;|k| (GeV/c);|r| (fm)",
        100, 0., 2., 100, 0., 20.);
    fHJacobiKvsR->SetDirectory(nullptr);

    fHWignerKvsR = std::make_unique<TH2D>(
        "hWignerKvsR_A4",
        "Wigner |k_{j}| vs |r_{j}|;|k| (GeV/c);|r| (fm)",
        100, 0., 2., 100, 0., 20.);
    fHWignerKvsR->SetDirectory(nullptr);
    for (int xbin = 1; xbin <= fHWignerKvsR->GetNbinsX(); ++xbin) {
        for (int ybin = 1; ybin <= fHWignerKvsR->GetNbinsY(); ++ybin) {
            const double k = fHWignerKvsR->GetXaxis()->GetBinCenter(xbin);
            const double r = fHWignerKvsR->GetYaxis()->GetBinCenter(ybin);
            const double w = fWigner->evaluate({TVector3(k, 0, 0)}, {TVector3(r, 0, 0)});
            fHWignerKvsR->SetBinContent(xbin, ybin, w);
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// run
// ─────────────────────────────────────────────────────────────────────────────

void IntegratorA4::run() {

    const long long nSamples   = fCfg.nSamples;
    const long long printEvery = std::max(1LL, nSamples / 100);
    const int       A          = fCfg.A;

    // SA for Li4 (all excited states)
    constexpr double SA = 3./16.;

    std::vector<TVector3> positions(A);
    std::vector<TVector3> momenta(A);

    auto hLi4PtNormalisation = (TH1F*)fHLi4Pt->Clone("hLi4PtNormalisation");
    hLi4PtNormalisation->Reset();

    for (long long iS = 0; iS < nSamples; ++iS) {
        if (iS % printEvery == 0)
            std::cout << " [IntegratorA4] Sample " << iS << " / " << nSamples
                      << "  (" << (100 * iS / nSamples) << "%)\n";

        // ── 1. Draw positions and momenta ─────────────────────────────────────
        drawSample(positions, momenta);
        
        // ── 2. Transform to Jacobi coordinates ───────────────────────────────
        auto r_jacobi = JacobiTransform::relative(positions); // A-1 vectors [fm]
        auto k_jacobi = JacobiTransform::relative(momenta);   // A-1 vectors [GeV/c]


        
        std::vector<TVector3> r_jacobiForIntegration, k_jacobiForIntegration;
        r_jacobiForIntegration.reserve(A-1);
        k_jacobiForIntegration.reserve(A-1);
        for (int j = 1; j < A; ++j) {
            r_jacobiForIntegration.push_back(r_jacobi[j]);
            k_jacobiForIntegration.push_back(k_jacobi[j]);
        }

        // ── 3. Evaluate A-body Wigner density ────────────────────────────────
        const double D = fWigner->evaluate(k_jacobiForIntegration, r_jacobiForIntegration);

        // ── 4. Compute Li4 pT and fill weighted histogram ────────────────────
        const double pT = li4pT(momenta);
        fHLi4Pt->Fill(pT, SA * D);
        hLi4PtNormalisation->Fill(pT, 1.);

        // ── 5. Fill diagnostic histograms (unweighted) ───────────────────────
        for (int i = 0; i < A; ++i)
            fHNucleonPt->Fill(std::sqrt(momenta[i].X()*momenta[i].X()
                                      + momenta[i].Y()*momenta[i].Y()));

        // Inter-nucleon distances (all pairs)
        for (int i = 0; i < A; ++i)
            for (int j = i + 1; j < A; ++j)
                fHSourceDist->Fill((positions[i] - positions[j]).Mag());

        // Jacobi k vs r for each relative coordinate
        for (int j = 0; j < A - 1; ++j)
            fHJacobiKvsR->Fill(k_jacobi[j].Mag(), r_jacobi[j].Mag());
    }

    // ── Normalise Li4 spectrum ────────────────────────────────────────────────
    for (int ib = 1; ib <= fHLi4Pt->GetNbinsX(); ++ib) {
        const double w = hLi4PtNormalisation->GetBinContent(ib);
        if (w > 0.)            fHLi4Pt->SetBinContent(ib, fHLi4Pt->GetBinContent(ib) / w);
        if (w > 0. && fHLi4Pt->GetBinError(ib) > 0.) fHLi4Pt->SetBinError  (ib, fHLi4Pt->GetBinError(ib)   / w);
        
        // Alternative normalisation using bin width:
        //const double w = fHLi4Pt->GetBinWidth(ib);
        //fHLi4Pt->SetBinContent(ib, fHLi4Pt->GetBinContent(ib) / w);
        //fHLi4Pt->SetBinError  (ib, fHLi4Pt->GetBinError(ib)   / w);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// drawSample
// ─────────────────────────────────────────────────────────────────────────────

void IntegratorA4::drawSample(std::vector<TVector3>& positions,
                               std::vector<TVector3>& momenta) const
{
    for (int i = 0; i < fCfg.A; ++i) {
        fMomSampler->sampleMomentum(momenta[i]);

        //const double mT = p.mT();
        //const double r0 = fSrcSize.r0();
        fSrcSize.samplePosition(positions[i]);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// li4pT
// ─────────────────────────────────────────────────────────────────────────────

double IntegratorA4::li4pT(const std::vector<TVector3>& momenta) {
    double px = 0., py = 0.;
    for (const auto& p : momenta) {
        px += p.X();
        py += p.Y();
    }
    return std::sqrt(px*px + py*py);
}

// ─────────────────────────────────────────────────────────────────────────────
// write
// ─────────────────────────────────────────────────────────────────────────────

void IntegratorA4::write() const {
    fHLi4Pt     ->Write();
    fHNucleonPt ->Write();
    fHSourceDist->Write();
    fHJacobiKvsR->Write();
    fHWignerKvsR->Write();
}