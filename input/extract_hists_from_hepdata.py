import ROOT
import sys
import array

ROOT.gROOT.SetBatch(True)


# ── file paths ────────────────────────────────────────────────────────────────
HE3_FILE      = "HEPData-He3-spectra-PbPb-5p02-GeV.root"
PROTON_FILE   = "HEPData-p-spectra-PbPb-5p02-GeV.root"
THERMAL_MACRO = "../include/ThermalModels.h"    # your AliPWGFunc-style blast-wave fitter

MASS_PROTON = 0.9382720813  # GeV/c^2

ROOT.gROOT.LoadMacro(THERMAL_MACRO + "+")
gThermalModels = ROOT.ThermalModels()
gThermalModels.SetVarType(ROOT.ThermalModels.kdNdpt)  # matches "(1/Nev) d^2N/dpTdy" (NOT 1/pT-weighted)


def extend_with_blastwave_fit(hist, mass, pt_min=0.0, low_bin_width=0.05,
                               beta0=0.65, T0=0.10, n0=0.9, fix_shape=False,
                               tf1_name="fBW", outfile=None):
    """
    Fit `hist` (a native (1/Nev) d2N/dpTdy density -- NOT yet multiplied by
    bin width) with a Boltzmann-Gibbs Blast-Wave (ThermalModels::GetBGBW)
    over its measured range, then return a NEW histogram extended down to
    pt_min: below the lowest measured point the BGBW fit is used (integrated
    per bin), at/above it the original measured bins are copied verbatim
    (content AND error -- so hist's bin errors should already hold whatever
    combined stat+syst uncertainty you want propagated).

    Physics caveat: a single-species (proton-only) BGBW fit is under-
    constrained -- T, beta and n are strongly correlated, and it's normally
    the pi/K/p MASS ordering in a combined fit that breaks that degeneracy.
    If you already have T, beta, n from a published combined fit at this
    centrality/energy, call with fix_shape=True and those values as
    beta0/T0/n0: only the overall normalisation is then fit to the proton
    data, which is far more robust than letting all three shape parameters
    float against one species alone.

    Returns (extended_hist, fBW, fitResult).
    """
    pt_data_min = hist.GetXaxis().GetXmin()
    pt_data_max = hist.GetXaxis().GetXmax()

    if pt_data_min <= pt_min:
        print(f"extend_with_blastwave_fit: {hist.GetName()} already starts at/below "
              f"pt_min={pt_min}, nothing to extrapolate.")
        return hist.Clone(hist.GetName() + "_extrap"), None, None

    fBW = gThermalModels.GetBGBW(mass, beta0, T0, n0, hist.GetMaximum(), tf1_name)
    fBW.SetRange(pt_data_min, 2.5)

    if fix_shape:
        fBW.FixParameter(1, beta0)  # beta  (transverse flow velocity, units of c)
        fBW.FixParameter(2, T0)     # T     (kinetic freeze-out temperature, GeV)
        fBW.FixParameter(3, n0)     # n     (velocity profile exponent)
    else:
        fBW.SetParLimits(1, 0.30, 0.90)
        fBW.SetParLimits(2, 0.05, 0.20)
        fBW.SetParLimits(3, 0.30, 2.50)

    fitResult = hist.Fit(fBW, "RSQM")  # R: use range, S: TFitResultPtr, Q: quiet, M: more robust minimization
    status = int(fitResult)
    print(f"BGBW fit to {hist.GetName()}: status={status}  "
          f"chi2/ndf = {fitResult.Chi2():.2f}/{fitResult.Ndf()} "
          f"= {fitResult.Chi2()/max(fitResult.Ndf(), 1):.2f}")
    print(f"  beta = {fBW.GetParameter(1):.3f} +- {fBW.GetParError(1):.3f}")
    print(f"  T    = {fBW.GetParameter(2):.3f} +- {fBW.GetParError(2):.3f} GeV")
    print(f"  n    = {fBW.GetParameter(3):.3f} +- {fBW.GetParError(3):.3f}")
    if status != 0:
        print(f"  WARNING: fit did not converge cleanly (status={status}) -- "
              f"check the QA plot before trusting the low-pT extrapolation.")

    # ── build the extended, variable-binning histogram ──────────────────────
    n_low_bins = max(1, int(round((pt_data_min - pt_min) / low_bin_width)))
    low_edges  = [pt_min + i * low_bin_width for i in range(n_low_bins)] + [pt_data_min]
    data_edges = [hist.GetXaxis().GetBinLowEdge(i) for i in range(1, hist.GetNbinsX() + 2)]
    all_edges  = array.array('d', sorted(set(low_edges + data_edges)))

    hExt = ROOT.TH1D(hist.GetName() + "_extrap", hist.GetTitle(),
                      len(all_edges) - 1, all_edges)

    # Diagnostic-only error for the extrapolated bins: MomentumSampler's
    # GetRandom() only ever reads bin CONTENT, never bin error, so this
    # choice has zero effect on the toy MC itself -- only on the QA plot.
    # (TF1::IntegralError needs the FULL n x n covariance matrix including
    # zeroed rows/columns for fixed parameters -- mass is always fixed via
    # FixParameter(0, mass) inside ThermalModels, so fitResult's covariance
    # matrix alone is the wrong shape to feed it directly; not worth the
    # extra bookkeeping for a number with no physics consequence here.)
    rel_err = (hist.GetBinError(1) / hist.GetBinContent(1)
               if hist.GetBinContent(1) > 0 else 0.2)

    for ib in range(1, hExt.GetNbinsX() + 1):
        lo, hi = hExt.GetXaxis().GetBinLowEdge(ib), hExt.GetXaxis().GetBinUpEdge(ib)
        width = hi - lo

        if hi <= pt_data_min + 1e-9:
            content = fBW.Integral(lo, hi) / width
            error   = content * rel_err
        else:
            origBin = hist.GetXaxis().FindBin(0.5 * (lo + hi))
            content = hist.GetBinContent(origBin)
            error   = hist.GetBinError(origBin)

        hExt.SetBinContent(ib, content)
        hExt.SetBinError(ib, error)

    hExt.GetXaxis().SetTitle(hist.GetXaxis().GetTitle())
    hExt.GetYaxis().SetTitle(hist.GetYaxis().GetTitle())
    hExt.SetDirectory(0)  # avoid ROOT ownership issues
    if outfile:
        outfile.mkdir("blastwave_fit")
        outfile.cd("blastwave_fit")
        fBW.Write()
        hExt.Write()
        hist.Write()
    return hExt, fBW, fitResult


def scale_by_bin_width(hist):
    """(1/Nev) dN/dpTdy density -> 'counts per bin', for TH1::GetRandom() in
    the C++ toy MC (MomentumSampler treats bin content as a raw PDF)."""
    for i in range(1, hist.GetNbinsX() + 1):
        w = hist.GetBinWidth(i)
        hist.SetBinContent(i, hist.GetBinContent(i) * w)
        hist.SetBinError(i, hist.GetBinError(i) * w)

if __name__ == "__main__":
    
    f_he3    = ROOT.TFile.Open(HE3_FILE)
    f_proton = ROOT.TFile.Open(PROTON_FILE)

    if not f_he3 or f_he3.IsZombie():
        sys.exit(f"Cannot open {HE3_FILE}")
    if not f_proton or f_proton.IsZombie():
        sys.exit(f"Cannot open {PROTON_FILE}")

    hHe3      = f_he3.Get("hHe30_10")      # central values
    hHe3Stat  = f_he3.Get("hHe3Stat0_10")  # statistical errors
    hHe3Syst  = f_he3.Get("hHe3Syst0_10")  # systematic errors

    if not hHe3:
        sys.exit("He3 histogram 'hHe30_10' not found — check the histogram name.")

    if hHe3Stat and hHe3Syst:
        for i in range(1, hHe3.GetNbinsX() + 1):
            bin_width = hHe3.GetBinWidth(i)
            err_stat = hHe3Stat.GetBinError(i)
            err_syst = hHe3Syst.GetBinError(i)
            hHe3.SetBinContent(i, hHe3.GetBinContent(i) * bin_width)  # convert to density
            hHe3.SetBinError(i, (err_stat**2 + err_syst**2)**0.5 * bin_width)  # combine errors in quadrature

    hHe3.SetName("hHe3_0_10")
    hHe3.SetTitle("He3 p_{T} spectrum, 0-10% centrality, #sqrt{s_{NN}} = 5.02 TeV")
    hHe3.GetXaxis().SetTitle("p_{T} (GeV/c)")
    hHe3.GetYaxis().SetTitle("(1/N_{ev}) d^{2}N / dp_{T}dy [(GeV/c)^{-1}]")

    # ── Protons: 0-10% centrality (merge y1=0-5% and y2=5-10%) ───────────────────
    # Table 5: y1 → 0-5%, y2 → 5-10%
    # Each Hist1D_yN carries the bin content; _e1=stat, _e2=syst, _e3=syst.uncorr.
    hP_05  = f_proton.Get("Table 5/Hist1D_y1")
    hP_510 = f_proton.Get("Table 5/Hist1D_y2")

    hP_05_stat  = f_proton.Get("Table 5/Hist1D_y1_e1")
    hP_510_stat = f_proton.Get("Table 5/Hist1D_y2_e1")

    hP_05_syst  = f_proton.Get("Table 5/Hist1D_y1_e2")
    hP_510_syst = f_proton.Get("Table 5/Hist1D_y2_e2")

    if not hP_05 or not hP_510:
        sys.exit("Proton histograms not found in 'Table 5/' — check directory/names.")

    # Average the two centrality bins (they cover equal width: 5% each)
    hProtonDensity = hP_05.Clone("hProton_0_10")
    hProtonDensity.Add(hP_510)
    hProtonDensity.Scale(0.5)   # arithmetic mean → representative 0-10% spectrum
    
    #hProtonNoReweight = hProton.Clone("hProton_0_10_no_reweight")  # for comparison, if needed

    # Transfer stat errors (averaged in quadrature / linearly — linear here)
    if hP_05_stat and hP_510_stat:
        for i in range(1, hProtonDensity.GetNbinsX() + 1):
            e1_stat = hP_05_stat.GetBinContent(i)
            e2_stat = hP_510_stat.GetBinContent(i)
            err_stat = 0.5 * (e1_stat + e2_stat)

            e1_syst = hP_05_syst.GetBinContent(i) if hP_05_syst else 0
            e2_syst = hP_510_syst.GetBinContent(i) if hP_510_syst else 0
            err_syst = 0.5 * (e1_syst + e2_syst)

            hProtonDensity.SetBinError(i, (err_stat**2 + err_syst**2)**0.5)

    #for hist in [hProton, hProtonNoReweight]:
    for hist in [hProtonDensity]:
        hist.SetTitle("Proton p_{T} spectrum, 0-10% centrality, #sqrt{s_{NN}} = 5.02 TeV")
        hist.GetXaxis().SetTitle("p_{T} (GeV/c)")
        hist.GetYaxis().SetTitle("(1/N_{ev}) d^{2}N / dp_{T}dy [(GeV/c)^{-1}]")
            
    OUT = "spectra_0_10_blastwave.root"
    fout = ROOT.TFile(OUT, "RECREATE")
        
    hProtonExtrap, fBW, fitResult = extend_with_blastwave_fit(
        hProtonDensity, MASS_PROTON, pt_min=0.0, low_bin_width=0.05,
        beta0=0.65, T0=0.10, n0=0.9, fix_shape=False, outfile=fout
    )

    # QA plot: raw data points + fit + hybrid extrapolated spectrum -- eyeball
    # this before trusting the extrapolation for anything downstream
    canvas = ROOT.TCanvas("cBlastWaveQA", "Blast-wave low-pT extrapolation", 800, 600)
    canvas.SetLogy()
    hProtonExtrap.SetLineColor(ROOT.kAzure + 1)
    hProtonExtrap.SetMarkerStyle(20)
    hProtonExtrap.SetMarkerSize(0.6)
    hProtonExtrap.Draw("E1")
    hProtonDensity.SetMarkerStyle(24)
    hProtonDensity.SetMarkerColor(ROOT.kBlack)
    hProtonDensity.SetLineColor(ROOT.kBlack)
    hProtonDensity.Draw("E1 SAME")
    if fBW:
        fBW.SetLineColor(ROOT.kRed)
        fBW.Draw("SAME")
    legend = ROOT.TLegend(0.55, 0.65, 0.88, 0.88)
    legend.AddEntry(hProtonDensity, "HEPData (measured)", "ep")
    legend.AddEntry(hProtonExtrap, "hybrid (fit + data)", "ep")
    if fBW:
        legend.AddEntry(fBW, "BGBW fit", "l")
    legend.Draw()
    canvas.SaveAs("blastwave_lowpt_extrapolation.pdf")

    hProton = hProtonExtrap
    hProton.SetName("hProton_0_10")
    hProton.SetTitle("Proton p_{T} spectrum, 0-10% centrality, #sqrt{s_{NN}} = 5.02 TeV "
                    f"(BGBW-extrapolated below #it{{p}}_{{T}} = {hProtonDensity.GetXaxis().GetXmin():.2f} GeV/#it{{c}})")

    hProtonNoReweight = hProton.Clone("hProton_0_10_no_reweight")  # density units, pre bin-width scaling

    # ── NEW: convert to 'counts per bin' for TH1::GetRandom() (previously done
    #         inline for the data-only histogram; now applied once, after
    #         extrapolation, over the full hybrid spectrum) ────────────────────
    scale_by_bin_width(hProton)

    # ── save to output file ───────────────────────────────────────────────────────
    
    fout.cd()
    hHe3.Write()
    hProton.Write()
    hProtonNoReweight.Write()
    fout.Close()
    print(f"Saved '{hHe3.GetName()}' and '{hProton.GetName()}' to {OUT}")

    f_he3.Close()
    f_proton.Close()