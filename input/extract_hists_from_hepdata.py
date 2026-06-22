import ROOT
import sys

# ── file paths ────────────────────────────────────────────────────────────────
HE3_FILE    = "HEPData-He3-spectra-PbPb-5p02-GeV.root"   # adjust to your actual filename
PROTON_FILE = "HEPData-p-spectra-PbPb-5p02-GeV.root"  # adjust to your actual filename

# ── open files ────────────────────────────────────────────────────────────────
f_he3    = ROOT.TFile.Open(HE3_FILE)
f_proton = ROOT.TFile.Open(PROTON_FILE)

if not f_he3 or f_he3.IsZombie():
    sys.exit(f"Cannot open {HE3_FILE}")
if not f_proton or f_proton.IsZombie():
    sys.exit(f"Cannot open {PROTON_FILE}")

# ── He3: 0-10% centrality ─────────────────────────────────────────────────────
# The file contains separate 0-10% histogram directly
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
hProton = hP_05.Clone("hProton_0_10")
hProton.Add(hP_510)
hProton.Scale(0.5)   # arithmetic mean → representative 0-10% spectrum
hProtonNoReweight = hProton.Clone("hProton_0_10_no_reweight")  # for comparison, if needed

# Transfer stat errors (averaged in quadrature / linearly — linear here)
if hP_05_stat and hP_510_stat:
    for i in range(1, hProton.GetNbinsX() + 1):
        bin_width = hProton.GetBinWidth(i)
        e1_stat = hP_05_stat.GetBinContent(i)
        e2_stat = hP_510_stat.GetBinContent(i)
        err_stat = 0.5 * (e1_stat + e2_stat)  # average stat error
        
        e1_syst = hP_05_syst.GetBinContent(i) if hP_05_syst else 0
        e2_syst = hP_510_syst.GetBinContent(i) if hP_510_syst else 0
        err_syst = 0.5 * (e1_syst + e2_syst)  # average syst error
        
        error = (err_stat**2 + err_syst**2)**0.5  # combine in quadrature
        hProton.SetBinContent(i, hProton.GetBinContent(i) * bin_width)  # convert to density
        hProton.SetBinError(i, error * bin_width)  # average stat error
        
        hProtonNoReweight.SetBinError(i, error)

for hist in [hProton, hProtonNoReweight]:
    hist.SetTitle("Proton p_{T} spectrum, 0-10% centrality, #sqrt{s_{NN}} = 5.02 TeV")
    hist.GetXaxis().SetTitle("p_{T} (GeV/c)")
    hist.GetYaxis().SetTitle("(1/N_{ev}) d^{2}N / dp_{T}dy [(GeV/c)^{-1}]")



# ── save to output file ───────────────────────────────────────────────────────
OUT = "spectra_0_10.root"
fout = ROOT.TFile(OUT, "RECREATE")
hHe3.Write()
hProton.Write()
hProtonNoReweight.Write()
fout.Close()
print(f"Saved '{hHe3.GetName()}' and '{hProton.GetName()}' to {OUT}")

f_he3.Close()
f_proton.Close()