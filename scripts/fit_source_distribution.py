import array

from ROOT import TFile, TF1, TH1D, TCanvas, TLine

from torchic.core.graph import load_graph
from torchic.core.histogram import load_hist
from torchic.utils.root import set_alice_global_style, init_legend, set_root_object
from torchic.utils.colors import get_color

def graph_to_hist(graph, name="h_fromgraph"):
    """Build a TH1D from a TGraphAsymmErrors, using the x-errors as bin edges
    and the y-errors (max of low/high) as bin errors."""
    n = graph.GetN()
 
    xs = graph.GetX()
    ys = graph.GetY()
 
    edges = []
    for i in range(n):
        exl = graph.GetErrorXlow(i)
        edges.append(xs[i] - exl)
    exh_last = graph.GetErrorXhigh(n - 1)
    edges.append(xs[n - 1] + exh_last)
 
    edges_arr = array.array('d', edges)
    h = TH1D(name, name, n, edges_arr)
 
    for i in range(n):
        h.SetBinContent(i + 1, ys[i])
        eyl = graph.GetErrorYlow(i)
        eyh = graph.GetErrorYhigh(i)
        h.SetBinError(i + 1, max(eyl, eyh))
 
    return h

def ratio_with_different_bins(h1, h2, name="h_ratio"):
    '''Skips the bins that are different in the two histograms. Returns a new histogram with the ratio.'''
    
    if h1.GetNbinsX() != h2.GetNbinsX():
        print("Warning: histograms have different number of bins. Skipping bins with different edges.")
    
    edges = []
    for i in range(1, min(h1.GetNbinsX(), h2.GetNbinsX()) + 1):
        if (abs(h1.GetBinLowEdge(i) - h2.GetBinLowEdge(i)) < 1e-9) and (abs(h1.GetBinLowEdge(i + 1) - h2.GetBinLowEdge(i + 1)) < 1e-9):
            edges.append(h1.GetBinLowEdge(i))
    edges.append(h1.GetBinLowEdge(min(h1.GetNbinsX(), h2.GetNbinsX()) + 1))
    edges_arr = array.array('d', edges)
    h_ratio = TH1D(name, name, len(edges) - 1, edges_arr)
    for i in range(1, min(h1.GetNbinsX(), h2.GetNbinsX()) + 1):
        if (abs(h1.GetBinLowEdge(i) - h2.GetBinLowEdge(i)) < 1e-9) and (abs(h1.GetBinLowEdge(i + 1) - h2.GetBinLowEdge(i + 1)) < 1e-9):
            ratio = h2.GetBinContent(i) / h1.GetBinContent(i) if h1.GetBinContent(i) != 0 else 0
            ratio_error = (ratio * ((h1.GetBinError(i) / h1.GetBinContent(i))**2 + (h2.GetBinError(i) / h2.GetBinContent(i))**2)**0.5 
                           if h1.GetBinContent(i) != 0 and h2.GetBinContent(i) != 0 else 0)
            h_ratio.SetBinContent(i, ratio)
            h_ratio.SetBinError(i, ratio_error)
    return h_ratio

def compare_nuclear_spectra(hepdata_file, coalescence_file, graph_name, coalescence_hist_name, n_events,
                            outfile, nucleus_name="deuteron"):
        
    g_hepdata = load_graph(hepdata_file, graph_name)
    h_hepdata = graph_to_hist(g_hepdata, "h_hepdata")
    h_coalescence = load_hist(coalescence_file, coalescence_hist_name)
    h_coalescence.Scale(1.0 / (n_events * h_coalescence.GetBinWidth(1))) # assuming uniform bin width, scale to 1/N_events * 1/bin_width
 
    c1 = TCanvas(f"comparison_{nucleus_name}", f"comparison_{nucleus_name}", 800, 600)
    set_root_object(h_hepdata, marker_style=20, marker_color=get_color(0), line_color=get_color(0))
    hframe = c1.DrawFrame(h_hepdata.GetXaxis().GetXmin(), 0, h_hepdata.GetXaxis().GetXmax(), 
                          1.5 * max(h_hepdata.GetMaximum(), h_coalescence.GetMaximum()),
                          "Comparison of nuclear spectra; #it{{p}}_{{T}} (GeV/#it{{c}});"
                           "1/#it{{N}}_{{events}} d^{{2}}#it{a{N}}/d#it{{p}}_{{T}}d#it{{y}}")
    h_hepdata.Draw("E SAME")
    
    set_root_object(h_coalescence, marker_style=21, marker_color=get_color(1), line_color=get_color(1))
    h_coalescence.Draw("E SAME")
 
    leg = init_legend(0.4, 0.2, 0.7, 0.4)
    leg.AddEntry(h_hepdata, "Run 2, Pb-Pb", "lep")
    leg.AddEntry(h_coalescence, "AV18", "lep")
    leg.Draw()
    
    h_ratio = ratio_with_different_bins(h_coalescence, h_hepdata, "h_ratio")
 
    c2 = TCanvas(f"ratio_{nucleus_name}", f"ratio_{nucleus_name}", 800, 600)
    h_ratio.SetTitle(f"AV18 / Run 2, Pb-Pb;#it{{p}}_{{T}} (GeV/#it{{c}}); Ratio")
    h_ratio.SetMinimum(0.0)
    h_ratio.SetMaximum(2.0)
 
    line = TLine(h_ratio.GetXaxis().GetXmin(), 1, h_ratio.GetXaxis().GetXmax(), 1)
    line.SetLineStyle(2)
    line.Draw()
    h_ratio.Draw("E")
 
    outfile.cd()
    c1.Write()
    c2.Write()

def fit_source_distribution(hist, f_source, output_file, init_radius):
    
    f_source.SetParameters(hist.Integral(), init_radius)
    
    hist.Fit(f_source, "RMS+")
    output_file.cd()
    hist.Write()
    f_source.Write()

if __name__ == "__main__":
    
    set_alice_global_style()
    
    output_file = TFile("../output/fit_source_d_av18_Rp_6p02_fm.root", "RECREATE")
    f_source = TF1("f_source", "[0]* x^2 *exp(-x^2/(2*[1]*[1]))", 0, 40)
    
    COALESCENCE_FILE = "../output/output_d_av18_Rp_6p02_fm.root"
    
    h_source_nucleons = load_hist(COALESCENCE_FILE, "hPositionNucleons_0")
    fit_source_distribution(h_source_nucleons, f_source, output_file, 6.02)
    
    h_source_nuclei = load_hist(COALESCENCE_FILE, "hPositionNucleus_0")
    fit_source_distribution(h_source_nuclei, f_source, output_file, 4.92)
    
    h_source_relative = load_hist(COALESCENCE_FILE, "hRelativePosition_0")
    fit_source_distribution(h_source_relative, f_source, output_file, 7.72)
    
    HEPDATA_FILE        = "/home/galucia/CoalescenceAfterburner/input/HEPData-ins2667337-v1-Deuteron_spectrum_in_0-5%_V0M_centrality_class.root"
    HEPDATA_GRAPH_PATH  = "Deuteron spectrum in 0-5% V0M centrality class/Graph1D_y1"
    HIST_NAME           = "hPtNucleus_0"
    N_EVENTS            = 100_000
    compare_nuclear_spectra(HEPDATA_FILE, COALESCENCE_FILE, HEPDATA_GRAPH_PATH, HIST_NAME, N_EVENTS,
                            output_file, nucleus_name="deuteron")
    
    HEPDATA_FILE        = "/home/galucia/CoalescenceAfterburner/input/HEPData-ins2667337-v1-Deuteron_spectrum_in_0-5%_V0M_centrality_class.root"
    HEPDATA_GRAPH_PATH  = "Deuteron spectrum in 0-5% V0M centrality class/Graph1D_y1"
    HIST_NAME           = "hPtNucleus_0"
    N_EVENTS            = 100_000
    COALESCENCE_FILE_GAUS = "../output/output_d_gaus_Rp_6p02_fm.root"
    compare_nuclear_spectra(HEPDATA_FILE, COALESCENCE_FILE_GAUS, HEPDATA_GRAPH_PATH, HIST_NAME, N_EVENTS,
                            output_file, nucleus_name="deuteron_gaus")
    
    # protons
    HEPDATA_FILE_P       = "/home/galucia/CoalescenceAfterburner/input/HEPData-p-spectra-PbPb-5p02-GeV.root"
    HEPDATA_GRAPH_PATH_P = "Table 5/Graph1D_y1" # note: this is 0-5%
    HIST_NAME_p          = "hPtProtonOneRapidityUnit_0"
    
    compare_nuclear_spectra(HEPDATA_FILE_P, COALESCENCE_FILE, HEPDATA_GRAPH_PATH_P, HIST_NAME_p, N_EVENTS,
                                output_file, nucleus_name="proton")
    
    output_file.Close()