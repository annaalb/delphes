import os,sys
from ROOT import TFile,TTree, TCanvas, TH1F
import numpy as np

root_file = TFile("/afs/cern.ch/user/a/aalbrech/Delphes/cards/CMS_PhaseII/resolution/dzres_tracker_MS_Pt000.root", "READ");
tcl_file = "/afs/cern.ch/user/a/aalbrech/Delphes/cards/CMS_PhaseII/resolution/trackResolutionCMSPhaseII.tcl";

canvas = root_file.Get("dzres_tracker_MS_Pt000")
hist1 = canvas.GetPrimitive("Resolution_tracker_dzres_tracker_MS_Pt_z_vs_eta1tracker_profile") # get name of hist
hist2 = canvas.GetPrimitive("Resolution_tracker_dzres_tracker_MS_Pt_z_vs_eta10tracker_profile") # get name of hist
hist3 = canvas.GetPrimitive("Resolution_tracker_dzres_tracker_MS_Pt_z_vs_eta100tracker_profile") # get name of hist

hists = [hist1, hist2, hist3]

ptstring1 = """  ) * ( pt > 0.0 && pt <= 8.0 ) * """
ptstring2 = """  ) * ( pt > 8.0 && pt <= 50.0 ) * """
ptstring3 = """  ) * ( pt > 50.0 ) * """

ptstrings = [ptstring1, ptstring2, ptstring3]

nbins = hist1.GetNbinsX()

f = open(tcl_file, 'w')
f.write("set DZResolutionFormula {" + "\n")

for i in range(len(ptstrings)):
    for k in range(nbins):
        min = hists[i].GetBinLowEdge(k+1)
        max = hists[i].GetBinLowEdge(k+1) + hist1.GetBinWidth(k+1)
        res = hists[i].GetBinContent(k+1) / 1000
        f.write(""" ( abs(eta) > """ +str(min)+ " && abs(eta) <= "+ str(max) + ptstrings[i] + str(res) + " +\ " + "\n")
        pass
    pass
f.write("}" + "\n")
f.close()
