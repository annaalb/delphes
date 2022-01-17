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
f.write("""set PResolutionFormula { 0.0 } """ + "\n"+ """set CtgThetaResolutionFormula { 0.0 } """ + "\n"+"""set PhiResolutionFormula { 0.0 } """ + "\n"+"""set D0ResolutionFormula { 0.0 } """ + "\n"+"""# taken from https://cms-tklayout.web.cern.ch/cms-tklayout/layouts-work/recent-layouts/OT800_IT700/errorstracker.html """ + "\n"+"""set DZResolutionFormula {""" + "\n")

for i in range(len(ptstrings)):
    for k in range(nbins):
        min = round((hists[i].GetBinLowEdge(k+1)), 4)
        max = round((hists[i].GetBinLowEdge(k+1) + hist1.GetBinWidth(k+1)), 4)
        res = round((hists[i].GetBinContent(k+1) / 1000), 4)
        f.write(""" ( abs(eta) > """ +str(min)+ " && abs(eta) <= "+ str(max) + ptstrings[i] + str(res) + " +\ " + "\n")
        pass
    pass
f.write("}" + "\n")
f.close()