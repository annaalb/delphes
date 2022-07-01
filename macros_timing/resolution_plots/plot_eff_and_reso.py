#!/usr/bin/env python
from __future__ import division
import subprocess
import glob
from os import system
import sys
from os import mkdir
from os.path import exists
from array import array
import math
import re
import time
import datetime
from ROOT import TFile,TCanvas,gROOT,gStyle,TLegend,TGraphAsymmErrors,kGreen,kOrange,kSpring,TF1,kAzure, TH2F,TH1F,gPad, TPaveText, TH1,kRed,SetOwnership, TMath, kBlue, kBlack, kFALSE, kTRUE,kSienna,TLatex, TMultiGraph,  TLine, kDashed
from collections import OrderedDict
import CMSPlotStyle
import plot_functions as pf


gROOT.SetBatch(kTRUE)

style = CMSPlotStyle.getStyle()
style.cd()
gROOT.SetStyle("CMS_Style")
gStyle.SetOptStat(0)
gStyle.SetPaintTextFormat("2.3f")


properties = {}
properties["eta0to1p2"] = [kBlue,1,"0 < |#eta| < 1.2"]
properties["eta1p2to2p5"] = [kOrange+2,3,"1.2 < |#eta| < 2.5"]
properties["eta2p5to4p0"] = [kGreen+2,5,"2.5 < |#eta| < 4.0"]

def plot_ChargedHadronTrackingEfficiency(xmax = 5):
    hist_to_draw = {}

    tmp_dict={}
    tmp_dict["eta0to1p2"]={}
    tmp_dict["eta1p2to2p5"]={}
    tmp_dict["eta2p5to4p0"]={}

    tmp_dict["eta0to1p2"]["kpt"] = 0.96
    tmp_dict["eta0to1p2"]["gpt"] = 0.97

    tmp_dict["eta1p2to2p5"]["kpt"] = 0.85
    tmp_dict["eta1p2to2p5"]["gpt"] = 0.87

    tmp_dict["eta2p5to4p0"]["kpt"] = 0.8
    tmp_dict["eta2p5to4p0"]["gpt"] = 0.82


    for etabin in sorted(tmp_dict):

        print(etabin)

        c = TCanvas()

        hist = TH1F(etabin,etabin,100,0,xmax)
        for xbin in range(hist.GetNbinsX()+2):
            eff = 0
            if hist.GetBinCenter(xbin) < 1:
                eff = tmp_dict[etabin]["kpt"] * hist.GetBinCenter(xbin)
            else:
                eff = tmp_dict[etabin]["gpt"]
            hist.SetBinContent(xbin,eff)

        hist.SetMarkerSize(2)
        hist.SetLineWidth(2)
        hist.GetXaxis().SetTitle("track p_{T} [GeV]")
        hist.GetYaxis().SetTitle("Efficiency")
        hist.GetYaxis().SetRangeUser(0,1.01)

        hist_to_draw[etabin] = hist

        #hist.Draw("HIST")

        #line = TLine (0,1,20,1)
        #line.SetLineColor(kRed)
        #line.SetLineStyle(kDashed)
        #line.Draw("same")


        #c.SaveAs("../../PLOTS/resolution/ChargedHadronTrackingEfficiency_"+etabin+".pdf")

    pf.plot(hist_to_draw,"plots/","ChargedHadronTrackingEfficiency", "track p_{T} [GeV]", "Efficiency", 0, 1.01, 0.5,0.2,0.9,0.5, properties, cmsextratext = "", blogo = False, btline = True)





##########
#
# MAIN
#
########


plot_ChargedHadronTrackingEfficiency()
