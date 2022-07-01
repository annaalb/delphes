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
properties["eta0to3p0"] = [kAzure+2,1,"0 < |#eta| < 3.0"]
properties["eta3p0to5p0"] = [kBlack,8,"3.0 < |#eta| < 5.0"]

def plot_Hcal_resolution(xmax = 100):
    hist_to_draw = {}

    tmp_dict={}
    tmp_dict["eta0to3p0"]={}
    tmp_dict["eta3p0to5p0"]={}

    tmp_dict["eta0to3p0"]["1"]=0.05
    tmp_dict["eta0to3p0"]["2"]=1.0

    tmp_dict["eta3p0to5p0"]["1"]=0.11
    tmp_dict["eta3p0to5p0"]["2"]=2.80


#(abs(eta) <= 1.5) * sqrt(energy^2*0.05^2 + energy*1.00^2) + \
#(abs(eta) > 1.5 && abs(eta) <= 3.0) * sqrt(energy^2*0.05^2 + energy*1.00^2) + \
#(abs(eta) > 3.0 && abs(eta) <= 5.0) * sqrt(energy^2*0.11^2 + energy*2.80^2)}

    for etabin in sorted(tmp_dict):

        print(etabin)

        c = TCanvas()

        hist = TH1F(etabin,etabin,100,0,xmax)
        for xbin in range(hist.GetNbinsX()+2):
            eff = 0
            energy = hist.GetBinCenter(xbin)
            eff = math.sqrt((energy**2*tmp_dict[etabin]["1"]**2) + ((energy*tmp_dict[etabin]["2"])**2))
            hist.SetBinContent(xbin,eff)

        hist_to_draw[etabin] = hist

    pf.plot(hist_to_draw,"../../PLOTS/resolution/","Hcal_resolution", "Energy [GeV]", "HCal Resolution [GeV]", 0, 1.01, 0.2,0.7,0.5,0.9, properties, cmsextratext = "", blogo = False,yscalefactor=0, btline = False)





##########
#
# MAIN
#
########


plot_Hcal_resolution()
