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
properties["eta0to1p5"] = [kBlue,1,"0 < |#eta| < 1.5"]
properties["eta2p15to1p75"] = [kGreen+2,5,"1.5 < |#eta| < 1.75"]
properties["eta1p75to2p15"] = [kOrange+2,3,"1.75 < |#eta| < 2.15"]
properties["eta2p15to3p0"] = [kAzure+2,7,"2.15 < |#eta| < 3.0"]
properties["eta3p0to5p0"] = [kBlack,8,"3.0 < |#eta| < 5.0"]

def plot_Ecal_resolution(xmax = 100):
    hist_to_draw = {}

    tmp_dict={}
    tmp_dict["eta0to1p5"]={}
    tmp_dict["eta2p15to1p75"]={}
    tmp_dict["eta1p75to2p15"]={}
    tmp_dict["eta2p15to3p0"]={}
    tmp_dict["eta3p0to5p0"]={}

    tmp_dict["eta0to1p5"]["1"]=0.009
    tmp_dict["eta0to1p5"]["2"]=0.12
    tmp_dict["eta0to1p5"]["3"]=0.45

    tmp_dict["eta2p15to1p75"]["1"]=0.006
    tmp_dict["eta2p15to1p75"]["2"]=0.20
    tmp_dict["eta2p15to1p75"]["3"]=0.0

    tmp_dict["eta1p75to2p15"]["1"]=0.007
    tmp_dict["eta1p75to2p15"]["2"]=0.21
    tmp_dict["eta1p75to2p15"]["3"]=0.0

    tmp_dict["eta2p15to3p0"]["1"]=0.008
    tmp_dict["eta2p15to3p0"]["2"]=0.24
    tmp_dict["eta2p15to3p0"]["3"]=0.0

    tmp_dict["eta3p0to5p0"]["1"]=0.08
    tmp_dict["eta3p0to5p0"]["2"]=1.98
    tmp_dict["eta3p0to5p0"]["3"]=0.0


# {  0.5*( 0.5*(abs(eta) <= 1.50)          * sqrt(energy^2*0.009^2 + energy*0.12^2 + 0.45^2) +
# (abs(eta) > 1.50 && abs(eta) <= 1.75) * sqrt(energy^2*0.006^2 + energy*0.20^2) + \
#(abs(eta) > 1.75 && abs(eta) <= 2.15) * sqrt(energy^2*0.007^2 + energy*0.21^2) + \
#(abs(eta) > 2.15 && abs(eta) <= 3.00) * sqrt(energy^2*0.008^2 + energy*0.24^2) + \
#(abs(eta) >= 3.0 && abs(eta) <= 5.0)  * sqrt(energy^2*0.08^2 + energy*1.98^2))}

    for etabin in sorted(tmp_dict):

        print(etabin)

        c = TCanvas()

        hist = TH1F(etabin,etabin,100,0,xmax)
        for xbin in range(hist.GetNbinsX()+2):
            eff = 0
            energy = hist.GetBinCenter(xbin)
            eff = math.sqrt((energy**2*tmp_dict[etabin]["1"]**2) + ((energy*tmp_dict[etabin]["2"])**2) + (energy*tmp_dict[etabin]["3"])**2)
            if etabin == "eta0to1p5":
                eff = eff*0.5
                pass
            hist.SetBinContent(xbin,eff)

        hist_to_draw[etabin] = hist

    pf.plot(hist_to_draw,"../../PLOTS/resolution/","Ecal_resolution", " Energy [GeV]", "ECal Resolution [GeV]", 0, 1.01, 0.2,0.6,0.5,0.9, properties, cmsextratext = "", blogo = False, yscalefactor=0, btline = False)




##########
#
# MAIN
#
########


plot_Ecal_resolution()
