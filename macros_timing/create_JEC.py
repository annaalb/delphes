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
from ROOT import TFile,TCanvas,gROOT,gStyle,TLegend,TGraphAsymmErrors,kGreen,kOrange,kSpring,TF1,kAzure, TH2F,TH1F,gPad, TPaveText, TH1,kRed,SetOwnership, TMath, kBlue, kBlack, kFALSE, kTRUE,kSienna,TLatex, Double_t, TMultiGraph
from collections import OrderedDict
#import CMSPlotStyle
#import plot_functions as pf
#from hist_properties import *

gROOT.SetBatch(kTRUE)

#style = CMSPlotStyle.getStyle()
#style.cd()
#gROOT.SetStyle("CMS_Style")
gStyle.SetOptStat(0)
gStyle.SetPaintTextFormat("2.3f")

PUPPI_folder = "/eos/user/a/aalbrech/workspace/delphes_timing/Plots/study_timing_cut/VBF/Raw_response/"

PUPPI_EtaMax0 = TFile(PUPPI_folder+"EtaMax0/results.root")
PUPPI_EtaMax3 = TFile(PUPPI_folder+"EtaMax3/results.root")
PUPPI_EtaMax4 = TFile(PUPPI_folder+"EtaMax4/results.root")

PUPPI_file_list = {}
PUPPI_file_list["PUPPI_EtaMax0"] = PUPPI_EtaMax0
PUPPI_file_list["PUPPI_EtaMax3"] = PUPPI_EtaMax3
PUPPI_file_list["PUPPI_EtaMax4"] = PUPPI_EtaMax4

def create_JEC(file_list, filename, CHSorPUPPI):

    etabins = ["eta0p0to1p3","eta1p3to2p0","eta2p0to3","eta3to4","eta4p0toInf"]
    ptbins = ["pt0to20","pt20to30","pt30to50","pt50to80","pt80to100","pt100toInf"]

    resolution_dict = {}

    for el in file_list:
        resolution_dict[el] = {}
        for eta in etabins:
            resolution_dict[el][eta] = {}
            for pt in ptbins:

                print(el)
                hist = file_list[el].Get(CHSorPUPPI+"_RecoJet_response_"+eta+"_"+pt)
                print(CHSorPUPPI+"_RecoJet_response"+eta+pt)
                gaussian_fit = TF1("gaussian_fit", "gaus", 0.2, 2.0);
                gaussian_fit.SetParameter(1, 1.0);
                gaussian_fit.SetParameter(2, 0.1);
                hist.Fit(gaussian_fit, "R");


                for i in range(0,2):
                    lower_bound = gaussian_fit.GetParameter(1) -1.5*gaussian_fit.GetParameter(2);
                    higher_bound = gaussian_fit.GetParameter(1)+1.5*gaussian_fit.GetParameter(2);

                    gaussian_fit = TF1("gaussian_fit", "gaus",lower_bound,higher_bound);
                    gaussian_fit.SetParameter(1, gaussian_fit.GetParameter(1));
                    gaussian_fit.SetParameter(2, gaussian_fit.GetParameter(2));

                    fit_result = hist.Fit(gaussian_fit,"R");

                resolution_gauss = 0
                if (gaussian_fit.GetParameter(2) !=0 and gaussian_fit.GetParameter(1) !=0):
                    resolution_gauss = gaussian_fit.GetParameter(2)/gaussian_fit.GetParameter(1)


                print("#mu = %.3f"%(gaussian_fit.GetParameter(1)))
                print("#sigma = %.3f"%(resolution_gauss))

                if (gaussian_fit.GetParameter(2) !=0 and gaussian_fit.GetParameter(1) !=0):
                    resolution_dict[el][eta][pt] = 1+((1-gaussian_fit.GetParameter(1))/gaussian_fit.GetParameter(1))
                else:
                    resolution_dict[el][eta][pt] = 1
                print(gaussian_fit.GetParameter(1) * resolution_dict[el][eta][pt])



    print(resolution_dict)



    print()
    print()

    for el in resolution_dict:
        f = open(filename+"JEC_"+el+".txt", "w")
        for eta in sorted(resolution_dict[el]):
            correction_string="JEC_"+el+"_"+eta +" = ["
            f.write("\n")
            etamin = eta
            etamax = eta
            etamin = etamin.replace("eta","").split("to")[0].replace("p",".")
            etamax = etamax.replace("eta","").split("to")[1].replace("p",".")
            if "Inf" in etamax: etamax = '10000'

            for pt in sorted(resolution_dict[el][eta]):


                correction_string+=str(resolution_dict[el][eta][pt])+","
                f.write(etamin+"\t"+etamax+"\t")

                ptmin = pt
                ptmax = pt
                ptmin = ptmin.replace("pt","").split("to")[0].replace("p",".")
                ptmax = ptmax.replace("pt","").split("to")[1].replace("p",".")
                if "Inf" in ptmax: ptmax = '10000'


                f.write(ptmin+"\t"+ptmax + "\t")
                f.write(str(resolution_dict[el][eta][pt])+"\n")

            correction_string[:-1]
            correction_string+="]"

            print(correction_string)


    f.close()


#create_JEC(PUPPI_file_list, "PUPPI_", "PUPPI")
create_JEC(PUPPI_file_list, "CHS_", "CHS")
