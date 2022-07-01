#!/usr/bin/env python
from ROOT import TFile,TCanvas,gROOT,gStyle,TLegend,TGraphAsymmErrors,kGreen,kOrange,kSpring,TF1,kAzure, TH2F,TH1F,gPad, TPaveText, TH1,kRed,SetOwnership, TMath, kBlue, kBlack, kFALSE, kTRUE,kSienna,TLatex, TMultiGraph, kMagenta, TPad, TLine, kDashed, kGray, TGaxis,kWhite
import CMSPlotStyle
import os

class bcolors:
    HEADER = '\033[95m'
    RED = '\033[91m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW  = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    BR = '\033[94m' #BLue
    PROCESS = '\033[93m' #Yellow
    CHANNEL = '\033[95m' #Magenta
    SYSTEMATIC = '\033[96m' #Cyan

gROOT.SetBatch(kTRUE)


gStyle.SetOptStat(0)
gStyle.SetOptFit(0)
gStyle.SetPaintTextFormat("2.3f")

yOffset = 1.2
xOffset = 1.2
ySize = 0.05
xSize =0.05
xLabelSize = 0.05
yLabelSize = 0.05
ndivisions = 505

colors = [kBlue,kBlack, kRed]
beps = False

def plot(hist_to_draw,folder,filename,xtitle,ytitle,ymin,ymax,legxmin,legymin,legxmax,legymax,properties,blogy=False, draw_opt = "HIST", extratext = "", b_draw_leg = True, xmin = 0, xmax =100, legheader = "",blogo = True,cmsextratext = "work in progress",yscalefactor = 0.3, btline = False):

    if blogy and ymin==0: ymin=0.0001
    if blogy and "weight" in filename: ymin = 0.000001
    if "4D" in folder and "dt_dzk0p1" in filename: xmax = 0.6
    c = TCanvas()
    leg=TLegend(legxmin,legymin,legxmax,legymax,"","brNDC")
    leg.SetBorderSize(0);
    leg.SetTextSize(0.035);
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    leg.SetLineColor(1);
    leg.SetTextFont(42);
    leg.SetHeader(legheader)

    ymaxi = 0
    ypos = .8
    extra_dict = {}
    firstkey = ''
    for key in sorted(hist_to_draw.keys()):
        if not ymaxi: firstkey = key
        if hist_to_draw[key].GetMaximum() > ymaxi: ymaxi = hist_to_draw[key].GetMaximum()
        if xmin ==0: xmin = hist_to_draw[key].GetBinLowEdge(hist_to_draw[key].FindFirstBinAbove(0)-2)

        if key in properties:
            hist_to_draw[key].SetLineColor(properties[key][0])
            hist_to_draw[key].SetMarkerColor(properties[key][0])
            hist_to_draw[key].SetLineStyle(properties[key][1])
        else:
            hist_to_draw[key].SetLineColor(properties["others"][0])
            hist_to_draw[key].SetMarkerColor(properties["others"][0])
            hist_to_draw[key].SetLineStyle(properties["others"][1])

        hist_to_draw[key].SetMarkerSize(2)
        hist_to_draw[key].SetLineWidth(2)
        hist_to_draw[key].GetXaxis().SetTitle(xtitle)
        hist_to_draw[key].GetXaxis().SetTitleSize(xSize)
        hist_to_draw[key].GetXaxis().SetTitleOffset(xOffset)
        hist_to_draw[key].GetXaxis().SetLabelSize(xLabelSize)
        hist_to_draw[key].GetYaxis().SetTitle(ytitle)
        hist_to_draw[key].GetYaxis().SetTitleSize(ySize)
        hist_to_draw[key].GetYaxis().SetTitleOffset(yOffset)
        hist_to_draw[key].GetYaxis().SetLabelSize(yLabelSize)

        hist_to_draw[key].GetXaxis().SetRangeUser(xmin,xmax)
        hist_to_draw[key].GetYaxis().SetNdivisions(ndivisions)

        hist_to_draw[key].Draw(draw_opt+"same")



        if key in properties:
            leg.AddEntry(hist_to_draw[key], properties[key][2],"l")
        else:
            leg.AddEntry(hist_to_draw[key], properties["others"][2],"lpe")

        if "JetResponse" in filename:
            gaussian_fit = TF1("gaussian_fit", "gaus", 0.3, 2.0);
            gaussian_fit.SetParameter(1, 1.0);
            gaussian_fit.SetParameter(2, 0.1);
            hist_to_draw[key].Fit(gaussian_fit, "R");


            for i in range(0,2):
                lower_bound = gaussian_fit.GetParameter(1) -1.5*gaussian_fit.GetParameter(2);
                higher_bound = gaussian_fit.GetParameter(1)+1.5*gaussian_fit.GetParameter(2);

                gaussian_fit = TF1("gaussian_fit", "gaus",lower_bound,higher_bound);
                gaussian_fit.SetParameter(1, gaussian_fit.GetParameter(1));
                gaussian_fit.SetParameter(2, gaussian_fit.GetParameter(2));

                fit_result = hist_to_draw[key].Fit(gaussian_fit,"R");

            if key in properties:
                gaussian_fit.SetLineColor(properties[key][0])
                gaussian_fit.SetLineStyle(1)
                gaussian_fit.SetLineWidth(1)

            gaussian_fit.Draw("same")
            resolution_gauss = 0
            resolution_gauss_error = 0
            if gaussian_fit.GetParameter(2) !=0:
                resolution_gauss = gaussian_fit.GetParameter(2)/gaussian_fit.GetParameter(1)
                resolution_gauss_error = (gaussian_fit.GetParError(2)/gaussian_fit.GetParameter(1)) + (gaussian_fit.GetParError(1) * (gaussian_fit.GetParameter(2)/(gaussian_fit.GetParameter(1)*gaussian_fit.GetParameter(1))))


            # extra_jer = TLatex(3.5, 24, "#splitline{#mu = %.2f"%(hist_to_draw[key].GetMean())+"}{#sigma = %.2f"%(hist_to_draw[key].GetRMS())+"}");
            extra_jer = TLatex(3.5, 24, "#splitline{#mu = %.3f #pm %.3f"%(gaussian_fit.GetParameter(1),gaussian_fit.GetParError(1))+"}{#sigma/#mu = %.3f #pm %.3f"%(resolution_gauss, resolution_gauss_error)+"}");
            extra_jer.SetNDC()
            extra_jer.SetX(0.65)
            extra_jer.SetY(ypos)
            extra_jer.SetTextFont(52)
            extra_jer.SetTextSize(0.035)
            if key in properties:
                extra_jer.SetTextColor(properties[key][0])
            extra_dict[key] = extra_jer
            ypos-=0.1

        if "eff" in filename:
            extra_eff = TLatex(3.5, 24, "#splitline{eff = %.3f "%(hist_to_draw[key].GetBinContent(1))+"}{purity = %.3f"%(hist_to_draw[key].GetBinContent(2))+"}");
            extra_eff.SetNDC()
            extra_eff.SetX(0.65)
            extra_eff.SetY(ypos)
            extra_eff.SetTextFont(52)
            extra_eff.SetTextSize(0.035)
            if key in properties:
                extra_eff.SetTextColor(properties[key][0])
            extra_dict[key] = extra_eff
            ypos-=0.1


    if firstkey in hist_to_draw:
        if ymaxi < ymax:         hist_to_draw[firstkey].GetYaxis().SetRangeUser(ymin,ymax)
        else: hist_to_draw[firstkey].GetYaxis().SetRangeUser(ymin,ymaxi+yscalefactor * ymaxi)


    if b_draw_leg: leg.Draw()

    extra_text = TLatex(3.5, 24, extratext);
    extra_text.SetNDC()
    extra_text.SetX(0.5)
    extra_text.SetY(0.5)
    extra_text.SetTextFont(52)
    extra_text.SetTextSize(0.035)
    extra_text.Draw()

    for key in extra_dict:
        extra_dict[key].Draw()

    CMSPlotStyle.extratext = cmsextratext
    if blogo:
        text = CMSPlotStyle.draw_cmstext("left", True,1.5)
        text[0].Draw()
        text[1].Draw()
    lumi = CMSPlotStyle.draw_lumi(True)
    lumi.Draw()

    if btline:
        line = TLine (0,1,hist_to_draw[key].GetBinLowEdge(hist_to_draw[key].GetNbinsX()+1),1)
        line.SetLineColor(kRed)
        line.SetLineStyle(kDashed)
        line.Draw("same")


    if blogy: c.SetLogy()

    if not os.path.exists(folder):
        os.makedirs(folder)

    if beps: c.Print(folder+filename+".eps")
    c.Print(folder+filename+".pdf")


def plot_fit(hist,fitfunction,folder,filename,xtitle,ytitle,ymin,ymax,legxmin,legymin,legxmax,legymax,blogy, draw_opt, markerstyle, extratext = "", extrax=0.5,extray=0.5,b_draw_leg = True, plot_stats=False, xmin = 0 ,xmax=100, mcsample = "Z+jets"):

    if plot_stats:
        gStyle.SetOptFit(1111);
        gStyle.SetStatX(0.9);
        gStyle.SetStatY(0.9);

    c = TCanvas()
    leg=TLegend(legxmin,legymin,legxmax,legymax,"","brNDC")
    leg.SetBorderSize(0);
    leg.SetTextSize(0.035);
    leg.SetFillColor(0);
    leg.SetLineColor(1);
    leg.SetTextFont(42);

    i=0
    hist.SetLineColor(colors[i])
    hist.SetMarkerColor(colors[i])
    hist.SetMarkerStyle(markerstyle)
    hist.GetXaxis().SetTitle(xtitle)
    hist.GetXaxis().SetTitleSize(xSize)
    hist.GetXaxis().SetTitleOffset(xOffset)
    hist.GetYaxis().SetTitle(ytitle)
    hist.GetYaxis().SetTitleSize(ySize)
    hist.GetYaxis().SetTitleOffset(yOffset)
    hist.GetYaxis().SetRangeUser(ymin,ymax)
    hist.GetXaxis().SetRangeUser(xmin,xmax)

    hist.Fit(fitfunction,"R")
    leg.AddEntry(hist,"MC, "+mcsample,"lpe")
    leg.AddEntry(fitfunction,"N_{vertices} = %.3f #mu + %.3f"%(fitfunction.GetParameter(1),fitfunction.GetParameter(0)),"l")

    hist.Draw(draw_opt+"same")

    extrainfo = CMSPlotStyle.draw_info(extratext,extrax,extray)
    extrainfo.Draw()

#    CMSPlotStyle.extratext = "#splitline{Simulation}{Supplementary}"
    CMSPlotStyle.extratext = "Simulation"
    text = CMSPlotStyle.draw_cmstext("left", True,1.2)
    text[0].Draw()
    text[1].Draw()
    lumi = CMSPlotStyle.draw_lumi(True)
    lumi.Draw()

    print("======================= %.2f"%(fitfunction.GetNDF()))
    info = CMSPlotStyle.draw_info("#chi^{2}/ndf = %.2f"%(fitfunction.GetChisquare()/fitfunction.GetNDF()), legxmin+0.27,legymin-0.05)
    info.Draw()

    if b_draw_leg: leg.Draw()
    if blogy: c.SetLogy()
    if beps: c.Print(folder+filename+".eps")
    c.Print(folder+filename+".pdf")



def plot_ratio(hist_to_draw,ratio_to_draw,folder,filename,xtitle,ytitle,ymin,ymax,legxmin,legymin,legxmax,legymax,blogy, draw_opt,leg_opt = "pe", extratext = "", ratiotitle = "Data / MC", xmin=0, xmax=200, ratioymin=0.5, ratioymax=1.5,extrax=0.9, extray=0.82, extratext2="", extrax2=0.9, extray2=0.82, markers = {1,1,1,1,1,1,1}, bratioleg = False, lumitext = "35.9",legheader = "",blogo=True, antix = 0.9, antiy = 0.87,errors = {},xtitlesize = 0.18,xtitleoffset = 1.2,b_anti = True):

    c = TCanvas()

    yplot = 0.69
    yratio = 0.25                                 #  y6 +-------------+
    y6 = 0.97                                     #     |     pad1    |
    y5 = y6-yplot                                 #  y5 |-------------|
    y4 = y5-yratio                                #     |     rp1     |
    y4 = 0.01
    x1 = 0.01                                     #  y4 +-------------+
    x2 = 0.95
    leftmargin = 0.15
    rightmargin = 0.045

    pad_top = TPad("pad", "Control Plots", x1, y5, x2, y6)
    pad_top.SetTopMargin(0.05);
    pad_top.SetBottomMargin(0.025);
    pad_top.SetLeftMargin(leftmargin);
    pad_top.SetRightMargin(rightmargin);
    pad_top.Draw();
    pad_ratio = TPad("rp", "Ratio", x1, y4, x2, y5);
    pad_ratio.SetTopMargin(0);
    pad_ratio.SetBottomMargin(0.47);
    pad_ratio.SetLeftMargin(leftmargin);
    pad_ratio.SetRightMargin(rightmargin);
    pad_ratio.Draw();

    leg=TLegend(legxmin,legymin,legxmax,legymax,"","brNDC")
    leg.SetBorderSize(0);
    leg.SetTextSize(0.06);
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    leg.SetLineColor(1);
    leg.SetTextFont(42);
    leg.SetHeader(legheader)
    c.GetFrame().SetFillColor(21);
    c.GetFrame().SetBorderSize(12);

    i=0
    for key in sorted(hist_to_draw.keys()):
        hist_to_draw[key].SetLineColor(colors[i])
        hist_to_draw[key].SetFillColor(colors[i])
        hist_to_draw[key].SetLineWidth(2)
        hist_to_draw[key].SetMarkerColor(colors[i])

        hist_to_draw[key].SetMarkerStyle(markers[i])
        hist_to_draw[key].SetMarkerSize(1.2)
        hist_to_draw[key].GetXaxis().SetLabelSize(0)
        hist_to_draw[key].GetYaxis().SetLabelSize(0.05)
        hist_to_draw[key].GetXaxis().SetTitle(xtitle)
        hist_to_draw[key].GetYaxis().SetTitle(ytitle)
        hist_to_draw[key].GetYaxis().SetTitleSize(0.09)
        hist_to_draw[key].GetYaxis().SetTitleOffset(0.8)
        hist_to_draw[key].GetYaxis().SetRangeUser(ymin,ymax)
        hist_to_draw[key].GetXaxis().SetRangeUser(xmin,xmax)
        hist_to_draw[key].GetXaxis().SetTitleSize(0)





        pad_top.cd()
        title = key[1:]
        if "ata" in key:
            hist_to_draw[key].Draw("E same")
            leg.AddEntry(hist_to_draw[key], title,"pe")
        else:
            hist_to_draw[key].Draw("HIST same")
            leg.AddEntry(hist_to_draw[key], title,leg_opt)



        i+=1

    if "staterror_top" in errors:
        errors["staterror_top"].SetFillColor(kGray)
        errors["staterror_top"].SetLineColor(kGray)
        errors["staterror_top"].SetLineWidth(0)
        errors["staterror_top"].SetFillStyle(3244)
        errors["staterror_top"].SetMarkerColor(kGray)
        errors["staterror_top"].Draw("E2 same")
        leg.AddEntry(errors["staterror_top"],"Stat. unc.", "F")




    hist_to_draw["1Data"].Draw("E same")
    i=0
    pad_ratio.cd()

    ###### stat error
    if "staterr" in errors:
        errors["staterr"].SetFillColor(kGray)
        errors["staterr"].SetLineColor(kGray)
        errors["staterr"].GetXaxis().SetTitle(xtitle)
        errors["staterr"].GetXaxis().SetLabelSize(0.15)
        errors["staterr"].GetXaxis().SetTitleSize(xtitlesize)
        errors["staterr"].GetXaxis().SetTitleOffset(xtitleoffset)
        errors["staterr"].GetXaxis().SetRangeUser(xmin,xmax)
        errors["staterr"].GetYaxis().SetTitle(ratiotitle)
        errors["staterr"].GetYaxis().SetLabelSize(0.15)
        errors["staterr"].GetYaxis().SetTitleSize(0.15)
        errors["staterr"].GetYaxis().SetTitleOffset(0.45)
        errors["staterr"].GetYaxis().SetNdivisions(503)
        errors["staterr"].GetYaxis().SetRangeUser(ratioymin,ratioymax)
        errors["staterr"].GetYaxis().SetNdivisions(ndivisions)

        errors["staterr"].Draw("E2")


    ratio_leg=TLegend(0.17,0.81,0.37,0.95,"","brNDC")
    ratio_leg.SetBorderSize(0);
    ratio_leg.SetTextSize(0.12);
    ratio_leg.SetFillColor(0);
    ratio_leg.SetLineColor(1);
    ratio_leg.SetTextFont(42);

    ratio2_leg=TLegend(0.45,0.81,0.7,0.95,"","brNDC")
    ratio2_leg.SetBorderSize(0);
    ratio2_leg.SetTextSize(0.12);
    ratio2_leg.SetFillColor(0);
    ratio2_leg.SetLineColor(1);
    ratio2_leg.SetTextFont(42);


    for key in sorted(ratio_to_draw.keys()):

        ratio_to_draw[key].SetLineColor(kBlack)
        ratio_to_draw[key].SetMarkerColor(kBlack)
        ratio_to_draw[key].SetMarkerStyle(20)


        ratio_to_draw[key].GetXaxis().SetTitle(xtitle)
        ratio_to_draw[key].GetXaxis().SetLabelSize(0.15)
        ratio_to_draw[key].GetXaxis().SetTitleSize(0.18)
        ratio_to_draw[key].GetXaxis().SetTitleOffset(1.2)
        ratio_to_draw[key].GetXaxis().SetRangeUser(xmin,xmax)
        ratio_to_draw[key].GetYaxis().SetTitle(ratiotitle)
        ratio_to_draw[key].GetYaxis().SetLabelSize(0.15)
        ratio_to_draw[key].GetYaxis().SetTitleSize(0.15)
        ratio_to_draw[key].GetYaxis().SetTitleOffset(0.45)
        ratio_to_draw[key].GetYaxis().SetNdivisions(503)
        ratio_to_draw[key].GetYaxis().SetRangeUser(ratioymin,ratioymax)
        ratio_to_draw[key].GetYaxis().SetNdivisions(ndivisions)

        ratio_to_draw[key].Draw(draw_opt+" same")
        i+=2


    if bratioleg: ratio_leg.Draw()
    if bratioleg: ratio2_leg.Draw()


    line = TLine (xmin,1,xmax,1)
    line.SetLineColor(kRed)
    line.SetLineStyle(kDashed)
    line.Draw("same")

    gPad.Update();
    gPad.RedrawAxis();

    pad_top.cd()

    leg.Draw()
    gPad.Update();
    gPad.RedrawAxis();

    if extratext not in "":
        text = CMSPlotStyle.draw_info(extratext, extrax,extray, factor = 1.5)

    if extratext2 not in "":
        text2 = CMSPlotStyle.draw_info(extratext2, extrax2,extray2, factor = 1.5)

    if blogo: cmslogo = CMSPlotStyle.draw_cmstext("left", False,2)


    lumi = CMSPlotStyle.draw_lumi(False,0.4,lumitext)
    lumi.Draw()

    if b_anti:
        anti = CMSPlotStyle.draw_info("Anti-k_{T}, R = 0.4",antix,antiy,factor=1.5)
        anti.Draw()

    if blogy: pad_top.SetLogy()
    if beps: c.Print(folder+filename+".eps")
    c.Print(folder+filename+".pdf")

    CMSPlotStyle.extratext = "Preliminary"
    cmslogo = CMSPlotStyle.draw_cmstext("left", True,2)
    c.Print(folder+filename+"_pre.pdf")
