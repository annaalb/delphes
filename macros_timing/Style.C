#include <iostream>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TFile.h>
#include <TLegend.h>
#include <TArray.h>
#include <TPad.h>
#include <TGraph.h>
#include <TGraphErrors.h>

using namespace std;

void SetStyle(){
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetLegendBorderSize(0);
}

void SetupGlobalStyle()
{
  // general appearance and style

  gROOT->SetStyle("Plain");
  //gStyle->SetOptStat(1100);
  //gStyle->SetOptStat(false);
  gStyle -> SetPadTickX(1);
  gStyle -> SetPadTickY(1);

  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(kWhite);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gStyle->SetGridColor(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);

  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(2);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameFillStyle(0);
  gStyle->SetFrameLineColor(1);
  gStyle->SetFrameLineStyle(1);
  gStyle->SetFrameLineWidth(2);

  gStyle->SetOptTitle(0);
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");

  gStyle->SetAxisColor(1, "XYZ");
  gStyle->SetStripDecimals(kTRUE);
  gStyle->SetTickLength(0.05, "XYZ");
  gStyle->SetNdivisions(510, "XYZ");

  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetCanvasDefH(700);
  gStyle->SetCanvasDefW(800);

  // Margins:
  // gStyle->SetPadTopMargin(    0.05  );
  // gStyle->SetPadBottomMargin( 0.16  );
  // gStyle->SetPadRightMargin(  0.05  );
  // gStyle->SetPadLeftMargin(   0.19  );

  gStyle->SetPadTopMargin(0.10);
  gStyle->SetPadRightMargin(0.10);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);

  gStyle->SetLegendBorderSize(0);

  gStyle->UseCurrentStyle();

}

TCanvas* SetupCanvas(TString name = "canvas")
{
  // optimised plots for including in theses or publications and documents

  Int_t CanWidth;
  Int_t CanHeight;
  CanWidth = 500;
  CanHeight = 400;

  // set up the canvas
  TCanvas* can = new TCanvas(name,"Control Plots", CanWidth, CanHeight);

  can->SetTopMargin(    0.05  );
  can->SetBottomMargin( 0.16  );
  can->SetRightMargin(  0.05  );
  can->SetLeftMargin(   0.16  );

  return can;
}

TCanvas* SetupCanvas2d(TString name = "canvas")
{
  // optimised plots for including in theses or publications and documents

  Int_t CanWidth;
  Int_t CanHeight;
  CanWidth = 500;
  CanHeight = 500;

  // set up the canvas
  TCanvas* can = new TCanvas(name,"Control Plots", CanWidth, CanHeight);

  can->SetTopMargin(    0.05  );
  can->SetBottomMargin( 0.14  );
  can->SetRightMargin(  0.05  );
  can->SetLeftMargin(   0.14  );

  return can;
}

void Cosmetics(TH1* hist)
{

    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetXaxis()->SetLabelOffset(0.008);
    hist->GetXaxis()->SetTickLength(0.04);
    hist->GetXaxis()->SetTitleSize(0.055);
    hist->GetXaxis()->SetTitleOffset(1.1);
    hist->GetXaxis()->SetNdivisions(510);

    hist->GetYaxis()->SetTitleOffset(1.3);
    hist->GetYaxis()->SetTitleSize(0.055);
    hist->GetYaxis()->SetLabelSize(0.045);
    hist->GetYaxis()->SetTickLength(0.04);
    hist->GetYaxis()->SetLabelOffset(0.011);
    hist->GetYaxis()->SetNdivisions(505);

}

void Cosmetics2d(TH2* hist)
{

    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetXaxis()->SetLabelOffset(0.008);
    hist->GetXaxis()->SetTickLength(0.03);
    hist->GetXaxis()->SetTitleSize(0.055);
    hist->GetXaxis()->SetTitleOffset(1.1);
    hist->GetXaxis()->SetNdivisions(510);

    hist->GetYaxis()->SetTitleOffset(1.1);
    hist->GetYaxis()->SetTitleSize(0.055);
    hist->GetYaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetTickLength(0.03);
    hist->GetYaxis()->SetLabelOffset(0.011);
    hist->GetYaxis()->SetNdivisions(510);

}

TLegend* SetupLegend(double shiftleft=0.)
{
  Double_t xl1=.21+shiftleft, yl1=0.65, xl2=xl1+0.3, yl2=yl1+.15;
  TLegend* leg = new TLegend(xl1,yl1,xl2,yl2, NULL,"brNDC");
  leg->SetFillColor(0);
  leg->SetLineColor(1);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.055);
  return leg;
}

TLegend* SetupLegend2(double shiftleft=0.)
{
  Double_t xl1=.21+shiftleft, yl1=0.62, xl2=xl1+0.3, yl2=yl1+.19;
  TLegend* leg = new TLegend(xl1,yl1,xl2,yl2, NULL,"brNDC");
  leg->SetFillColor(0);
  leg->SetLineColor(1);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.045);
  leg->SetEntrySeparation(0.15);
  return leg;
}


TLegend* SetupEffLegend()
{
  Double_t xl1=.21, yl1=0.75, xl2=xl1+0.55, yl2=yl1+.16;
  TLegend* leg = new TLegend(xl1,yl1,xl2,yl2, NULL,"brNDC");
  leg-> SetNColumns(2);
  leg->SetFillColor(0);
  leg->SetLineColor(1);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.045);
  return leg;
}

void Simulation_label()
{
  Int_t c_LGray = TColor::GetColor( "#e6e6e6" );

  TLatex *text1 = new TLatex(3.5, 24, "Delphes Simulation");
  text1->SetNDC();
  text1->SetTextAlign(13);
  text1->SetX(0.72);
  text1->SetTextFont(62);
  text1->SetTextSize(0.075);
  text1->SetY(0.892);
  text1->SetTextColor(c_LGray);
  text1->Draw();
}
