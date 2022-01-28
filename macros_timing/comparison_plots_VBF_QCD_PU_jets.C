#include "TH1.h"
#include "TSystem.h"
#include <vector>

#include "Style.C"
//------------------------------------------------------------------------------

void comparison_plots_VBF_QCD_PU_jets()
{
  //SetStyle();
  SetupGlobalStyle();
  // read from results.root
  TFile *f = new TFile("PLOTS/10kevents/macro_jets_snowmass_21_01/results.root");
  //TFile *f_QCD = new TFile("PLOTS/10kevents/QCD/macro_jets_snowmass/results.root");
  TFile *f_QCD = new TFile("PLOTS/10kevents/QCD/macro_QCD_jets_12_01/results.root");

  TCanvas* canvas = new TCanvas("c", "c", 600, 600);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);
  TH1F* track_VBF[3][2];
  TH1F* track_PU[3][2];

  TPaveStats *stats;
  TPaveStats *stats1;
  TPaveStats *stats2;

  TString name[3] = {"CHS", "PUPPI", "GenJet"};
  TString category[3] = {"VBF", "b", "PU"};
  TString IsPU[3] = {"", "_signal", "_PU"};
  TString ptbin[5] = {"", "_pt20", "_20pt30", "_30pt50", "_pt50"};
  TString eta[6] = {"", "_eta01p3", "_eta1p32", "_eta23", "_eta34", "_eta4"};

  TString variable[2] = {"_jet_pt", "_jet_eta"};
  TString variablenames[2] = {"_jet_pt", "_jet_eta"};

  // TString variable[7] = {"_jet_time_minus_PV", "_jet_pt_weighted_time_minus_PV", "_track_DZ", "_track_T", "_jet_n", "_jet_n_PU_tracks", "_jet_n_signal_tracks"};
  // TString variablenames[7] = {"jet_time_minus_PV [ns]", "jet_pt_weighted_time_minus_PV [ns]", "track dZ [m]", "track dT [ns]", "Number of PU jets", "Number of PU tracks", "Number of signal tracks"};

  Int_t k = 0;
  Int_t e = 0;
  // get the histograms
  for (size_t i = 0; i < 2; i++) { // CHS or PUPPI
    for (size_t l = 0; l < 2; l++) {  // variables
      //for (size_t l = 0; l < 7; l++) {  // variables
        //for (size_t e = 0; e < 6; e++) { // eta bins
          //if (l>2) { for (size_t k = 0; k < 3; k++) { // split tracks into signal / PU
    track_VBF[i][k] = ((TH1F*)f->Get(category[2]+"_"+name[i]+IsPU[k]+variable[l]+eta[e]));
    track_PU[i][k] = ((TH1F*)f_QCD->Get(category[2]+"_"+name[i]+IsPU[k]+variable[l]+eta[e]));

    // VBF jets
    track_VBF[i][k]->SetLineColor(kRed);
    track_VBF[i][k]->Scale(1/track_VBF[i][k]->Integral());
    track_VBF[i][k]->Draw("HIST");

    // set stat boxes
    canvas->Update();
    stats = (TPaveStats*)track_VBF[i][k]->GetListOfFunctions()->FindObject("stats");
    stats->SetName("h1stats");
    stats->SetY1NDC(.75);
    stats->SetY2NDC(.85);
    stats->SetX1NDC(.6);
    stats->SetX2NDC(.85);
    stats->SetTextColor(2);
    stats->SetTextSize(0.02);
    stats->SetBorderSize(0);

    // PU jet
    track_PU[i][k]->SetLineColor(kBlack);
    track_PU[i][k]->Scale(1/track_PU[i][k]->Integral());
    track_PU[i][k]->Draw(" HIST same");

    // set stat boxes
    canvas->Update();
    stats1 = (TPaveStats*)track_PU[i][k]->GetListOfFunctions()->FindObject("stats");
    stats1->SetName("h1stats1");
    stats1->SetY1NDC(.65);
    stats1->SetY2NDC(.75);
    stats1->SetX1NDC(.6);
    stats1->SetX2NDC(.85);
    stats1->SetTextColor(1);
    stats1->SetTextSize(0.02);
    stats1->SetBorderSize(0);


   stats1->Draw();
   stats->Draw();

   canvas->Update();
   canvas->Modified();

   Double_t ymax = track_VBF[i][k]->GetMaximum() > track_PU[i][k]->GetMaximum() ? track_VBF[i][k]->GetMaximum() : track_PU[i][k]->GetMaximum();
   Double_t ymin = 0;
   track_VBF[i][k]->GetYaxis()->SetRangeUser(ymin,ymax *1.8);

   track_VBF[i][k]->GetXaxis()->SetTitle(variablenames[l]+eta[e]);

   //if(l==2) track_VBF[i][k]->GetXaxis()->SetNdivisions(504);
  // if(l==0 || l==1) track_VBF[i][k]->GetXaxis()->SetRangeUser(-0.4,0.4);
   track_VBF[i][k]->GetYaxis()->SetTitleOffset(2);
   track_VBF[i][k]->GetYaxis()->SetTitle("fraction of jets");

   auto legend = new TLegend(0.25, 0.7, 0.55, 0.85);
   legend->SetTextSize(.03);
   legend->SetFillStyle(0);
   legend->AddEntry(track_VBF[0][k],"VBF sample", "l");
   legend->AddEntry(track_PU[0][k],"QCD sample", "l");

   legend->Draw();

  // save output
  TString outdir = "PLOTS/10kevents/comparison_VBF_QCD_PU_jets/";
  canvas->SaveAs(outdir  + name[i]+IsPU[k]+variable[l]+eta[e]+".pdf");
  canvas->SaveAs(outdir  + name[i]+IsPU[k]+variable[l]+eta[e]+ ".root");

//} // end if track variable
//} // end signal / PU (for tracks)
//}// end eta bins
} // end variables
} // end jet loop
} // end void
