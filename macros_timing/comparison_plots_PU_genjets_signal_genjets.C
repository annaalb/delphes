#include "TH1.h"
#include "TSystem.h"
#include <vector>

#include "Style.C"
//------------------------------------------------------------------------------

void comparison_plots_PU_genjets_signal_genjets()
{
  //SetStyle();
  SetupGlobalStyle();
  // read from results.root
  TFile *f = new TFile("PLOTS/10kevents/macro_jets_snowmass_signalTrackRejection/results.root"); // change here

  TCanvas* canvas = new TCanvas("c", "c", 600, 600);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);
  TH1F* track_VBF[4][2];
  TH1F* track_PU[4][2];

  TPaveStats *stats;
  TPaveStats *stats1;
  TPaveStats *stats2;

  TString name[4] = {"CHS", "PUPPI", "GenJet", "PUGenJet"};
  TString category[3] = {"VBF", "b", "PU"};
  TString IsPU[3] = {"_signal", "_PU", ""};
  TString eta[4] = {"", "_eta02", "_eta23", "_eta34"};

  TString variable[4] = {"_jet_time_minus_PV", "_jet_pt_weighted_time_minus_PV", "_track_DZ", "_track_T"};


  Int_t k = 2;
  Int_t i = 0;

  // get the histograms
  //for (size_t i = 0; i < 3; i++) { // CHS or PUPPI, genjet or PUgenjet // change here
    //for (size_t k = 0; k < 3; k++) {
      for (size_t l = 0; l < 4; l++) {  // variables
        for (size_t e = 0; e < 4; e++) { // eta bins

    track_VBF[i][k] = ((TH1F*)f->Get(category[0]+"_"+name[2]+IsPU[k]+variable[l]+eta[e]));
    track_PU[i][k] = ((TH1F*)f->Get(category[2]+"_"+name[3]+IsPU[k]+variable[l]+eta[e]));

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

    // set stat boxes
    canvas->Update();
    stats1 = (TPaveStats*)track_PU[i][k]->GetListOfFunctions()->FindObject("stats");
    stats1->SetName("h1stats1");
    stats1->SetY1NDC(.65);
    stats1->SetY2NDC(.75);
    stats1->SetX1NDC(.6);
    stats1->SetX2NDC(.85);
    stats1->SetTextColor(4);
    stats1->SetTextSize(0.02);
    stats1->SetBorderSize(0);

    // PU jet
    track_PU[i][k]->SetLineColor(kBlue);
    track_PU[i][k]->Scale(1/track_PU[i][k]->Integral());
    track_PU[i][k]->Draw(" HIST same");

   stats1->Draw();
   stats->Draw();

   canvas->Update();
   canvas->Modified();

   Double_t ymax = track_VBF[i][k]->GetMaximum() > track_PU[i][k]->GetMaximum() ? track_VBF[i][k]->GetMaximum() : track_PU[i][k]->GetMaximum();
   Double_t ymin = 0;
   track_VBF[i][k]->GetYaxis()->SetRangeUser(ymin,ymax *1.8);

   track_VBF[i][k]->GetXaxis()->SetTitle(variable[l]+eta[e]+" [ns]");
   //track_VBF[i][k]->GetXaxis()->SetTitle("Track dz [m]");
  // track_VBF[i][k]->GetXaxis()->SetNdivisions(504);
   //track_VBF[i][k]->GetXaxis()->SetRangeUser(0,20);
   track_VBF[i][k]->GetYaxis()->SetTitleOffset(2);
   track_VBF[i][k]->GetYaxis()->SetTitle("fraction of jets");

  // legend
  // auto legend = new TLegend(0.25, 0.7, 0.55, 0.85);
  // legend->SetTextSize(.04);
  // legend->SetFillStyle(0);
  // //legend->SetHeader(xx);
  // legend->AddEntry(track_VBF[0][k],"signal", "l");
  // legend->AddEntry(track_PU[0][k],"PU", "l");

  //legend->Draw();
  // save output
  TString outdir = "PLOTS/10kevents/comparison_plots/"; // change here
  canvas->SaveAs(outdir  +"genjets"+variable[l]+eta[e]+".pdf");
  canvas->SaveAs(outdir  +"genjets"+variable[l]+eta[e]+ ".root");
} // end eta bins
  }
//  }
//}
}
