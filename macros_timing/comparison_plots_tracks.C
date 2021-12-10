#include "TH1.h"
#include "TSystem.h"
#include <vector>

#include "Style.C"
//------------------------------------------------------------------------------

void comparison_plots_tracks()
{
  //SetStyle();
  SetupGlobalStyle();
  // read from results.root
  TFile *f = new TFile("PLOTS/10kevents/track_control_plots/results.root");

  TCanvas* canvas = new TCanvas("c", "c", 600, 600);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);
  TH1F* track_signal[6];
  TH1F* track_PU[6];

  TPaveStats *stats;
  TPaveStats *stats2;

  TString name[6] = {"eflowTrack", "RecoPUTrack", "track_merger", "smeared_track", "CHSeflow", "PUPPIeflow"};
  TString IsPU[3] = {"_signal", "_PU", "_inclusive"};
  // get the histograms
  for (size_t i = 0; i < 6; i++) {
    track_signal[i] = ((TH1F*)f->Get("dt_"+name[i]+IsPU[0]));
    track_PU[i] = ((TH1F*)f->Get("dt_"+name[i]+IsPU[1]));

    track_signal[i]->SetLineColor(kRed);
    //track_signal[i]->Scale(1/track_signal[i]->Integral());
    track_PU[i]->SetLineColor(kBlue);
    //track_PU[i]->Scale(1/track_PU[i]->Integral());

  //track_signal[i]->Draw("H");
  //track_PU[i]->Draw("H same");

  // set stat boxes
  canvas->Update();
  stats = (TPaveStats*)track_signal[i]->GetListOfFunctions()->FindObject("stats");
  stats->SetName("h1stats");
  stats->SetY1NDC(.75);
  stats->SetY2NDC(.85);
  stats->SetX1NDC(.65);
  stats->SetX2NDC(.85);
  stats->SetTextColor(2);
  stats->SetTextSize(0.03);
  stats->SetBorderSize(0);

  canvas->Update();
  stats2 = (TPaveStats*)track_PU[i]->GetListOfFunctions()->FindObject("stats");
   stats2->SetName("h1stats2");
   stats2->SetY1NDC(.65);
   stats2->SetY2NDC(.75);
   stats2->SetX1NDC(.65);
   stats2->SetX2NDC(.85);
   stats2->SetTextColor(4);
   stats2->SetTextSize(0.03);
   stats2->SetBorderSize(0);

   THStack *hs = new THStack("hs","Stacked 1D histograms");

   hs->Add(track_signal[i]);
   hs->Add(track_PU[i]);
   hs->Draw("nostack");

   hs->GetXaxis()->SetTitle("Track dt [ns]");
   hs->GetYaxis()->SetTitleOffset(2);
   hs->GetYaxis()->SetTitle("number of tracks");

  // legend
  auto legend = new TLegend(0.25, 0.7, 0.55, 0.85);
  legend->SetTextSize(.04);
  legend->SetFillStyle(0);
  //legend->SetHeader(xx);
  legend->AddEntry(track_signal[0],"signal", "l");
  legend->AddEntry(track_PU[0],"PU", "l");

  //legend->Draw();
  // save output
  TString outdir = "PLOTS/10kevents/track_control_plots/comparison/";
  canvas->SaveAs(outdir +"dt_" + name[i]+ ".eps");
  canvas->SaveAs(outdir +"dt_" + name[i]+ ".root");

  }
}
