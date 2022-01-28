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
  TFile *f = new TFile("PLOTS/10kevents/track_vertex_association/scenario_b/results.root");

  TCanvas* canvas = new TCanvas("c", "c", 600, 600);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);
  TH1F* track_signal[2];
  TH1F* track_PU[2];

  TPaveStats *stats;
  TPaveStats *stats2;

  //TString name[6] = {"eflowTrack", "RecoPUTrack", "track_merger", "smeared_track", "CHSeflow", "PUPPIeflow"};
  TString IsPU[3] = {"_signal", "_PU", "_inclusive"};
  TString variable[2] = {"track_dz_smeared_all", "track_dt_smeared_all"};
  TString variablenames[2] = { "track dZ [m]", "track dT [ns]"};

  // get the histograms
  //for (size_t i = 0; i < 6; i++) {
    for (size_t l = 0; l < 2; l++) {
    track_signal[l] = ((TH1F*)f->Get(variable[l]+IsPU[0]));
    track_PU[l] = ((TH1F*)f->Get(variable[l]+IsPU[1]));

    track_signal[l]->SetLineColor(kRed);
  //  track_signal[l]->Scale(1/track_signal[l]->Integral());
    track_PU[l]->SetLineColor(kBlue);
    //track_PU[l]->Scale((1/track_PU[l]->Integral()));

  track_signal[l]->Draw("H");
  track_PU[l]->Draw("H same");

  // set stat boxes
  canvas->Update();
  stats = (TPaveStats*)track_signal[l]->GetListOfFunctions()->FindObject("stats");
  stats->SetName("h1stats");
  stats->SetY1NDC(.75);
  stats->SetY2NDC(.85);
  stats->SetX1NDC(.6);
  stats->SetX2NDC(.85);
  stats->SetTextColor(2);
  stats->SetTextSize(0.03);
  stats->SetBorderSize(0);

  canvas->Update();
  stats2 = (TPaveStats*)track_PU[l]->GetListOfFunctions()->FindObject("stats");
   stats2->SetName("h1stats2");
   stats2->SetY1NDC(.65);
   stats2->SetY2NDC(.75);
   stats2->SetX1NDC(.6);
   stats2->SetX2NDC(.85);
   stats2->SetTextColor(4);
   stats2->SetTextSize(0.03);
   stats2->SetBorderSize(0);

   Double_t ymax = track_signal[l]->GetMaximum() > track_PU[l]->GetMaximum() ? track_signal[l]->GetMaximum() : track_PU[l]->GetMaximum();
   Double_t ymin = 0;
   track_signal[l]->GetYaxis()->SetRangeUser(ymin,ymax *1.8);

   if(l==0) track_signal[l]->GetXaxis()->SetNdivisions(506);
   //if(l==0) track_signal[l]->GetXaxis()->SetRangeUser(-0.1,0.1);
   if(l==1) track_signal[l]->GetXaxis()->SetRangeUser(-0.5,0.5);

   track_signal[l]->GetXaxis()->SetTitle(variablenames[l]);
   track_signal[l]->GetYaxis()->SetTitleOffset(2);
   track_signal[l]->GetYaxis()->SetTitle("number of tracks");

  // legend
  auto legend = new TLegend(0.25, 0.7, 0.55, 0.85);
  legend->SetTextSize(.04);
  legend->SetFillStyle(0);
  //legend->SetHeader(xx);
  legend->AddEntry(track_signal[l],"signal", "l");
  legend->AddEntry(track_PU[l],"PU", "l");

  //legend->Draw();
  // save output
  TString outdir = "PLOTS/10kevents/track_vertex_association/";
  canvas->SaveAs(outdir +variable[l]+ ".pdf");
  canvas->SaveAs(outdir +variable[l]+ ".root");

} // end variables
//} // end track branches
}
