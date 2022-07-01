#include "TH1.h"
#include "TSystem.h"
#include <vector>

#include "Style.C"
//------------------------------------------------------------------------------

void dz_comparison()
{
  SetupGlobalStyle();
  // read from results.root
  TFile *f[3];

  f[0] = new TFile("PLOTS/study_timing_cut/VBF/Debug_signalRejection/results.root");
  f[1] = new TFile("PLOTS/study_timing_cut/VBF/Debug_nosignalRejection/results.root");

  TCanvas* canvas = new TCanvas("c", "c", 600, 600);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);
  TH1F* track_signal[2][6][4];
  int integral[2][6][4];

  TPaveStats *stats;
  TPaveStats *stats2;

  TString name[2] = {"Signal_rej", "No_signal_rej"};
  TString IsPU[3] = {"_signal", "_PU", "_inclusive"};
  TString etabin[6] = {"", "_eta01p3", "_eta1p32", "_eta23", "_eta34", "_eta4"};
  TString ptbin[4] = {"_pt0to1","_pt1to2","_pt2to10","_pt10toInf"};

  for (size_t k = 0; k < 6; k++) { // eta bins
    for (size_t p = 0; p < 4; p++) { // pt bins
      for (size_t l = 0; l < 2; l++) {
        track_signal[l][k][p] = ((TH1F*)f[l]->Get(name[l]+"_dz"+IsPU[0]+etabin[k]+ptbin[p]));
        cout << name[l]+"_dz"+IsPU[0]+etabin[k]+"_"+ptbin[p] << endl;
        integral[l][k][p] = track_signal[l][k][p]->Integral(26,75) / track_signal[l][k][p]->Integral() * 100;
      }

  track_signal[0][k][p]->SetLineStyle(33);

  track_signal[1][k][p]->SetLineColor(30);
  track_signal[1][k][p]->SetLineStyle(22);

  TPaveStats *stats;
  TPaveStats *stats2;

  // create a multigraph and draw it
  THStack  *mg  = new THStack();
  mg->Add(track_signal[0][k][p]);
  mg->Add(track_signal[1][k][p]);
  mg->Draw("nostack");
  Double_t ymax = track_signal[0][k][p]->GetMaximum() > track_signal[1][k][p]->GetMaximum() ? track_signal[0][k][p]->GetMaximum() : track_signal[1][k][p]->GetMaximum();
  Double_t ymin = 0;
  //mg->GetYaxis()->SetRangeUser(ymin,ymax *1.8);
  mg->SetMaximum(ymax *1.8);
  mg->GetXaxis()->SetNdivisions(506);

  mg->GetXaxis()->SetTitle("Track dz [m]");
  mg->GetYaxis()->SetTitle("Number of tracks");
  canvas->Modified();
  mg->Draw("nostack");

  // set stat boxes
  canvas->Update();
  stats = (TPaveStats*)track_signal[0][k][p]->GetListOfFunctions()->FindObject("stats");
  stats->SetName("h1stats");
  stats->SetY1NDC(.75);
  stats->SetY2NDC(.85);
  stats->SetX1NDC(.6);
  stats->SetX2NDC(.85);
  stats->SetTextSize(0.03);
  stats->SetBorderSize(0);

  canvas->Update();
  stats2 = (TPaveStats*)track_signal[1][k][p]->GetListOfFunctions()->FindObject("stats");
   stats2->SetName("h1stats2");
   stats2->SetY1NDC(.65);
   stats2->SetY2NDC(.75);
   stats2->SetX1NDC(.6);
   stats2->SetX2NDC(.85);
   stats2->SetTextColor(30);
   stats2->SetTextSize(0.03);
   stats2->SetBorderSize(0);

  // legend
  auto legend = new TLegend(0.25, 0.7, 0.55, 0.85);
  legend->SetTextSize(.04);
  legend->SetFillStyle(0);
  //legend->SetHeader(xx);
  legend->AddEntry(track_signal[0][k][p],"Signal rej. ", "l");
  legend->AddEntry(track_signal[1][k][p],"No signal rej. ", "l");

  TPaveText *comment = new TPaveText(0.25, 0.5, 0.55, 0.7, "brNDC");

  comment->SetTextSize(.04);
//  comment->SetTextFont(kExRootFont);
  comment->SetTextAlign(12);
  comment->SetFillStyle(0);
  comment->SetBorderSize(0);
  const char* str = to_string(integral[0][k][p]).c_str();
  const char* str1 = to_string(integral[1][k][p]).c_str();
  comment->AddText("Integral (-dz, dz)");
  comment->AddText(str);
  comment->AddText(str1); ((TText*)comment->GetListOfLines()->Last())->SetTextColor(30);

  legend->Draw();
  comment->Draw();
  // save output
  TString outdir = "PLOTS/study_timing_cut/VBF/Debug_signalRejection/";
  canvas->SaveAs(outdir + "dz"+etabin[k]+ptbin[p]+".pdf");
  //canvas->SaveAs(outdir +variable[l]+ ".root");

  }
  }
}
