#include "TH1.h"
#include "TSystem.h"
#include <vector>

#include "Style.C"
//------------------------------------------------------------------------------

void jet_eta_pt()
{
  //SetStyle();
  SetupGlobalStyle();
  TFile *f[3];
  TFile *f_QCD[3];
  // read from results.root
  f[0] = new TFile("PLOTS/study_timing_cut/VBF/macro_jets_snowmass_EtaMax0_correct/results.root");
  f[1] = new TFile("PLOTS/study_timing_cut/VBF/macro_jets_snowmass_EtaMax3_correct/results.root");
  f[2] = new TFile("PLOTS/study_timing_cut/VBF/macro_jets_snowmass_EtaMax4_correct/results.root");

  f_QCD[0] = new TFile("PLOTS/study_timing_cut/QCD_100kPU/macro_QCD_jets_EtaMax0/results.root");
  f_QCD[1] = new TFile("PLOTS/study_timing_cut/QCD_100kPU/macro_QCD_jets_EtaMax3/results.root");
  f_QCD[2] = new TFile("PLOTS/study_timing_cut/QCD_100kPU/macro_QCD_jets_EtaMax4/results.root");

  TCanvas* canvas = new TCanvas("c", "c", 600, 600);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);

  TH1 *hist_eta[3][3][2]; // eta max, CHS or PUPPI, eta bins

  TString category[3] = {"VBF", "b", "PU"};
  TString name[3] = {"CHS", "PUPPI", "GenJet"};

  // get the histograms
  for (size_t i = 0; i < 2; i++) { // CHS or PUPPI
      for (size_t k = 0; k < 3; k++) {
        if (k<2) {
          hist_eta[0][i][k] = ((TH1F*)f[0]->Get(category[k]+"_"+name[i]+"_jet_eta"));
          hist_eta[1][i][k] = ((TH1F*)f[1]->Get(category[k]+"_"+name[i]+"_jet_eta"));
          hist_eta[2][i][k] = ((TH1F*)f[2]->Get(category[k]+"_"+name[i]+"_jet_eta"));
        }

        if (k==2) {
          hist_eta[0][i][k] = ((TH1F*)f_QCD[0]->Get(category[k]+"_"+name[i]+"_jet_eta"));
          hist_eta[1][i][k] = ((TH1F*)f_QCD[1]->Get(category[k]+"_"+name[i]+"_jet_eta"));
          hist_eta[2][i][k] = ((TH1F*)f_QCD[2]->Get(category[k]+"_"+name[i]+"_jet_eta"));
        }

        // rebin eta plots
        hist_eta[0][i][k]->Rebin(2);
        hist_eta[1][i][k]->Rebin(2);
        hist_eta[2][i][k]->Rebin(2);

      hist_eta[0][i][k]->SetLineColor(1);
      hist_eta[1][i][k]->SetLineColor(30);
      hist_eta[2][i][k]->SetLineColor(46);

    // create a multigraph and draw it
    THStack  *mg  = new THStack();
    mg->Add(hist_eta[0][i][k]);
    mg->Add(hist_eta[1][i][k]);
    mg->Add(hist_eta[2][i][k]);
    mg->Draw("nostack");

    mg->GetXaxis()->SetTitle("#eta");
    mg->GetYaxis()->SetTitle("Number of jets");

    Double_t ymax = hist_eta[0][i][k]->GetMaximum() > hist_eta[1][i][k]->GetMaximum() ? hist_eta[0][i][k]->GetMaximum() : hist_eta[1][i][k]->GetMaximum();
    Double_t ymin = hist_eta[0][i][k]->GetMinimum() < hist_eta[1][i][k]->GetMinimum() ? hist_eta[0][i][k]->GetMinimum() : hist_eta[1][i][k]->GetMinimum();
  //  mg->GetYaxis()->SetRangeUser(0,ymax *2);
    mg->SetMaximum(ymax * 1.8);
  //  mg->SetMinimum(0.);
    //mg->GetXaxis()->SetRangeUser(0,1);


  TPaveStats *stats;
  TPaveStats *stats2;
  TPaveStats *stats3;

  canvas->Update();
  stats = (TPaveStats*)hist_eta[0][i][k]->GetListOfFunctions()->FindObject("stats");
  stats->SetName("h1stats");
  stats->SetY1NDC(.8);
  stats->SetY2NDC(.9);
  stats->SetX1NDC(.7);
  stats->SetX2NDC(.9);
  //stats->SetTextColor(2);
  stats->SetTextSize(0.03);
  stats->SetBorderSize(0);
  stats->SetOptStat(1100);

  canvas->Update();
  stats2 = (TPaveStats*)hist_eta[1][i][k]->GetListOfFunctions()->FindObject("stats");
   stats2->SetName("h1stats2");
   stats2->SetY1NDC(.7);
   stats2->SetY2NDC(.8);
   stats2->SetX1NDC(.7);
   stats2->SetX2NDC(.9);
   stats2->SetTextColor(30);
   stats2->SetTextSize(0.03);
   stats2->SetBorderSize(0);
   stats2->SetOptStat(1);

   canvas->Update();
   stats3 = (TPaveStats*)hist_eta[2][i][k]->GetListOfFunctions()->FindObject("stats");
    stats3->SetName("h1stats3");
    stats3->SetY1NDC(.6);
    stats3->SetY2NDC(.7);
    stats3->SetX1NDC(.7);
    stats3->SetX2NDC(.9);
    stats3->SetTextColor(46);
    stats3->SetTextSize(0.03);
    stats3->SetBorderSize(0);
    stats3->SetOptStat(1110);

    canvas->Modified();
    mg->Draw("nostack");

    relPosX    = 0.045;
    CMS_lumi(canvas, 0, 10);

    // // legend
    //auto legend = new TLegend(0.25,0.7,0.85,0.9); // x1, y1,x2, y2
    auto legend = new TLegend(0.25,0.65,0.85,0.85); // x1, y1,x2, y2
    legend->SetTextSize(.03);
    legend->SetFillStyle(0);
    legend->SetHeader("AK4 "+ name[i]+", "+category[k]+" jets , p_{T}>20 GeV");
    legend->AddEntry(hist_eta[0][i][k],"no timing", "l");
    legend->AddEntry(hist_eta[1][i][k],"timing for |#eta|<3", "l");
    legend->AddEntry(hist_eta[2][i][k],"timing for |#eta|<4", "l");

    legend->Draw();

    // // save output
    TString outdir = "PLOTS/study_timing_cut/eta_comparison/";
    canvas->SaveAs(outdir  + category[k]+"_"+name[i]+"_eta_100kPU_correctSignal.pdf");
    //canvas->SaveAs(outdir  + category[k]+"_"+name[i] +"_eta.root");
    }// end categories
  } // end jet loop


} // end void
