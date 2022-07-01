#include "TH1.h"
#include "TSystem.h"
#include <vector>

#include "Style.C"
//------------------------------------------------------------------------------

void create_MTD_plot()
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

  TGraphAsymmErrors *jet_rate[3][3][2];

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
        hist_eta[0][i][k]->Rebin(10);
        hist_eta[1][i][k]->Rebin(10);
        hist_eta[2][i][k]->Rebin(10);

        jet_rate[0][i][k] = new TGraphAsymmErrors(hist_eta[0][i][k], hist_eta[0][i][k], "pois");
        jet_rate[1][i][k] = new TGraphAsymmErrors(hist_eta[1][i][k], hist_eta[0][i][k], "pois");
        jet_rate[2][i][k] = new TGraphAsymmErrors(hist_eta[2][i][k], hist_eta[0][i][k], "pois");

      jet_rate[0][i][k]->SetMarkerStyle(33);
      jet_rate[0][i][k]->SetLineColor(1);

      jet_rate[1][i][k]->SetLineColor(30);
      jet_rate[1][i][k]->SetMarkerColor(30);
      jet_rate[1][i][k]->SetMarkerStyle(22);

      jet_rate[2][i][k]->SetLineColor(46);
      jet_rate[2][i][k]->SetMarkerColor(46);
      jet_rate[2][i][k]->SetMarkerStyle(20);

    // create a multigraph and draw it
    TMultiGraph  *mg  = new TMultiGraph();
    mg->Add(jet_rate[0][i][k]);
    mg->Add(jet_rate[1][i][k]);
    mg->Add(jet_rate[2][i][k]);
    mg->GetXaxis()->SetTitle("#eta");
    mg->GetYaxis()->SetTitle("Jet rate");

    Double_t ymax = jet_rate[0][i][k]->GetMaximum() > jet_rate[1][i][k]->GetMaximum() ? jet_rate[0][i][k]->GetMaximum() : jet_rate[1][i][k]->GetMaximum();
    Double_t ymin = jet_rate[0][i][k]->GetMinimum() < jet_rate[1][i][k]->GetMinimum() ? jet_rate[0][i][k]->GetMinimum() : jet_rate[1][i][k]->GetMinimum();
  //  mg->GetYaxis()->SetRangeUser(0,ymax *2);
    mg->SetMaximum(1.8);
    mg->SetMinimum(0.);

    canvas->Modified();
    mg->Draw("AP");

    // CMS in the upper left corner (outside the plot)
    CMS_lumi(canvas);

    // // legend
  //  auto legend = new TLegend(0.25,0.25,0.85,0.45); // x1, x2, y1, y2
    auto legend = new TLegend(0.25,0.7,0.85,0.9); // x1, y1,x2, y2
    legend->SetTextSize(.03);
    legend->SetFillStyle(0);
    legend->SetHeader("AK4 "+ name[i]+", "+category[k]+" jets, p_{T}>20 GeV");
    legend->AddEntry(jet_rate[0][i][k],"no timing", "l");
    legend->AddEntry(jet_rate[1][i][k],"timing for |#eta|<3", "l");
    legend->AddEntry(jet_rate[2][i][k],"timing for |#eta|<4", "l");

    legend->Draw();

    // // save output
    TString outdir = "PLOTS/study_timing_cut/jet_rates/";
    canvas->SaveAs(outdir  + category[k]+"_"+name[i]+"_jet_rate_correct.pdf");
    canvas->SaveAs(outdir  + category[k]+"_"+name[i] +"_jet_rate_correct.root");
    }// end categories
  } // end jet loop


} // end void
