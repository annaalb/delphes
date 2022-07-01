#include "TH1.h"
#include "TSystem.h"
#include <vector>

#include "Style.C"
//------------------------------------------------------------------------------

void general_plotting()
{
  //SetStyle();
  SetupGlobalStyle();
  TFile *f[3];
  TFile *f_QCD[3];
  // read from results.root
  //f[0] = new TFile("PLOTS/study_timing_cut/VBF/macro_jets_snowmass_EtaMax0/results.root");
  f[0] = new TFile("PLOTS/study_timing_cut/VBF/macro_jets_snowmass_EtaMax0_noSignalCut/results.root");

  TCanvas* canvas = new TCanvas("c", "c", 600, 600);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);

  TH1 *hist; // eta max, CHS or PUPPI, eta bins

  TString category[3] = {"VBF", "b", "PU"};
  TString name[3] = {"CHS", "PUPPI", "GenJet"};
  TString variables[3] = {"_jet_deltaR", "_jet_eta", "_jet_pt"};

  // get the histograms
  for (size_t i = 0; i < 1; i++) { // CHS or PUPPI
      for (size_t k = 0; k < 1; k++) {
        for (size_t l = 0; l < 1; l++) {
          //hist = ((TH1F*)f[0]->Get(category[k]+"_"+name[i]+variables[l]));
          hist = ((TH1F*)f[0]->Get("genparticle_VBF_eta"));

          hist->Draw("H");

          //hist->GetXaxis()->SetTitle(variables[l]);
          hist->GetXaxis()->SetTitle("#eta");

        //  hist->GetYaxis()->SetTitle("Number of jets");
          hist->GetYaxis()->SetTitle("Number of particles");

          Double_t ymax = hist->GetMaximum();
          hist->SetMaximum(ymax * 1.8);

          TPaveStats *stats;
          canvas->Update();
          stats = (TPaveStats*)hist->GetListOfFunctions()->FindObject("stats");
          stats->SetName("h1stats");
          stats->SetY1NDC(.8);
          stats->SetY2NDC(.9);
          stats->SetX1NDC(.6);
          stats->SetX2NDC(.95);
          stats->SetTextSize(0.03);
          stats->SetBorderSize(0);
          stats->SetOptStat(1100);

          canvas->Modified();

          relPosX    = 0.045;
          CMS_lumi(canvas, 0, 10);

          TPaveText *comment = new TPaveText(0.25, 0.7, 0.55, 0.8, "brNDC");
          comment->SetTextSize(.03);
          comment->SetTextAlign(12);
          comment->SetFillStyle(0);
          comment->SetBorderSize(0);
          //comment->AddText("AK4 "+ name[i]+", "+category[k]+" jets , p_{T}>20 GeV, no timing");
          comment->AddText("VBF generator signal particles, p_{T}>20 GeV");
          comment->Draw();

          // // save output
          TString outdir = "PLOTS/study_timing_cut/eta_comparison/";
          //canvas->SaveAs(outdir  + category[k]+"_"+name[i]+variables[l]+"_noSignalCut.pdf");
          canvas->SaveAs(outdir  + "genparticle_VBF_eta.pdf");

        }
    }// end categories
  } // end jet loop
} // end void
