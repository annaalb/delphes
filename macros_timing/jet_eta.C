#include "TH1.h"
#include "TSystem.h"
#include <vector>

#include "Style.C"
//------------------------------------------------------------------------------

void jet_eta()
{
  //SetStyle();
  SetupGlobalStyle();
  TFile *f[3];
  // read from results.root
  f[0] = new TFile("PLOTS/study_timing_cut/VBF/EtaMax0/results.root");
  f[1] = new TFile("PLOTS/study_timing_cut/VBF/EtaMax3/results.root");
  f[2] = new TFile("PLOTS/study_timing_cut/VBF/EtaMax4/results.root");

  TCanvas* canvas = new TCanvas("c", "c", 600, 600);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);

  TH1 *hist_all[3][2][6]; // eta max, CHS or PUPPI, eta bins, pt bins

  TString name[3] = {"CHS", "PUPPI", "GenJet"};
  TString eta[6] = {"", "_eta01p3", "_eta1p32", "_eta23", "_eta34", "_eta4"};
  TString eta_names[6] = {"inclusive", "0 < #eta < 1.3", "1.3 < #eta < 2", "2 < #eta < 3", "3 < #eta < 4", "#eta > 4"};
  TString ptbin[7] = {"", "_pt20", "_20pt30", "_30pt50", "_50pt80", "_80pt100", "_pt100"};
  TString pt_names[7] = {"inclusive", "p_{T, gen}<20", "20<p_{T, gen}<30", "30<p_{T, gen}<50", "50<p_{T, gen}<80", "80<p_{T, gen}<100", "p_{T, gen}>100"};

  // get the histograms
  for (size_t i = 0; i < 2; i++) { // CHS or PUPPI
    for (size_t k = 0; k < 6; k++) { // eta bins

      for (size_t l = 0; l < 3; l++) {
        // Hists for efficiency
      //  hist_all[l][i][k] = ((TH1F*)f[l]->Get("GenJet_pt_all"+eta[k]));
        //hist_all[l][i][k] = ((TH1F*)f[l]->Get("GenJet_pt_"+name[i]+"_matched"+eta[k]));
        // Hists for purity
      //  hist_all[l][i][k] = ((TH1F*)f[l]->Get(name[i]+"_RecoJet_pt_all"+eta[k]));
        hist_all[l][i][k] = ((TH1F*)f[l]->Get(name[i]+"_RecoJet_pt_matched"+eta[k]));

      }

      hist_all[0][i][k]->SetLineStyle(33);

      hist_all[1][i][k]->SetLineColor(30);
      hist_all[1][i][k]->SetLineStyle(22);

      hist_all[2][i][k]->SetLineColor(46);
      hist_all[2][i][k]->SetLineStyle(20);

      TPaveStats *stats;
      TPaveStats *stats2;
      TPaveStats *stats3;

    // create a multigraph and draw it
    THStack  *mg  = new THStack();
    mg->Add(hist_all[0][i][k]);
    mg->Add(hist_all[1][i][k]);
    mg->Add(hist_all[2][i][k]);
    mg->Draw("nostack");
    Double_t ymax = hist_all[0][i][k]->GetMaximum() > hist_all[1][i][k]->GetMaximum() ? hist_all[0][i][k]->GetMaximum() : hist_all[1][i][k]->GetMaximum();
    Double_t ymin = 0;
    //mg->GetYaxis()->SetRangeUser(ymin,ymax *1.8);
    mg->SetMaximum(ymax *1.5);

    mg->GetXaxis()->SetTitle("P_{T}");
    mg->GetYaxis()->SetTitle("Number of jets");
    canvas->Modified();
    mg->Draw("nostack");

    canvas->Update();
    stats = (TPaveStats*)hist_all[0][i][k]->GetListOfFunctions()->FindObject("stats");
    stats->SetName("h1stats");
    stats->SetY1NDC(.75);
    stats->SetY2NDC(.85);
    stats->SetX1NDC(.65);
    stats->SetX2NDC(.85);
    //stats->SetTextColor(2);
    stats->SetTextSize(0.03);
    stats->SetBorderSize(0);
    stats->SetOptStat(1100);

    canvas->Update();
    stats2 = (TPaveStats*)hist_all[1][i][k]->GetListOfFunctions()->FindObject("stats");
     stats2->SetName("h1stats2");
     stats2->SetY1NDC(.65);
     stats2->SetY2NDC(.75);
     stats2->SetX1NDC(.65);
     stats2->SetX2NDC(.85);
     stats2->SetTextColor(30);
     stats2->SetTextSize(0.03);
     stats2->SetBorderSize(0);
     stats2->SetOptStat(1);

     canvas->Update();
     stats3 = (TPaveStats*)hist_all[2][i][k]->GetListOfFunctions()->FindObject("stats");
      stats3->SetName("h1stats3");
      stats3->SetY1NDC(.55);
      stats3->SetY2NDC(.65);
      stats3->SetX1NDC(.65);
      stats3->SetX2NDC(.85);
      stats3->SetTextColor(46);
      stats3->SetTextSize(0.03);
      stats3->SetBorderSize(0);
      stats3->SetOptStat(1110);

      canvas->Modified();
      mg->Draw("nostack");

    // // legend
    auto legend = new TLegend(0.25,0.65,0.85,0.85); // x1, y1,x2, y2
    legend->SetTextSize(.03);
    legend->SetFillStyle(0);
    legend->SetHeader(name[i]+", "+eta_names[k]);
    legend->AddEntry(hist_all[0][i][k],"no timing", "l");
    legend->AddEntry(hist_all[1][i][k],"timing for |#eta|<3", "l");
    legend->AddEntry(hist_all[2][i][k],"timing for |#eta|<4", "l");

    legend->Draw();

    // // save output
    TString outdir = "PLOTS/study_timing_cut/jet_efficiency_purity/";
    canvas->SaveAs(outdir  + name[i]+eta[k]+"_RecoJet_matched_pt.pdf");
    canvas->SaveAs(outdir  + name[i]+eta[k]+ "_RecoJet_matched_pt.root");

    }// end eta bins
  } // end jet loop


} // end void
