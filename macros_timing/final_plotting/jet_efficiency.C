#include "TH1.h"
#include "TSystem.h"
#include <vector>

#include "Style.C"
//------------------------------------------------------------------------------

void jet_efficiency()
{
  //SetStyle();
  SetupGlobalStyle();
  TFile *f[3];
  // read from results.root
  f[0] = new TFile("PLOTS/study_timing_cut/VBF/Corrected_response/EtaMax0/results.root");
  f[1] = new TFile("PLOTS/study_timing_cut/VBF/Corrected_response/EtaMax3/results.root");
  f[2] = new TFile("PLOTS/study_timing_cut/VBF/Corrected_response/EtaMax4/results.root");

  TCanvas* canvas = new TCanvas("c", "c", 600, 600);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);

  TH1 *hist_all[3][2][6]; // eta max, CHS or PUPPI, eta bins
  TH1 *hist_matched[3][2][6];
  TGraphAsymmErrors *eff[3][2][6];

  TString name[3] = {"CHS", "PUPPI", "GenJet"};
  TString ptbin[7] = {"pt0to20","pt20to30","pt30to50","pt50to80","pt80to100","pt100toInf"};
  TString eta[6] = {"eta0p0to1p3","eta1p3to2p0","eta2p0to3","eta3to4","eta4p0toInf"};
  TString eta_names[6] = {"0 < #eta < 1.3", "1.3 < #eta < 2", "2 < #eta < 3", "3 < #eta < 4", "#eta > 4"};
  Double_t ptbins[8] = {0,10,20,30,50,80,100,200};

  // get the histograms
  for (size_t i = 0; i < 2; i++) { // CHS or PUPPI
    for (size_t k = 0; k < 5; k++) { // eta bins

      for (size_t l = 0; l < 3; l++) {
        hist_all[l][i][k] = ((TH1F*)f[l]->Get("GenJet_pt_all_"+eta[k]));
        hist_matched[l][i][k] = ((TH1F*)f[l]->Get("GenJet_pt_"+name[i]+"_matched_"+eta[k]));

        // rebin pt plots
        hist_all[l][i][k] = hist_all[l][i][k]->Rebin(8, hist_all[l][i][k]->GetTitle(), ptbins);
        hist_matched[l][i][k] = hist_matched[l][i][k]->Rebin(8, hist_all[l][i][k]->GetTitle(), ptbins);

        eff[l][i][k] = new TGraphAsymmErrors(hist_matched[l][i][k], hist_all[l][i][k], "cl=0.683 b(1,1) mode");
      }

      eff[0][i][k]->SetMarkerStyle(33);

      eff[1][i][k]->SetLineColor(30);
      eff[1][i][k]->SetMarkerColor(30);
      eff[1][i][k]->SetMarkerStyle(22);

      eff[2][i][k]->SetLineColor(46);
      eff[2][i][k]->SetMarkerColor(46);
      eff[2][i][k]->SetMarkerStyle(20);

    // create a multigraph and draw it
    TMultiGraph  *mg  = new TMultiGraph();
    mg->Add(eff[0][i][k]);
    mg->Add(eff[1][i][k]);
    mg->Add(eff[2][i][k]);
    mg->GetXaxis()->SetTitle("Genjet p_{T}");
    mg->GetYaxis()->SetTitle("Efficiency");
    mg->Draw("AP");

    //gPad->SetLogy();

    // // legend
    auto legend = new TLegend(0.3,0.25,0.85,0.45); // x1, x2, y1, y2
    legend->SetTextSize(.03);
    legend->SetFillStyle(0);
    legend->SetHeader(name[i]+", p_{T, reco}>10 GeV, p_{T, gen}>20 GeV, "+eta_names[k]);
    legend->AddEntry(eff[0][i][k],"no timing", "l");
    legend->AddEntry(eff[1][i][k],"timing for |#eta|<3", "l");
    legend->AddEntry(eff[2][i][k],"timing for |#eta|<4", "l");

    legend->Draw();

    // // save output
    TString outdir = "PLOTS/study_timing_cut/Corrected_eff_purity/";
    canvas->SaveAs(outdir  + name[i]+eta[k]+"_efficiency.pdf");
  //  canvas->SaveAs(outdir  + name[i]+eta[k]+ "_efficiency.root");
    }// end eta bins
  } // end jet loop


} // end void
