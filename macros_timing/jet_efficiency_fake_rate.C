#include "TH1.h"
#include "TSystem.h"
#include <vector>

#include "Style.C"
//------------------------------------------------------------------------------

void jet_efficiency_fake_rate()
{
  //SetStyle();
  SetupGlobalStyle();
  // read from results.root
  TFile *f = new TFile("PLOTS/10kevents/VBF_sample/macro_jets_snowmass_dz_smeared_JES_applied/results.root");
  TFile *f_QCD = new TFile("PLOTS/10kevents/QCD/macro_QCD_jets_dz_smeared_JES_applied/results.root");

  TCanvas* canvas = new TCanvas("c", "c", 600, 600);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);
  TH1F* hist_VBF[3][6];
  TH1F* hist_b[3][6];
  TH1F* hist_PU[3][6];

  Int_t n_VBF_jets[3][6]; // CHS, PUPPI genjets ; eta bins
  Int_t n_b_jets[3][6]; // CHS, PUPPI genjets ; eta bins
  Int_t n_PU_jets[3][6]; // CHS, PUPPI genjets ; eta bins

  TString name[3] = {"CHS", "PUPPI", "GenJet"};
  TString category[3] = {"VBF", "b", "PU"};
  TString eta[6] = {"", "_eta01p3", "_eta1p32", "_eta23", "_eta34", "_eta4"};

  // get the histograms
  for (size_t i = 0; i < 2; i++) { // CHS or PUPPI
    for (size_t k = 0; k < 6; k++) { // eta bins

    hist_VBF[i][k] = ((TH1F*)f->Get(category[0]+"_"+name[i]+"_jet_time_minus_PV"+eta[k]));
    hist_b[i][k] = ((TH1F*)f->Get(category[1]+"_"+name[i]+"_jet_time_minus_PV"+eta[k]));
    hist_PU[i][k] = ((TH1F*)f_QCD->Get(category[2]+"_"+name[i]+"_jet_time_minus_PV"+eta[k]));

    n_VBF_jets[i][k] = hist_VBF[i][k]->Integral();

    cout << n_VBF_jets[i][k] << endl;
    }// end eta bins
  } // end jet loop



  //  Double_t ymax = hist_VBF[i][k]->GetMaximum() > hist_PU[i][k]->GetMaximum() ? hist_VBF[i][k]->GetMaximum() : hist_PU[i][k]->GetMaximum();
  //  Double_t ymin = 0;
  //  hist_VBF[i][k]->GetYaxis()->SetRangeUser(ymin,ymax *1.8);
  //
  //  hist_VBF[i][k]->GetXaxis()->SetTitle(variablenames[l]+eta[e]);
  //
  // // if(l==2) hist_VBF[i][k]->GetXaxis()->SetNdivisions(504);
  //  if(l==0 || l==1) hist_VBF[i][k]->GetXaxis()->SetRangeUser(-0.4,0.4);
  //  hist_VBF[i][k]->GetYaxis()->SetTitleOffset(2);
  //  hist_VBF[i][k]->GetYaxis()->SetTitle("fraction of jets");
  //
  // // legend
  // // auto legend = new TLegend(0.25, 0.7, 0.55, 0.85);
  // // legend->SetTextSize(.04);
  // // legend->SetFillStyle(0);
  // // //legend->SetHeader(xx);
  // // legend->AddEntry(hist_VBF[0][k],"signal", "l");
  // // legend->AddEntry(hist_PU[0][k],"PU", "l");
  //
  // //legend->Draw();
  // // save output
  // TString outdir = "PLOTS/10kevents/jet_efficiency_fake_rate/";
  // canvas->SaveAs(outdir  + name[i]+variable[l]+".pdf");
  // canvas->SaveAs(outdir  + name[i]+variable[l]+ ".root");


} // end void
