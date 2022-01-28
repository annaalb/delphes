#include "TH1.h"
#include "TSystem.h"
#include <vector>
#include "TGraphErrors.h"

#include "Style.C"
//------------------------------------------------------------------------------

void number_of_tracks()
{
  //SetStyle();
  SetupGlobalStyle();
  // read from results.root
  TFile *f = new TFile("PLOTS/10kevents/VBF_sample/macro_jets_snowmass_dz_smeared_JES_applied/results.root");
  //TFile *f_QCD = new TFile("PLOTS/10kevents/QCD/macro_jets_snowmass/results.root");
  TFile *f_QCD = new TFile("PLOTS/10kevents/QCD/macro_QCD_jets_dz_smeared_JES_applied/results.root");

  TCanvas* canvas = new TCanvas("c", "c", 600, 600);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);
  TH1F* track_VBF[3];
  TH1F* track_PU[3];
  TH1F* track_VBF_signal[3];
  TH1F* track_PU_signal[3];

  TString name[3] = {"CHS", "PUPPI", "GenJet"};
  TString category[3] = {"VBF", "b", "PU"};
  TString eta[6] = {"", "_eta01p3", "_eta1p32", "_eta23", "_eta34", "_eta4"};
  Double_t etamiddle[6] = {-1, 0, 1.3, 2, 3, 4};

  TString variable[3] = {"_jet_n_tracks", "_jet_n_signal_tracks", "_jet_n_PU_tracks"};

  TGraphErrors* _jet_n_PU_tracks[2];
  TGraphErrors* _jet_n_signal_tracks[2];

  Double_t n_signal_tracks[6];
  Double_t n_PU_tracks[6];

  Int_t k = 0;
  // get the histograms
  for (size_t i = 0; i < 1; i++) { // CHS or PUPPI
    cout << "------------------" << name[i] << "---------------" << endl;

        for (size_t e = 0; e < 6; e++) { // eta bins
          cout << "------------------" << eta[e] << "---------------" << endl;
    track_VBF[i] = ((TH1F*)f->Get(category[0]+"_"+name[i]+variable[3]+eta[e]));
    //track_PU[i] = ((TH1F*)f_QCD->Get(category[2]+"_"+name[i]+variable[3]+eta[e]));

    track_VBF_signal[i] = ((TH1F*)f->Get(category[0]+"_"+name[i]+variable[2]+eta[e]));
    //track_PU_signal[i] = ((TH1F*)f_QCD->Get(category[2]+"_"+name[i]+variable[2]+eta[e]));

    // cout << "VBF "<< track_VBF[i]->GetMean() << endl;
    // cout << "PU "<< track_PU[i]->GetMean() << endl;

    n_signal_tracks[e]=track_VBF_signal[i]->GetMean();
    n_PU_tracks[e]=track_VBF[i]->GetMean();


}// end eta bins

_jet_n_signal_tracks[i] = new TGraphErrors(6, etamiddle, n_signal_tracks[k]);
_jet_n_PU_tracks[i] = new TGraphErrors(6, etamiddle, n_PU_tracks[k]);

_jet_n_PU_tracks[i]->Draw();

TString outdir = "PLOTS/10kevents/comparison_VBF_QCD_dz_smearing_JES_applied/";
canvas->SaveAs(outdir  + name[i]+"_n_tracks.pdf");
canvas->SaveAs(outdir  + name[i]+"_n_tracks.root");
} // end jet loop
} // end void
