#include "TH1.h"
#include "TSystem.h"
#include <vector>

#include "Style.C"
//------------------------------------------------------------------------------

void comparison_plots_VBF_QCD()
{
  //SetStyle();
  SetupGlobalStyle();
  // read from results.root
  TFile *f_Eta0 = new TFile("PLOTS/10kevents/VBF_sample/macro_jets_snowmass_dz_smeared_JES_applied/results.root");
  TFile *f_Eta0_QCD = new TFile("PLOTS/10kevents/QCD/macro_QCD_jets_dz_smeared_JES_applied/results.root");
  TFile *f_Eta3 = new TFile("PLOTS/study_timing_cut/VBF/macro_jets_snowmass_EtaMax3/results.root");
  TFile *f_Eta3_QCD = new TFile("PLOTS/study_timing_cut/QCD/macro_QCD_jets_EtaMax3/results.root");
  TFile *f_Eta4 = new TFile("PLOTS/study_timing_cut/VBF/macro_jets_snowmass_EtaMax3/results.root");
  TFile *f_Eta4_QCD = new TFile("PLOTS/study_timing_cut/QCD/macro_QCD_jets_EtaMax3/results.root");

  TCanvas* canvas = new TCanvas("c", "c", 600, 600);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);
  TH1F* track_VBF[3][2];
  TH1F* track_b[3][2];
  TH1F* track_PU[3][2];

  TPaveStats *stats;
  TPaveStats *stats1;
  TPaveStats *stats2;

  TString name[3] = {"CHS", "PUPPI", "GenJet"};
  TString category[3] = {"VBF", "b", "PU"};

  // get the histograms
    track_VBF[i][k] = ((TH1F*)f->Get(category[0]+"_"+name[i]+IsPU[k]+variable[l]+eta[e]));
    track_b[i][k] = ((TH1F*)f->Get(category[1]+"_"+name[i]+IsPU[k]+variable[l]+eta[e]));
    track_PU[i][k] = ((TH1F*)f_QCD->Get(category[2]+"_"+name[i]+IsPU[k]+variable[l]+eta[e]));

    // VBF jets
    track_VBF[i][k]->SetLineColor(kRed);
    track_VBF[i][k]->Scale(1/track_VBF[i][k]->Integral());
    track_VBF[i][k]->Draw("HIST");

    // set stat boxes
    canvas->Update();
    stats = (TPaveStats*)track_VBF[i][k]->GetListOfFunctions()->FindObject("stats");
    stats->SetName("h1stats");
    stats->SetY1NDC(.75);
    stats->SetY2NDC(.85);
    stats->SetX1NDC(.6);
    stats->SetX2NDC(.85);
    stats->SetTextColor(2);
    stats->SetTextSize(0.02);
    stats->SetBorderSize(0);

    // b jets
    track_b[i][k]->SetLineColor(kBlack);
    track_b[i][k]->Scale(1/track_b[i][k]->Integral());
    track_b[i][k]->Draw("HIST same");

    // set stat boxes
    canvas->Update();
    stats1 = (TPaveStats*)track_b[i][k]->GetListOfFunctions()->FindObject("stats");
    stats1->SetName("h1stats1");
    stats1->SetY1NDC(.65);
    stats1->SetY2NDC(.75);
    stats1->SetX1NDC(.6);
    stats1->SetX2NDC(.85);
    stats1->SetTextColor(1);
    stats1->SetTextSize(0.02);
    stats1->SetBorderSize(0);

    // PU jet
    track_PU[i][k]->SetLineColor(kBlue);
    track_PU[i][k]->Scale(1/track_PU[i][k]->Integral());
    track_PU[i][k]->Draw(" HIST same");

    canvas->Update();
    stats2 = (TPaveStats*)track_PU[i][k]->GetListOfFunctions()->FindObject("stats");
   stats2->SetName("h1stats2");
   stats2->SetY1NDC(.55);
   stats2->SetY2NDC(.65);
   stats2->SetX1NDC(.6);
   stats2->SetX2NDC(.85);
   stats2->SetTextColor(4);
   stats2->SetTextSize(0.02);
   stats2->SetBorderSize(0);

   stats1->Draw();
   stats2->Draw();
   stats->Draw();

   canvas->Update();
   canvas->Modified();

   Double_t ymax = track_VBF[i][k]->GetMaximum() > track_PU[i][k]->GetMaximum() ? track_VBF[i][k]->GetMaximum() : track_PU[i][k]->GetMaximum();
   Double_t ymin = 0;
   track_VBF[i][k]->GetYaxis()->SetRangeUser(ymin,ymax *1.8);

   track_VBF[i][k]->GetXaxis()->SetTitle(variablenames[l]+eta[e]);

  // if(l==2) track_VBF[i][k]->GetXaxis()->SetNdivisions(504);
   if(l==0 || l==1) track_VBF[i][k]->GetXaxis()->SetRangeUser(-0.4,0.4);
   track_VBF[i][k]->GetYaxis()->SetTitleOffset(2);
   track_VBF[i][k]->GetYaxis()->SetTitle("fraction of jets");

  // legend
  // auto legend = new TLegend(0.25, 0.7, 0.55, 0.85);
  // legend->SetTextSize(.04);
  // legend->SetFillStyle(0);
  // //legend->SetHeader(xx);
  // legend->AddEntry(track_VBF[0][k],"signal", "l");
  // legend->AddEntry(track_PU[0][k],"PU", "l");

  //legend->Draw();
  // save output
  TString outdir = "PLOTS/10kevents/comparison_VBF_QCD_dz_smearing_JES_applied/";
  canvas->SaveAs(outdir  + name[i]+IsPU[k]+variable[l]+eta[e]+".pdf");
  canvas->SaveAs(outdir  + name[i]+IsPU[k]+variable[l]+eta[e]+ ".root");

//} // end if track variable
//} // end signal / PU (for tracks)
}// end eta bins
} // end variables
} // end jet loop
} // end void
