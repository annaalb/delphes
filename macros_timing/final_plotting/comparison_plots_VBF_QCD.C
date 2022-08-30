#include "TH1.h"
#include "TSystem.h"
#include <vector>

#include "Style.C"
//------------------------------------------------------------------------------
// run with
// root -l -b -q comparison_plots_VBF_QCD.C
void comparison_plots_VBF_QCD()
{
  //SetStyle();
  SetupGlobalStyle();
  TFile *f[3];
  TFile *f_QCD[3];
  // read from results.root
  f[0] = new TFile("PLOTS/study_timing_cut/VBF/macro_jets_snowmass_EtaMax0_correct/results.root");
  f[1] = new TFile("PLOTS/study_timing_cut/VBF/macro_jets_snowmass_EtaMax3_correct/results.root");
  f[2] = new TFile("PLOTS/study_timing_cut/VBF/macro_jets_snowmass_EtaMax4_correct/results.root");

  // f_QCD[0] = new TFile("PLOTS/study_timing_cut/QCD/macro_QCD_jets_EtaMax0/results.root");
  // f_QCD[1] = new TFile("PLOTS/study_timing_cut/QCD/macro_QCD_jets_EtaMax3/results.root");
  // f_QCD[2] = new TFile("PLOTS/study_timing_cut/QCD/macro_QCD_jets_EtaMax4/results.root");

  f_QCD[0] = new TFile("PLOTS/study_timing_cut/QCD_100kPU/macro_QCD_jets_track_eta_EtaMax0/results.root");
  f_QCD[1] = new TFile("PLOTS/study_timing_cut/QCD_100kPU/macro_QCD_jets_track_eta_EtaMax3/results.root");
  f_QCD[2] = new TFile("PLOTS/study_timing_cut/QCD_100kPU/macro_QCD_jets_track_eta_EtaMax4/results.root");

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
  TString IsPU[3] = {"", "_signal", "_PU"};
  TString eta[6] = {"", "_eta01p3", "_eta1p32", "_eta23", "_eta34", "_eta4"};
  TString etaname[6] = {"", "|#eta_{jet}| < 1.3", "1.3<|#eta_{jet}|<2", "2<|#eta_{jet}|<3", "3<|#eta_{jet}|<4", "|#eta_{jet}|>4"};

  TString timing_scenario[3] = {"_no_timing", "_timing_etaMax3", "_timing_etaMax4"};
  TString timing_scenario_name[3] = {"No timing", "Timing for |#eta|<3", "Timing for |#eta|<4"};

  TString variable[7] = {"_jet_time_minus_PV", "_jet_pt_weighted_time_minus_PV",  "_track_DZ", "_track_T", "_track_eta", "_jet_pt", "_jet_eta" };
  TString variablenames[7] = {"t (jet - LV) [ns]", "t_{p_{T} weighted} (jet - LV) [ns]" ,"track dz [m]", "track dt [ns]", "Track #eta", "Jet p_{T}", "Jet #eta"};

  Int_t k = 0;
  Int_t e = 0;
//  Int_t l = 3;

  // get the histograms
  for (size_t s = 0; s < 3; s++) { // timing_scenario
  for (size_t i = 1; i < 2; i++) { // CHS or PUPPI
      for (size_t l = 5; l < 7; l++) {  // variables
        //for (size_t e = 0; e < 5; e++) { // eta bins
        //  for (size_t k = 0; k < 3; k++) { // split tracks into signal / PU
    track_VBF[i][k] = ((TH1F*)f[s]->Get(category[0]+"_"+name[i]+IsPU[k]+variable[l]+eta[e]));
    track_b[i][k] = ((TH1F*)f[s]->Get(category[1]+"_"+name[i]+IsPU[k]+variable[l]+eta[e]));
    track_PU[i][k] = ((TH1F*)f_QCD[s]->Get(category[2]+"_"+name[i]+IsPU[k]+variable[l]+eta[e]));

    cout << "Get hists...." << endl;
    cout << category[0]+"_"+name[i]+IsPU[k]+variable[l]+eta[e] << endl;
    cout << category[1]+"_"+name[i]+IsPU[k]+variable[l]+eta[e] << endl;
    cout << category[2]+"_"+name[i]+IsPU[k]+variable[l]+eta[e] << endl;

    // VBF jets
    track_VBF[i][k]->SetLineColor(kRed);
    track_VBF[i][k]->Scale(1/track_VBF[i][k]->Integral());
    track_VBF[i][k]->Draw("HIST");

    // set stat boxes
    canvas->Update();
    stats = (TPaveStats*)track_VBF[i][k]->GetListOfFunctions()->FindObject("stats");
    stats->SetName("h1stats");
    stats->SetY1NDC(.8);
    stats->SetY2NDC(.9);
    stats->SetX1NDC(.65);
    stats->SetX2NDC(.9);
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
    stats1->SetY1NDC(.7);
    stats1->SetY2NDC(.8);
    stats1->SetX1NDC(.65);
    stats1->SetX2NDC(.9);
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
   stats2->SetY1NDC(.6);
   stats2->SetY2NDC(.7);
   stats2->SetX1NDC(.65);
   stats2->SetX2NDC(.9);
   stats2->SetTextColor(4);
   stats2->SetTextSize(0.02);
   stats2->SetBorderSize(0);
   //
    stats1->Draw();
    stats2->Draw();
    stats->Draw();

   canvas->Update();
   canvas->Modified();

   Double_t ymax = track_VBF[i][k]->GetMaximum() > track_PU[i][k]->GetMaximum() ? track_VBF[i][k]->GetMaximum() : track_PU[i][k]->GetMaximum();
  // Double_t ymax = track_VBF[i][k]->GetMaximum();

   Double_t ymin = 0;
   track_VBF[i][k]->GetYaxis()->SetRangeUser(ymin,ymax *1.8);

   track_VBF[i][k]->GetXaxis()->SetTitle(variablenames[l]);

  // if(l==2) track_VBF[i][k]->GetXaxis()->SetNdivisions(504);
   if(l==0 || l==1) track_VBF[i][k]->GetXaxis()->SetRangeUser(-0.4,0.4);
   track_VBF[i][k]->GetYaxis()->SetTitleOffset(2);
   track_VBF[i][k]->GetYaxis()->SetTitle("fraction of jets");

   relPosX    = 0.045;
   CMS_lumi(canvas, 0, 10);
   // // legend
   // auto legend = new TLegend(0.25,0.75,0.85,0.85); // x1, y1,x2, y2
   // legend->SetTextSize(.03);
   // legend->SetFillStyle(0);
   // legend->SetHeader(timing_scenario_name[s]);
   //legend->Draw();

   TPaveText *comment = new TPaveText(0.25, 0.75, 0.85, 0.85, "brNDC");
   comment->SetTextSize(.03);
   comment->SetTextAlign(12);
   comment->SetFillStyle(0);
   comment->SetBorderSize(0);
   comment->AddText(timing_scenario_name[s]);
   comment->AddText(etaname[e]);
   comment->Draw();

  // save output
  TString outdir = "PLOTS/study_timing_cut/comparison_plots_correct/";
  canvas->SaveAs(outdir  + name[i]+IsPU[k]+variable[l]+eta[e]+timing_scenario[s]+".pdf");
  //canvas->SaveAs(outdir  + name[i]+IsPU[k]+variable[l]+eta[e]+timing_scenario[s]+".root");

//}// end eta bins
} // end variables or signal / PU
} // end jet loop
} // end timing_scenario
} // end void
