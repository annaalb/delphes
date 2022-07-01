#include "TH1.h"
#include "TSystem.h"
#include <vector>

#include "Style.C"
//------------------------------------------------------------------------------

void reco_jet_response()
{
  //SetStyle();
  SetupGlobalStyle();
  TFile *f[3];
  // read from results.root
  // f[0] = new TFile("PLOTS/study_timing_cut/VBF/Corrected_response/EtaMax0/results.root");
  // f[1] = new TFile("PLOTS/study_timing_cut/VBF/Corrected_response/EtaMax3/results.root");
  // f[2] = new TFile("PLOTS/study_timing_cut/VBF/Corrected_response/EtaMax4/results.root");

  //f[0] = new TFile("PLOTS/study_timing_cut/VBF/EtaMax0_origin/results.root");
  f[0] = new TFile("PLOTS/study_timing_cut/VBF/EtaMax0_noSignalRejection/results.root");
  f[1] = new TFile("PLOTS/study_timing_cut/VBF/EtaMax0_origin_noSignalCut/results.root");
  //f[2] = new TFile("PLOTS/study_timing_cut/VBF/EtaMax0_origin_SignalCut/results.root");
  f[2] = new TFile("PLOTS/study_timing_cut/VBF/Raw_response/EtaMax0/results.root");


  TCanvas* canvas = new TCanvas("c", "c", 600, 600);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);

  TH1 *hist_all[3][2][6][7]; // eta max, CHS or PUPPI, eta bins, pt bins

  TString name[3] = {"PUPPI", "CHS", "GenJet"};
  TString ptbin[7] = {"pt0to20","pt20to30","pt30to50","pt50to80","pt80to100","pt100toInf"};
  TString eta[6] = {"eta0p0to1p3","eta1p3to2p0","eta2p0to3","eta3to4","eta4p0toInf"};
  TString eta_names[6] = {"0 < #eta < 1.3", "1.3 < #eta < 2", "2 < #eta < 3", "3 < #eta < 4", "#eta > 4"};
  TString pt_names[7] = {"p_{T, gen}<20", "20<p_{T, gen}<30", "30<p_{T, gen}<50", "50<p_{T, gen}<80", "80<p_{T, gen}<100", "p_{T, gen}>100"};

  // get the histograms
  for (size_t i = 0; i < 1; i++) { // CHS or PUPPI
    for (size_t k = 0; k < 5; k++) { // eta bins
      for (size_t p = 0; p < 6; p++) { // pt bins

      for (size_t l = 0; l < 2; l++) {
        hist_all[l][i][k][p] = ((TH1F*)f[l]->Get(name[i]+"_RecoJet_response_"+eta[k]+"_"+ptbin[p]));
        cout << name[i]+"_RecoJet_response_"+eta[k]+"_"+ptbin[p] << endl;
        hist_all[l][i][k][p]->Rebin(2);
        // if (l==2) {
        //   hist_all[l][i][k][p]->Scale(1/10); // norm plots to 1000 events
        // }
      }

      hist_all[0][i][k][p]->SetLineStyle(33);

      hist_all[1][i][k][p]->SetLineColor(30);
      hist_all[1][i][k][p]->SetLineStyle(22);

      // hist_all[2][i][k][p]->SetLineColor(46);
      // hist_all[2][i][k][p]->SetLineStyle(20);

      TPaveStats *stats;
      TPaveStats *stats2;
      TPaveStats *stats3;

    // create a multigraph and draw it
    THStack  *mg  = new THStack();
    mg->Add(hist_all[0][i][k][p]);
    mg->Add(hist_all[1][i][k][p]);
  //  mg->Add(hist_all[2][i][k][p]);
    mg->Draw("nostack");
    Double_t ymax = hist_all[0][i][k][p]->GetMaximum() > hist_all[1][i][k][p]->GetMaximum() ? hist_all[0][i][k][p]->GetMaximum() : hist_all[1][i][k][p]->GetMaximum();
    Double_t ymin = 0;
    //mg->GetYaxis()->SetRangeUser(ymin,ymax *1.8);
    mg->SetMaximum(ymax *1.8);

    mg->GetXaxis()->SetTitle("P_{T, reco} / P_{T, gen}");
    mg->GetYaxis()->SetTitle("Number of jets");
  //  mg->GetYaxis()->SetTitle("Fraction of jets");
    canvas->Modified();
    mg->Draw("nostack");

    canvas->Update();
    stats = (TPaveStats*)hist_all[0][i][k][p]->GetListOfFunctions()->FindObject("stats");
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
    stats2 = (TPaveStats*)hist_all[1][i][k][p]->GetListOfFunctions()->FindObject("stats");
     stats2->SetName("h1stats2");
     stats2->SetY1NDC(.65);
     stats2->SetY2NDC(.75);
     stats2->SetX1NDC(.65);
     stats2->SetX2NDC(.85);
     stats2->SetTextColor(30);
     stats2->SetTextSize(0.03);
     stats2->SetBorderSize(0);
     stats2->SetOptStat(1);

     // canvas->Update();
     // stats3 = (TPaveStats*)hist_all[2][i][k][p]->GetListOfFunctions()->FindObject("stats");
     //  stats3->SetName("h1stats3");
     //  stats3->SetY1NDC(.55);
     //  stats3->SetY2NDC(.65);
     //  stats3->SetX1NDC(.65);
     //  stats3->SetX2NDC(.85);
     //  stats3->SetTextColor(46);
     //  stats3->SetTextSize(0.03);
     //  stats3->SetBorderSize(0);
     //  stats3->SetOptStat(1110);

      canvas->Modified();
      mg->Draw("nostack");

    // // legend
    auto legend = new TLegend(0.25,0.65,0.85,0.85); // x1, y1,x2, y2
    legend->SetTextSize(.03);
    legend->SetFillStyle(0);
    legend->SetHeader(name[i]+", "+pt_names[p]+", "+eta_names[k]);
    // legend->AddEntry(hist_all[0][i][k][p],"no timing", "l");
    // legend->AddEntry(hist_all[1][i][k][p],"timing for |#eta|<3", "l");
    // legend->AddEntry(hist_all[2][i][k][p],"timing for |#eta|<4", "l");
    legend->AddEntry(hist_all[0][i][k][p],"Our card + no Sig. Rej.", "l");
    legend->AddEntry(hist_all[1][i][k][p],"Snowmass", "l");
  //  legend->AddEntry(hist_all[2][i][k][p],"Our card (Sig. Rej.)", "l");

    legend->Draw();

    // // save output
    TString outdir = "PLOTS/study_timing_cut/Debug_response/";
    canvas->SaveAs(outdir  + name[i]+eta[k]+ptbin[p]+"_response.pdf");
    //canvas->SaveAs(outdir  + name[i]+eta[k]+ptbin[p]+ "_response.root");
    }
    }// end eta bins
  } // end jet loop


} // end void
