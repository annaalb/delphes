/*
Simple macro showing how to access branches from the delphes output root file,
loop over events, store histograms in a root file and print them as image files.

root -l examples/macro_jets.C'("delphes_output.root")'
*/

#include "TH1.h"
#include "TSystem.h"
#include <vector>

#include "Style.C"
#include "helper_functions.C"

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#endif

//------------------------------------------------------------------------------

struct MyPlots
{
  TH1* fGenJetPT[6]; // 6 eta bins
  TH1* fmatchedGenJetPT[2][6]; // CHS or PUPPI

  TH1* fRecoJetPT[2][6]; // CHS and PUPPI reco jets; 6 eta bins
  TH1* fmatchedRecoJetPT[2][6];

  TH1* fmatchedRecoJetResponse[2][6][7];
};

//------------------------------------------------------------------------------

class ExRootResult;
class ExRootTreeReader;

//------------------------------------------------------------------------------

void BookHistogramsBasic(ExRootResult *result, MyPlots *plots)
{
  TString name[3] = {"PUPPI", "CHS", "GenJet"};
  TString ptbin[7] = {"pt0toInf","pt0to20","pt20to30","pt30to50","pt50to80","pt80to100","pt100toInf"};
  TString etabin[6] = {"eta0p0toInf", "eta0p0to1p3","eta1p3to2p0","eta2p0to3","eta3to4","eta4p0toInf"};

  for (size_t e = 0; e < 6; e++) {
  plots->fGenJetPT[e] = result->AddHist1D(
    "GenJet_pt_all_"+etabin[e], "all jet P_{T}",
    " GenJet P_{T}, GeV/c", "number of jets", 300, 0.0, 300.0);

  for (size_t i = 0; i < 2; i++) {
    plots->fmatchedGenJetPT[i][e] = result->AddHist1D(
      "GenJet_pt_"+name[i]+"_matched_"+etabin[e], "matched GenJet P_{T}",
      " GenJet ("+name[i]+" matched) P_{T}, GeV/c", "number of jets", 300, 0.0, 300.0);
    plots->fRecoJetPT[i][e] = result->AddHist1D(
      name[i]+"_RecoJet_pt_all_"+etabin[e], "all jet P_{T}",
      name[i]+" RecoJet P_{T}, GeV/c", "number of jets", 300, 0.0, 300.0);
    plots->fmatchedRecoJetPT[i][e] = result->AddHist1D(
      name[i]+"_RecoJet_pt_matched_"+etabin[e], "matched RecoJet P_{T}",
      name[i]+" RecoJet (matched) P_{T}, GeV/c", "number of jets", 300, 0.0, 300.0);

      for (size_t p = 0; p < 7; p++) {
        plots->fmatchedRecoJetResponse[i][e][p] = result->AddHist1D(
          name[i]+"_RecoJet_response_"+etabin[e]+"_"+ptbin[p], " RecoJet response",
          name[i]+" P_{T, reco}/P_{T, gen}", "number of jets", 200, 0.0, 2.0);
      }

    }
  }

}

  //------------------------------------------------------------------------------

  void AnalyseEvents(ExRootTreeReader *treeReader, MyPlots *plots, string etamax)
  {

  //  TClonesArray *branchCHSeflow = treeReader->UseBranch("ParticleFlowCandidateCHS");
    TClonesArray *branchPuppiParticle = treeReader->UseBranch("ParticleFlowCandidate");

    TClonesArray *branchParticle = treeReader->UseBranch("Particle"); // for identification of VBF and b quarks
  //  TClonesArray *branchfilteredParticle = treeReader->UseBranch("filteredParticle"); // input to genjets

  //  TClonesArray *branchMergerParticle = treeReader->UseBranch("mergerParticle"); // input to tracks

  //  TClonesArray *branchJetCHS = treeReader->UseBranch("JetCHS");
    TClonesArray *branchJetPUPPI = treeReader->UseBranch("JetPUPPI");

    TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");

    TClonesArray *branchVtx = treeReader->UseBranch("Vertex");

    Int_t nbranches = 2;

    //TClonesArray *branchJets[3] = {branchJetPUPPI, branchJetCHS, branchGenJet};
    TClonesArray *branchJets[2] = {branchJetPUPPI, branchGenJet};

    Long64_t allEntries = treeReader->GetEntries();

    cout << "** Chain contains " << allEntries << " events" << endl;

    GenParticle *genparticle;
    Jet *genjet;

    Jet *matched_to_genjet = nullptr;
    Bool_t matched;

    ParticleFlowCandidate *pf;

    TObject* object;

    GenParticle* particle;
    GenParticle* constituent;
    Track* track;
    Tower* tower;
    Vertex* vtx;

    Long64_t entry;

    Bool_t debug=false;
    Bool_t do_JEC = false;

    string jecfile_CHS = "/afs/cern.ch/user/a/aalbrech/Delphes/macros_timing/CHS_JEC_PUPPI_"+etamax+".txt";
    string jecfile_PUPPI = "/afs/cern.ch/user/a/aalbrech/Delphes/macros_timing/PUPPI_JEC_PUPPI_"+etamax+".txt";

    string jecfile;

    // Loop over all events
    for(entry = 0; entry < allEntries; ++entry)
    {
      cout << "Process event " << entry+1 << " of total " << allEntries << endl;
      if(debug){cout << "-------- begin event ---------"<< endl;}

      // Load selected branches with data from specified event
      treeReader->ReadEntry(entry);
      if(debug){cout << "-------- begin analyse genparticles ---------"<< endl;}

      // efficiency plots
      if(branchJets[nbranches-1]->GetEntriesFast() > 0) // genjets
      {
        Jet* jet;
        matched = false;
        for (size_t k = 0; k < branchJets[nbranches-1]->GetEntriesFast(); k++) {
          jet = (Jet*) branchJets[nbranches-1]->At(k);
          Int_t eta_index;
          if(jet->PT>20){
          if (abs(jet->Eta) < 1.3) {eta_index=1;};
          if (abs(jet->Eta) > 1.3 && abs(jet->Eta) < 2) {eta_index=2;};
          if (2 < abs(jet->Eta) && abs(jet->Eta) < 3) {eta_index=3;};
          if (3 < abs(jet->Eta) && abs(jet->Eta) < 4) {eta_index=4;};
          if (abs(jet->Eta) > 4) {eta_index=5;};

          plots->fGenJetPT[0]->Fill(jet->PT);
          plots->fGenJetPT[eta_index]->Fill(jet->PT);
          // now match to reco jet ( CHS or PUPPI )
          for (size_t m = 0; m < 2; m++) {
            if (m==0) {
              jecfile = jecfile_CHS;
            }
            else{
              jecfile = jecfile_PUPPI;
            }
            // apply SF on reco jets (stored in jecfile)
            matched = matched_to_recojet(branchJets[m], jet, 0.4, 10, branchJets[nbranches-1], do_JEC, jecfile); // is there a close recojet?
            if (matched) {
              plots->fmatchedGenJetPT[m][0]->Fill(jet->PT);
              plots->fmatchedGenJetPT[m][eta_index]->Fill(jet->PT);
            }
          }
        }
        }
      } // end genjets

    // --------Analyse reco jets for purity -- 0 CHS, 1 PUPPI ---------
    for (Int_t m = 0; m < nbranches-1; m++) {
      if (m==0) {
        jecfile = jecfile_CHS;
      }
      else{
        jecfile = jecfile_PUPPI;
      }
      if(branchJets[m]->GetEntriesFast() > 0) // reco
      {
        Jet* jet;
        matched = false;
        for (size_t k = 0; k < branchJets[m]->GetEntriesFast(); k++) {
          jet = (Jet*) branchJets[m]->At(k);
          Int_t eta_index;
          Int_t pt_index;

          // get jet pt, SF, store corrected_pt
          Double_t reco_pt = jet->PT;
          Double_t sf;
          if (do_JEC) {
            sf = get_SF(jet, branchJets[nbranches-1], jecfile);
            reco_pt = reco_pt*sf;
          }

          if(reco_pt>20){
          if (abs(jet->Eta) < 1.3) {eta_index=1;};
          if (abs(jet->Eta) > 1.3 && abs(jet->Eta) < 2) {eta_index=2;};
          if (2 < abs(jet->Eta) && abs(jet->Eta) < 3) {eta_index=3;};
          if (3 < abs(jet->Eta) && abs(jet->Eta) < 4) {eta_index=4;};
          if (abs(jet->Eta) > 4) {eta_index=5;};

          plots->fRecoJetPT[m][0]->Fill(reco_pt);
          plots->fRecoJetPT[m][eta_index]->Fill(reco_pt);
          // now match to gen jet
          matched = matched_to_jet(branchJets[nbranches-1], jet, 0.4, 10); // is there a close genjet?
            if (matched) {
              genjet = get_closest_jet(branchJets[nbranches-1], jet, 0.4);

              plots->fmatchedRecoJetPT[m][0]->Fill(reco_pt);
              plots->fmatchedRecoJetPT[m][eta_index]->Fill(reco_pt);

              if (abs(genjet->Eta) < 1.3) {eta_index=1;};
              if (abs(genjet->Eta) > 1.3 && abs(genjet->Eta) < 2) {eta_index=2;};
              if (2 < abs(genjet->Eta) && abs(genjet->Eta) < 3) {eta_index=3;};
              if (3 < abs(genjet->Eta) && abs(genjet->Eta) < 4) {eta_index=4;};
              if (abs(genjet->Eta) > 4) {eta_index=5;};

              if (abs(genjet->PT) < 20) {pt_index=1;};
              if (abs(genjet->PT) > 20 && abs(genjet->PT) < 30) {pt_index=2;};
              if (abs(genjet->PT) > 30 && abs(genjet->PT) < 50) {pt_index=3;};
              if (abs(genjet->PT) > 50 && abs(genjet->PT) < 80) {pt_index=4;};
              if (abs(genjet->PT) > 80 && abs(genjet->PT) < 100) {pt_index=5;};
              if (abs(genjet->PT) > 100) {pt_index=6;};

              if (debug) {cout << " Fill response plot: reco pt = " << reco_pt << endl;}

              plots->fmatchedRecoJetResponse[m][0][0]->Fill(reco_pt/genjet->PT);
              plots->fmatchedRecoJetResponse[m][0][pt_index]->Fill(reco_pt/genjet->PT);
              plots->fmatchedRecoJetResponse[m][eta_index][0]->Fill(reco_pt/genjet->PT);
              plots->fmatchedRecoJetResponse[m][eta_index][pt_index]->Fill(reco_pt/genjet->PT);
            }
        }
        }
      } // end recojet

    } // end loop over jet branches


    }// end loop over entries

  }//end void

  //------------------------------------------------------------------------------

  void PrintHistograms(ExRootResult *result, MyPlots *plots)
  {
    result->Print("png");
  }

  //------------------------------------------------------------------------------

// run with:
// root -l -b -q macros_timing/macro_efficiency_purity.C'("ROOTOUTPUT/VBF/EtaMax3/CMS_PhaseII_Snowmass_200PU_VBF_10k_ptmin10_EtaMax3.root", "PLOTS/study_timing_cut/VBF/EtaMax3")'
  void macro_efficiency_purity(const char *inputFile, const char *outdir, string etamax="EtaMax0")
  {
    gSystem->Load("libDelphes");

    TChain *chain = new TChain("Delphes");
    chain->Add(inputFile);

    ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
    ExRootResult *result = new ExRootResult();

    MyPlots *plots = new MyPlots;
cout << "Book hists "<< endl;
    BookHistogramsBasic(result, plots);

    SetupGlobalStyle();
    //Simulation_label();
cout << "Analyse event "<< endl;
    AnalyseEvents(treeReader, plots, etamax);
    gSystem->cd(outdir);
    cout << "Print hists "<< endl;

    PrintHistograms(result, plots);

    result->Write("results.root");

    cout << "** Exiting..." << endl;

    delete plots;
    delete result;
    delete treeReader;
    delete chain;
  }
