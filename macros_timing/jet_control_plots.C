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
  TH1 *fJetPT[4];
  TH1 *fJetEta[4];
  TH1 *fJetTrackT[4][2];
  TH1 *fJetTrackZ[4][2];
};

//------------------------------------------------------------------------------

class ExRootResult;
class ExRootTreeReader;

//------------------------------------------------------------------------------

void BookHistogramsBasic(ExRootResult *result, MyPlots *plots)
{
  TString branchJetNames[3] = {"branchJetCHS", "branchJetPUPPI", "branchGenJet"};
  TString name2[2] = {"_signal", "_PU"};

  // book more histograms
  for (size_t i = 0; i < 3; i++) {
    plots->fJetPT[i] = result->AddHist1D(
      "jet_pt_"+branchJetNames[i], "all jet P_{T}",
      "jet P_{T}, GeV/c", "number of jets", 100, 0.0, 200.0);

    plots->fJetEta[i] = result->AddHist1D(
      "jet_eta_"+branchJetNames[i], "all jet #eta",
      "jet #eta", "number of jets", 100, -5.0, 5.0);


      for (size_t k = 0; k < 2; k++) {
        plots->fJetTrackT[i][k] = result->AddHist1D(
          "track_dT_"+branchJetNames[i]+name2[k], "track T",
          "track dT [ns]", "number of tracks", 100, -1, 1);
        plots->fJetTrackZ[i][k] = result->AddHist1D(
          "track_dz_"+branchJetNames[i]+name2[k], "track Z ",
          "track dZ[m]", "number of tracks", 100, -0.1, 0.1);
              }

  }

  }

  //------------------------------------------------------------------------------

  void AnalyseEvents(ExRootTreeReader *treeReader, MyPlots *plots)
  {
    //TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
    //TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
    //TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");
    //TClonesArray *branchEFlowCHS = treeReader->UseBranch("CHSeflow");
    TClonesArray *branchEFlowCHS2 = treeReader->UseBranch("ParticleFlowCandidateCHS");

    TClonesArray *branchPuppiParticle = treeReader->UseBranch("ParticleFlowCandidate");
    //TClonesArray *branchPuppiTrack = treeReader->UseBranch("PuppiTrack");
    //TClonesArray *branchPuppiTower = treeReader->UseBranch("PuppiTower");

    TClonesArray *branchJetCHS = treeReader->UseBranch("JetCHS");
    TClonesArray *branchJetPUPPI = treeReader->UseBranch("JetPUPPI");
    TClonesArray *branchJetConstituentsCHS = treeReader->UseBranch("JetConstituentsCHS"); // jet constituents

    TClonesArray *branchParticle = treeReader->UseBranch("Particle"); // for identification of VBF and b quarks
    TClonesArray *branchfilteredParticle = treeReader->UseBranch("filteredParticle"); // input to genjets
    TClonesArray *branchMergerParticle = treeReader->UseBranch("mergerParticle"); // input to tracks

    TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
    TClonesArray *branchVtx = treeReader->UseBranch("Vertex");

    TClonesArray *branchJets[3] = {branchJetCHS, branchJetPUPPI, branchGenJet};

    Long64_t allEntries = treeReader->GetEntries();

    cout << "** Chain contains " << allEntries << " events" << endl;

    Long64_t entry;

    TObject* object;
    Track* track;
    ParticleFlowCandidate* pf;

    Vertex* vtx;
    Double_t Pvtx_T, Pvtx_Z;

    Int_t i;
    // Loop over all events
    for(entry = 0; entry < allEntries; ++entry)
    {

      // Load selected branches with data from specified event
      treeReader->ReadEntry(entry);

      // -------------- Analyse vertex ---------------------
      // PV time and Z
      for(Int_t i = 0; i < branchVtx->GetEntriesFast(); ++i)
      {
        // Take vtx
        vtx = (Vertex*) branchVtx->At(i);
        if (vtx->Index == 0) {
          Pvtx_Z = vtx->Z;
          Pvtx_T = vtx->T;
          }
      }

      // --------Analyse jets -- 0 CHS, 1 PUPPI,  2 GenJets ---------
      for (Int_t m = 0; m < 3; m++) {
        if(branchJets[m]->GetEntriesFast() > 0)
        {
          Jet* jet;
          for (size_t k = 0; k < branchJets[m]->GetEntriesFast(); k++) {
            jet = (Jet*) branchJets[m]->At(k);
            //if(jet->PT>30){
             plots->fJetPT[m]->Fill(jet->PT);
             plots->fJetEta[m]->Fill(jet->Eta);

             // Loop over all jet's constituents
             for(size_t j = 0; j < jet->Constituents.GetEntriesFast(); ++j)
             {
               object = jet->Constituents.At(j);
               // Check if the constituent is accessible
               if(object == 0)
               {
                 continue;
               }
               if(object->IsA() == ParticleFlowCandidate::Class())
               {
                 pf = (ParticleFlowCandidate*) object;
                 if (pf->Charge != 0) {
                 GenParticle *p = static_cast<GenParticle*>(pf->Particles.At(0));
                 if (p->IsPU == 0) { // signal
                   plots->fJetTrackT[m][0]->Fill((pf->T - Pvtx_T) * 1000000000);
                   plots->fJetTrackZ[m][0]->Fill((pf->Z - Pvtx_Z) / 1000);
                   //plots->fJetTrackZ[m][0]->Fill((pf->dZ) / 1000);

                  }
                  else if (p->IsPU == 1) { // PU tracks
                    plots->fJetTrackT[m][1]->Fill((pf->T - Pvtx_T) * 1000000000);
                    plots->fJetTrackZ[m][1]->Fill((pf->Z - Pvtx_Z) / 1000);
                  }

                  }
                }
              } // end loop over constituents
            //}
          }
        }
      } // loop over jet branches


    }// end loop over entries
  }//end void

  //------------------------------------------------------------------------------

  void PrintHistograms(ExRootResult *result, MyPlots *plots)
  {
    result->Print("png");
  }

  //------------------------------------------------------------------------------

  void jet_control_plots(const char *inputFile)
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
    AnalyseEvents(treeReader, plots);
    gSystem->cd("Plots/final/jet_control_plots/");
    cout << "Print hists "<< endl;

    PrintHistograms(result, plots);

    result->Write("results.root");

    cout << "** Exiting..." << endl;

    delete plots;
    delete result;
    delete treeReader;
    delete chain;
  }
