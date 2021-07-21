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
  TH1 *fPUgenjetPT[6];
  TH1 *fPUgenjetRho;

  TH1 *fPUgenjetETA[3];
  TH1 *fPUgenjetNCharged[3];
  TH1 *fPUgenjetNNeutrals[3];

  TH1 *fPUgenjetConstPT;
  TH1 *fPUgenjetConstETA;

  TH1 *fPUgenjetCHPT[3];
  TH1 *fPUgenjetCHEta[3];
  TH1 *fPUgenjetCHdT[3];
  TH1 *fPUgenjetCHdZ[3];

  TH1 *fPUgenparticlePT;
  TH1 *fPUgenparticleETA;

  TH1 *fconstituentsPT;
  TH1 *fconstituentsETA;

  TH1 *fbmatchedjetDeltaR;
  TH1 *fVBFmatchedjetDeltaR;

  TH1 *fvtxDeltaT;
  TH1 *fPvtxT;

};

//------------------------------------------------------------------------------

class ExRootResult;
class ExRootTreeReader;

//------------------------------------------------------------------------------

void BookHistogramsBasic(ExRootResult *result, MyPlots *plots)
{
  TPaveText *comment[6];

  TString eta_min[6] = {"0", "1.5", "3", "4", "0", "0"};
  TString eta_max[6] = {"1.5", "3", "4", "5", "5", "10"};

    // genjets
    for (Int_t k = 0; k < 6; k++) {
      TString hname = TString::Format("PU_genjet_pt_bin_%i",k);

      plots->fPUgenjetPT[k] = result->AddHist1D(
        hname, "PU genjet P_{T}",
        "PU genjet P_{T}, GeV", "number of genjets", 50, 0.0, 200.0);

      comment[k] = result->AddComment(0.70, 0.60, 0.95, 0.70);
      comment[k]->AddText(eta_min[k] + " < abs(#eta) <= "+eta_max[k]);

      result->Attach(plots->fPUgenjetPT[k], comment[k]);
    }

    plots->fPUgenjetRho = result->AddHist1D(
      "PU_rho", "PU genjet #rho",
      " #rho [GeV]", "number of events", 50, 0.0, 200.0);

    // genjets constituents
  plots->fPUgenjetConstPT = result->AddHist1D(
    "PU_genjet_constituents_pt", "PU genjet constituents P_{T}",
    "PU genjet constituents P_{T}, GeV", "number of constituents", 50, 0.0, 10.0);
  plots->fPUgenjetConstETA = result->AddHist1D(
    "PU_genjet_constituents_eta", "PU constituents #eta",
    "PU genjet constituents #eta ", "number of constituents", 100, -5.0, 5.0);


    TString name[3] = {"", "_b_matched", "_VBF_matched"};
    for (size_t i = 0; i < 3; i++) {
      plots->fPUgenjetETA[i] = result->AddHist1D(
        "PU_genjet"+name[i]+"_eta", "PU genjet #eta",
        "PU genjet #eta ", "number of genjets", 100, -8.0, 8.0);
        plots->fPUgenjetNCharged[i] = result->AddHist1D(
          "PU_genjet"+name[i]+"_NCharged", "PU genjet NCharged",
          "number of NCharged", "number of jets", 100, 0.0, 200.0);
          plots->fPUgenjetNNeutrals[i] = result->AddHist1D(
            "PU_genjet"+name[i]+"_NNeutrals", "PU genjet NNeutrals",
            "number of NNeutrals", "number of jets", 100, 0.0, 200.0);
      // PU genjets Charged hadrons
      plots->fPUgenjetCHPT[i] = result->AddHist1D(
        "PU_genjet"+name[i]+"_CH_pt", "PU genjet CH P_{T}",
        "PU genjet CH P_{T}, GeV", "number of particles", 50, 0.0, 10.0);
      plots->fPUgenjetCHEta[i] = result->AddHist1D(
        "PU_genjet"+name[i]+"_CH_eta", "PU CH #eta",
        "PU genjet CH #eta ", "number of particles", 100, -5.0, 5.0);
      plots->fPUgenjetCHdT[i] = result->AddHist1D(
        "PU_genjet"+name[i]+"_CH_dT", "CH T - PV T",
        "PU CH T - PV T [ns]", "number of particles", 100, -1, 1);
      plots->fPUgenjetCHdZ[i] = result->AddHist1D(
        "PU_genjet"+name[i]+"_CH_dZ", "CH Z - PV Z",
        "PU CH Z - PV Z [m]", "number of particles", 100, -0.25, 0.25);
    }


  plots->fbmatchedjetDeltaR = result->AddHist1D(
    "b_matched_PU_genjet_deltaR", "b matched jet #Delta R",
    "b matched jet #Delta R", "number of jets", 50, 0.0, 5.0);

  plots->fVBFmatchedjetDeltaR = result->AddHist1D(
    "VBF_matched_PU_genjet_deltaR", "VBF matched jet #Delta R",
    "VBF matched jet #Delta R", "number of jets", 50, 0.0, 5.0);

    plots->fvtxDeltaT = result->AddHist1D(
      "vtx_dT", "vtx T - PV T",
      "vtx T - PV T [ns]", "number of vtx", 100, -1, 1);

      plots->fPvtxT = result->AddHist1D(
        "P_vtx_T", " PV T",
        "PV T [ns]", "number of vtx", 100, -1, 1);
  }

  //------------------------------------------------------------------------------

  void AnalyseEvents(ExRootTreeReader *treeReader, MyPlots *plots)
  {
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    TClonesArray *branchGenJetPU = treeReader->UseBranch("GenJetPU");
    TClonesArray *branchGenJetPURho = treeReader->UseBranch("Rho");
    TClonesArray *branchParticlePU = treeReader->UseBranch("filteredPUParticle");

    TClonesArray *branchVtx = treeReader->UseBranch("Vertex");

    Long64_t allEntries = treeReader->GetEntries();

    cout << "** Chain contains " << allEntries << " events" << endl;

    GenParticle *genparticle;
    Jet *genjet;
    Jet *jet;

    TObject* object;
    GenParticle* particle;
    Candidate* candidate;

    Vertex* vtx;
    Double_t Pvtx_T, Pvtx_Z;

    Long64_t entry;

    Rho *rho;

    Double_t n_eta_min[5] = {0, 1.5, 3, 4, 0};
    Double_t n_eta_max[5] = {1.5, 3, 4, 5, 5};

    Int_t i;
    // Loop over all events
    for(entry = 0; entry < allEntries; ++entry)
    {
      //cout << entry << endl;
      // Load selected branches with data from specified event
      treeReader->ReadEntry(entry);

      int nstableParts = 0;
      int nPUParts = 0;
      std::vector<GenParticle*> genparts;

      //------Analyse GenParticle---------
      if(branchParticle->GetEntriesFast() > 0)
      {
        for (size_t k = 0; k < branchParticle->GetEntriesFast(); k++) {
          genparticle = (GenParticle*) branchParticle->At(k);
          // select b quarks from higgs decay
          if (is_b(genparticle) && abs(genparticle->IsPU)==0){
              genparts.push_back(genparticle);
          }
          // select VBF jets
          if (is_VBF(genparticle) && abs(genparticle->IsPU)==0){
              genparts.push_back(genparticle);
          }
        } // end for genparticle branch
      }// end if genparticle branch not empty

      // PV time and Z
      for(Int_t i = 0; i < branchVtx->GetEntriesFast(); ++i)
      {
        // Take vtx
        vtx = (Vertex*) branchVtx->At(i);
        if (vtx->Index == 0) {
          Pvtx_Z = vtx->Z;
          Pvtx_T = vtx->T;
          plots->fPvtxT->Fill(( Pvtx_T) * 1000000000);
          }
        plots->fvtxDeltaT->Fill((vtx->T - Pvtx_T) * 1000000000);

      }
    //---------Analyse GenJet-------
      if(branchGenJetPU->GetEntriesFast() > 0)
      {
        for (size_t k = 0; k < branchGenJetPU->GetEntriesFast(); k++) {
          genjet = (Jet*) branchGenJetPU->At(k);
          if (genjet->PT > 30) {
          for (size_t l = 0; l < genjet->Constituents.GetEntriesFast(); l++) {
            object = genjet->Constituents.At(l);
            // Check if the constituent is accessible
            if(object == 0) continue;
            if(object->IsA() == GenParticle::Class())
            { particle = (GenParticle*) object;
              plots->fPUgenjetConstPT->Fill(particle->PT);
              plots->fPUgenjetConstETA->Fill(particle->Eta);
              // cout << "eta " << particle->Eta << endl;
              // cout << "pt " << particle->PT << endl;
              // cout << "Charge " << particle->Charge << endl;
              if (particle->Charge != 0) {
                plots->fPUgenjetCHPT[0]->Fill(particle->PT);
                plots->fPUgenjetCHEta[0]->Fill(particle->Eta);
                plots->fPUgenjetCHdT[0]->Fill((particle->T  - Pvtx_T)* 1000000000);
                plots->fPUgenjetCHdZ[0]->Fill((particle->Z  - Pvtx_Z)/1000);
              }
            }
          }
          for (size_t k = 0; k < 5; k++) {
            if(n_eta_min[k] < abs(genjet->Eta)&& abs(genjet->Eta)<=n_eta_max[k]){
              plots->fPUgenjetPT[k]->Fill(genjet->PT);
            }
          }
          plots->fPUgenjetPT[5]->Fill(genjet->PT);

          plots->fPUgenjetETA[0]->Fill(genjet->Eta);
          plots->fPUgenjetNCharged[0]->Fill(genjet->NCharged);
          plots->fPUgenjetNNeutrals[0]->Fill(genjet->NNeutrals);
          // match to b and VBF particles
          particle = get_closest_particle(genjet, genparts, 0.4);
          if (is_b(particle)) {
              //put into list of b jets
              plots->fbmatchedjetDeltaR->Fill(get_distance(genjet, particle));
              plots->fPUgenjetETA[1]->Fill(genjet->Eta);
              plots->fPUgenjetNCharged[1]->Fill(genjet->NCharged);
              plots->fPUgenjetNNeutrals[1]->Fill(genjet->NNeutrals);
              for (size_t l = 0; l < genjet->Constituents.GetEntriesFast(); l++) {
                object = genjet->Constituents.At(l);
                // Check if the constituent is accessible
                if(object == 0) continue;
                if(object->IsA() == GenParticle::Class())
                { particle = (GenParticle*) object;
                  if (particle->Charge != 0) {
                    plots->fPUgenjetCHPT[1]->Fill(particle->PT);
                    plots->fPUgenjetCHEta[1]->Fill(particle->Eta);
                    plots->fPUgenjetCHdT[1]->Fill((particle->T  - Pvtx_T)* 1000000000);
                    plots->fPUgenjetCHdZ[1]->Fill((particle->Z  - Pvtx_Z)/1000);
                  }
                }
              }
            }
          else if(is_VBF(particle)){
            // put into list of VBF jets
            plots->fVBFmatchedjetDeltaR->Fill(get_distance(genjet, particle));
            plots->fPUgenjetETA[2]->Fill(genjet->Eta);
            plots->fPUgenjetNCharged[2]->Fill(genjet->NCharged);
            plots->fPUgenjetNNeutrals[2]->Fill(genjet->NNeutrals);
            for (size_t l = 0; l < genjet->Constituents.GetEntriesFast(); l++) {
              object = genjet->Constituents.At(l);
              // Check if the constituent is accessible
              if(object == 0) continue;
              if(object->IsA() == GenParticle::Class())
              { particle = (GenParticle*) object;
                if (particle->Charge != 0) {
                  plots->fPUgenjetCHPT[2]->Fill(particle->PT);
                  plots->fPUgenjetCHEta[2]->Fill(particle->Eta);
                  plots->fPUgenjetCHdT[2]->Fill((particle->T  - Pvtx_T)* 1000000000);
                  plots->fPUgenjetCHdZ[2]->Fill((particle->Z  - Pvtx_Z)/1000);
                }
              }
            }
          }

          // TODO check vertices of PU genjets (real PU jets from one vertex?)
        } // end pt cut 30 GeV
        }// end PU genjets
      }

      //---------Analyse Rho-------
        if(branchGenJetPURho->GetEntriesFast() > 0)
        {
          for (size_t k = 0; k < branchGenJetPURho->GetEntriesFast(); k++) {
            rho = (Rho*) branchGenJetPURho->At(k);
            plots->fPUgenjetRho->Fill(rho->Rho);
            }
          }

    }// end loop over entries
  }//end void

  //------------------------------------------------------------------------------

  void PrintHistograms(ExRootResult *result, MyPlots *plots)
  {
    result->Print("png");
  }

  //------------------------------------------------------------------------------

  void PU_gen_jets(const char *inputFile)
  {
    gSystem->Load("libDelphes");
    gROOT->ForceStyle();

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
    gSystem->cd("Plots/Only_PU/Constituents/");
    cout << "Print hists "<< endl;

    PrintHistograms(result, plots);

    result->Write("results.root");

    cout << "** Exiting..." << endl;

    delete plots;
    delete result;
    delete treeReader;
    delete chain;
  }
