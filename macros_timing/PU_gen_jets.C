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
  TH1 *fPUgenjetETA;
  TH1 *fPUgenjetNCharged;
  TH1 *fPUgenjetNNeutrals;

  TH1 *fPUgenjetConstPT;
  TH1 *fPUgenjetConstETA;

  TH1 *fPUgenjetCHPT;
  TH1 *fPUgenjetCHEta;
  TH1 *fPUgenjetCHT;
  TH1 *fPUgenjetCHZ;
  TH1 *fPUgenjetCHdT;
  TH1 *fPUgenjetCHdZ;
  TH1 *fPUgenjetCHPosT;

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

  plots->fPUgenjetETA = result->AddHist1D(
    "PU_genjet_eta", "PU genjet #eta",
    "PU genjet #eta ", "number of genjets", 100, -10.0, 10.0);

    plots->fPUgenjetNCharged = result->AddHist1D(
      "PU_genjet_NCharged", "PU genjet NCharged",
      "number of NCharged", "number of jets", 100, 0.0, 200.0);
      plots->fPUgenjetNNeutrals = result->AddHist1D(
        "PU_genjet_NNeutrals", "PU genjet NNeutrals",
        "number of NNeutrals", "number of jets", 100, 0.0, 200.0);
    // genjets constituents
  plots->fPUgenjetConstPT = result->AddHist1D(
    "PU_genjet_constituents_pt", "PU genjet constituents P_{T}",
    "PU genjet constituents P_{T}, GeV", "number of constituents", 50, 0.0, 10.0);
  plots->fPUgenjetConstETA = result->AddHist1D(
    "PU_genjet_constituents_eta", "PU constituents #eta",
    "PU genjet constituents #eta ", "number of constituents", 100, -5.0, 5.0);
    // PU genjets Carged hadrons
    plots->fPUgenjetCHPT = result->AddHist1D(
      "PU_genjet_CH_pt", "PU genjet CH P_{T}",
      "PU genjet CH P_{T}, GeV", "number of particles", 50, 0.0, 10.0);
    plots->fPUgenjetCHEta = result->AddHist1D(
      "PU_genjet_CH_eta", "PU CH #eta",
      "PU genjet CH #eta ", "number of particles", 100, -5.0, 5.0);
    plots->fPUgenjetCHT = result->AddHist1D(
      "PU_genjet_CH_T", "CH T - PV T",
      "PU CH T  [ns]", "number of particles", 100, -1, 1);
    plots->fPUgenjetCHZ = result->AddHist1D(
      "PU_genjet_CH_Z", "track Z - vtx Z",
      "PU CH Z  [m]", "number of particles", 100, -0.25, 0.25);

      plots->fPUgenjetCHdT = result->AddHist1D(
        "PU_genjet_CH_dT", "CH T - PV T",
        "PU CH T - PV T [ns]", "number of particles", 100, -1, 1);
      plots->fPUgenjetCHdZ = result->AddHist1D(
        "PU_genjet_CH_dZ", "track Z - vtx Z",
        "PU CH Z - PV Z [m]", "number of particles", 100, -0.25, 0.25);

      plots->fPUgenjetCHPosT = result->AddHist1D(
        "PU_genjet_CH_Pos_dT", "CH T - PV T",
        "PU CH Position.T - PV T [ns]", "number of particles", 100, -1, 1);
    // genparticles
  plots->fPUgenparticlePT = result->AddHist1D(
    "PU_genparticle_pt", "PU genparticle P_{T}",
    "PU genparticle P_{T}, GeV", "number of genparticle", 50, 0.0, 10.0);
  plots->fPUgenparticleETA = result->AddHist1D(
    "PU_genparticle_eta", "PU genparticle #eta",
    "PU genparticle #eta ", "number of genparticle", 100, -5.0, 5.0);

    // genparticles
  plots->fconstituentsPT = result->AddHist1D(
    "PU_constituents_pt", "PU constituents P_{T}",
    "PU constituents P_{T}, GeV", "number of constituents", 50, 0.0, 10.0);
  plots->fconstituentsETA = result->AddHist1D(
    "PU_constituents_eta", "PU constituents #eta",
    "PU constituents #eta ", "number of constituents", 100, -5.0, 5.0);

  plots->fbmatchedjetDeltaR = result->AddHist1D(
    "b_matched_PU_genjet_deltaR", "b matched jet #Delta R",
    "b matched jet #Delta R", "number of jets", 50, 0.0, 5.0);

  plots->fVBFmatchedjetDeltaR = result->AddHist1D(
    "VBF_matched_PU_genjet_deltaR", "VBF matched jet #Delta R",
    "VBF matched jet #Delta R", "number of jets", 50, 0.0, 5.0);

    // plots->fdeltaT[i] = result->AddHist1D(
    //   "PU_track_T_minus_vtx_T_", "track T - vtx T",
    //   "PU track T - PV T [ns]", "number of tracks", 100, -1, 1);
    // plots->fdeltaZ[i] = result->AddHist1D(
    //   "PU_track_Z_minus_vtx_Z_", "track Z - vtx Z",
    //   "PU track Z - PV Z [m]", "number of tracks", 100, -0.25, 0.25);

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
    // TClonesArray *branchJet = treeReader->UseBranch("Jet");
    // TClonesArray *branchJetPUPPI = treeReader->UseBranch("JetPUPPI");
    // TClonesArray *branchJetAK8 = treeReader->UseBranch("JetAK8");
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    TClonesArray *branchGenJetPU = treeReader->UseBranch("GenJetPU");
    TClonesArray *branchGenJetPURho = treeReader->UseBranch("Rho");
    //TClonesArray *branchParticlePU = treeReader->UseBranch("PUParticle");
    TClonesArray *branchConstituentsPU = treeReader->UseBranch("ConstituentsPU");
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
    for(entry = 0; entry < 1; ++entry)
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
          for (size_t l = 0; l < genjet->Constituents.GetEntriesFast(); l++) {
            object = genjet->Constituents.At(l);
            // Check if the constituent is accessible
            if(object == 0) continue;
            if(object->IsA() == GenParticle::Class())
            {  particle = (GenParticle*) object;
              // plots->fPUgenjetConstPT->Fill(particle->PT);
              // plots->fPUgenjetConstETA->Fill(particle->Eta);
              // cout << "eta " << particle->Eta << endl;
              // cout << "pt " << particle->PT << endl;
              // cout << "Charge " << particle->Charge << endl;
              if (particle->Charge != 0) {
                // plots->fPUgenjetCHPT->Fill(particle->PT);
                // plots->fPUgenjetCHEta->Fill(particle->Eta);
                // plots->fPUgenjetCHT->Fill((particle->T  )* 1000000000);
                // plots->fPUgenjetCHZ->Fill((particle->Z  )/1000);
                // cout << "entry " << entry << endl;
                // //
                // cout << "PV time: " << Pvtx_T*1000000000 << endl;
                // cout << "particle time: " << particle->T*1000000000 << endl;

              }
              // if (particle->IsPU == 0) {
              //   ++nstableParts;
              // }
              // if (particle->IsPU == 1) {
              //   ++nPUParts;
              // }
            }
            if(object->IsA() == Candidate::Class())
            {  candidate = (Candidate*) object;
              if (candidate->Charge != 0) {
                plots->fPUgenjetCHPosT->Fill((candidate->Position.T() - Pvtx_T )* 1000000000);
              }
            }
            // get tracks
          }
          //if (genjet->PT > 20){
          for (size_t k = 0; k < 5; k++) {
            if(n_eta_min[k] < abs(genjet->Eta)&& abs(genjet->Eta)<=n_eta_max[k]){
              plots->fPUgenjetPT[k]->Fill(genjet->PT);
            }
          }
          plots->fPUgenjetPT[5]->Fill(genjet->PT);
          plots->fPUgenjetETA->Fill(genjet->Eta);
          plots->fPUgenjetNCharged->Fill(genjet->NCharged);
          plots->fPUgenjetNNeutrals->Fill(genjet->NNeutrals);
          //}
          particle = get_closest_particle(genjet, genparts, 5);
          if (is_b(particle)) {
              //put into list of b jets
              plots->fbmatchedjetDeltaR->Fill(get_distance(genjet, particle));
            }
          else if(is_VBF(particle)){
            // put into list of VBF jets
            plots->fVBFmatchedjetDeltaR->Fill(get_distance(genjet, particle));
          }
        }
      }

      Int_t ncharged=0;
      Int_t nneutral=0;
      for (size_t l = 0; l < branchConstituentsPU->GetEntriesFast(); l++) {
        object = branchConstituentsPU->At(l);
        // Check if the constituent is accessible
        if(object == 0) continue;
        if(object->IsA() == GenParticle::Class())
        {  particle = (GenParticle*) object;
          plots->fPUgenjetConstPT->Fill(particle->PT);
          plots->fPUgenjetConstETA->Fill(particle->Eta);
          // cout << "eta " << particle->Eta << endl;
          // cout << "pt " << particle->PT << endl;
          // cout << "Charge " << particle->Charge << endl;
          if (particle->Charge != 0 ) {
            plots->fPUgenjetCHPT->Fill(particle->PT);
            plots->fPUgenjetCHEta->Fill(particle->Eta);
            plots->fPUgenjetCHT->Fill((particle->T )* 1000000000);
            plots->fPUgenjetCHZ->Fill((particle->Z )/1000);
            plots->fPUgenjetCHdT->Fill((particle->T  - Pvtx_T)* 1000000000);
            plots->fPUgenjetCHdZ->Fill((particle->Z  - Pvtx_Z)/1000);
            //cout << "CH PID " << particle->PID << endl;
            ++ncharged;
            // cout << "entry " << entry << endl;
            // //
            // cout << "PV time: " << Pvtx_T*1000000000 << endl;
            // cout << "particle time: " << particle->T*1000000000 << endl;

          }
          if (particle->Charge == 0 && particle->PID != 22) {
            cout << "NH PID " << particle->PID << endl;
            ++nneutral;
          }
          if (particle->IsPU == 0) {
            ++nstableParts;
          }
          if (particle->IsPU == 1) {
            ++nPUParts;
          }
        }

      }
      //---------Analyse Rho-------
        if(branchGenJetPURho->GetEntriesFast() > 0)
        {
          for (size_t k = 0; k < branchGenJetPURho->GetEntriesFast(); k++) {
            rho = (Rho*) branchGenJetPURho->At(k);
            plots->fPUgenjetRho->Fill(rho->Rho);
            }
          }

        cout << "Number of stable Parts: "<< nstableParts << endl;
       cout << "Number of PU Parts: "<< nPUParts << endl;
       cout << "Charged hadrons " <<ncharged << endl;
       cout << "Neutral hadrons " <<nneutral << endl;

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
    gSystem->cd("Plots/test_plots/");
    cout << "Print hists "<< endl;

    PrintHistograms(result, plots);

    result->Write("results.root");

    cout << "** Exiting..." << endl;

    delete plots;
    delete result;
    delete treeReader;
    delete chain;
  }
