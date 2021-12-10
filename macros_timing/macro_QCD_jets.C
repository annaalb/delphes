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
  TH1 *fgenparticlePT;
  TH1 *fgenparticleEta;

  TH1 *fvtxZ;

  // plots for all genjet collections (3)
  TH1 *fjetDeltaR[4];
  TH1 *fjetPT[4];
  TH1 *fjetEta[4];
  TH1 *fjetDeltaEta[4];
  TH1 *fjetDeltaPhi[4];

// jet time ( jet branches, eta bins)
  TH1 *fjettimeminusPV[4][4];
  TH1 *fjetarrivaltimeminusPV[4][4];
  TH1 *fjetptweightedtimeminusPV[4][4];
  TH1 *fjetRweightedtimeminusPV[4][4];

  // jets -> tracks Plots (jet branch (CHS and PUPPI, genjet), signal or PU, eta)
  TH1 *fjetTrackPT[4][3];
  TH1 *fjetTrackEta[4][3];
  // in eta bins
  TH1 *fjetTrackdT[4][3][4];
  TH1 *fjetTrackDZ[4][3][4];
  TH1 *fjetTrackTOuter[4][3][4];
};

//------------------------------------------------------------------------------

class ExRootResult;
class ExRootTreeReader;

//------------------------------------------------------------------------------

void BookHistogramsBasic(ExRootResult *result, MyPlots *plots)
{
  TLegend *legend;

  THStack *stack[7][4];

  THStack *stack_PU_1;
  THStack *stack_PU_2;
  THStack *stack_PU_3;
  THStack *stack_PU_4;

  THStack *stack_track_dz_signal[2];
  THStack *stack_track_dz_PU[2];
  THStack *stack_track_dt_signal[2];
  THStack *stack_track_dt_PU[2];

// PV control plots
    plots->fvtxZ = result->AddHist1D(
      "Pvtx_Z", " PV Z",
      " PV Z [m]", "number of PV", 100, -0.1, 0.1);
  //------genparticles------
  //pt
  plots->fgenparticlePT = result->AddHist1D(
    "genparticle_pt", "genparticle P_{T}",
    "genparticle P_{T}, GeV/c", "number of genparticles", 100, 0.0, 200.0);
  //eta
  plots->fgenparticleEta = result->AddHist1D(
    "genparticle_eta", "genparticle #eta",
    "genparticle #eta", "number of genparticles", 100, -5.0, 5.0);

  //-----------------------------------------------------------------
  //---------------VBF matched jets----------
  //------------------------------------------------------------------
  TString name[3] = {"PF", "PUPPI", "GenJet"};
  TString name2[3] = {"","_signal", "_PU"};
  TString eta[4] = {"", "_eta02", "_eta23", "_eta34"};

  for (size_t i = 0; i < 1; i++) { // jet collections
    // average initial time
    for (size_t e = 0; e < 4; e++) { // eta regions
    plots->fjettimeminusPV[i][e] = result->AddHist1D(
       name[i]+"_jet_time_minus_PV"+eta[e],  " jet time",
       " jet time - PV [ns]", "number of jets", 100, -1.0, 1.0);
    // average arrival time
    plots->fjetarrivaltimeminusPV[i][e] = result->AddHist1D(
       name[i]+"_jet_arrivaltime_minus_PV"+eta[e],  " jet time",
       " jet arrival time - PV [ns]", "number of jets", 100, -1.0, 20.0);
    // pt weighted initial time
    plots->fjetptweightedtimeminusPV[i][e] = result->AddHist1D(
       name[i]+"_jet_pt_weighted_time_minus_PV"+eta[e],  " jet time",
       " jet pt weighted time - PV [ns]", "number of jets", 100, -1.0, 1.0);
      // R weighted initial time
      plots->fjetRweightedtimeminusPV[i][e] = result->AddHist1D(
         name[i]+"_jet_R_weighted_time_minus_PV"+eta[e],  " jet time",
         " jet R weighted time - PV [ns]", "number of jets", 100, -1.0, 1.0);
    }

    // pt
    plots->fjetPT[i] = result->AddHist1D(
       name[i]+"_jet_pt",  "  jet P_{T}",
       " jet P_{T}, GeV/c", "number of jets", 100, 0.0, 200.0);
    // eta
    plots->fjetEta[i] = result->AddHist1D(
       name[i]+"_jet_eta",  " jet #eta",
       " jet #eta", "number of jets", 100, -5.0, 5.0);

    // track timing Plots
    //if (i<2 ) { // skip for GenJets and b matched
    for (size_t k = 0; k < 3; k++) { // split into signal and PU
      plots->fjetTrackPT[i][k] = result->AddHist1D(
         name[i]+name2[k]+"_track_pt",  " track P_{T}",
         " track P_{T}, GeV/c", "number of tracks", 100, 0.0, 10.0);
      plots->fjetTrackEta[i][k] = result->AddHist1D(
         name[i]+name2[k]+"_track_eta",  " track #eta",
         " track #eta", "number of tracks", 100, -5.0, 5.0);
      for (size_t e = 0; e < 4; e++) { // eta regions
      plots->fjetTrackdT[i][k][e] = result->AddHist1D(
         name[i]+name2[k]+"_track_T"+eta[e],  " track T",
         " track dT [ns]", "number of tracks", 100, -1, 1);
      plots->fjetTrackDZ[i][k][e] = result->AddHist1D(
         name[i]+name2[k]+"_track_DZ"+eta[e],  " track DZ",
         " track dZ [m]", "number of tracks", 100, -0.01, 0.01);
      plots->fjetTrackTOuter[i][k][e] = result->AddHist1D(
         name[i]+name2[k]+"_track_TOuter"+eta[e],  " track TOuter",
         " track TOuter [ns]", "number of tracks", 100, 0.0, 25.0);
      } // end eta regions
      } // end signal or PU
    //}
  } // end jet collections

  }

  //------------------------------------------------------------------------------

  void AnalyseEvents(ExRootTreeReader *treeReader, MyPlots *plots)
  {

    TClonesArray *branchCHSeflow = treeReader->UseBranch("ParticleFlowCandidateCHS");
    TClonesArray *branchPuppiParticle = treeReader->UseBranch("ParticleFlowCandidate");
    TClonesArray *branchPFeflow = treeReader->UseBranch("ParticleFlowCandidatePF");

    TClonesArray *branchParticle = treeReader->UseBranch("Particle"); // for identification of VBF and b quarks
    TClonesArray *branchfilteredParticle = treeReader->UseBranch("filteredParticle"); // input to genjets

    TClonesArray *branchMergerParticle = treeReader->UseBranch("mergerParticle"); // input to tracks

    TClonesArray *branchJetCHS = treeReader->UseBranch("JetCHS");
    TClonesArray *branchJetPUPPI = treeReader->UseBranch("JetPUPPI");
    TClonesArray *branchJetPF = treeReader->UseBranch("FastJetFinderPF");

    TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");

    TClonesArray *branchVtx = treeReader->UseBranch("Vertex");

    Int_t nbranches = 3;

    TClonesArray *branchJets[3] = {branchJetPF, branchJetPUPPI, branchGenJet};

    Long64_t allEntries = treeReader->GetEntries();

    cout << "** Chain contains " << allEntries << " events" << endl;

    GenParticle *genparticle;
    Jet *genjet;

    ParticleFlowCandidate *pf;

    TObject* object;
    TObject* object2;

    GenParticle* particle;
    GenParticle* constituent;
    Track* track;
    Tower* tower;
    Vertex* vtx;
    Double_t Pvtx_T, Pvtx_Z;
    Double_t TOF, deltaT, deltaZ, TOuter;

    Int_t n_jets[nbranches];
    Double_t fakerate;

    Double_t t_jet[nbranches], pt_t_jet[nbranches], dt_jet[nbranches], arrt_jet[nbranches], darrt_jet[nbranches], pt_sum[nbranches];
    Double_t R, R_sum[nbranches], R_t_jet[nbranches];
    Int_t n_tracks[nbranches], n_signal_tracks[nbranches], n_PU_tracks[nbranches];
    Double_t mistag_rate[nbranches];

    Long64_t entry;

    Int_t i;

    Bool_t debug=false;

    // Loop over all events
    for(entry = 0; entry < allEntries; ++entry)
    {
      cout << "Process event " << entry+1 << " of total " << allEntries << endl;
      if(debug){cout << "-------- begin event ---------"<< endl;}
      std::vector<GenParticle*> genparts;

      // Load selected branches with data from specified event
      treeReader->ReadEntry(entry);
      if(debug){cout << "-------- begin analyse genparticles ---------"<< endl;}
      //------Analyse GenParticle---------
      if(branchParticle->GetEntriesFast() > 0)
      {
        for (size_t k = 0; k < branchParticle->GetEntriesFast(); k++) {
          genparticle = (GenParticle*) branchParticle->At(k);
           plots->fgenparticlePT->Fill(genparticle->PT);
           plots->fgenparticleEta->Fill(genparticle->Eta);

        } // end for genparticle branch
      }// end if genparticle branch not empty
      if(debug){cout << "-------- begin analyse vertices ---------"<< endl;}
      // -------------- Analyse vertex ---------------------
      // PV time and Z
      for(Int_t i = 0; i < branchVtx->GetEntriesFast(); ++i)
      {
        // Take vtx
        vtx = (Vertex*) branchVtx->At(i);
        if (vtx->Index == 0) {
          Pvtx_Z = vtx->Z;
          Pvtx_T = vtx->T;
          plots->fvtxZ->Fill(Pvtx_Z/1000);
          }
      }
      if(debug){cout << "-------- begin analyse jets -> put them into categories ---------"<< endl;}
    // --------Analyse jets -- 0 CHS, 1 PUPPI,  2 GenJets ---------
    for (Int_t m = 0; m < 1; m++) {
      // get jets
      if(branchJets[m]->GetEntriesFast() > 0)
      {
        Jet* jet;
        for (size_t k = 0; k < branchJets[m]->GetEntriesFast(); k++) {
          jet = (Jet*) branchJets[m]->At(k);
          if(jet->PT>30){
           plots->fjetPT[m]->Fill(jet->PT);
           plots->fjetEta[m]->Fill(jet->Eta);

        ++n_jets [m];

        t_jet [m] = 0;
        n_tracks [m]=0;
        n_signal_tracks [m]=0;
        n_PU_tracks [m]=0;
        mistag_rate [m]=0;

        pt_sum [m]=0;
        pt_t_jet [m] = 0; // pt weighted time
        R = 0;
        R_sum [m] = 0;
        R_t_jet [m] = 0;
        dt_jet [m] = 0;
        arrt_jet [m] = 0;
        darrt_jet [m] = 0;
        // Loop over all jet's constituents
        if(debug){cout << "-------- fill plots: loop over jet constituents ---------"<< endl;}

        for(size_t j = 0; j < jet->Constituents.GetEntriesFast(); ++j)
        {
          object = jet->Constituents.At(j);
          // Check if the constituent is accessible
          if(object == 0){continue;}
          // --------particle candidates (GenJets)----------------
          else if(object->IsA() == GenParticle::Class())
          {
            if(debug){cout << "-------- fill plots: jet constituent is genparticle ---------"<< endl;}
            particle = (GenParticle*) object;
            deltaT = (particle->T - Pvtx_T) * 1000000000;
            deltaZ = (particle->Z - Pvtx_Z) / 1000;
            if (particle->Charge != 0) { // charged particles
              if(debug){cout << "-------- charged genparticle ---------"<< endl;}
              plots->fjetTrackPT[m][0]->Fill(particle->PT);
              plots->fjetTrackEta[m][0]->Fill(particle->Eta);
              // in eta bins (first inclusive)
              plots->fjetTrackdT[m][0][0]->Fill(deltaT);
              plots->fjetTrackDZ[m][0][0]->Fill(deltaZ);
              if (2 > abs(jet->Eta) ) {
                plots->fjetTrackdT[m][0][1]->Fill(deltaT);
                plots->fjetTrackDZ[m][0][1]->Fill(deltaZ);
              }
              if (2 < abs(jet->Eta) && abs(jet->Eta) < 3) {
                plots->fjetTrackdT[m][0][2]->Fill(deltaT);
                plots->fjetTrackDZ[m][0][2]->Fill(deltaZ);
              }
              if (3 < abs(jet->Eta) && abs(jet->Eta) < 4) {
                plots->fjetTrackdT[m][0][3]->Fill(deltaT);
                plots->fjetTrackDZ[m][0][3]->Fill(deltaZ);
              }

              if (particle->IsPU == 0) { // signal genjets
                ++n_signal_tracks [m];
                if(debug) {cout << "--- process Signal genjets ---" << endl;};
              }
              else if (particle->IsPU == 1) { // PU genjets
                ++n_PU_tracks [m];
                if(debug) {cout << "--- process PU genjets ---" << endl;};
              }
              // for the jet time
              t_jet [m] += (particle->T ); // mean jet time
              pt_t_jet [m] += (particle->T * particle->PT); // pt weighted time
              pt_sum [m] += (particle->PT); // summed pt of tracks
              R = get_distance(jet, particle);
              R_sum [m] += R;
              R_t_jet [m] += (particle->T * R); // dr weighted time

              dt_jet [m] += (particle->T - Pvtx_T);

              ++n_tracks [m];
            } // end charged pf
          } // --------END particle candidates----------------
          // --------pf candidates----------------
          else if(object->IsA() == ParticleFlowCandidate::Class())
          {
            if(debug){cout << "-------- fill plots: jet constituent is pf candidate ---------"<< endl;}

            pf = (ParticleFlowCandidate*) object;
            deltaT = (pf->T - Pvtx_T) * 1000000000;
            deltaZ = (pf->Z - Pvtx_Z) / 1000;
            TOuter = pf->TOuter * 1000000000;
            TOF = TOuter - (pf->T * 1000000000);
            if (pf->Charge != 0) { // charged pf -> tracks
              GenParticle *p = static_cast<GenParticle*>(pf->Particles.At(0));
              plots->fjetTrackPT[m][0]->Fill(pf->PT);
              plots->fjetTrackEta[m][0]->Fill(pf->Eta);

              // in eta bins (first inclusive)
              plots->fjetTrackdT[m][0][0]->Fill(deltaT);
              plots->fjetTrackDZ[m][0][0]->Fill(deltaZ);
              plots->fjetTrackTOuter[m][0][0]->Fill(TOuter);
              if (2 > abs(jet->Eta) ) {
                plots->fjetTrackdT[m][0][1]->Fill(deltaT);
                plots->fjetTrackDZ[m][0][1]->Fill(deltaZ);
                plots->fjetTrackTOuter[m][0][1]->Fill(TOuter);
              }
              if (2 < abs(jet->Eta) && abs(jet->Eta) < 3) {
                plots->fjetTrackdT[m][0][2]->Fill(deltaT);
                plots->fjetTrackDZ[m][0][2]->Fill(deltaZ);
                plots->fjetTrackTOuter[m][0][2]->Fill(TOuter);
              }
              if (3 < abs(jet->Eta) && abs(jet->Eta) < 4) {
                plots->fjetTrackdT[m][0][3]->Fill(deltaT);
                plots->fjetTrackDZ[m][0][3]->Fill(deltaZ);
                plots->fjetTrackTOuter[m][0][3]->Fill(TOuter);
              }

              if (p->IsPU == 0) { // signal tracks
                ++n_signal_tracks [m];
                plots->fjetTrackPT[m][1]->Fill(pf->PT);
                plots->fjetTrackEta[m][1]->Fill(pf->Eta);

                // in eta bins (first inclusive)
                plots->fjetTrackdT[m][1][0]->Fill(deltaT);
                plots->fjetTrackDZ[m][1][0]->Fill(deltaZ);
                plots->fjetTrackTOuter[m][1][0]->Fill(TOuter);
                if (2 > abs(jet->Eta) ) {
                  plots->fjetTrackdT[m][1][1]->Fill(deltaT);
                  plots->fjetTrackDZ[m][1][1]->Fill(deltaZ);
                  plots->fjetTrackTOuter[m][1][1]->Fill(TOuter);
                }
                if (2 < abs(jet->Eta) && abs(jet->Eta) < 3) {
                  plots->fjetTrackdT[m][1][2]->Fill(deltaT);
                  plots->fjetTrackDZ[m][1][2]->Fill(deltaZ);
                  plots->fjetTrackTOuter[m][1][2]->Fill(TOuter);
                }
                if (3 < abs(jet->Eta) && abs(jet->Eta) < 4) {
                  plots->fjetTrackdT[m][1][3]->Fill(deltaT);
                  plots->fjetTrackDZ[m][1][3]->Fill(deltaZ);
                  plots->fjetTrackTOuter[m][1][3]->Fill(TOuter);
                }
              }
              else if (p->IsPU == 1) { // PU tracks
                ++n_PU_tracks [m];
                plots->fjetTrackPT[m][2]->Fill(pf->PT);
                plots->fjetTrackEta[m][2]->Fill(pf->Eta);

                // in eta bins (first inclusive)
                plots->fjetTrackdT[m][2][0]->Fill(deltaT);
                plots->fjetTrackDZ[m][2][0]->Fill(deltaZ);
                plots->fjetTrackTOuter[m][2][0]->Fill(TOuter);
                if (2 > abs(jet->Eta) ) {
                  plots->fjetTrackdT[m][2][1]->Fill(deltaT);
                  plots->fjetTrackDZ[m][2][1]->Fill(deltaZ);
                  plots->fjetTrackTOuter[m][2][1]->Fill(TOuter);
                }
                if (2 < abs(jet->Eta) && abs(jet->Eta) < 3) {
                  plots->fjetTrackdT[m][2][2]->Fill(deltaT);
                  plots->fjetTrackDZ[m][2][2]->Fill(deltaZ);
                  plots->fjetTrackTOuter[m][2][2]->Fill(TOuter);
                }
                if (3 < abs(jet->Eta) && abs(jet->Eta) < 4) {
                  plots->fjetTrackdT[m][2][3]->Fill(deltaT);
                  plots->fjetTrackDZ[m][2][3]->Fill(deltaZ);
                  plots->fjetTrackTOuter[m][2][3]->Fill(TOuter);
                }
              }
              // for the jet time
              t_jet [m] += (pf->T ); // mean jet time
              pt_t_jet [m] += (pf->T * pf->PT); // pt weighted time
              pt_sum [m] += (pf->PT); // summed pt of tracks
              R = get_distance(jet, pf);
              R_sum [m] += 1/R;
              R_t_jet [m] += (pf->T * 1/R); // dr weighted time

              dt_jet [m] += (pf->T - Pvtx_T);
              arrt_jet [m] += (pf->TOuter);
              darrt_jet [m] += (pf->TOuter - Pvtx_T);

              ++n_tracks [m];
            } // end charged pf
          } // --------END pf candidates----------------
        } // END JET CONSTITUENTS

        if(debug){cout << "-------- fill plots: finished jet constituents ---------"<< endl;}

        if (n_tracks [m]!=0) {
          // signal efficiency and mistag rate of jets
          mistag_rate [m]=n_PU_tracks [m]/n_tracks [m];

        plots->fjetarrivaltimeminusPV [m][0]->Fill((darrt_jet [m]/n_tracks [m]) * 1000000000); // average jet time from arrival time to PV
        plots->fjettimeminusPV [m][0]->Fill((t_jet [m]/n_tracks [m] - Pvtx_T) * 1000000000 ); // average jet time to PV
        plots->fjetptweightedtimeminusPV [m][0]->Fill((pt_t_jet [m]/pt_sum [m] - Pvtx_T) * 1000000000 ); // pt weighted jet time minus P
        plots->fjetRweightedtimeminusPV [m][0]->Fill((R_t_jet [m]/R_sum [m] - Pvtx_T) * 1000000000 ); // r weighted jet time minus PV
        if (abs(jet->Eta) < 2) {
          plots->fjetarrivaltimeminusPV [m][1]->Fill((darrt_jet [m]/n_tracks [m]) * 1000000000); // average jet time from arrival time to PV
          plots->fjettimeminusPV [m][1]->Fill((t_jet [m]/n_tracks [m] - Pvtx_T) * 1000000000 ); // average jet time to PV
          plots->fjetptweightedtimeminusPV [m][1]->Fill((pt_t_jet [m]/pt_sum [m] - Pvtx_T) * 1000000000 ); // pt weighted jet time minus PV
          plots->fjetRweightedtimeminusPV [m][1]->Fill((R_t_jet [m]/R_sum [m] - Pvtx_T) * 1000000000 ); // r weighted jet time minus PV
        }
        if (2 < abs(jet->Eta) && abs(jet->Eta) < 3) {
          plots->fjetarrivaltimeminusPV [m][2]->Fill((darrt_jet [m]/n_tracks [m]) * 1000000000); // average jet time from arrival time to PV
          plots->fjettimeminusPV [m][2]->Fill((t_jet [m]/n_tracks [m] - Pvtx_T) * 1000000000 ); // average jet time to PV
          plots->fjetptweightedtimeminusPV [m][2]->Fill((pt_t_jet [m]/pt_sum [m] - Pvtx_T) * 1000000000 ); // pt weighted jet time minus PV
          plots->fjetRweightedtimeminusPV [m][2]->Fill((R_t_jet [m]/R_sum [m] - Pvtx_T) * 1000000000 ); // r weighted jet time minus PV
        }
        if (3 < abs(jet->Eta) && abs(jet->Eta) < 4) {
          plots->fjetarrivaltimeminusPV [m][3]->Fill((darrt_jet [m]/n_tracks [m]) * 1000000000); // average jet time from arrival time to PV
          plots->fjettimeminusPV [m][3]->Fill((t_jet [m]/n_tracks [m] - Pvtx_T) * 1000000000 ); // average jet time to PV
          plots->fjetptweightedtimeminusPV [m][3]->Fill((pt_t_jet [m]/pt_sum [m] - Pvtx_T) * 1000000000 ); // pt weighted jet time minus PV
          plots->fjetRweightedtimeminusPV [m][3]->Fill((R_t_jet [m]/R_sum [m] - Pvtx_T) * 1000000000 ); // r weighted jet time minus PV
        }
        }

      } // end loop over jets
    }
    } // end loop over jet branches
  }

    }// end loop over entries

  }//end void

  //------------------------------------------------------------------------------

  void PrintHistograms(ExRootResult *result, MyPlots *plots)
  {
    result->Print("png");
  }

  //------------------------------------------------------------------------------

  void macro_QCD_jets(const char *inputFile)
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
    gSystem->cd("PLOTS/10kevents/QCD/macro_QCD_jets/");
    cout << "Print hists "<< endl;

    PrintHistograms(result, plots);

    result->Write("results.root");

    cout << "** Exiting..." << endl;

    delete plots;
    delete result;
    delete treeReader;
    delete chain;
  }
