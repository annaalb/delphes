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

  // plots for all genjet collections (4)
  TH1 *fjetPT[5];
  TH1 *fjetEta[5];

  TH1 *fsignaljetPT[2];
  TH1 *fsignaljetEta[2];
  TH1 *fjetResponse[2][5][6];

  TH1 *fjetDeltaR[2];

  TH1 *fPUjetPT[2];
  TH1 *fPUjetEta[2];
  TH1 *fnPUjet[2][5][6];

  TH1 *fPUjetnSignalTracks[4][5][6];
  TH1 *fPUjetnPUTracks[4][5][6];

  TH1 *fPUjet_matched_PT;
  TH1 *fPUjet_matched_Eta;
  TH1 *fPUjet_not_matched_PT;
  TH1 *fPUjet_not_matched_Eta;

// jet time ( jet branches, pt bins, eta bins)
  TH1 *fjettimeminusPV[4][5][6];
  TH1 *fjetarrivaltimeminusPV[4][5][6];
  TH1 *fjetptweightedtimeminusPV[4][5][6];
  TH1 *fjetRweightedtimeminusPV[4][5][6];

  // jets -> tracks Plots (jet branch (CHS and PUPPI, genjet), signal or PU)
  TH1 *fjetTrackPT[4][3];
  TH1 *fjetTrackEta[4][3];
  // in eta bins
  TH1 *fjetTrackdT[4][3][6];
  TH1 *fjetTrackDZ[4][3][6];
  TH1 *fjetTrackTOuter[4][3][6];
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
  TString jet_branch[4] = {"CHS", "PUPPI", "GenJet", "PUGenJet"};
  TString name[4] = {"PU_CHS", "PU_PUPPI", "GenJet","PUGenJet"};
  TString name2[3] = {"","_signal", "_PU"};
  TString ptbin[5] = {"", "_pt20", "_20pt30", "_30pt50", "_pt50"};
  TString etabin[6] = {"", "_eta01p3", "_eta1p32", "_eta23", "_eta34", "_eta4"};

  // CHS, PUPPI, Genjets, PU genjets
  for (size_t i = 0; i < 4; i++) { // jet collections
    // pt
    plots->fjetPT[i] = result->AddHist1D(
       jet_branch[i]+"_jet_pt",  "  jet P_{T}",
       " jet P_{T}, GeV/c", "number of jets", 100, 0.0, 200.0);
    // eta
    plots->fjetEta[i] = result->AddHist1D(
       jet_branch[i]+"_jet_eta",  " jet #eta",
       " jet #eta", "number of jets", 100, -5.0, 5.0);

       if (i==3) { // PU genjets
       for (size_t p = 0; p < 5; p++) { // pt bins
         for (size_t e = 0; e < 6; e++) { // eta bins
           // # signal and PU tracks (or charged particles for genjets)
           plots->fPUjetnSignalTracks[i][p][e] = result->AddHist1D(
              name[i]+"_jet_n_signal_tracks"+ptbin[p]+etabin[e],  " jet signal tracks",
              " Number of signal tracks", "number of jets", 100, 0.0, 30.0);
           plots->fPUjetnPUTracks[i][p][e] = result->AddHist1D(
             name[i]+"_jet_n_PU_tracks"+ptbin[p]+etabin[e],  " jet PU tracks",
             " Number of PU tracks", "number of jets", 100, 0.0, 100.0);
           }
         }
       }

       if (i<2) { // reco jets only

         // pt
         plots->fsignaljetPT[i] = result->AddHist1D(
            jet_branch[i]+"_signal_jet_pt",  "  jet P_{T}",
            " jet P_{T}, GeV/c", "number of jets", 100, 0.0, 200.0);
         // eta
         plots->fsignaljetEta[i] = result->AddHist1D(
            jet_branch[i]+"_signal_jet_eta",  " jet #eta",
            " jet #eta", "number of jets", 100, -5.0, 5.0);
            // jet response
         for (size_t p = 0; p < 5; p++) { // pt bins
           for (size_t e = 0; e < 6; e++) { // eta bins
             plots->fjetResponse[i][p][e] = result->AddHist1D(
             jet_branch[i]+"_response"+ptbin[p]+etabin[e], "jet P_{T} / genjet P_{T}",
             "jet P_{T} / genjet P_{T}", "number of genparticles", 100, 0.0, 2.0);

             plots->fnPUjet[i][p][e] = result->AddHist1D(
                name[i]+"_jet_n"+ptbin[p]+etabin[e],  " Number of PU jets",
                " number of jets ", "", 100, 0.0, 10.0);

                plots->fPUjetnSignalTracks[i][p][e] = result->AddHist1D(
                   name[i]+"_jet_n_signal_tracks"+ptbin[p]+etabin[e],  " jet signal tracks",
                   " Number of signal tracks", "number of jets", 100, 0.0, 20.0);
                plots->fPUjetnPUTracks[i][p][e] = result->AddHist1D(
                  name[i]+"_jet_n_PU_tracks"+ptbin[p]+etabin[e],  " jet PU tracks",
                  " Number of PU tracks", "number of jets", 100, 0.0, 20.0);

               plots->fjettimeminusPV[i][p][e] = result->AddHist1D(
                  name[i]+"_jet_time_minus_PV"+ptbin[p]+etabin[e],  " jet time",
                  " jet time - PV [ns]", "number of jets", 100, -1.0, 1.0);
               // average arrival time
               plots->fjetarrivaltimeminusPV[i][p][e] = result->AddHist1D(
                  name[i]+"_jet_arrivaltime_minus_PV"+ptbin[p]+etabin[e],  " jet time",
                  " jet arrival time - PV [ns]", "number of jets", 100, -1.0, 20.0);
               // pt weighted initial time
               plots->fjetptweightedtimeminusPV[i][p][e] = result->AddHist1D(
                  name[i]+"_jet_pt_weighted_time_minus_PV"+ptbin[p]+etabin[e],  " jet time",
                  " jet pt weighted time - PV [ns]", "number of jets", 100, -1.0, 1.0);
                 // R weighted initial time
                 plots->fjetRweightedtimeminusPV[i][p][e] = result->AddHist1D(
                    name[i]+"_jet_R_weighted_time_minus_PV"+ptbin[p]+etabin[e],  " jet time",
                    " jet R weighted time - PV [ns]", "number of jets", 100, -1.0, 1.0);
           }
          }

          // ------ PU jets -----
          plots->fjetDeltaR[i] = result->AddHist1D(
             name[i]+"_jet_deltaR",  " jet dR",
             " dR (recojet - genjet)", "number of jets", 200, 0.0, 2.0);
         // pt
         plots->fPUjetPT[i] = result->AddHist1D(
            name[i]+"_jet_pt",  "  jet P_{T}",
            " jet P_{T}, GeV/c", "number of jets", 100, 0.0, 200.0);
         // eta
         plots->fPUjetEta[i] = result->AddHist1D(
            name[i]+"_jet_eta",  " jet #eta",
            " jet #eta", "number of jets", 100, -5.0, 5.0);

          // track timing Plots
          for (size_t k = 0; k < 3; k++) { // split into signal and PU
            plots->fjetTrackPT[i][k] = result->AddHist1D(
               name[i]+name2[k]+"_track_pt",  " track P_{T}",
               " track P_{T}, GeV/c", "number of tracks", 100, 0.0, 10.0);
            plots->fjetTrackEta[i][k] = result->AddHist1D(
               name[i]+name2[k]+"_track_eta",  " track #eta",
               " track #eta", "number of tracks", 100, -5.0, 5.0);
            for (size_t e = 0; e < 6; e++) { // eta regions
            plots->fjetTrackdT[i][k][e] = result->AddHist1D(
               name[i]+name2[k]+"_track_T"+etabin[e],  " track T",
               " track dT [ns]", "number of tracks", 100, -1, 1);
            plots->fjetTrackDZ[i][k][e] = result->AddHist1D(
               name[i]+name2[k]+"_track_DZ"+etabin[e],  " track DZ",
               " track dZ [m]", "number of tracks", 100, -0.01, 0.01);
            plots->fjetTrackTOuter[i][k][e] = result->AddHist1D(
               name[i]+name2[k]+"_track_TOuter"+etabin[e],  " track TOuter",
               " track TOuter [ns]", "number of tracks", 100, 0.0, 25.0);
            } // end eta regions
            } // end signal or PU

       } // end reco jets

   }

   // some checks for CHS jets (match to PUPPI jets)
   // pt
   plots->fPUjet_matched_PT = result->AddHist1D(
      "PU_CHS_jet_matched_pt",  "  jet P_{T}",
      " jet P_{T}, GeV/c", "number of jets", 100, 0.0, 200.0);
   // eta
   plots->fPUjet_matched_Eta = result->AddHist1D(
      "PU_CHS_jet_matched_eta",  " jet #eta",
      " jet #eta", "number of jets", 100, -5.0, 5.0);
      // pt
      plots->fPUjet_not_matched_PT = result->AddHist1D(
         "PU_CHS_jet_not_matched_pt",  "  jet P_{T}",
         " jet P_{T}, GeV/c", "number of jets", 100, 0.0, 200.0);
      // eta
      plots->fPUjet_not_matched_Eta = result->AddHist1D(
         "PU_CHS_jet_not_matched_eta",  " jet #eta",
         " jet #eta", "number of jets", 100, -5.0, 5.0);

  }

  //------------------------------------------------------------------------------

  void AnalyseEvents(ExRootTreeReader *treeReader, MyPlots *plots)
  {

    TClonesArray *branchCHSeflow = treeReader->UseBranch("ParticleFlowCandidateCHS");
    TClonesArray *branchPuppiParticle = treeReader->UseBranch("ParticleFlowCandidate");
    TClonesArray *branchPFeflow = treeReader->UseBranch("ParticleFlowCandidatePF");

    TClonesArray *branchParticle = treeReader->UseBranch("Particle"); // for identification of VBF and b quarks
    TClonesArray *branchfilteredParticle = treeReader->UseBranch("filteredParticle"); // input to genjets
    TClonesArray *branchfilteredPUParticle = treeReader->UseBranch("filteredPUParticle"); // input to PU genjets

    TClonesArray *branchMergerParticle = treeReader->UseBranch("mergerParticle"); // input to tracks

    TClonesArray *branchJetCHS = treeReader->UseBranch("JetCHS");
    TClonesArray *branchJetPUPPI = treeReader->UseBranch("JetPUPPI");
    TClonesArray *branchJetPF = treeReader->UseBranch("FastJetFinderPF");

    TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
    TClonesArray *branchPUGenJet = treeReader->UseBranch("GenJetPU");

    TClonesArray *branchVtx = treeReader->UseBranch("Vertex");

    Int_t nbranches = 4;

    TClonesArray *branchJets[4] = {branchJetCHS, branchJetPUPPI, branchGenJet, branchPUGenJet};

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
    Int_t nPUjet[nbranches][5][6]; // CHS and PUPPI, pt bins, eta bins
    Double_t fakerate;

    Double_t t_jet[nbranches], pt_t_jet[nbranches], dt_jet[nbranches], arrt_jet[nbranches], darrt_jet[nbranches], pt_sum[nbranches];
    Double_t R, R_sum[nbranches], R_t_jet[nbranches];
    Int_t n_tracks[nbranches], n_signal_tracks[nbranches], n_PU_tracks[nbranches];

    Long64_t entry;

    Int_t i;

    Bool_t debug=false;

    // Loop over all events n= 0 -> 0-1000, n=1 -> 1000-2000, ... n=9
    //for(entry = allEntries*(n*1000); entry < allEntries*((n+1)*1000); ++entry)
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
    // --------Analyse jets -- 0 CHS, 1 PUPPI,  2 GenJets, 3 PU genjets ---------
    for (Int_t m = 0; m < nbranches; m++) {
      // get jets
      if(branchJets[m]->GetEntriesFast() > 0)
      {
        Jet* jet;
        for (size_t k = 0; k < branchJets[m]->GetEntriesFast(); k++) {
          jet = (Jet*) branchJets[m]->At(k);
          plots->fjetPT[m]->Fill(jet->PT);
          if (jet->PT > 20) {
            plots->fjetEta[m]->Fill(jet->Eta);
          }
         }
       }
     }

// for each event set n PU jets to zero
for (size_t i = 0; i < 5; i++) {
  for (size_t j = 0; j < 6; j++) {
    nPUjet[0][i][j]=0; // CHS
    nPUjet[1][i][j]=0; // PUPPI
  }
}

    // -------- now only PU genjets -------
    // get jets
    if(branchJets[3]->GetEntriesFast() > 0)
    {
      Jet* jet;
      for (size_t k = 0; k < branchJets[3]->GetEntriesFast(); k++) {
        jet = (Jet*) branchJets[3]->At(k);

        Int_t eta_index;
        Int_t pt_index;
        if (abs(jet->Eta) < 1.3) {eta_index=1;};
        if (abs(jet->Eta) > 1.3 && abs(jet->Eta) < 2) {eta_index=2;};
        if (2 < abs(jet->Eta) && abs(jet->Eta) < 3) {eta_index=3;};
        if (3 < abs(jet->Eta) && abs(jet->Eta) < 4) {eta_index=4;};
        if (abs(jet->Eta) > 4) {eta_index=5;};

        if (abs(jet->PT) < 20) {pt_index=1;};
        if (abs(jet->PT) > 20 && abs(jet->PT) < 30) {pt_index=2;};
        if (abs(jet->PT) > 30 && abs(jet->PT) < 50) {pt_index=3;};
        if (abs(jet->PT) > 50) {pt_index=4;};

        if (jet->PT > 20) {

          n_tracks [3]=0;
          n_signal_tracks [3]=0;
          n_PU_tracks [3]=0;

          for(size_t j = 0; j < jet->Constituents.GetEntriesFast(); ++j)
          {
            object = jet->Constituents.At(j);
            // Check if the constituent is accessible
            if(object == 0){continue;}
            // --------pf candidates----------------
            else if(object->IsA() == GenParticle::Class())
            {
              ++n_tracks[3];
              particle = (GenParticle*) object;
              if (particle->Charge == 1 && particle->IsPU == 1) {
                ++n_PU_tracks[3];
              }
              else if (particle->Charge == 1 && particle->IsPU == 0) {
                ++n_signal_tracks[3];
              }

            }

          }
          if (n_tracks [3]!=0) {
            plots->fPUjetnPUTracks[3][0][0]->Fill(n_PU_tracks[3]);
            plots->fPUjetnSignalTracks[3][0][0]->Fill(n_signal_tracks[3]);
          }
        }
      }
    }


    // -------- now only reco jets --------
    for (Int_t m = 0; m < 2; m++) {
      // get jets
      if(branchJets[m]->GetEntriesFast() > 0)
      {
        Jet* jet;
        Jet* close_genjet;
        for (size_t k = 0; k < branchJets[m]->GetEntriesFast(); k++) {
          jet = (Jet*) branchJets[m]->At(k);
          if (jet->PT > 20) {
          // match to genjets
          close_genjet = get_closest_jet(branchGenJet, jet, 10);
          plots->fjetDeltaR[m]->Fill(get_distance(jet, close_genjet));

          Bool_t signal_matched;
          signal_matched = matched_to_jet(branchGenJet, jet, 0.2, 0); // is there a close genjet?
          if (signal_matched && get_distance(jet, close_genjet)<0.2 && close_genjet->PT>20) { // if there is a genjet with dr > 0.2 and pt > 20
            plots->fsignaljetEta[m]->Fill(jet->Eta);
            plots->fsignaljetPT[m]->Fill(jet->PT);
            // plot jet response in pt bins ( 0-20, 20, 30 , 50 of genjet)
            //and in eta bins ( 0-1.3, 1.3-2.0, 2.0-3.0 , 3-4, >4 of genjet)
            Int_t eta_index;
            Int_t pt_index;
            if (abs(close_genjet->Eta) < 1.3) {eta_index=1;};
            if (abs(close_genjet->Eta) > 1.3 && abs(close_genjet->Eta) < 2) {eta_index=2;};
            if (2 < abs(close_genjet->Eta) && abs(close_genjet->Eta) < 3) {eta_index=3;};
            if (3 < abs(close_genjet->Eta) && abs(close_genjet->Eta) < 4) {eta_index=4;};
            if (abs(close_genjet->Eta) > 4) {eta_index=5;};

            if (abs(close_genjet->PT) < 20) {pt_index=1;};
            if (abs(close_genjet->PT) > 20 && abs(close_genjet->PT) < 30) {pt_index=2;};
            if (abs(close_genjet->PT) > 30 && abs(close_genjet->PT) < 50) {pt_index=3;};
            if (abs(close_genjet->PT) > 50) {pt_index=4;};

            plots->fjetResponse[m][0][0]->Fill(jet->PT/close_genjet->PT);
            plots->fjetResponse[m][0][eta_index]->Fill(jet->PT/close_genjet->PT);
            plots->fjetResponse[m][pt_index][0]->Fill(jet->PT/close_genjet->PT);
            plots->fjetResponse[m][pt_index][eta_index]->Fill(jet->PT/close_genjet->PT);

          } // end matched

          Bool_t matched = true;
          if(!signal_matched) {matched = matched_to_jet(branchGenJet, jet, 0.6, 20);} // is there a close genjet with at least 20 GeV ?
          if (!matched) { // if there is no genjet
           plots->fPUjetPT[m]->Fill(jet->PT);
           plots->fPUjetEta[m]->Fill(jet->Eta);

           Int_t eta_index;
           Int_t pt_index;
           if (abs(jet->Eta) < 1.3) {eta_index=1;};
           if (abs(jet->Eta) > 1.3 && abs(jet->Eta) < 2) {eta_index=2;};
           if (2 < abs(jet->Eta) && abs(jet->Eta) < 3) {eta_index=3;};
           if (3 < abs(jet->Eta) && abs(jet->Eta) < 4) {eta_index=4;};
           if (abs(jet->Eta) > 4) {eta_index=5;};

           if (abs(jet->PT) < 20) {pt_index=1;};
           if (abs(jet->PT) > 20 && abs(jet->PT) < 30) {pt_index=2;};
           if (abs(jet->PT) > 30 && abs(jet->PT) < 50) {pt_index=3;};
           if (abs(jet->PT) > 50) {pt_index=4;};

           //calculate n PU jet per eta bin here
           ++nPUjet[m][0][0];
           ++nPUjet[m][0][eta_index];
           ++nPUjet[m][pt_index][0];
           ++nPUjet[m][pt_index][eta_index];

           // some checks on CHS jets
           Bool_t matched_to_PUPPI_jet;
           matched_to_PUPPI_jet = matched_to_jet(branchJetPUPPI, jet, 0.6, 0); // is there a close PUPPI jet with pt > 30?
           if (m==0) {
             // find closest PUPPI jet
             if (matched_to_PUPPI_jet) {
               plots->fPUjet_matched_PT->Fill(jet->PT);
               plots->fPUjet_matched_Eta->Fill(jet->Eta);
              }
             else if (!matched_to_PUPPI_jet) {
               plots->fPUjet_not_matched_PT->Fill(jet->PT);
               plots->fPUjet_not_matched_Eta->Fill(jet->Eta);
              }
           }

        ++n_jets [m];

        t_jet [m] = 0;
        n_tracks [m]=0;
        n_signal_tracks [m]=0;
        n_PU_tracks [m]=0;

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

              plots->fjetTrackdT[m][0][eta_index]->Fill(deltaT);
              plots->fjetTrackDZ[m][0][eta_index]->Fill(deltaZ);
              plots->fjetTrackTOuter[m][0][eta_index]->Fill(TOuter);

              if (p->IsPU == 0) { // signal tracks
                ++n_signal_tracks [m];
                plots->fjetTrackPT[m][1]->Fill(pf->PT);
                plots->fjetTrackEta[m][1]->Fill(pf->Eta);

                // in eta bins (first inclusive)
                plots->fjetTrackdT[m][1][0]->Fill(deltaT);
                plots->fjetTrackDZ[m][1][0]->Fill(deltaZ);
                plots->fjetTrackTOuter[m][1][0]->Fill(TOuter);

                  plots->fjetTrackdT[m][1][eta_index]->Fill(deltaT);
                  plots->fjetTrackDZ[m][1][eta_index]->Fill(deltaZ);
                  plots->fjetTrackTOuter[m][1][eta_index]->Fill(TOuter);

              }
              else if (p->IsPU == 1) { // PU tracks
                ++n_PU_tracks [m];
                plots->fjetTrackPT[m][2]->Fill(pf->PT);
                plots->fjetTrackEta[m][2]->Fill(pf->Eta);

                // in eta bins (first inclusive)
                plots->fjetTrackdT[m][2][0]->Fill(deltaT);
                plots->fjetTrackDZ[m][2][0]->Fill(deltaZ);
                plots->fjetTrackTOuter[m][2][0]->Fill(TOuter);

                  plots->fjetTrackdT[m][2][eta_index]->Fill(deltaT);
                  plots->fjetTrackDZ[m][2][eta_index]->Fill(deltaZ);
                  plots->fjetTrackTOuter[m][2][eta_index]->Fill(TOuter);

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
          plots->fPUjetnPUTracks[m][0][0]->Fill(n_PU_tracks[m]);
          plots->fPUjetnSignalTracks[m][0][0]->Fill(n_signal_tracks[m]);

        plots->fjetarrivaltimeminusPV [m][0][0]->Fill((darrt_jet [m]/n_tracks [m]) * 1000000000); // average jet time from arrival time to PV
        plots->fjettimeminusPV [m][0][0]->Fill((t_jet [m]/n_tracks [m] - Pvtx_T) * 1000000000 ); // average jet time to PV
        plots->fjetptweightedtimeminusPV [m][0][0]->Fill((pt_t_jet [m]/pt_sum [m] - Pvtx_T) * 1000000000 ); // pt weighted jet time minus P
        plots->fjetRweightedtimeminusPV [m][0][0]->Fill((R_t_jet [m]/R_sum [m] - Pvtx_T) * 1000000000 ); // r weighted jet time minus PV

        Int_t eta_index;
        Int_t pt_index;
        if (abs(jet->Eta) < 1.3) {eta_index=1;};
        if (abs(jet->Eta) > 1.3 && abs(jet->Eta) < 2) {eta_index=2;};
        if (2 < abs(jet->Eta) && abs(jet->Eta) < 3) {eta_index=3;};
        if (3 < abs(jet->Eta) && abs(jet->Eta) < 4) {eta_index=4;};
        if (abs(jet->Eta) > 4) {eta_index=5;};

        if (abs(jet->PT) < 20) {pt_index=1;};
        if (abs(jet->PT) > 20 && abs(jet->PT) < 30) {pt_index=2;};
        if (abs(jet->PT) > 30 && abs(jet->PT) < 50) {pt_index=3;};
        if (abs(jet->PT) > 50) {pt_index=4;};

        plots->fPUjetnPUTracks[m][0][eta_index]->Fill(n_PU_tracks[m]);
        plots->fPUjetnSignalTracks[m][0][eta_index]->Fill(n_signal_tracks[m]);

        plots->fjetarrivaltimeminusPV [m][0][eta_index]->Fill((darrt_jet [m]/n_tracks [m]) * 1000000000); // average jet time from arrival time to PV
        plots->fjettimeminusPV [m][0][eta_index]->Fill((t_jet [m]/n_tracks [m] - Pvtx_T) * 1000000000 ); // average jet time to PV
        plots->fjetptweightedtimeminusPV [m][0][eta_index]->Fill((pt_t_jet [m]/pt_sum [m] - Pvtx_T) * 1000000000 ); // pt weighted jet time minus PV
        plots->fjetRweightedtimeminusPV [m][0][eta_index]->Fill((R_t_jet [m]/R_sum [m] - Pvtx_T) * 1000000000 ); // r weighted jet time minus PV

        plots->fPUjetnPUTracks[m][pt_index][0]->Fill(n_PU_tracks[m]);
        plots->fPUjetnSignalTracks[m][pt_index][0]->Fill(n_signal_tracks[m]);

        plots->fjetarrivaltimeminusPV [m][pt_index][0]->Fill((darrt_jet [m]/n_tracks [m]) * 1000000000); // average jet time from arrival time to PV
        plots->fjettimeminusPV [m][pt_index][0]->Fill((t_jet [m]/n_tracks [m] - Pvtx_T) * 1000000000 ); // average jet time to PV
        plots->fjetptweightedtimeminusPV [m][pt_index][0]->Fill((pt_t_jet [m]/pt_sum [m] - Pvtx_T) * 1000000000 ); // pt weighted jet time minus PV
        plots->fjetRweightedtimeminusPV [m][pt_index][0]->Fill((R_t_jet [m]/R_sum [m] - Pvtx_T) * 1000000000 ); // r weighted jet time minus PV

        plots->fPUjetnPUTracks[m][pt_index][eta_index]->Fill(n_PU_tracks[m]);
        plots->fPUjetnSignalTracks[m][pt_index][eta_index]->Fill(n_signal_tracks[m]);

        plots->fjetarrivaltimeminusPV [m][pt_index][eta_index]->Fill((darrt_jet [m]/n_tracks [m]) * 1000000000); // average jet time from arrival time to PV
        plots->fjettimeminusPV [m][pt_index][eta_index]->Fill((t_jet [m]/n_tracks [m] - Pvtx_T) * 1000000000 ); // average jet time to PV
        plots->fjetptweightedtimeminusPV [m][pt_index][eta_index]->Fill((pt_t_jet [m]/pt_sum [m] - Pvtx_T) * 1000000000 ); // pt weighted jet time minus PV
        plots->fjetRweightedtimeminusPV [m][pt_index][eta_index]->Fill((R_t_jet [m]/R_sum [m] - Pvtx_T) * 1000000000 ); // r weighted jet time minus PV

        }

      } // jet pt cut
      } // end loop over jets
    }
    } // end loop over jet branches
  } // end reco jets

// Number of PU jets for CHS and PUPPI
for (size_t m = 0; m < 2; m++) {
  for (size_t i = 0; i < 5; i++) { // pt bins
    for (size_t j = 0; j < 6; j++) { // eta bins
      plots->fnPUjet[m][i][j]->Fill(nPUjet[m][i][j]);
    }
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
    gSystem->cd("PLOTS/10kevents/QCD/macro_QCD_jets_dz_smeared_JES_applied/");
    cout << "Print hists "<< endl;

    PrintHistograms(result, plots);

    result->Write("results.root");

    cout << "** Exiting..." << endl;

    delete plots;
    delete result;
    delete treeReader;
    delete chain;
  }
