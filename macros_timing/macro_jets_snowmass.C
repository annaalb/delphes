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

  TH1 *fgenparticlePT[4];
  TH1 *fgenparticleEta[4];

  TH1 *fvtxZ;

  // plots for all genjet collections (3) for VBF and b matched jets and not matched (3)
  TH1 *fmatchedjetDeltaR[4][4];
  TH1 *fmatchedjetPT[4][4];
  TH1 *fmatchedjetEta[4][4];
  TH1 *fmatchedjetDeltaEta[4][4];
  TH1 *fmatchedjetDeltaPhi[4][4];

// jet time ( jet branches, categories, eta bins)
  TH1 *fmatchedjetdt[4][4][6];
  TH1 *fmatchedjettimeminusPV[4][4][6];
  TH1 *fmatchedjetarrivaltimeminusPV[4][4][6];
  TH1 *fmatchedjetptweightedtimeminusPV[4][4][6];
  TH1 *fmatchedjetRweightedtimeminusPV[4][4][6];

  // jets -> tracks Plots (jet branch (CHS and PUPPI, genjet), signal or PU, category)
  TH1 *fmatchedTrackPT[4][3][4];
  TH1 *fmatchedTrackEta[4][3][4];
  // in eta bins
  TH1 *fmatchedTrackdT[4][3][4][6];
  TH1 *fmatchedTrackDZ[4][3][4][6];
  TH1 *fmatchedTrackTOuter[4][3][4][6];

  TH1 *fmatchedTrackZ[4][3][4];
  TH1 *fmatchedTrackTOF[4][3][4];
  TH1 *fmatchedTrackPartT[4][3][3];

  TH1 *fVBFmatchedGenPartPT[3];
  TH1 *fVBFmatchedGenPartEta[3];
  TH1 *fVBFmatchedGenPartT[3];

  TH1 *fnogenjetPT[3];
  TH1 *fnogenjetEta[3];
  TH1 *fmatchedgenjetPT[3];
  TH1 *fmatchedgenjetEta[3];

  // number of PU jets
  TH1 *fnPUjet[2][5][6];
  // number of tracks in jets (jet branch, category, pt, eta)
  TH1 *fPUjetnSignalTracks[3][3][5][6];
  TH1 *fPUjetnPUTracks[3][3][5][6];
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
  plots->fgenparticlePT[0] = result->AddHist1D(
    "genparticle_b_pt", "genparticle b P_{T}",
    "genparticle b P_{T}, GeV/c", "number of genparticles",100, 0.0, 200.0);
  plots->fgenparticlePT[1] = result->AddHist1D(
    "genparticle_VBF_pt", "genparticle VBF P_{T}",
    "genparticle VBF P_{T}, GeV/c", "number of genparticles", 100, 0.0, 200.0);
  plots->fgenparticlePT[2] = result->AddHist1D(
    "genparticle_pt", "genparticle P_{T}",
    "genparticle P_{T}, GeV/c", "number of genparticles", 100, 0.0, 200.0);
  //eta
  plots->fgenparticleEta[0] = result->AddHist1D(
    "genparticle_b_eta", "genparticle b #eta",
    "genparticle b #eta", "number of genparticles", 100, -5.0, 5.0);
  plots->fgenparticleEta[1] = result->AddHist1D(
    "genparticle_VBF_eta", "genparticle VBF #eta",
    "genparticle VBF #eta", "number of genparticles", 100, -5.0, 5.0);
  plots->fgenparticleEta[2] = result->AddHist1D(
    "genparticle_eta", "genparticle #eta",
    "genparticle #eta", "number of genparticles", 100, -5.0, 5.0);

  //-----------------------------------------------------------------
  //---------------VBF matched jets----------
  //------------------------------------------------------------------
  TString name[3] = {"CHS", "PUPPI", "GenJet"};
  TString name2[3] = {"","_signal", "_PU"};
  TString category[3] = {"VBF", "b", "PU"};
  TString eta[4] = {"", "_eta02", "_eta23", "_eta34"};
  TString ptbin[5] = {"", "_pt20", "_20pt30", "_30pt50", "_pt50"};
  TString etabin[6] = {"", "_eta01p3", "_eta1p32", "_eta23", "_eta34", "_eta4"};

  for (size_t m = 0; m < 3; m++) {
    plots->fJetPT[m] = result->AddHist1D(
      "jet_"+name[m]+"_pt_all", "all jet P_{T}",
      " jet P_{T}, GeV/c", "number of jets", 100, 0.0, 200.0);

    plots->fJetEta[m] = result->AddHist1D(
      "jet_"+name[m]+"_eta_all", "all jet #eta",
      " jet #eta", "number of jets", 100, -5.0, 5.0);
  }

  for (size_t b = 0; b < 3; b++) { // begin VBF or b matched
  stack[0][b] = result->AddHistStack(category[b]+"_jet_pt", category[b]+" matched jets P_{T}");
  stack[1][b] = result->AddHistStack(category[b]+"_jet_eta", category[b]+" matched jets #eta");
  stack[2][b] = result->AddHistStack(category[b]+"_jet_deltaR", category[b]+" matched jets #Delta R");
  stack[3][b] = result->AddHistStack(category[b]+"_jet_deltaPhi", category[b]+" matched jets #Delta #Phi");
  stack[4][b] = result->AddHistStack(category[b]+"_jet_deltaEta", category[b]+" matched jets #Delta #Eta");

  // TODO pt bins
  Int_t p = 0;
  for (size_t i = 0; i < 3; i++) { // jet collections (all 3)
    // pt
    plots->fmatchedjetPT[i][b] = result->AddHist1D(
      category[b]+"_"+name[i]+"_jet_pt", category[b]+" matched jet P_{T}",
      category[b]+" matched jet P_{T}, GeV/c", "number of jets", 100, 0.0, 200.0);
    plots->fmatchedjetPT[i][b]->SetLineColor(i+1);
    stack[0][b]->Add(plots->fmatchedjetPT[i][b]);
    // eta
    plots->fmatchedjetEta[i][b] = result->AddHist1D(
      category[b]+"_"+name[i]+"_jet_eta", category[b]+" matched jet #eta",
      category[b]+" matched jet #eta", "number of jets", 100, -5.0, 5.0);
    plots->fmatchedjetEta[i][b]->SetLineColor(i+1);
    stack[1][b]->Add(plots->fmatchedjetEta[i][b]);
    // delta R
    plots->fmatchedjetDeltaR[i][b] = result->AddHist1D(
      category[b]+"_"+name[i]+"_jet_deltaR", category[b]+" matched jet #Delta R",
      category[b]+" matched jet #Delta R", "number of jets", 50, 0.0, 7.0);
    plots->fmatchedjetDeltaR[i][b]->SetLineColor(i+1);
    stack[2][b]->Add(plots->fmatchedjetDeltaR[i][b]);

    // delta phi
    plots->fmatchedjetDeltaPhi[i][b] = result->AddHist1D(
      category[b]+"_"+name[i]+"_jet_deltaPhi", category[b]+" matched jet #Delta #phi",
      category[b]+" matched jet #Delta #phi", "number of jets", 50, 0.0, 1.0);
    plots->fmatchedjetDeltaPhi[i][b]->SetLineColor(i+1);
    stack[3][b]->Add(plots->fmatchedjetDeltaPhi[i][b]);
    // // delta eta
    plots->fmatchedjetDeltaEta[i][b] = result->AddHist1D(
      category[b]+"_"+name[i]+"_jet_deltaEta", category[b]+" matched jet #Delta #eta",
      category[b]+" matched jet #Delta #eta", "number of jets", 50, 0.0, 1.0);
    plots->fmatchedjetDeltaEta[i][b]->SetLineColor(i+1);
    stack[4][b]->Add(plots->fmatchedjetDeltaEta[i][b]);

    for (size_t k = 0; k < 3; k++) { // split into signal and PU
      plots->fmatchedTrackPT[i][k][b] = result->AddHist1D(
        category[b]+"_"+name[i]+name2[k]+"_track_pt", category[b]+" matched track P_{T}",
        category[b]+" matched track P_{T}, GeV/c", "number of tracks", 100, 0.0, 10.0);
      plots->fmatchedTrackEta[i][k][b] = result->AddHist1D(
        category[b]+"_"+name[i]+name2[k]+"_track_eta", category[b]+" matched track #eta",
        category[b]+" matched track #eta", "number of tracks", 100, -5.0, 5.0);
      for (size_t e = 0; e < 6; e++) { // eta regions
      plots->fmatchedTrackdT[i][k][b][e] = result->AddHist1D(
        category[b]+"_"+name[i]+name2[k]+"_track_T"+etabin[e], category[b]+" matched track T",
        category[b]+" matched track dT [ns]", "number of tracks", 100, -1, 1);
      plots->fmatchedTrackDZ[i][k][b][e] = result->AddHist1D(
        category[b]+"_"+name[i]+name2[k]+"_track_DZ"+etabin[e], category[b]+" matched track DZ",
        category[b]+" matched track dZ [m]", "number of tracks", 100, -0.01, 0.01);
      plots->fmatchedTrackTOuter[i][k][b][e] = result->AddHist1D(
        category[b]+"_"+name[i]+name2[k]+"_track_TOuter"+etabin[e], category[b]+" matched track TOuter",
        category[b]+" matched track TOuter [ns]", "number of tracks", 100, 0.0, 25.0);
      } // end eta regions
      plots->fmatchedTrackZ[i][k][b] = result->AddHist1D(
        category[b]+"_"+name[i]+name2[k]+"_track_Z", category[b]+" matched track Z",
        category[b]+" matched track Z [m]", "number of tracks", 100, -0.1, 0.1);
      plots->fmatchedTrackTOF[i][k][b] = result->AddHist1D(
        category[b]+"_"+name[i]+name2[k]+"_track_TOF", category[b]+" matched track TOF",
        category[b]+" matched track TOF [ns]", "number of tracks", 100, 0.0, 25.0);
      } // end signal or PU

      for (size_t e = 0; e < 6; e++) { // eta regions
        for (size_t p = 0; p < 5; p++) { // pt bins
      // # signal and PU tracks
      plots->fPUjetnSignalTracks[i][b][p][e] = result->AddHist1D(
         category[b]+"_"+name[i]+"_jet_n_signal_tracks"+ptbin[p]+etabin[e],  " jet signal tracks",
         " Number of signal tracks", "number of jets", 100, 0.0, 20.0);
      plots->fPUjetnPUTracks[i][b][p][e] = result->AddHist1D(
        category[b]+"_"+name[i]+"_jet_n_PU_tracks"+ptbin[p]+etabin[e],  " jet PU tracks",
        " Number of PU tracks", "number of jets", 100, 0.0, 20.0);
      }

        plots->fmatchedjetdt[i][b][e] = result->AddHist1D(
          category[b]+"_"+name[i]+"_jet_dt_summed"+etabin[e], category[b]+" matched jet time",
          category[b]+" jet dt sum [ns]", "number of jets", 100, -1.0, 1.0);
        plots->fmatchedjettimeminusPV[i][b][e] = result->AddHist1D(
          category[b]+"_"+name[i]+"_jet_time_minus_PV"+etabin[e], category[b]+" matched jet time",
          category[b]+" jet time - PV [ns]", "number of jets", 100, -1.0, 1.0);
        // average arrival time
        plots->fmatchedjetarrivaltimeminusPV[i][b][e] = result->AddHist1D(
          category[b]+"_"+name[i]+"_jet_arrivaltime_minus_PV"+etabin[e], category[b]+" matched jet time",
          category[b]+" jet arrival time - PV [ns]", "number of jets", 100, -1.0, 20.0);
        // pt weighted initial time
        plots->fmatchedjetptweightedtimeminusPV[i][b][e] = result->AddHist1D(
          category[b]+"_"+name[i]+"_jet_pt_weighted_time_minus_PV"+etabin[e], category[b]+" matched jet time",
          category[b]+" jet pt weighted time - PV [ns]", "number of jets", 100, -1.0, 1.0);
          // R weighted initial time
          plots->fmatchedjetRweightedtimeminusPV[i][b][e] = result->AddHist1D(
            category[b]+"_"+name[i]+"_jet_R_weighted_time_minus_PV"+etabin[e], category[b]+" matched jet time",
            category[b]+" jet R weighted time - PV [ns]", "number of jets", 100, -1.0, 1.0);
      }
    // track timing Plots
    if (i<2 ) { // skip for GenJets
    // average initial time

    for (size_t p = 0; p < 5; p++) { // pt bins
      for (size_t e = 0; e < 6; e++) { // eta regions
          // number of PU jets (only for PU jets )
          if (b==2) {
            plots->fnPUjet[i][p][e] = result->AddHist1D(
               category[2]+"_"+name[i]+"_jet_n"+ptbin[p]+etabin[e],  " Number of PU jets",
               " number of jets ", "", 100, 0.0, 10.0);
          }
        }
    }

    } // end only reco jets
  } // end jet collections
} // end VBF or b or not matched

TString charge[3] = {"", "0", "1"};

for (size_t i = 0; i < 3; i++) {
  plots->fVBFmatchedGenPartPT[i] = result->AddHist1D(
    "VBF_GenPart_pT"+charge[i], "VBF matched GenParts pT"+charge[i],
    "VBF matched GenParts P_{T}, GeV/c", "number of particles", 100, 0.0, 10.0);
    plots->fVBFmatchedGenPartEta[i] = result->AddHist1D(
      "VBF_GenPart_eta"+charge[i], "VBF matched GenParts #eta"+charge[i],
      "VBF matched GenParts #eta", "number of particles", 100, -5, 5);
      plots->fVBFmatchedGenPartT[i] = result->AddHist1D(
        "VBF_GenPart_T"+charge[i], "VBF matched GenParts T"+charge[i],
        "VBF matched GenParts dT [ns]", "number of particles", 100, -0.1, 0.1);
      }

  // book legend for stack of 3 histograms
  legend = result->AddLegend(0.25, 0.86, 0.45, 0.98);
  legend->AddEntry(plots->fmatchedjetPT[0][0], "CHS jet", "l");
  legend->AddEntry(plots->fmatchedjetPT[1][0], "PUPPI jet", "l");
  legend->AddEntry(plots->fmatchedjetPT[2][0], "Genjet", "l");

  // attach legend to stacks (legend will be printed over stack in .eps file)
  for (size_t p = 0; p < 5; p++) {
    result->Attach(stack[p][0], legend);
    result->Attach(stack[p][1], legend);
    result->Attach(stack[p][2], legend);
    result->Attach(stack[p][3], legend);
  }


 // plots for unmatched jets (no b or VBF ) that are matched to genjets
  for (size_t g = 0; g < 3; g++) {
    plots->fnogenjetPT[g] = result->AddHist1D(
      "no_gen_jet_pt_"+name[g], "jet P_{T}"+name[g],
      "jet P_{T}, GeV/c", "number of jets", 100, 0.0, 200.0);

    plots->fnogenjetEta[g] = result->AddHist1D(
      "no_gen_jet_eta_"+name[g], "jet #eta"+name[g],
      "jet #eta", "number of jets", 100, -5.0, 5.0);

      plots->fmatchedgenjetPT[g] = result->AddHist1D(
        "matched_gen_jet_pt_"+name[g], "jet P_{T}"+name[g],
        "jet P_{T}, GeV/c", "number of jets", 100, 0.0, 200.0);

      plots->fmatchedgenjetEta[g] = result->AddHist1D(
        "matched_gen_jet_eta_"+name[g], "jet #eta"+name[g],
        "jet #eta", "number of jets", 100, -5.0, 5.0);
     }

  }

  //------------------------------------------------------------------------------

  void AnalyseEvents(ExRootTreeReader *treeReader, MyPlots *plots)
  {

    TClonesArray *branchCHSeflow = treeReader->UseBranch("ParticleFlowCandidateCHS");
    TClonesArray *branchPuppiParticle = treeReader->UseBranch("ParticleFlowCandidate");

    TClonesArray *branchParticle = treeReader->UseBranch("Particle"); // for identification of VBF and b quarks
    TClonesArray *branchfilteredParticle = treeReader->UseBranch("filteredParticle"); // input to genjets

    TClonesArray *branchMergerParticle = treeReader->UseBranch("mergerParticle"); // input to tracks

    TClonesArray *branchJetCHS = treeReader->UseBranch("JetCHS");
    TClonesArray *branchJetPUPPI = treeReader->UseBranch("JetPUPPI");

    TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");

    TClonesArray *branchVtx = treeReader->UseBranch("Vertex");

    Int_t nbranches = 3;

    TClonesArray *branchJets[3] = {branchJetCHS, branchJetPUPPI, branchGenJet};

    Long64_t allEntries = treeReader->GetEntries();

    cout << "** Chain contains " << allEntries << " events" << endl;

    GenParticle *genparticle;
    Jet *genjet;
    Jet *VBF_matched_jet;
    GenParticle *VBF_genpart;
    Jet *b_matched_jet;
    GenParticle *b_genpart;
    Jet *PU_matched_jet;
    GenParticle *PU_genpart;

    Jet *matched_to_genjet = nullptr;

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

    Double_t matching_radius_large = 1;
    Double_t matching_radius = 0.4;
    Double_t gen_matching_radius = 0.2;

    Int_t n_jets[nbranches][3];
    Int_t nPUjet[2][5][6]; // CHS and PUPPI, pt bins, eta bins
    Double_t fakerate;

    Double_t t_jet[nbranches][3], pt_t_jet[nbranches][3], dt_jet[nbranches][3], arrt_jet[nbranches][3], darrt_jet[nbranches][3], pt_sum[nbranches][3];
    Double_t R, R_sum[nbranches][3], R_t_jet[nbranches][3];
    Int_t n_tracks[nbranches][3], n_signal_tracks[nbranches][3], n_PU_tracks[nbranches][3];

    Long64_t entry;

    Int_t i;

    // TODO write a module that stores these jets in a tree
    std::vector<Jet*> matched_jets[nbranches][3]; // 4 branches and 3 categories

    Bool_t debug=false;

    // Loop over all events
    for(entry = 0; entry < 1000; ++entry)
    {
      cout << "Process event " << entry+1 << " of total " << allEntries << endl;
      if(debug){cout << "-------- begin event ---------"<< endl;}
      std::vector<GenParticle*> VBF_genparts;
      std::vector<GenParticle*> b_genparts;
      std::vector<GenParticle*> PU_genparts;
      std::vector<GenParticle*> genparts;

      // Load selected branches with data from specified event
      treeReader->ReadEntry(entry);
      if(debug){cout << "-------- begin analyse genparticles ---------"<< endl;}
      //------Analyse GenParticle---------
      if(branchParticle->GetEntriesFast() > 0)
      {
        for (size_t k = 0; k < branchParticle->GetEntriesFast(); k++) {
          genparticle = (GenParticle*) branchParticle->At(k);

          // select b quarks from higgs decay
          if (is_b(genparticle) && abs(genparticle->IsPU)==0){
              b_genparts.push_back(genparticle);
              genparts.push_back(genparticle);
                plots->fgenparticlePT[0]->Fill(genparticle->PT);
                plots->fgenparticleEta[0]->Fill(genparticle->Eta);
          }
          // select VBF particles
          if (is_VBF(genparticle) && abs(genparticle->IsPU)==0){
              VBF_genparts.push_back(genparticle);
              genparts.push_back(genparticle);
                plots->fgenparticlePT[1]->Fill(genparticle->PT);
                plots->fgenparticleEta[1]->Fill(genparticle->Eta);
          }
           plots->fgenparticlePT[2]->Fill(genparticle->PT);
           plots->fgenparticleEta[2]->Fill(genparticle->Eta);

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

      // for each event set n PU jets to zero
      for (size_t i = 0; i < 5; i++) {
        for (size_t j = 0; j < 6; j++) {
          nPUjet[0][i][j]=0; // CHS
          nPUjet[1][i][j]=0; // PUPPI
        }
      }
    // --------Analyse jets -- 0 CHS, 1 PUPPI,  2 GenJets ---------
    for (Int_t m = 0; m < nbranches; m++) {
      // get jets
      if(branchJets[m]->GetEntriesFast() > 0)
      {
        Jet* jet;
        for (size_t k = 0; k < branchJets[m]->GetEntriesFast(); k++) {
          jet = (Jet*) branchJets[m]->At(k);
          plots->fJetPT[m]->Fill(jet->PT);
          if(jet->PT>20){ // TODO stay with 20 GeV??
           plots->fJetEta[m]->Fill(jet->Eta);

          GenParticle* particle = get_closest_particle(jet, genparts, 0.4);
            //------ b matched ------
            if (is_b(particle) && get_distance(jet, particle)<0.2 && particle->PT > 20) { // b matched if dr < 0.2
              //put into list of b jets [1]
              matched_jets[m][1].push_back(jet);
              plots->fmatchedjetDeltaR[m][1]->Fill(get_distance(jet, particle));
              plots->fmatchedjetDeltaPhi[m][1]->Fill(abs(jet->Phi - particle->Phi));
              plots->fmatchedjetDeltaEta[m][1]->Fill(abs(jet->Eta - particle->Eta));
            }
            //----- VBF matched ------
            else if(is_VBF(particle) && particle->PT > 20){
              matched_jets[m][0].push_back(jet);
              // put into list of VBF jets [0]
              plots->fmatchedjetDeltaR[m][0]->Fill(get_distance(jet, particle));
              plots->fmatchedjetDeltaPhi[m][0]->Fill(abs(jet->Phi - particle->Phi));
              plots->fmatchedjetDeltaEta[m][0]->Fill(abs(jet->Eta - particle->Eta));

              // Loop over all jet's constituents
              for(size_t j = 0; j < jet->Constituents.GetEntriesFast(); ++j)
              {
                object = jet->Constituents.At(j);
                // Check if the constituent is accessible
                if(object == 0){continue;}
                if(object->IsA() == GenParticle::Class())
                {
                  particle = (GenParticle*) object;
                  plots->fVBFmatchedGenPartPT[0]->Fill((particle->PT));
                  plots->fVBFmatchedGenPartEta[0]->Fill((particle->Eta));
                  plots->fVBFmatchedGenPartT[0]->Fill((particle->T - Pvtx_T) * 1000000000);
                  if (particle->Charge == 0) {
                   plots->fVBFmatchedGenPartPT[1]->Fill((particle->PT));
                   plots->fVBFmatchedGenPartEta[1]->Fill((particle->Eta));
                   plots->fVBFmatchedGenPartT[1]->Fill((particle->T - Pvtx_T) * 1000000000);
                  }
                  else {
                   plots->fVBFmatchedGenPartPT[2]->Fill((particle->PT));
                   plots->fVBFmatchedGenPartEta[2]->Fill((particle->Eta));
                   plots->fVBFmatchedGenPartT[2]->Fill((particle->T - Pvtx_T) * 1000000000);
                  }
                }
              } // END JET CONSTITUENTS
            } // end is VBF matched
            // ----- not matched -----
            else{
              // match "not matched jets" to genjets
              matched_to_genjet = get_closest_jet(branchGenJet, jet, 10);
              // create dr plot
              plots->fmatchedjetDeltaR[m][2]->Fill(get_distance(jet, matched_to_genjet));
             Bool_t matched;
             // is there a genjet closer than 0.1?
             matched = matched_to_jet(branchGenJet, jet, 0.1, 20); // is there a close genjet?
             if (matched && matched_to_genjet->PT > 20 && get_distance(jet, matched_to_genjet)<0.1) { // if there is a genjet with pt > 20 GeV closer than 0.1
               plots->fmatchedgenjetPT[m]->Fill(jet->PT);
               plots->fmatchedgenjetEta[m]->Fill(jet->Eta);
             }
             // ----- "PU jets" -------
            else if (!matched && get_distance(jet, matched_to_genjet)>0.6) { // if there is no genjet
              matched_jets[m][2].push_back(jet);


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
             }
           }//end else

          }
        }
      }

    } // end loop over jet branches

    if(debug){cout << "-------- begin analyse jets -> fill plots ---------"<< endl;}
    // Fill hists for all three jet branches in all categories
    for (size_t k = 0; k < 3; k++) { // for each category
      for (size_t m = 0; m < nbranches; m++) { // CHS, PUPPI, genjets
        for (size_t i = 0; i < matched_jets[m][k].size(); i++) {

        Jet* jet = matched_jets[m][k].at(i);
        ++n_jets[m][k];
        plots->fmatchedjetPT[m][k]->Fill(jet->PT);
        plots->fmatchedjetEta[m][k]->Fill(jet->Eta);

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

        t_jet[m][k] = 0;
        n_tracks[m][k]=0;
        n_signal_tracks[m][k]=0;
        n_PU_tracks[m][k]=0;

        pt_sum[m][k]=0;
        pt_t_jet[m][k] = 0; // pt weighted time
        R = 0;
        R_sum[m][k] = 0;
        R_t_jet[m][k] = 0;
        dt_jet[m][k] = 0;
        arrt_jet[m][k] = 0;
        darrt_jet[m][k] = 0;
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
              plots->fmatchedTrackPT[m][0][k]->Fill(particle->PT);
              plots->fmatchedTrackEta[m][0][k]->Fill(particle->Eta);
              plots->fmatchedTrackZ[m][0][k]->Fill((particle->Z) / 1000);
              // in eta bins (first inclusive)
              plots->fmatchedTrackdT[m][0][k][0]->Fill(deltaT);
              plots->fmatchedTrackDZ[m][0][k][0]->Fill(deltaZ);

                plots->fmatchedTrackdT[m][0][k][eta_index]->Fill(deltaT);
                plots->fmatchedTrackDZ[m][0][k][eta_index]->Fill(deltaZ);


              if (particle->IsPU == 0) { // signal genjets
                ++n_signal_tracks[m][k];
                if(debug) {cout << "--- process Signal genjets ---" << endl;};
              }
              else if (particle->IsPU == 1) { // PU genjets
                ++n_PU_tracks[m][k];
                if(debug) {cout << "--- process PU genjets ---" << endl;};
              }
              // for the jet time
              t_jet[m][k] += (particle->T ); // mean jet time
              pt_t_jet[m][k] += (particle->T * particle->PT); // pt weighted time
              pt_sum[m][k] += (particle->PT); // summed pt of tracks
              R = get_distance(jet, particle);
              R_sum[m][k] += R;
              R_t_jet[m][k] += (particle->T * R); // dr weighted time

              dt_jet[m][k] += (particle->T - Pvtx_T);

              ++n_tracks[m][k];
            } // end charged particle
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
              plots->fmatchedTrackPT[m][0][k]->Fill(pf->PT);
              plots->fmatchedTrackEta[m][0][k]->Fill(pf->Eta);
              plots->fmatchedTrackZ[m][0][k]->Fill((pf->Z) / 1000);
              plots->fmatchedTrackTOF[m][0][k]->Fill(TOF);

              // in eta bins (first inclusive)
              plots->fmatchedTrackdT[m][0][k][0]->Fill(deltaT);
              plots->fmatchedTrackDZ[m][0][k][0]->Fill(deltaZ);
              plots->fmatchedTrackTOuter[m][0][k][0]->Fill(TOuter);

                plots->fmatchedTrackdT[m][0][k][eta_index]->Fill(deltaT);
                plots->fmatchedTrackDZ[m][0][k][eta_index]->Fill(deltaZ);
                plots->fmatchedTrackTOuter[m][0][k][eta_index]->Fill(TOuter);

              if (p->IsPU == 0) { // signal tracks
                ++n_signal_tracks[m][k];
                plots->fmatchedTrackPT[m][1][k]->Fill(pf->PT);
                plots->fmatchedTrackEta[m][1][k]->Fill(pf->Eta);
                plots->fmatchedTrackZ[m][1][k]->Fill((pf->Z) / 1000);
                plots->fmatchedTrackTOF[m][1][k]->Fill(TOF);

                // in eta bins (first inclusive)
                plots->fmatchedTrackdT[m][1][k][0]->Fill(deltaT);
                plots->fmatchedTrackDZ[m][1][k][0]->Fill(deltaZ);
                plots->fmatchedTrackTOuter[m][1][k][0]->Fill(TOuter);

                  plots->fmatchedTrackdT[m][1][k][eta_index]->Fill(deltaT);
                  plots->fmatchedTrackDZ[m][1][k][eta_index]->Fill(deltaZ);
                  plots->fmatchedTrackTOuter[m][1][k][eta_index]->Fill(TOuter);

              }
              else if (p->IsPU == 1) { // PU tracks
                ++n_PU_tracks[m][k];
                plots->fmatchedTrackPT[m][2][k]->Fill(pf->PT);
                plots->fmatchedTrackEta[m][2][k]->Fill(pf->Eta);
                plots->fmatchedTrackZ[m][2][k]->Fill((pf->Z) / 1000);
                plots->fmatchedTrackTOF[m][2][k]->Fill(TOF);

                // in eta bins (first inclusive)
                plots->fmatchedTrackdT[m][2][k][0]->Fill(deltaT);
                plots->fmatchedTrackDZ[m][2][k][0]->Fill(deltaZ);
                plots->fmatchedTrackTOuter[m][2][k][0]->Fill(TOuter);

                  plots->fmatchedTrackdT[m][2][k][eta_index]->Fill(deltaT);
                  plots->fmatchedTrackDZ[m][2][k][eta_index]->Fill(deltaZ);
                  plots->fmatchedTrackTOuter[m][2][k][eta_index]->Fill(TOuter);
              }
              // for the jet time
              t_jet[m][k] += (pf->T ); // mean jet time
              pt_t_jet[m][k] += (pf->T * pf->PT); // pt weighted time
              pt_sum[m][k] += (pf->PT); // summed pt of tracks
              R = get_distance(jet, pf);
              R_sum[m][k] += R;
              R_t_jet[m][k] += (pf->T * R); // dr weighted time

              dt_jet[m][k] += (pf->T - Pvtx_T);
              arrt_jet[m][k] += (pf->TOuter);
              darrt_jet[m][k] += (pf->TOuter - Pvtx_T);

              ++n_tracks[m][k];
            } // end charged pf
          } // --------END pf candidates----------------
        } // END JET CONSTITUENTS

        if(debug){cout << "-------- fill plots: finished jet constituents ---------"<< endl;}

        if (n_tracks[m][k]!=0) {
          plots->fPUjetnPUTracks[m][k][0][0]->Fill(n_PU_tracks[m][k]);
          plots->fPUjetnSignalTracks[m][k][0][0]->Fill(n_signal_tracks[m][k]);

          plots->fPUjetnPUTracks[m][k][0][eta_index]->Fill(n_PU_tracks[m][k]);
          plots->fPUjetnSignalTracks[m][k][0][eta_index]->Fill(n_signal_tracks[m][k]);

          plots->fPUjetnPUTracks[m][k][pt_index][0]->Fill(n_PU_tracks[m][k]);
          plots->fPUjetnSignalTracks[m][k][pt_index][0]->Fill(n_signal_tracks[m][k]);

          plots->fPUjetnPUTracks[m][k][pt_index][eta_index]->Fill(n_PU_tracks[m][k]);
          plots->fPUjetnSignalTracks[m][k][pt_index][eta_index]->Fill(n_signal_tracks[m][k]);

        plots->fmatchedjetarrivaltimeminusPV[m][k][0]->Fill((darrt_jet[m][k]/n_tracks[m][k]) * 1000000000); // average jet time from arrival time to PV
        plots->fmatchedjetdt[m][k][0]->Fill((dt_jet[m][k]/n_tracks[m][k]) * 1000000000); // // average dt
        plots->fmatchedjettimeminusPV[m][k][0]->Fill((t_jet[m][k]/n_tracks[m][k] - Pvtx_T) * 1000000000 ); // average jet time to PV
        plots->fmatchedjetptweightedtimeminusPV[m][k][0]->Fill((pt_t_jet[m][k]/pt_sum[m][k] - Pvtx_T) * 1000000000 ); // pt weighted jet time minus P
        plots->fmatchedjetRweightedtimeminusPV[m][k][0]->Fill((R_t_jet[m][k]/R_sum[m][k] - Pvtx_T) * 1000000000 ); // r weighted jet time minus PV

          plots->fmatchedjetarrivaltimeminusPV[m][k][eta_index]->Fill((darrt_jet[m][k]/n_tracks[m][k]) * 1000000000); // average jet time from arrival time to PV
          plots->fmatchedjetdt[m][k][eta_index]->Fill((dt_jet[m][k]/n_tracks[m][k]) * 1000000000); // // average dt
          plots->fmatchedjettimeminusPV[m][k][eta_index]->Fill((t_jet[m][k]/n_tracks[m][k] - Pvtx_T) * 1000000000 ); // average jet time to PV
          plots->fmatchedjetptweightedtimeminusPV[m][k][eta_index]->Fill((pt_t_jet[m][k]/pt_sum[m][k] - Pvtx_T) * 1000000000 ); // pt weighted jet time minus PV
          plots->fmatchedjetRweightedtimeminusPV[m][k][eta_index]->Fill((R_t_jet[m][k]/R_sum[m][k] - Pvtx_T) * 1000000000 ); // r weighted jet time minus PV
        }
      } // end loop over jets
    } // end loop over jet branches
  } // end loop over categories

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

  void macro_jets_snowmass(const char *inputFile)
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
    gSystem->cd("PLOTS/10kevents/VBF_sample/macro_jets_snowmass_dz_smeared_JES_applied/");
    cout << "Print hists "<< endl;

    PrintHistograms(result, plots);

    result->Write("results.root");

    cout << "** Exiting..." << endl;

    delete plots;
    delete result;
    delete treeReader;
    delete chain;
  }
