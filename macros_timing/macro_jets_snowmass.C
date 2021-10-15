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
  TH1 *fJetPT[3];
  TH1 *fJetEta[3];

  TH1 *fgenparticlePT[4];
  TH1 *fgenparticleEta[4];

  TH1 *fvtxZ;

  // plots for all genjet collections (3) for VBF and b matched jets and not matched (3)
  TH1 *fmatchedjetDeltaR[3][3];
  TH1 *fmatchedjetPT[3][3];
  TH1 *fmatchedjetEta[3][3];
  TH1 *fmatchedjetDeltaEta[3][3];
  TH1 *fmatchedjetDeltaPhi[3][3];

  // jets -> tracks Plots
  TH1 *fmatchedTrackPT[2][3][3];
  TH1 *fmatchedTrackEta[2][3][3];
  TH1 *fmatchedTrackT[2][3][3];
  TH1 *fmatchedTrackDZ[2][3][3];
  TH1 *fmatchedTrackZ[2][3][3];

  TH1 *fmatchedTrackPartT[2][3][3];

  TH1 *fVBFmatchedGenPartPT[3];
  TH1 *fVBFmatchedGenPartEta[3];
  TH1 *fVBFmatchedGenPartT[3];

  TH1 *fnogenjetPT[3];
  TH1 *fnogenjetEta[3];
  TH1 *fmatchedgenjetPT[3];
  TH1 *fmatchedgenjetEta[3];
};

//------------------------------------------------------------------------------

class ExRootResult;
class ExRootTreeReader;

//------------------------------------------------------------------------------

void BookHistogramsBasic(ExRootResult *result, MyPlots *plots)
{
  TLegend *legend;

  THStack *stack[7][2];

  THStack *stack_PU_1;
  THStack *stack_PU_2;
  THStack *stack_PU_3;
  THStack *stack_PU_4;

  // book more histograms
  plots->fJetPT[0] = result->AddHist1D(
    "jet_pt_all", "all jet P_{T}",
    "jet P_{T}, GeV/c", "number of jets", 100, 0.0, 200.0);

  plots->fJetEta[0] = result->AddHist1D(
    "jet_eta_all", "all jet #eta",
    "jet #eta", "number of jets", 100, -5.0, 5.0);

  plots->fJetPT[1] = result->AddHist1D(
    "jet_PUPPI_pt_all", "all PUPPI jet P_{T}",
    "PUPPI jet P_{T}, GeV/c", "number of jets", 100, 0.0, 200.0);

  plots->fJetEta[1] = result->AddHist1D(
    "jet_PUPPI_eta_all", "all PUPPI jet #eta",
    "PUPPI jet #eta", "number of jets", 100, -5.0, 5.0);

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

    // genjets
  plots->fJetPT[2] = result->AddHist1D(
    "genjet_pt", "genjet P_{T}",
    "genjet P_{T}, GeV", "number of genjets", 50, 0.0, 200.0);
  plots->fJetEta[2] = result->AddHist1D(
    "genjet_eta", "genjet #eta",
    "genjet #eta ", "number of genjets", 100, -5.0, 5.0);
  //-----------------------------------------------------------------
  //---------------VBF matched jets----------
  //------------------------------------------------------------------
  TString name[3] = {"CHS", "PUPPI", "GenJet"};
  TString name2[3] = {"","_signal", "_PU"};
  TString category[3] = {"VBF", "b", "not"};

  for (size_t b = 0; b < 3; b++) { // begin VBF or b matched
  stack[0][b] = result->AddHistStack(category[b]+"_matched_jet_pt", category[b]+" matched jets P_{T}");
  stack[1][b] = result->AddHistStack(category[b]+"_matched_jet_eta", category[b]+" matched jets #eta");
  stack[2][b] = result->AddHistStack(category[b]+"_matched_jet_deltaR", category[b]+" matched jets #Delta R");
  stack[3][b] = result->AddHistStack(category[b]+"_matched_jet_deltaPhi", category[b]+" matched jets #Delta #Phi");
  stack[4][b] = result->AddHistStack(category[b]+"_matched_jet_deltaEta", category[b]+" matched jets #Delta #Eta");

  for (size_t i = 0; i < 3; i++) { // jet collections
    // pt
    plots->fmatchedjetPT[i][b] = result->AddHist1D(
      category[b]+"_matched_"+name[i]+"_jet_pt", category[b]+" matched jet P_{T}",
      category[b]+" matched jet P_{T}, GeV/c", "number of jets", 100, 0.0, 200.0);
    plots->fmatchedjetPT[i][b]->SetLineColor(i+1);
    stack[0][b]->Add(plots->fmatchedjetPT[i][b]);
    // eta
    plots->fmatchedjetEta[i][b] = result->AddHist1D(
      category[b]+"_matched_"+name[i]+"_jet_eta", category[b]+" matched jet #eta",
      category[b]+" matched jet #eta", "number of jets", 100, -5.0, 5.0);
    plots->fmatchedjetEta[i][b]->SetLineColor(i+1);
    stack[1][b]->Add(plots->fmatchedjetEta[i][b]);
    // delta R
    plots->fmatchedjetDeltaR[i][b] = result->AddHist1D(
      category[b]+"_matched_"+name[i]+"_jet_deltaR", category[b]+" matched jet #Delta R",
      category[b]+" matched jet #Delta R", "number of jets", 50, 0.0, 1.0);
    plots->fmatchedjetDeltaR[i][b]->SetLineColor(i+1);
    stack[2][b]->Add(plots->fmatchedjetDeltaR[i][b]);

    // delta phi
    plots->fmatchedjetDeltaPhi[i][b] = result->AddHist1D(
      category[b]+"_matched_"+name[i]+"_jet_deltaPhi", category[b]+" matched jet #Delta #phi",
      category[b]+" matched jet #Delta #phi", "number of jets", 50, 0.0, 1.0);
    plots->fmatchedjetDeltaPhi[i][b]->SetLineColor(i+1);
    stack[3][b]->Add(plots->fmatchedjetDeltaPhi[i][b]);
    // // delta eta
    plots->fmatchedjetDeltaEta[i][b] = result->AddHist1D(
      category[b]+"_matched_"+name[i]+"_jet_deltaEta", category[b]+" matched jet #Delta #eta",
      category[b]+" matched jet #Delta #eta", "number of jets", 50, 0.0, 1.0);
    plots->fmatchedjetDeltaEta[i][b]->SetLineColor(i+1);
    stack[4][b]->Add(plots->fmatchedjetDeltaEta[i][b]);

    // track timing Plots
    if (i<2 && (b==0 || b==2)) { // skip for GenJets and b matched
    for (size_t k = 0; k < 3; k++) { // split into signal and PU
      plots->fmatchedTrackPT[i][k][b] = result->AddHist1D(
        category[b]+"_matched_"+name[i]+name2[k]+"_track_pt", category[b]+" matched track P_{T}",
        category[b]+" matched track P_{T}, GeV/c", "number of tracks", 100, 0.0, 10.0);
        plots->fmatchedTrackEta[i][k][b] = result->AddHist1D(
          category[b]+"_matched_"+name[i]+name2[k]+"_track_eta", category[b]+" matched track #eta",
          category[b]+" matched track #eta", "number of tracks", 100, -5.0, 5.0);
      plots->fmatchedTrackT[i][k][b] = result->AddHist1D(
        category[b]+"_matched_"+name[i]+name2[k]+"_track_T", category[b]+" matched track T",
        category[b]+" matched track dT [ns]", "number of tracks", 100, -1, 1);
        plots->fmatchedTrackPartT[i][k][b] = result->AddHist1D(
          category[b]+"_matched_"+name[i]+name2[k]+"_track_p_T", "track->p T",
          category[b]+" matched track->p dT [ns]", "number of particles", 100, -1, 1);
          plots->fmatchedTrackDZ[i][k][b] = result->AddHist1D(
            category[b]+"_matched_"+name[i]+name2[k]+"_track_DZ", category[b]+" matched track DZ",
            category[b]+" matched track dZ [m]", "number of tracks", 100, -0.01, 0.01);

            plots->fmatchedTrackZ[i][k][b] = result->AddHist1D(
              category[b]+"_matched_"+name[i]+name2[k]+"_track_Z", category[b]+" matched track Z",
              category[b]+" matched track Z [m]", "number of tracks", 100, -0.1, 0.1);

        } // end signal or PU
      }
  } // end jet collections
} // end VBF or b or not matched

TString charge[3] = {"", "0", "1"};

for (size_t i = 0; i < 3; i++) {
  plots->fVBFmatchedGenPartPT[i] = result->AddHist1D(
    "VBF_matched_GenPart_pT"+charge[i], "VBF matched GenParts pT"+charge[i],
    "VBF matched GenParts P_{T}, GeV/c", "number of particles", 100, 0.0, 10.0);
    plots->fVBFmatchedGenPartEta[i] = result->AddHist1D(
      "VBF_matched_GenPart_eta"+charge[i], "VBF matched GenParts #eta"+charge[i],
      "VBF matched GenParts #eta", "number of particles", 100, -5, 5);
      plots->fVBFmatchedGenPartT[i] = result->AddHist1D(
        "VBF_matched_GenPart_T"+charge[i], "VBF matched GenParts T"+charge[i],
        "VBF matched GenParts dT [ns]", "number of particles", 100, -1, 1);
      }

  // book legend for stack of 3 histograms
  legend = result->AddLegend(0.25, 0.86, 0.45, 0.98);
  legend->AddEntry(plots->fmatchedjetPT[0][0], "PUPPI jet", "l");
  legend->AddEntry(plots->fmatchedjetPT[1][0], "Genjet", "l");

  // attach legend to stacks (legend will be printed over stack in .eps file)
  for (size_t p = 0; p < 5; p++) {
    result->Attach(stack[p][0], legend);
    result->Attach(stack[p][1], legend);
  }


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

    Double_t matching_radius_large = 1;
    Double_t matching_radius = 0.4;
    Double_t gen_matching_radius = 0.2;

    Long64_t entry;

    Int_t i;
    // Loop over all events
    for(entry = 0; entry < allEntries; ++entry)
    {
      //cout << "-------- begin event ---------"<< endl;
      std::vector<GenParticle*> VBF_genparts;
      std::vector<GenParticle*> b_genparts;
      std::vector<GenParticle*> PU_genparts;
      std::vector<GenParticle*> genparts;


      // Load selected branches with data from specified event
      treeReader->ReadEntry(entry);

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
          // select VBF jets
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

    // --------Analyse jets -- 0 CHS, 1 PUPPI,  2 GenJets ---------
    for (Int_t m = 0; m < 3; m++) {
    //  cout <<"................................"<< branchJets[m]->GetName() << endl;
      if(branchJets[m]->GetEntriesFast() > 0)
      {
        Jet* jet;
        for (size_t k = 0; k < branchJets[m]->GetEntriesFast(); k++) {
          jet = (Jet*) branchJets[m]->At(k);
          if(jet->PT>30){
           plots->fJetPT[m]->Fill(jet->PT);
           plots->fJetEta[m]->Fill(jet->Eta);

          GenParticle* particle = get_closest_particle(jet, genparts, 0.4);
          if (is_b(particle)) {
              //put into list of b jets [1]
              plots->fmatchedjetDeltaR[m][1]->Fill(get_distance(jet, particle));
              plots->fmatchedjetPT[m][1]->Fill(jet->PT);
              plots->fmatchedjetEta[m][1]->Fill(jet->Eta);
              plots->fmatchedjetDeltaPhi[m][1]->Fill(abs(jet->Phi - particle->Phi));
              plots->fmatchedjetDeltaEta[m][1]->Fill(abs(jet->Eta - particle->Eta));
            }
          else if(is_VBF(particle)){
            // put into list of VBF jets [0]
            plots->fmatchedjetDeltaR[m][0]->Fill(get_distance(jet, particle));
            plots->fmatchedjetPT[m][0]->Fill(jet->PT);
            plots->fmatchedjetEta[m][0]->Fill(jet->Eta);
            plots->fmatchedjetDeltaPhi[m][0]->Fill(abs(jet->Phi - particle->Phi));
            plots->fmatchedjetDeltaEta[m][0]->Fill(abs(jet->Eta - particle->Eta));

            // Loop over all jet's constituents
            for(size_t j = 0; j < jet->Constituents.GetEntriesFast(); ++j)
            {
              object = jet->Constituents.At(j);
              // Check if the constituent is accessible
              if(object == 0)
              {
                continue;
              }

              if(object->IsA() == GenParticle::Class())
              {
                // todo split into charged and neutral particles
                particle = (GenParticle*) object;
                 plots->fVBFmatchedGenPartPT[0]->Fill((particle->PT));
                 plots->fVBFmatchedGenPartEta[0]->Fill((particle->Eta));
                 plots->fVBFmatchedGenPartT[0]->Fill((particle->T - Pvtx_T) * 1000000000);
                 if (particle->Charge == 0) {
                   plots->fVBFmatchedGenPartPT[1]->Fill((particle->PT));
                   plots->fVBFmatchedGenPartEta[1]->Fill((particle->Eta));
                   plots->fVBFmatchedGenPartT[1]->Fill((particle->T - Pvtx_T) * 1000000000);
                 }
                 else { // only else?
                   plots->fVBFmatchedGenPartPT[2]->Fill((particle->PT));
                   plots->fVBFmatchedGenPartEta[2]->Fill((particle->Eta));
                   plots->fVBFmatchedGenPartT[2]->Fill((particle->T - Pvtx_T) * 1000000000);
                 }
              }
              else if(object->IsA() == Track::Class())
              {
                if (track->Charge != 0) {

                track = (Track*) object;
                GenParticle *p = static_cast<GenParticle*>(track->Particle.GetObject());
                plots->fmatchedTrackPT[m][0][0]->Fill(track->PT);
                plots->fmatchedTrackEta[m][0][0]->Fill(track->Eta);
                plots->fmatchedTrackT[m][0][0]->Fill((track->T - Pvtx_T) * 1000000000);
                plots->fmatchedTrackDZ[m][0][0]->Fill((track->Z - Pvtx_Z) / 1000);
                plots->fmatchedTrackZ[m][0][0]->Fill((track->Z) / 1000);

                if (p->IsPU == 0) { // signal tracks
                  plots->fmatchedTrackPT[m][1][0]->Fill(track->PT);
                  plots->fmatchedTrackEta[m][1][0]->Fill(track->Eta);
                  plots->fmatchedTrackT[m][1][0]->Fill((track->T - Pvtx_T) * 1000000000);
                  plots->fmatchedTrackDZ[m][1][0]->Fill((track->Z - Pvtx_Z) / 1000);
                  plots->fmatchedTrackZ[m][1][0]->Fill((track->Z) / 1000);

                }
                else if (p->IsPU == 1) { // PU tracks
                  plots->fmatchedTrackPT[m][2][0]->Fill(track->PT);
                  plots->fmatchedTrackEta[m][2][0]->Fill(track->Eta);
                  plots->fmatchedTrackT[m][2][0]->Fill((track->T - Pvtx_T) * 1000000000);
                  plots->fmatchedTrackDZ[m][2][0]->Fill((track->Z - Pvtx_Z) / 1000);
                  plots->fmatchedTrackZ[m][2][0]->Fill((track->Z) / 1000);

                }
              }
              }
              else if(object->IsA() == Tower::Class())
              {
                tower = (Tower*) object;
                   GenParticle *p = static_cast<GenParticle*>(tower->Particles.At(0)); // reference is always the same particle

              }
              // pf candidates
              else if(object->IsA() == ParticleFlowCandidate::Class())
              {
                pf = (ParticleFlowCandidate*) object;
                if (pf->Charge != 0) {
                GenParticle *p = static_cast<GenParticle*>(pf->Particles.At(0));
                plots->fmatchedTrackPT[m][0][0]->Fill(pf->PT);
                plots->fmatchedTrackEta[m][0][0]->Fill(pf->Eta);
                plots->fmatchedTrackT[m][0][0]->Fill((pf->T - Pvtx_T) * 1000000000);
                plots->fmatchedTrackDZ[m][0][0]->Fill((pf->Z - Pvtx_Z) / 1000);
                plots->fmatchedTrackZ[m][0][0]->Fill((pf->Z) / 1000);

                if (p->IsPU == 0) { // signal tracks
                  plots->fmatchedTrackPT[m][1][0]->Fill(pf->PT);
                  plots->fmatchedTrackEta[m][1][0]->Fill(pf->Eta);
                  plots->fmatchedTrackT[m][1][0]->Fill((pf->T - Pvtx_T) * 1000000000);
                  plots->fmatchedTrackDZ[m][1][0]->Fill((pf->Z - Pvtx_Z) / 1000);
                  plots->fmatchedTrackZ[m][1][0]->Fill((pf->Z) / 1000);

                }
                else if (p->IsPU == 1) { // PU tracks
                  plots->fmatchedTrackPT[m][2][0]->Fill(pf->PT);
                  plots->fmatchedTrackEta[m][2][0]->Fill(pf->Eta);
                  plots->fmatchedTrackT[m][2][0]->Fill((pf->T - Pvtx_T) * 1000000000);
                  plots->fmatchedTrackDZ[m][2][0]->Fill((pf->Z - Pvtx_Z) / 1000);
                  plots->fmatchedTrackZ[m][2][0]->Fill((pf->Z) / 1000);

                }
              }
              }
            } // END JET CONSTITUENTS

          } // end is VBF matched
          else{
            // rest of the jets
            plots->fmatchedjetPT[m][2]->Fill(jet->PT);
            plots->fmatchedjetEta[m][2]->Fill(jet->Eta);
            for (size_t k = 0; k < jet->Constituents.GetEntriesFast(); k++) {
              object = jet->Constituents.At(k);
              if (object == 0) continue;
              if (object->IsA() == GenParticle::Class()) {
                constituent = (GenParticle*) jet->Particles.At(k);
                if (constituent->Charge != 0) {
                  //TODO
                }
              }
              else if(object->IsA() == ParticleFlowCandidate::Class())
              {
                pf = (ParticleFlowCandidate*) object;
                GenParticle *p = static_cast<GenParticle*>(pf->Particles.At(0));
                if (pf->Charge != 0) {
                plots->fmatchedTrackPT[m][0][2]->Fill(pf->PT);
                plots->fmatchedTrackEta[m][0][2]->Fill(pf->Eta);
                plots->fmatchedTrackT[m][0][2]->Fill((pf->T - Pvtx_T) * 1000000000);
                plots->fmatchedTrackDZ[m][0][2]->Fill((pf->Z - Pvtx_Z) / 1000);

                if (p->IsPU == 0) { // signal tracks
                  plots->fmatchedTrackPT[m][1][2]->Fill(pf->PT);
                  plots->fmatchedTrackEta[m][1][2]->Fill(pf->Eta);
                  plots->fmatchedTrackT[m][1][2]->Fill((pf->T - Pvtx_T) * 1000000000);
                  plots->fmatchedTrackDZ[m][1][2]->Fill((pf->Z - Pvtx_Z) / 1000);
                }
                else if (p->IsPU == 1) { // PU tracks
                  plots->fmatchedTrackPT[m][2][2]->Fill(pf->PT);
                  plots->fmatchedTrackEta[m][2][2]->Fill(pf->Eta);
                  plots->fmatchedTrackT[m][2][2]->Fill((pf->T - Pvtx_T) * 1000000000);
                  plots->fmatchedTrackDZ[m][2][2]->Fill((pf->Z - Pvtx_Z) / 1000);
                }
              }
              }
            }
            // match not matched jets to genjets
             matched_to_genjet = get_closest_jet(branchGenJet, jet, 0.1);
             Bool_t matched;
             matched = matched_to_jet(branchGenJet, jet, 0.1); // is there a close genjet?
             if (matched) { // if there is a genjet
               plots->fmatchedgenjetPT[m]->Fill(jet->PT);
               plots->fmatchedgenjetEta[m]->Fill(jet->Eta);
             }
             if (!matched) { // if there is no genjet
               plots->fnogenjetPT[m]->Fill(jet->PT);
               plots->fnogenjetEta[m]->Fill(jet->Eta);
             }
          }//end else

          }
        }
      }
    } // end loop over jet branches

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
    gSystem->cd("Plots/final/macro_jets_snowmass/");
    cout << "Print hists "<< endl;

    PrintHistograms(result, plots);

    result->Write("results.root");

    cout << "** Exiting..." << endl;

    delete plots;
    delete result;
    delete treeReader;
    delete chain;
  }
