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

  TH1 *fJetAK8PT;
  TH1 *fJetAK8Eta;

  TH1 *fgenparticlePT[4];
  TH1 *fgenparticleEta[4];

  // plots for all genjet collections (3) for VBF and b matched jets (2)
  TH1 *fmatchedjetDeltaR[3][2];
  TH1 *fmatchedjetPT[3][2];
  TH1 *fmatchedjetEta[3][2];
  TH1 *fmatchedjetDeltaEta[3][2];
  TH1 *fmatchedjetDeltaPhi[3][2];

  TH1 *fVBFmatchedjetCHPT[3];
  TH1 *fVBFmatchedjetCHEta[3];
  TH1 *fVBFmatchedjetCHT[3];
  TH1 *fVBFmatchedjetCHZ[3];

  TH1 *fnotmatchedjetPT[3];
  TH1 *fnotmatchedjetEta[3];
  TH1 *fnotmatchedjetCHT[3];
  TH1 *fnotmatchedjetCHZ[3];

  // VBF tracks Plots
  TH1 *fVBFmatchedTrackPT[2][3];
  TH1 *fVBFmatchedTrackEta[2][3];
  TH1 *fVBFmatchedTrackT[2][3];
  TH1 *fVBFmatchedTrackPartT[2][3];

  TH1 *fVBFmatchedGenPartPT[3];
  TH1 *fVBFmatchedGenPartEta[3];
  TH1 *fVBFmatchedGenPartT[3];

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

  plots->fJetAK8PT = result->AddHist1D(
    "jet_AK8_pt_all", "all AK8 jet P_{T}",
    "AK8 jet P_{T}, GeV/c", "number of jets", 100, 0.0, 200.0);

  plots->fJetAK8Eta = result->AddHist1D(
    "jet_AK8_eta_all", "all AK8 jet #eta",
    "AK8 jet #eta", "number of jets", 100, -5.0, 5.0);

  plots->fJetPT[1] = result->AddHist1D(
    "jet_PUPPI_pt_all", "all PUPPI jet P_{T}",
    "PUPPI jet P_{T}, GeV/c", "number of jets", 100, 0.0, 200.0);

  plots->fJetEta[1] = result->AddHist1D(
    "jet_PUPPI_eta_all", "all PUPPI jet #eta",
    "PUPPI jet #eta", "number of jets", 100, -5.0, 5.0);
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
  TString name[2] = {"PUPPI", "GenJet"};
  TString name2[3] = {"","_signal", "_PU"};
  TString category[2] = {"VBF", "b"};

  for (size_t b = 0; b < 2; b++) { // begin VBF or b matched
  stack[0][b] = result->AddHistStack(category[b]+"_matched_jet_pt", category[b]+" matched jets P_{T}");
  stack[1][b] = result->AddHistStack(category[b]+"_matched_jet_eta", category[b]+" matched jets #eta");
  stack[2][b] = result->AddHistStack(category[b]+"_matched_jet_deltaR", category[b]+" matched jets #Delta R");
  stack[3][b] = result->AddHistStack(category[b]+"_matched_jet_deltaPhi", category[b]+" matched jets #Delta #Phi");
  stack[4][b] = result->AddHistStack(category[b]+"_matched_jet_deltaEta", category[b]+" matched jets #Delta #Eta");

  for (size_t i = 0; i < 2; i++) { // jet collections
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
    if (i<2 && b==0) { // skip for GenJets and b matched
    for (size_t k = 0; k < 3; k++) { // split into signal and PU
      plots->fVBFmatchedTrackPT[i][k] = result->AddHist1D(
        "VBF_matched_"+name[i]+name2[k]+"_track_pt", "VBF matched track P_{T}",
        "VBF matched track P_{T}, GeV/c", "number of tracks", 100, 0.0, 10.0);
        plots->fVBFmatchedTrackEta[i][k] = result->AddHist1D(
          "VBF_matched_"+name[i]+name2[k]+"_track_eta", "VBF matched track #eta",
          "VBF matched track #eta", "number of tracks", 100, -5.0, 5.0);
      plots->fVBFmatchedTrackT[i][k] = result->AddHist1D(
        "VBF_matched_"+name[i]+name2[k]+"_track_T", "VBF matched track T",
        "VBF matched track dT [ns]", "number of tracks", 100, -1, 1);
        plots->fVBFmatchedTrackPartT[i][k] = result->AddHist1D(
          "VBF_matched_"+name[i]+name2[k]+"_track_p_T", "track->p T",
          "VBF matched track->p dT [ns]", "number of particles", 100, -1, 1);
        } // end signal or PU
      }
  } // end jet collections
} // end VBF or b

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

  //----------not matched------------------------
  // book 3 histograms for PT of matched jets, PUPPIjets and genjet
  plots->fnotmatchedjetPT[0] = result->AddHist1D(
    "not_matched_AK4_jet_pt", "not matched AK4 jet P_{T}",
    "not matched jet P_{T}, GeV/c", "number of jets", 100, 0.0, 200.0);
  plots->fnotmatchedjetPT[1] = result->AddHist1D(
    "not_matched_PUPPI_jet_pt", "not matched PUPPI jet P_{T}",
    "not matched jet P_{T}, GeV/c", "number of jets", 100, 0.0, 200.0);
  plots->fnotmatchedjetPT[2] = result->AddHist1D(
    "not_matched_genjet_pt", "not matched Genjet P_{T}",
    "not matched jet P_{T}, GeV/c", "number of jets", 100, 0.0, 200.0);

  plots->fnotmatchedjetPT[0]->SetLineColor(kRed);
  plots->fnotmatchedjetPT[1]->SetLineColor(kBlue);
  plots->fnotmatchedjetPT[2]->SetLineColor(kGreen+3);

  // book 1 stack of 2 histograms
  stack_PU_1 = result->AddHistStack("not_matched_jet_pt", "not matched jets P_{T}");
  stack_PU_1->Add(plots->fnotmatchedjetPT[0]);
  stack_PU_1->Add(plots->fnotmatchedjetPT[1]);
  stack_PU_1->Add(plots->fnotmatchedjetPT[2]);

  // attach legend to stack (legend will be printed over stack in .eps file)
  result->Attach(stack_PU_1, legend);

  plots->fnotmatchedjetEta[0] = result->AddHist1D(
    "not_matched_AK4_jet_eta", "not matched AK4 jet #eta",
    "not matched jet #eta", "number of jets", 100, -5.0, 5.0);
  plots->fnotmatchedjetEta[1] = result->AddHist1D(
    "not_matched_PUPPI_jet_eta", "not matched PUPPI jet #eta",
    "not matched jet #eta", "number of jets", 100, -5.0, 5.0);
  plots->fnotmatchedjetEta[2] = result->AddHist1D(
    "not_matched_genjet_eta", "not matched genjet #eta",
    "not matched jet #eta", "number of jets", 100, -5.0, 5.0);

  plots->fnotmatchedjetEta[0]->SetLineColor(kRed);
  plots->fnotmatchedjetEta[1]->SetLineColor(kBlue);
  plots->fnotmatchedjetEta[2]->SetLineColor(kGreen+3);

  // book 1 stack of 2 histograms
  stack_PU_2 = result->AddHistStack("not_matched_jet_eta", "not matched jets #eta");
  stack_PU_2->Add(plots->fnotmatchedjetEta[0]);
  stack_PU_2->Add(plots->fnotmatchedjetEta[1]);
  stack_PU_2->Add(plots->fnotmatchedjetEta[2]);

  // attach legend to stack (legend will be printed over stack in .eps file)
  result->Attach(stack_PU_2, legend);

  }

  //------------------------------------------------------------------------------

  void AnalyseEvents(ExRootTreeReader *treeReader, MyPlots *plots)
  {
    TClonesArray *branchJetPUPPI = treeReader->UseBranch("JetPUPPI");

    TClonesArray *branchParticle = treeReader->UseBranch("Particle"); // for identification of VBF and b quarks
    TClonesArray *branchfilteredParticle = treeReader->UseBranch("filteredParticle"); // input to genjets
    TClonesArray *branchMergerParticle = treeReader->UseBranch("mergerParticle"); // input to tracks

    TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
    TClonesArray *branchVtx = treeReader->UseBranch("Vertex");

    TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
    TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
    TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");

    TClonesArray *branchPuppiParticle = treeReader->UseBranch("ParticleFlowCandidate");

    TClonesArray *branchJets[2] = {branchJetPUPPI, branchGenJet};

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
          }
      }

    // --------Analyse jets -- 1 PUPPI,  2 GenJets ---------
    for (Int_t m = 0; m < 2; m++) {
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
                track = (Track*) object;
                GenParticle *p = static_cast<GenParticle*>(track->Particle.GetObject());
                plots->fVBFmatchedTrackPT[m][0]->Fill(track->PT);
                plots->fVBFmatchedTrackEta[m][0]->Fill(track->Eta);
                plots->fVBFmatchedTrackT[m][0]->Fill((track->T - Pvtx_T) * 1000000000);

                if (p->IsPU == 0) { // signal tracks
                  plots->fVBFmatchedTrackPT[m][1]->Fill(track->PT);
                  plots->fVBFmatchedTrackEta[m][1]->Fill(track->Eta);
                  plots->fVBFmatchedTrackT[m][1]->Fill((track->T - Pvtx_T) * 1000000000);
                }
                else if (p->IsPU == 1) { // PU tracks
                  plots->fVBFmatchedTrackPT[m][2]->Fill(track->PT);
                  plots->fVBFmatchedTrackEta[m][2]->Fill(track->Eta);
                  plots->fVBFmatchedTrackT[m][2]->Fill((track->T - Pvtx_T) * 1000000000);
                }
              }
              else if(object->IsA() == Tower::Class())
              {
                tower = (Tower*) object;
                   GenParticle *p = static_cast<GenParticle*>(tower->Particles.At(0)); // reference is always the same particle

              }
            } // END JET CONSTITUENTS

          } // end is VBF matched
          else{
            // rest of the jets
            plots->fnotmatchedjetPT[m]->Fill(jet->PT);
            plots->fnotmatchedjetEta[m]->Fill(jet->Eta);
            for (size_t k = 0; k < jet->Particles.GetEntriesFast(); k++) {
              object = jet->Particles.At(k);
              if (object == 0) continue;
              if (object->IsA() == GenParticle::Class()) {
                constituent = (GenParticle*) jet->Particles.At(k);
                if (constituent->Charge != 0) {
                  // TODO fill tracks and genparticles properly
                }
              }
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
    gSystem->cd("Plots/macro_jets_snowmass/");
    cout << "Print hists "<< endl;

    PrintHistograms(result, plots);

    result->Write("results.root");

    cout << "** Exiting..." << endl;

    delete plots;
    delete result;
    delete treeReader;
    delete chain;
  }
