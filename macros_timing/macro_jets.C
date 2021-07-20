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

  TH1 *fVBFmatchedjetDeltaR[3];
  TH1 *fVBFmatchedjetPT[3];
  TH1 *fVBFmatchedjetEta[3];
  TH1 *fVBFmatchedjetCHPT[3];
  TH1 *fVBFmatchedjetCHEta[3];
  TH1 *fVBFmatchedjetCHT[3];
  TH1 *fVBFmatchedjetCHZ[3];

  TH1 *fbmatchedjetDeltaR[3];
  TH1 *fbmatchedjetPT[3];
  TH1 *fbmatchedjetEta[3];

  TH1 *fnotmatchedjetPT[3];
  TH1 *fnotmatchedjetEta[3];
  TH1 *fnotmatchedjetCHT[3];
  TH1 *fnotmatchedjetCHZ[3];

  // VBF tracks Plots
  TH1 *fVBFmatchedTrackPT[3];
  TH1 *fVBFmatchedTrackT[3];
  TH1 *fVBFmatchedTrackPartT[3];

  TH1 *fVBFmatchedGenPartT;

};

//------------------------------------------------------------------------------

class ExRootResult;
class ExRootTreeReader;

//------------------------------------------------------------------------------

void BookHistogramsBasic(ExRootResult *result, MyPlots *plots)
{
  TLegend *legend;

  THStack *stack_1;
  THStack *stack_2;
  THStack *stack_3;
  THStack *stack_4;
  THStack *stack_5;
  THStack *stack_6;
  THStack *stack_7;

  THStack *stack_b_1;
  THStack *stack_b_2;
  THStack *stack_b_3;

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
  // plots->fgenparticlePT[3] = result->AddHist1D(
  //   "genparticle_PU_pt", "genparticle PU P_{T}",
  //   "genparticle PU P_{T}, GeV/c", "number of genparticles", 100, 0.0, 200.0);
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
  // plots->fgenparticleEta[3] = result->AddHist1D(
  //   "genparticle_PU_eta", "genparticle PU #eta",
  //   "genparticle PU #eta", "number of genparticles", 100, -5.0, 5.0);

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

  for (size_t i = 0; i < 3; i++) {
    // pt
    plots->fVBFmatchedjetPT[i] = result->AddHist1D(
      "VBF_matched_"+name[i]+"_jet_pt", "VBF matched jet P_{T}",
      "VBF matched jet P_{T}, GeV/c", "number of jets", 100, 0.0, 200.0);
    plots->fVBFmatchedjetPT[i]->SetLineColor(i+1);
    // eta
    plots->fVBFmatchedjetEta[i] = result->AddHist1D(
      "VBF_matched_"+name[i]+"_jet_eta", "VBF matched jet #eta",
      "VBF matched jet #eta", "number of jets", 100, -5.0, 5.0);
    plots->fVBFmatchedjetEta[i]->SetLineColor(i+1);
    // delta R
    plots->fVBFmatchedjetDeltaR[i] = result->AddHist1D(
      "VBF_matched_"+name[i]+"_jet_deltaR", "VBF matched jet #Delta R",
      "VBF matched jet #Delta R", "number of jets", 50, 0.0, 1.0);
    plots->fVBFmatchedjetDeltaR[i]->SetLineColor(i+1);
    // VBF CH pt
    plots->fVBFmatchedjetCHPT[i] = result->AddHist1D(
      "VBF_"+name[i]+"_CH_pt", "VBF CH P_{T}",
      "VBF CH P_{T}, GeV/c", "number of particles", 100, 0.0, 40.0);
    plots->fVBFmatchedjetCHPT[i]->SetLineColor(i+1);
    // VBF CH eta
    plots->fVBFmatchedjetCHEta[i] = result->AddHist1D(
      "VBF_"+name[i]+"_CH_eta", "VBF CH #eta",
      "VBF CH #eta", "number of particles", 100, -5.0, 5.0);
    plots->fVBFmatchedjetCHEta[i]->SetLineColor(i+1);
    // VBF CH Time
    plots->fVBFmatchedjetCHT[i] = result->AddHist1D(
      "VBF_"+name[i]+"_CH_vtx_T", "VBF CH vtx T",
      "VBF CH vtx T - PV T [ns]", "number of particles", 100, -0.5, 0.5);
    plots->fVBFmatchedjetCHT[i]->SetLineColor(i+1);
    // VBF CH Z
    plots->fVBFmatchedjetCHZ[i] = result->AddHist1D(
      "VBF_"+name[i]+"_CH_vtx_Z", "VBF CH vtx Z",
      "VBF CH vtx Z - PV Z [m]", "number of particles", 100, -0.25, 0.25);
    plots->fVBFmatchedjetCHZ[i]->SetLineColor(i+1);
    // timing Plots
    plots->fVBFmatchedTrackPT[i] = result->AddHist1D(
      "VBF_matched_"+name[i]+"_track_pt", "VBF matched track P_{T}",
      "VBF matched track P_{T}, GeV/c", "number of tracks", 100, 0.0, 20.0);
    plots->fVBFmatchedTrackT[i] = result->AddHist1D(
      "VBF_matched_"+name[i]+"_track_T", "VBF matched track T",
      "VBF matched track dT [ns]", "number of tracks", 100, -1, 1);
      plots->fVBFmatchedTrackPartT[i] = result->AddHist1D(
        "VBF_matched_"+name[i]+"_track_p_T", "track->p T",
        "VBF matched track->p dT [ns]", "number of particles", 100, -1, 1);
  }

  plots->fVBFmatchedGenPartT = result->AddHist1D(
    "VBF_matched_GenPart_T", "VBF matched GenParts T",
    "VBF matched GenParts dT [ns]", "number of particles", 100, -1, 1);

  // book 1 stack of 2 histograms
  stack_1 = result->AddHistStack("VBF_matched_jet_pt", "VBF matched jets P_{T}");
  stack_1->Add(plots->fVBFmatchedjetPT[0]);
  stack_1->Add(plots->fVBFmatchedjetPT[1]);
  stack_1->Add(plots->fVBFmatchedjetPT[2]);
  // book legend for stack of 3 histograms
  legend = result->AddLegend(0.25, 0.86, 0.45, 0.98);
  legend->AddEntry(plots->fVBFmatchedjetPT[0], "CHS jet", "l");
  legend->AddEntry(plots->fVBFmatchedjetPT[1], "PUPPI jet", "l");
  legend->AddEntry(plots->fVBFmatchedjetPT[2], "Genjet", "l");

  // attach legend to stack (legend will be printed over stack in .eps file)
  result->Attach(stack_1, legend);

  // book 1 stack of 2 histograms
  stack_2 = result->AddHistStack("VBF_matched_jet_eta", "VBF matched jets #eta");
  stack_2->Add(plots->fVBFmatchedjetEta[0]);
  stack_2->Add(plots->fVBFmatchedjetEta[1]);
  stack_2->Add(plots->fVBFmatchedjetEta[2]);

  // attach legend to stack (legend will be printed over stack in .eps file)
  result->Attach(stack_2, legend);

  //  book 1 stack of 2 histograms
  stack_3 = result->AddHistStack("VBF_matched_jet_deltaR", "VBF matched jets #Delta R");
  stack_3->Add(plots->fVBFmatchedjetDeltaR[0]);
  stack_3->Add(plots->fVBFmatchedjetDeltaR[1]);
  stack_3->Add(plots->fVBFmatchedjetDeltaR[2]);

  // attach legend to stack (legend will be printed over stack in .eps file)
  result->Attach(stack_3, legend);

  // book 1 stack of 2 histograms
  stack_4 = result->AddHistStack("VBF_CH_pt", "VBF CH P_{T}");
  stack_4->Add(plots->fVBFmatchedjetCHPT[0]);
  stack_4->Add(plots->fVBFmatchedjetCHPT[1]);
  stack_4->Add(plots->fVBFmatchedjetCHPT[2]);

  // attach legend to stack (legend will be printed over stack in .eps file)
  result->Attach(stack_4, legend);

  // book 1 stack of 2 histograms
  stack_5 = result->AddHistStack("VBF_CH_eta", "VBF CH #eta");
  stack_5->Add(plots->fVBFmatchedjetCHEta[0]);
  stack_5->Add(plots->fVBFmatchedjetCHEta[1]);
  stack_5->Add(plots->fVBFmatchedjetCHEta[2]);

  // attach legend to stack (legend will be printed over stack in .eps file)
  result->Attach(stack_5, legend);

  // book 1 stack of 2 histograms
  stack_6 = result->AddHistStack("VBF_CH_vtx_T", "VBF CH vtx T");
  stack_6->Add(plots->fVBFmatchedjetCHT[0]);
  stack_6->Add(plots->fVBFmatchedjetCHT[1]);
  stack_6->Add(plots->fVBFmatchedjetCHT[2]);

  // attach legend to stack (legend will be printed over stack in .eps file)
  result->Attach(stack_6, legend);

  // book 1 stack of 2 histograms
  stack_7 = result->AddHistStack("VBF_CH_vtx_Z", "VBF CH vtx Z");
  stack_7->Add(plots->fVBFmatchedjetCHZ[0]);
  stack_7->Add(plots->fVBFmatchedjetCHZ[1]);
  stack_7->Add(plots->fVBFmatchedjetCHZ[2]);

  // attach legend to stack (legend will be printed over stack in .eps file)
  result->Attach(stack_7, legend);

// TODO include b and not matched into loop, if possible also stacks
  //-----------b -matched----------
  // book 3 histograms for PT of matched jets, PUPPIjets and genjet
  plots->fbmatchedjetPT[0] = result->AddHist1D(
    "b_matched_AK4_jet_pt", "b matched AK4 jet P_{T}",
    "b matched jet P_{T}, GeV/c", "number of jets", 100, 0.0, 200.0);
  plots->fbmatchedjetPT[1] = result->AddHist1D(
    "b_matched_PUPPI_jet_pt", "b matched PUPPI jet P_{T}",
    "b matched jet P_{T}, GeV/c", "number of jets", 100, 0.0, 200.0);
  plots->fbmatchedjetPT[2] = result->AddHist1D(
    "b_matched_genjet_pt", "b matched Genjet P_{T}",
    "b matched jet P_{T}, GeV/c", "number of jets", 100, 0.0, 200.0);

  plots->fbmatchedjetPT[0]->SetLineColor(kRed);
  plots->fbmatchedjetPT[1]->SetLineColor(kBlue);
  plots->fbmatchedjetPT[2]->SetLineColor(kGreen+3);

  // // book 1 stack of 2 histograms
  stack_b_1 = result->AddHistStack("b_matched_jet_pt", "b matched jets P_{T}");
  stack_b_1->Add(plots->fbmatchedjetPT[0]);
  stack_b_1->Add(plots->fbmatchedjetPT[1]);
  stack_b_1->Add(plots->fbmatchedjetPT[2]);

  // attach legend to stack (legend will be printed over stack in .eps file)
  result->Attach(stack_b_1, legend);

  plots->fbmatchedjetEta[0] = result->AddHist1D(
    "b_matched_AK4_jet_eta", "b matched AK4 jet #eta",
    "b matched jet #eta", "number of jets", 100, -5.0, 5.0);
  plots->fbmatchedjetEta[1] = result->AddHist1D(
    "b_matched_PUPPI_jet_eta", "b matched PUPPI jet #eta",
    "b matched jet #eta", "number of jets", 100, -5.0, 5.0);
  plots->fbmatchedjetEta[2] = result->AddHist1D(
    "b_matched_genjet_eta", "b matched genjet #eta",
    "b matched jet #eta", "number of jets", 100, -5.0, 5.0);

  plots->fbmatchedjetEta[0]->SetLineColor(kRed);
  plots->fbmatchedjetEta[1]->SetLineColor(kBlue);
  plots->fbmatchedjetEta[2]->SetLineColor(kGreen+3);

  // book 1 stack of 2 histograms
  stack_b_2 = result->AddHistStack("b_matched_jet_eta", "b matched jets #eta");
  stack_b_2->Add(plots->fbmatchedjetEta[0]);
  stack_b_2->Add(plots->fbmatchedjetEta[1]);
  stack_b_2->Add(plots->fbmatchedjetEta[2]);

  // attach legend to stack (legend will be printed over stack in .eps file)
  result->Attach(stack_b_2, legend);

  plots->fbmatchedjetDeltaR[0] = result->AddHist1D(
    "b_matched_AK4_jet_deltaR", "b matched AK4 jet #Delta R",
    "b matched jet #Delta R", "number of jets", 50, 0.0, 1.0);
  plots->fbmatchedjetDeltaR[1] = result->AddHist1D(
    "b_matched_PUPPI_jet_deltaR", "b matched PUPPI jet #Delta R",
    "b matched jet #Delta R", "number of jets", 50, 0.0, 1.0);
  plots->fbmatchedjetDeltaR[2] = result->AddHist1D(
    "b_matched_genjet_deltaR", "b matched genjet #Delta R",
    "b matched jet #Delta R", "number of jets", 50, 0.0, 1.0);

  plots->fbmatchedjetDeltaR[0]->SetLineColor(kRed);
  plots->fbmatchedjetDeltaR[1]->SetLineColor(kBlue);
  plots->fbmatchedjetDeltaR[2]->SetLineColor(kGreen+3);

  //  book 1 stack of 2 histograms
  stack_b_3 = result->AddHistStack("b_matched_jet_deltaR", "b matched jets #Delta R");
  stack_b_3->Add(plots->fbmatchedjetDeltaR[0]);
  stack_b_3->Add(plots->fbmatchedjetDeltaR[1]);
  stack_b_3->Add(plots->fbmatchedjetDeltaR[2]);

  // attach legend to stack (legend will be printed over stack in .eps file)
  result->Attach(stack_b_3, legend);
  //----------PU matched------------------------
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

  plots->fnotmatchedjetCHT[0] = result->AddHist1D(
    "not_matched_AK4_CH_vtx_T", "not_matched CH vtx T",
    "not_matched CH vtx T - PV T [ns]", "number of particles", 100, -0.5, 0.5);
  plots->fnotmatchedjetCHT[1] = result->AddHist1D(
    "not_matched_PUPPI_CH_vtx_T", "not_matched CH vtx T",
    "not_matched CH vtx T - PV T [ns]", "number of particles", 100, -0.5, 0.5);
  plots->fnotmatchedjetCHT[2] = result->AddHist1D(
    "not_matched_genjet_CH_vtx_T", "not_matched CH vtx T",
    "not_matched CH vtx T - PV T [ns]", "number of particles", 100, -0.5, 0.5);

  plots->fnotmatchedjetCHT[0]->SetLineColor(kRed);
  plots->fnotmatchedjetCHT[1]->SetLineColor(kBlue);
  plots->fnotmatchedjetCHT[2]->SetLineColor(kGreen+3);

  // book 1 stack of 2 histograms
  stack_PU_3 = result->AddHistStack("not_matched_CH_vtx_T", "not_matched CH vtx T");
  stack_PU_3->Add(plots->fnotmatchedjetCHT[0]);
  stack_PU_3->Add(plots->fnotmatchedjetCHT[1]);
  stack_PU_3->Add(plots->fnotmatchedjetCHT[2]);

  // attach legend to stack (legend will be printed over stack in .eps file)
  result->Attach(stack_PU_3, legend);

  plots->fnotmatchedjetCHZ[0] = result->AddHist1D(
    "not_matched_AK4_CH_vtx_Z", "not_matched CH vtx Z",
    "not_matched CH vtx Z - PV Z [m]", "number of particles", 100, -0.25, 0.25);
  plots->fnotmatchedjetCHZ[1] = result->AddHist1D(
    "not_matched_PUPPI_CH_vtx_Z", "not_matched CH vtx Z",
    "not_matched CH vtx Z - PV Z [m]", "number of particles", 100, -0.25, 0.25);
  plots->fnotmatchedjetCHZ[2] = result->AddHist1D(
    "not_matched_genjet_CH_vtx_Z", "not_matched CH vtx Z",
    "not_matched CH vtx Z - PV Z [m]", "number of particles", 100, -0.25, 0.25);

  plots->fnotmatchedjetCHZ[0]->SetLineColor(kRed);
  plots->fnotmatchedjetCHZ[1]->SetLineColor(kBlue);
  plots->fnotmatchedjetCHZ[2]->SetLineColor(kGreen+3);

  // book 1 stack of 2 histograms
  stack_PU_4 = result->AddHistStack("not_matched_CH_vtx_Z", "not_matched CH vtx Z");
  stack_PU_4->Add(plots->fnotmatchedjetCHZ[0]);
  stack_PU_4->Add(plots->fnotmatchedjetCHZ[1]);
  stack_PU_4->Add(plots->fnotmatchedjetCHZ[2]);

  // attach legend to stack (legend will be printed over stack in .eps file)
  result->Attach(stack_PU_4, legend);
  }

  //------------------------------------------------------------------------------

  void AnalyseEvents(ExRootTreeReader *treeReader, MyPlots *plots)
  {
    TClonesArray *branchJet = treeReader->UseBranch("Jet");
    TClonesArray *branchJetPUPPI = treeReader->UseBranch("JetPUPPI");
    TClonesArray *branchJetAK8 = treeReader->UseBranch("JetAK8");

    TClonesArray *branchParticle = treeReader->UseBranch("Particle"); // for identification of VBF and b quarks
    TClonesArray *branchfilteredParticle = treeReader->UseBranch("filteredParticle"); // input to genjets
    TClonesArray *branchMergerParticle = treeReader->UseBranch("mergerSignalParticle"); // input to tracks

    TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
    TClonesArray *branchVtx = treeReader->UseBranch("Vertex");

    TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
    TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
    TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");

    TClonesArray *branchJets[3] = {branchJet, branchJetPUPPI, branchGenJet};

    Long64_t allEntries = treeReader->GetEntries();

    cout << "** Chain contains " << allEntries << " events" << endl;

    Jet *jetAK8;

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
      cout << "-------- begin event ---------"<< endl;
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

          //select genparticles from PU
          // if (abs(genparticle->IsPU)==1) {
          //   PU_genparts.push_back(genparticle);
          //    plots->fgenparticlePT[3]->Fill(genparticle->PT);
          //    plots->fgenparticleEta[3]->Fill(genparticle->Eta);
          // }
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

    // --------Analyse jets -- 0 AK4 jets, 1 PUPPI,  2 GenJets ---------
    for (Int_t m = 0; m < 3; m++) {
      cout <<"................................"<< branchJets[m]->GetName() << "  , "<< branchJets[m]->GetEntriesFast() << endl;
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
              //put into list of b jets
              plots->fbmatchedjetDeltaR[m]->Fill(get_distance(jet, particle));
              plots->fbmatchedjetPT[m]->Fill(jet->PT);
              plots->fbmatchedjetEta[m]->Fill(jet->Eta);
            }
          else if(is_VBF(particle)){
            // put into list of VBF jets
            plots->fVBFmatchedjetDeltaR[m]->Fill(get_distance(jet, particle));
            plots->fVBFmatchedjetPT[m]->Fill(jet->PT);
            plots->fVBFmatchedjetEta[m]->Fill(jet->Eta);
            // plots->fVBFmatchedjetNCharged[m]->Fill(jet->NCharged);
            // plots->fVBFmatchedjetNNeutrals[m]->Fill(jet->NNeutrals);

             cout<<"Looping over jet constituents. Jet pt: "<<jet->PT<<", eta: "<<jet->Eta<<", phi: "<<jet->Phi<<endl;
             cout << "Constituents Size "<< jet->Constituents.GetEntriesFast() << endl;
            // cout << " charged " << jet->NCharged << endl;
            // cout << " neutral " << jet->NNeutrals << endl;
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
                particle = (GenParticle*) object;
                // cout << "    GenPart pt: " << particle->PT << ", eta: " << particle->Eta << ", phi: " << particle->Phi << endl;
                 //cout << particle->IsPU << endl;
                 plots->fVBFmatchedGenPartT->Fill((particle->T - Pvtx_T) * 1000000000);
              }
              else if(object->IsA() == Track::Class())
              {
                track = (Track*) object;
                GenParticle *p = static_cast<GenParticle*>(track->Particle.GetObject());

                //plots->ftrackPT_VBF->Fill(track->PT);
                plots->fVBFmatchedTrackPT[m]->Fill(track->PT);
                plots->fVBFmatchedTrackT[m]->Fill((track->T - Pvtx_T) * 1000000000);
                plots->fVBFmatchedTrackPartT[m]->Fill((p->T - Pvtx_T) * 1000000000);

                  //cout << "    Track pt: " << track->PT << ", eta: " << track->Eta << ", phi: " << track->Phi << endl;
                 // cout << " Track Time " << track->T << endl;
                 // cout << "Particle PID " <<p->PID << endl;
                  //cout << "Particle IsPU " <<p->IsPU << endl;
                 // cout << "Particle Time " <<p->T << endl;
              }
              else if(object->IsA() == Tower::Class())
              {
                tower = (Tower*) object;
                  //cout << "    Tower pt: " << tower->ET << ", eta: " << tower->Eta << ", phi: " << tower->Phi << endl;
                   GenParticle *p = static_cast<GenParticle*>(tower->Particles.At(0)); // reference is always the same particle
                    // cout << "Tower Particle PID" << p->PID << endl;
                    // cout << "Tower Particle Time" << p->T << endl;

              }
            } // END JET CONSTITUENTS

              // TODO this gives always the same PID -> check again
                // for(Int_t j = 0; j < jet->Particles.GetEntriesFast(); ++j)
                // {
                //   GenParticle *p = static_cast<GenParticle*>(jet->Particles.At(j));
                //   cout<<"Particle PID "<< ( particle)->PID <<endl;
                // }
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
            plots->fnotmatchedjetCHT[m]->Fill((constituent->T - Pvtx_T )* 1000000000);
            plots->fnotmatchedjetCHZ[m]->Fill((constituent->Z - Pvtx_Z )/1000);
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

  void macro_jets(const char *inputFile)
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
