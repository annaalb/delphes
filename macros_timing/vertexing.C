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
  // general vertex variables
  TH1* fvtxX;
  TH1* fvtxY;
  TH1* fvtxZ;
  TH1* fvtxT;

  TH1* fvtxErrorT; // vertex position error (t component)
  TH1* fvtxErrorX; // vertex position error (x component)
  TH1* fvtxErrorY; // vertex position error (y component)
  TH1* fvtxErrorZ; // vertex position error (z component)

  TH1* fvtxIndex; // vertex index (Int_t)

  TH1* fvtxSigma; // vertex position (z component) error
  TH1* fvtxSumPT2; // sum pt^2 of tracks attached to the vertex
  TH1* fvtxGenSumPT2; // sum pt^2 of gen tracks attached to the vertex

  TH1* fvtxGenDeltaZ; // distance in z to closest generated vertex
  TH1* fvtxBTVSumPT2; // sum pt^2 of tracks attached to the secondary vertex

  // vertex from particle
  TH1* fgenparticleX[2];
  TH1* fgenparticleY[2];
  TH1* fgenparticleZ[2];
  TH1* fgenparticleT[2];

};

//------------------------------------------------------------------------------

class ExRootResult;
class ExRootTreeReader;

//------------------------------------------------------------------------------

void BookHistogramsBasic(ExRootResult *result, MyPlots *plots)
{
  plots->fvtxX = result->AddHist1D(
    "vtx_x", "vtx X",
    "vtx X [mm]", "number of vertices", 100, -2.0, 2.0);
  plots->fvtxY = result->AddHist1D(
    "vtx_y", "vtx Y",
    "vtx Y [mm]", "number of vertices", 100, -2.0, 2.0);
  plots->fvtxZ = result->AddHist1D(
    "vtx_z", "vtx Z",
    "vtx Z [m]", "number of vertices", 100, -0.25, 0.25);
  plots->fvtxT = result->AddHist1D(
    "vtx_T", "vtx T",
    "vtx T [ns]", "number of vertices", 100, -1, 1);

  plots->fvtxErrorT = result->AddHist1D(
    "vtx_error_t", "vtx error t",
    "vtx error T [ns]", "number of vertices", 100, 0, 0.01); // vertex position error (t component)
  plots->fvtxErrorX = result->AddHist1D(
    "vtx_error_x", "vtx error X",
    "vtx error X [mm]", "number of vertices", 100, 0.0, 0.01); // vertex position error (x component)
  plots->fvtxErrorY = result->AddHist1D(
    "vtx_error_y", "vtx error y",
    "vtx error Y [mm]", "number of vertices", 100, 0.0, 0.01); // vertex position error (y component)
  plots->fvtxErrorZ = result->AddHist1D(
    "vtx_error_z", "vtx error z",
    "vtx error Z [mm]", "number of vertices", 100, 0, 0.01); // vertex position error (z component)

  plots->fvtxIndex = result->AddHist1D(
    "vtx_index", "vtx index",
    "vtx index", "number of vertices", 250, 0.0, 250.0); // vertex index (Int_t)

  plots->fvtxSigma = result->AddHist1D(
    "vtx_sigma", "vtx sigma",
    "vtx sigma [mm]", "number of vertices", 100, -1.0, 1.0); // vertex position (z component) error
  plots->fvtxSumPT2 = result->AddHist1D(
    "vtx_sumpt2", "vtx sumpt2",
    "vtx sumpt2 [?]", "number of vertices", 100, 0.0, 5.0); // sum pt^2 of tracks attached to the vertex
  plots->fvtxGenSumPT2 = result->AddHist1D(
    "vtx_gensumpt2", "vtx GenSumPT2",
    "vtx GenSumPT2 [?]", "number of vertices", 100, 0.0, 5.0); // sum pt^2 of gen tracks attached to the vertex

  plots->fvtxGenDeltaZ = result->AddHist1D(
    "vtx_gendeltaZ", "vtx GenDeltaZ",
    "vtx GenDeltaZ [mm]", "number of vertices", 100, 0.0, 1.0); // distance in z to closest generated vertex
  plots->fvtxBTVSumPT2 = result->AddHist1D(
    "vtx_BTVSumPT2", "vtx BTVSumPT2",
    "vtx BTVSumPT2 [?]", "number of vertices", 100, -4.0, 4.0); // sum pt^2 of tracks attached to the secondary vertex

    TString name[2] = {"PV", "PU"};

for (size_t i = 0; i < 2; i++) {
  plots->fgenparticleX[i] = result->AddHist1D(
    name[i]+"_x", name[i]+"genparticle X",
     name[i]+" X [cm]", "number of vertices", 100, -4.0, 4.0);
  plots->fgenparticleY[i] = result->AddHist1D(
    name[i]+"_y", name[i]+"genparticle Y",
    name[i]+" Y [cm]", "number of vertices", 100, -4.0, 4.0);
  plots->fgenparticleZ[i] = result->AddHist1D(
    name[i]+"_z", name[i]+"genparticle Z",
    name[i]+" Z [m]", "number of vertices", 100, -0.25, 0.25);
  plots->fgenparticleT[i] = result->AddHist1D(
    name[i]+"_T", name[i]+"genparticle T",
    name[i]+" T [ns]", "number of vertices", 100, -1, 1);
  }

  }

  //------------------------------------------------------------------------------

  void AnalyseEvents(ExRootTreeReader *treeReader, MyPlots *plots)
  {
    Long64_t allEntries = treeReader->GetEntries();

    cout << "** Chain contains " << allEntries << " events" << endl;

    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    //TClonesArray *branchPUParticle = treeReader->UseBranch("PUParticle");
    //TClonesArray *branchGenvtx = treeReader->UseBranch("GenVertex");
    TClonesArray *branchVtx = treeReader->UseBranch("Vertex");
    //TClonesArray *branchTrack = treeReader->UseBranch("Track");

    Vertex *vtx, *recovtx;
    Track *track;
    GenParticle *genparticle;

    Long64_t entry;

    // Loop over all events
    for(entry = 0; entry < 10; ++entry)
    {
      // Load selected branches with data from specified event
      treeReader->ReadEntry(entry);

      // plot verizes, If event contains at least 1 vtx
      for(Int_t i = 0; i < branchVtx->GetEntriesFast(); ++i)
      {
        // Take vtx
        vtx = (Vertex*) branchVtx->At(i);

        plots->fvtxX->Fill(vtx->X);
        plots->fvtxY->Fill(vtx->Y);
        plots->fvtxZ->Fill(vtx->Z/1000); // to have it in m
        plots->fvtxT->Fill(vtx->T * 1000000000); // to have it in ns

        plots->fvtxErrorT->Fill(vtx->ErrorT* 1000000000); // vertex position error (t component)
        plots->fvtxErrorX->Fill(vtx->ErrorX); // vertex position error (x component)
        plots->fvtxErrorY->Fill(vtx->ErrorY); // vertex position error (y component)
        plots->fvtxErrorZ->Fill(vtx->ErrorZ); // vertex position error (z component)

        plots->fvtxIndex->Fill(vtx->Index); // vertex index (Int_t)

        plots->fvtxSigma->Fill(vtx->Sigma); // vertex position (z component) error
        plots->fvtxSumPT2->Fill(vtx->SumPT2); // sum pt^2 of tracks attached to the vertex
        plots->fvtxGenSumPT2->Fill(vtx->GenSumPT2); // sum pt^2 of gen tracks attached to the vertex

        plots->fvtxGenDeltaZ->Fill(vtx->GenDeltaZ); // distance in z to closest generated vertex
        plots->fvtxBTVSumPT2->Fill(vtx->BTVSumPT2); // sum pt^2 of tracks attached to the secondary vertex
      }
      // plot vertices of genparticles
      if(branchParticle->GetEntriesFast() > 0){
      for(Int_t i = 0; i < branchParticle->GetEntriesFast(); ++i)
      {
        // Take vtx
        genparticle = (GenParticle*) branchParticle->At(i);

        if (abs(genparticle->IsPU)==0) {
        plots->fgenparticleX[0]->Fill(genparticle->X);
        plots->fgenparticleY[0]->Fill(genparticle->Y);
        plots->fgenparticleZ[0]->Fill(genparticle->Z/ 1000); // to have it in m
        plots->fgenparticleT[0]->Fill(genparticle->T * 1000000000); // to have it in ns
        }
        else{
          plots->fgenparticleX[1]->Fill(genparticle->X);
          plots->fgenparticleY[1]->Fill(genparticle->Y);
          plots->fgenparticleZ[1]->Fill(genparticle->Z / 1000); // to have it in m
          plots->fgenparticleT[1]->Fill(genparticle->T * 1000000000); // to have it in ns
        }
      }
    }

    // plot vertices of genparticles
  //   if(branchPUParticle->GetEntriesFast() > 0){
  //   for(Int_t i = 0; i < branchPUParticle->GetEntriesFast(); ++i)
  //   {
  //     // Take vtx
  //     genparticle = (GenParticle*) branchPUParticle->At(i);
  //
  //     if (abs(genparticle->IsPU)==0) {
  //     plots->fgenparticleX[0]->Fill(genparticle->X);
  //     plots->fgenparticleY[0]->Fill(genparticle->Y);
  //     plots->fgenparticleZ[0]->Fill(genparticle->Z / 1000); // to have it in m
  //     plots->fgenparticleT[0]->Fill(genparticle->T * 1000000000); // to have it in ns
  //     }
  //     else{
  //       plots->fgenparticleX[1]->Fill(genparticle->X);
  //       plots->fgenparticleY[1]->Fill(genparticle->Y);
  //       plots->fgenparticleZ[1]->Fill(genparticle->Z / 1000); // to have it in m
  //       plots->fgenparticleT[1]->Fill(genparticle->T * 1000000000); // to have it in ns
  //     }
  //   }
  // }

    }// end loop over entries
  }//end void

  //------------------------------------------------------------------------------

  void PrintHistograms(ExRootResult *result, MyPlots *plots)
  {
    result->Print("png");
  }

  //------------------------------------------------------------------------------

  void vertexing(const char *inputFile)
  {
    gSystem->Load("libDelphes");

    TChain *chain = new TChain("Delphes");
    chain->Add(inputFile);

    ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
    ExRootResult *result = new ExRootResult();

    MyPlots *plots = new MyPlots;
    BookHistogramsBasic(result, plots);

    SetupGlobalStyle();
    AnalyseEvents(treeReader, plots);
    gSystem->cd("Plots/TimeSmearing/");

    PrintHistograms(result, plots);

    result->Write("results.root");

    cout << "** Exiting..." << endl;

    delete plots;
    delete result;
    delete treeReader;
    delete chain;
  }
