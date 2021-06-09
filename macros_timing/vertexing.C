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
  TH1* fvtxX;
  TH1* fvtxY;
  TH1* fvtxZ;
  TH1* fvtxT;
  TH1* fvtxGenDeltaZ;

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
    "vtx X [cm]", "number of vertices", 100, -4.0, 4.0);
  plots->fvtxY = result->AddHist1D(
    "vtx_y", "vtx Y",
    "vtx Y [cm]", "number of vertices", 100, -4.0, 4.0);
  plots->fvtxZ = result->AddHist1D(
    "vtx_z", "vtx Z",
    "vtx Z [cm]", "number of vertices", 100, -2.0, 2.0);
  plots->fvtxT = result->AddHist1D(
    "vtx_T", "vtx T",
    "vtx T [ns]", "number of vertices", 100, -1, 1);
    plots->fvtxGenDeltaZ = result->AddHist1D(
      "vtx_gen_delta_Z", "genvtx Z - vtx Z",
      "genvtx Z - vtx Z [mm]", "number of vertices", 100, -0.00000001, 0.00000001);

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
    name[i]+" Z [cm]", "number of vertices", 100, -2.0, 2.0);
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

    TClonesArray *branchParticle = treeReader->UseBranch("PUParticle");
    //TClonesArray *branchGenvtx = treeReader->UseBranch("GenVertex");
    TClonesArray *branchVtx = treeReader->UseBranch("Vertex");
    //TClonesArray *branchTrack = treeReader->UseBranch("Track");

    Vertex *vtx, *recovtx;
    Track *track;
    GenParticle *genparticle;

    Long64_t entry;

    // Loop over all events
    for(entry = 0; entry < 1; ++entry)
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
        plots->fvtxZ->Fill(vtx->Z / 100); // to have it in cm
        plots->fvtxT->Fill(vtx->T * 1000000000); // to have it in ns
        plots->fvtxGenDeltaZ->Fill(vtx->GenDeltaZ);
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
        plots->fgenparticleZ[0]->Fill(genparticle->Z / 100); // to have it in cm
        plots->fgenparticleT[0]->Fill(genparticle->T * 1000000000); // to have it in ns
        }
        else{
          plots->fgenparticleX[1]->Fill(genparticle->X);
          plots->fgenparticleY[1]->Fill(genparticle->Y);
          plots->fgenparticleZ[1]->Fill(genparticle->Z / 100); // to have it in cm
          plots->fgenparticleT[1]->Fill(genparticle->T * 1000000000); // to have it in ns
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
    gSystem->cd("Plots/vertexing/one_event/");

    PrintHistograms(result, plots);

    result->Write("results.root");

    cout << "** Exiting..." << endl;

    delete plots;
    delete result;
    delete treeReader;
    delete chain;
  }
