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
  TH1* fvtxZ;
  // general track variables
  TH1* ftrackPT[2];
  TH1* ftrackEta[2];

  TH1* fdeltaT[2][2];
  TH1* fdeltaZ[2][2];

};

//------------------------------------------------------------------------------

class ExRootResult;
class ExRootTreeReader;

//------------------------------------------------------------------------------

void BookHistogramsBasic(ExRootResult *result, MyPlots *plots)
{

  plots->fvtxZ = result->AddHist1D(
    "vtx_z", "vertex Z",
    "vertex Z [m]", "number of vertices", 100, -0.25, 0.25);

  TString name[2] = {"PUPPITrack", "PUPPIParticle"};

  for (size_t i = 0; i < 2; i++) { // Tracks or PUPPI particles
  plots->ftrackPT[i] = result->AddHist1D(
    "track_pt_"+name[i], "track pt",
    "track pt [GeV]", "number of tracks", 100, 0.0, 10.0);
  plots->ftrackEta[i] = result->AddHist1D(
    "track_eta_"+name[i], "track eta",
    "track eta ", "number of tracks", 100, -5.0, 5.0);

        TString IsPU[2] = {"_signal", "_PU"};

        for (size_t k = 0; k < 2; k++) { // signal or PU tracks
      plots->fdeltaT[i][k] = result->AddHist1D(
        "track_dt_"+name[i]+IsPU[k], "track T - PV T",
        "track T - PV T [ns]", "number of tracks", 100, -1, 1);
      plots->fdeltaZ[i][k] = result->AddHist1D(
        "track_dz_"+name[i]+IsPU[k], "track Z - vtx Z",
        "track Z - PV Z [m]", "number of tracks", 100, -0.1, 0.1);
      }
    }
  }

  //------------------------------------------------------------------------------

  void AnalyseEvents(ExRootTreeReader *treeReader, MyPlots *plots)
  {
    Long64_t allEntries = treeReader->GetEntries();

    cout << "** Chain contains " << allEntries << " events" << endl;

    TClonesArray *branchParticle = treeReader->UseBranch("mergerParticle"); // particles are input array to tracks
    TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");

    TClonesArray *branchPUPPITrack = treeReader->UseBranch("PuppiTrack");
    TClonesArray *branchPUPPIParticle = treeReader->UseBranch("ParticleFlowCandidate");

    TClonesArray *branchVtx = treeReader->UseBranch("Vertex");


    Track *track;
    Candidate *candidate;
    GenParticle *genparticle;
    Vertex *vtx;

    Long64_t entry;

    Double_t trackT;
    Double_t vtxT;
    Double_t deltaT;

    Double_t trackZ;
    Double_t vtxZ;
    Double_t deltaZ;

    // Loop over all events
    for(entry = 0; entry < allEntries; ++entry)
    {
      // Load selected branches with data from specified event
      treeReader->ReadEntry(entry);

      vtxT = 0;
      // get PV
      for (size_t v = 0; v < branchVtx->GetEntriesFast(); v++) {
        vtx = (Vertex*) branchVtx->At(v);
        plots->fvtxZ->Fill(vtx->Z / 1000);
        if (0 == vtx->Index) {
         vtxT = vtx->T;
         vtxZ = vtx->Z;
        }
      }

      // get tracks
        for(Int_t k = 0; k < branchPUPPITrack->GetEntriesFast(); ++k)
        {
          // Take track
          track = (Track*) branchPUPPITrack->At(k);
          Track *t = static_cast<Track*>(branchPUPPITrack->At(k));
          GenParticle *p = static_cast<GenParticle*>(t->Particle.GetObject());

          plots->ftrackPT[0]->Fill(t->PT);
          plots->ftrackEta[0]->Fill(t->Eta);
          trackT = t->T;
          trackZ = t->Z;

           deltaT = (trackT - vtxT)* 1000000000;// to have it in ns
           deltaZ = (trackZ - vtxZ)/ 1000;// to have it in m
           // Get track and associated particle

           // signal timediff
           if (p->IsPU == 0) {
             plots->fdeltaT[0][0]->Fill(deltaT );
             plots->fdeltaZ[0][0]->Fill(trackZ /1000 );
           }
           // PU timediff
           else if (p->IsPU == 1) {
             plots->fdeltaT[0][1]->Fill(deltaT );
             plots->fdeltaZ[0][1]->Fill(trackZ /1000 );
            }

        } // end loop over tracks

        for(Int_t k = 0; k < branchPUPPIParticle->GetEntriesFast(); ++k)
        {
          // Take track
          ParticleFlowCandidate *t = static_cast<ParticleFlowCandidate*>(branchPUPPIParticle->At(k));

          plots->ftrackPT[1]->Fill(t->PT);
          plots->ftrackEta[1]->Fill(t->Eta);
          trackT = t->T;
          trackZ = t->Z;

           deltaT = (trackT - vtxT)* 1000000000;// to have it in ns
           deltaZ = (trackZ - vtxZ)/ 1000;// to have it in m

           plots->fdeltaT[1][1]->Fill(deltaT );
           plots->fdeltaZ[1][1]->Fill(deltaZ );

        } // end loop over tracks


    }// end loop over entries

  }//end void

  //------------------------------------------------------------------------------

  void PrintHistograms(ExRootResult *result, MyPlots *plots)
  {
    result->Print("png");
  }

  //------------------------------------------------------------------------------

  void tracks_PUPPI(const char *inputFile)
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
    gSystem->cd("Plots/10kevents/tracks_PUPPI/");

    PrintHistograms(result, plots);

    result->Write("results.root");

    cout << "** Exiting..." << endl;

    delete plots;
    delete result;
    delete treeReader;
    delete chain;
  }
