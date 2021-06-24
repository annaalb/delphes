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
  // general track variables
  TH1* ftrackPT[2];
  TH1* ftrackEta[2];
  TH1* ftrackMass[2];

  TH1* ftrackX[2];
  TH1* ftrackY[2];
  TH1* ftrackZ[2];
  TH1* ftrackT[2];
  TH1* fdeltaT[2];
  TH1* fdeltaZ[2];
  TH1* fTres[2];

  TH1* ftrackErrorT[2]; // track position error (t component)
  // TH1* ftrackErrorX; // track position error (x component)
  // TH1* ftrackErrorY; // track position error (y component)
  // TH1* ftrackErrorZ; // track position error (z component)

};

//------------------------------------------------------------------------------

class ExRootResult;
class ExRootTreeReader;

//------------------------------------------------------------------------------

void BookHistogramsBasic(ExRootResult *result, MyPlots *plots)
{
  TString name[2] = {"no_smearing", "TimeSmearing"};

  for (size_t i = 0; i < 2; i++) {
  plots->ftrackPT[i] = result->AddHist1D(
    "track_pt_"+name[i], "track pt",
    "track pt [GeV]", "number of tracks", 100, 0.0, 10.0);
  plots->ftrackEta[i] = result->AddHist1D(
    "track_eta_"+name[i], "track eta",
    "track eta ", "number of tracks", 100, -5.0, 5.0);
  plots->ftrackMass[i] = result->AddHist1D(
    "track_mass_"+name[i], "track mass",
    "track mass [GeV]", "number of tracks", 100, 0.0, 1.0);

  plots->ftrackX[i] = result->AddHist1D(
    "track_x_"+name[i], "track X",
    "track X [mm]", "number of tracks", 100, -2.0, 2.0);
  plots->ftrackY[i] = result->AddHist1D(
    "track_y_"+name[i], "track Y",
    "track Y [mm]", "number of tracks", 100, -2.0, 2.0);
  plots->ftrackZ[i] = result->AddHist1D(
    "track_z_"+name[i], "track Z",
    "track Z [m]", "number of tracks", 100, -0.25, 0.25);
  plots->ftrackT[i] = result->AddHist1D(
    "track_T_"+name[i], "track T",
    "track T [ns]", "number of tracks", 100, -1, 1);

  plots->fdeltaT[i] = result->AddHist1D(
    "track_T_minus_vtx_T_"+name[i], "track T - vtx T",
    "track T - PV T [ns]", "number of tracks", 100, -1, 1);
    plots->fTres[i] = result->AddHist1D(
      "track_T_res_"+name[i], "track T - vtx T / Tres",
      "track T - PV T / T_{res} [ns]", "", 100, -1, 1);
  plots->fdeltaZ[i] = result->AddHist1D(
    "track_Z_minus_vtx_Z_"+name[i], "track Z - vtx Z",
    "track Z - PV Z [m]", "number of tracks", 100, -0.25, 0.25);

  plots->ftrackErrorT[i] = result->AddHist1D(
    "track_error_t_"+name[i], "track error t",
    "track error T [ns]", "number of tracks", 100, 0, 0.01); // track position error (t component)
  // plots->ftrackErrorX = result->AddHist1D(
  //   "track_error_x", "track error X",
  //   "track error X [mm]", "number of tracks", 100, 0.0, 0.01); // track position error (x component)
  // plots->ftrackErrorY = result->AddHist1D(
  //   "track_error_y", "track error y",
  //   "track error Y [mm]", "number of tracks", 100, 0.0, 0.01); // track position error (y component)
  // plots->ftrackErrorZ = result->AddHist1D(
  //   "track_error_z", "track error z",
  //   "track error Z [mm]", "number of tracks", 100, 0, 0.01); // track position error (z component)
}

  }

  //------------------------------------------------------------------------------

  void AnalyseEvents(ExRootTreeReader *treeReader, MyPlots *plots)
  {
    Long64_t allEntries = treeReader->GetEntries();

    cout << "** Chain contains " << allEntries << " events" << endl;

    TClonesArray *branchTrack = treeReader->UseBranch("Track");
    TClonesArray *branchTrackSmeared = treeReader->UseBranch("TrackSmeared");
    TClonesArray *branchVtx = treeReader->UseBranch("Vertex");

    TClonesArray *branchTracks[2] = {branchTrack, branchTrackSmeared};

    Track *track;
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

      for (size_t i = 0; i < 2; i++) {
        TClonesArray* branch = branchTracks[i];
        for(Int_t k = 0; k < branch->GetEntriesFast(); ++k)
        {
          // Take track
          track = (Track*) branch->At(k);

          plots->ftrackPT[i]->Fill(track->PT);
          plots->ftrackEta[i]->Fill(track->Eta);
          plots->ftrackMass[i]->Fill(track->Mass);

          plots->ftrackX[i]->Fill(track->X);
          plots->ftrackY[i]->Fill(track->Y);
          plots->ftrackZ[i]->Fill(track->Z/1000); // to have it in m
          plots->ftrackT[i]->Fill(track->T * 1000000000); // to have it in ns

          plots->ftrackErrorT[i]->Fill(track->ErrorT* 1000000000); // track position error (t component)

           trackT = track->T;
           trackZ = track->Z;
           vtxT = 0;
           // cout << "track-> VertexIndex" << track->VertexIndex << endl;
           // cout << "track T" << trackT << endl;
           // cout << "track Z" << trackZ << endl;
          for (size_t v = 0; v < branchVtx->GetEntriesFast(); v++) {
            vtx = (Vertex*) branchVtx->At(v);

            if (0 == vtx->Index) {
              vtxT = vtx->T;
              vtxZ = vtx->Z;
              // cout << "vtx" << vtx->Index << endl;
              // cout << "vtx T" << vtxT << endl;
              // cout << "vtx Z" << vtxZ << endl;

            }
          }
           deltaT = trackT - vtxT;
           deltaZ = trackZ - vtxZ;
          plots->fdeltaT[i]->Fill(deltaT * 1000000000); // to have it in ns
          plots->fdeltaZ[i]->Fill(deltaZ / 1000); // to have it in m
          plots->fTres[i]->Fill((deltaT * 1000000000)/0.03); // to have it in ns


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

  void tracks(const char *inputFile)
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
    gSystem->cd("Plots/TimeSmearing/tracks");

    PrintHistograms(result, plots);

    result->Write("results.root");

    cout << "** Exiting..." << endl;

    delete plots;
    delete result;
    delete treeReader;
    delete chain;
  }
