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
  TH1* fgenparticleT;
  TH1* fvtxZ;
  // general track variables
  TH1* ftrackPT[6];
  TH1* ftrackEta[6];

  TH1* ftrackZ[6];
  TH1* ftrackT[6];
  TH1* ftrackTOuter[6];

  TH1* fdeltaT[6][3];
  TH1* fdeltaZ[6][3];

};

//------------------------------------------------------------------------------

class ExRootResult;
class ExRootTreeReader;

//------------------------------------------------------------------------------

void BookHistogramsBasic(ExRootResult *result, MyPlots *plots)
{
  THStack *stack_PU;
  THStack *stack_signal;
  TLegend *legend;
  legend = result->AddLegend(0.25, 0.86, 0.45, 0.98);

  plots->fgenparticleT = result->AddHist1D(
    "genparticle_T", "genparticle T",
    "genparticle T [ns]", "number of genparticles", 100, -1, 1);

  plots->fvtxZ = result->AddHist1D(
    "vtx_z", "vertex Z",
    "vertex Z [m]", "number of vertices", 100, -0.25, 0.25);

  TString name[6] = {"eflowTrack", "RecoPUTrack", "track_merger", "smeared_track", "CHSeflow", "PUPPIeflow"};
  TString IsPU[3] = {"_signal", "_PU", "_inclusive"};

  for (size_t i = 0; i < 6; i++) { // different track and eflow branches
  plots->ftrackPT[i] = result->AddHist1D(
    "pt_"+name[i], "track pt",
    "track pt [GeV]", "number of tracks", 100, 0.0, 10.0);
  plots->ftrackEta[i] = result->AddHist1D(
    "eta_"+name[i], "track eta",
    "track eta ", "number of tracks", 100, -5.0, 5.0);

  plots->ftrackT[i] = result->AddHist1D(
    "T_"+name[i], "track T",
    "track T [ns]", "number of tracks", 100, -1, 1);
    plots->ftrackTOuter[i] = result->AddHist1D(
      "TOuter_"+name[i], "track TOuter",
      "track TOuter [ns]", "number of tracks", 100, 0.0, 30);

      plots->ftrackZ[i] = result->AddHist1D(
        "z_"+name[i], "track Z",
        "track Z [m]", "number of tracks", 100, -0.1, 0.1);

        for (size_t k = 0; k < 3; k++) { // signal or PU tracks or inclusive
      plots->fdeltaT[i][k] = result->AddHist1D(
        "dt_"+name[i]+IsPU[k], "track T - PV T",
        "track T - PV T [ns]", "number of tracks", 100, -1, 1);
      plots->fdeltaZ[i][k] = result->AddHist1D(
        "dz_"+name[i]+IsPU[k], "track Z - vtx Z",
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
    // track collections
    TClonesArray *branchTrack = treeReader->UseBranch("Track");
    TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
    TClonesArray *branchRecoPUTrack = treeReader->UseBranch("RecoPuTrack");
    TClonesArray *branchTrackSmeared = treeReader->UseBranch("TimeSmearedTrack");

    TClonesArray *branchCHSeflow = treeReader->UseBranch("ParticleFlowCandidateCHS");
    TClonesArray *branchPUPPIeflow = treeReader->UseBranch("ParticleFlowCandidate");

    TClonesArray *branchVtx = treeReader->UseBranch("Vertex");

    TClonesArray *branchTracks[4] = {branchEFlowTrack, branchRecoPUTrack, branchTrack, branchTrackSmeared};
    TClonesArray *branchEflow[2] = {branchCHSeflow, branchPUPPIeflow};

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
    for(entry = 0; entry < 10; ++entry)
    {
      // Load selected branches with data from specified event
      treeReader->ReadEntry(entry);

      // check particles
      for (size_t v = 0; v < branchParticle->GetEntriesFast(); v++) {
        genparticle = (GenParticle*) branchParticle->At(v);
        plots->fgenparticleT->Fill(genparticle->T * 1000000000); // to have it in ns
      }

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
      for (size_t i = 0; i < 4; i++) {
        TClonesArray* branch = branchTracks[i];
        for(Int_t k = 0; k < branch->GetEntriesFast(); ++k)
        {
          // Take track
          track = (Track*) branch->At(k);

          plots->ftrackPT[i]->Fill(track->PT);
          plots->ftrackEta[i]->Fill(track->Eta);
          trackT = track->T;
          trackZ = track->Z;
          plots->ftrackZ[i]->Fill(trackZ / 1000); // in m
          plots->ftrackT[i]->Fill(trackT * 1000000000); // to have it in ns
          plots->ftrackTOuter[i]->Fill(track->TOuter * 1000000000); // to have it in ns

           deltaT = (trackT - vtxT)* 1000000000;// to have it in ns
           deltaZ = (trackZ - vtxZ)/ 1000;// to have it in m

           plots->fdeltaT[i][2]->Fill(deltaT );
           plots->fdeltaZ[i][2]->Fill(deltaZ );
           // Get track and associated particle
           Track *t = static_cast<Track*>(branch->At(k));
           GenParticle *p = static_cast<GenParticle*>(t->Particle.GetObject());
           // signal timediff
           if (p->IsPU == 0) {
             plots->fdeltaT[i][0]->Fill(deltaT );
             plots->fdeltaZ[i][0]->Fill(deltaZ );
             //plots->ftrackZ[i][0]->Fill(trackZ / 1000); // in m

           }
           // PU timediff
           else if (p->IsPU == 1) {
             plots->fdeltaT[i][1]->Fill(deltaT );
             plots->fdeltaZ[i][1]->Fill(deltaZ );
             //plots->ftrackZ[i][1]->Fill(trackZ / 1000); // in m

            }

        } // end loop over tracks

      } // end loop over track branches

      for (size_t i = 4; i < 6; i++) {
        TClonesArray* branch = branchEflow[i-4];
        for(Int_t k = 0; k < branch->GetEntriesFast(); ++k)
        {
        // Take track
        ParticleFlowCandidate *t = static_cast<ParticleFlowCandidate*>(branch->At(k));
        GenParticle *p = static_cast<GenParticle*>(t->Particles.At(0));

        if (t->Charge != 0) {
          plots->ftrackPT[i]->Fill(t->PT);
          plots->ftrackEta[i]->Fill(t->Eta);
          trackT = t->T;
          trackZ = t->Z;
          plots->ftrackZ[i]->Fill(trackZ / 1000); // in m
          plots->ftrackT[i]->Fill(trackT * 1000000000); // to have it in ns
          plots->ftrackTOuter[i]->Fill(track->TOuter * 1000000000); // to have it in ns

           deltaT = (trackT - vtxT)* 1000000000;// to have it in ns
           deltaZ = (trackZ - vtxZ)/ 1000;// to have it in m

           plots->fdeltaT[i][2]->Fill(deltaT );
           plots->fdeltaZ[i][2]->Fill(deltaZ );

           // signal timediff
           if (p->IsPU == 0) {
             plots->fdeltaT[i][0]->Fill(deltaT );
             plots->fdeltaZ[i][0]->Fill(deltaZ );

           }
           // PU timediff
           else if (p->IsPU == 1) {
             plots->fdeltaT[i][1]->Fill(deltaT );
             plots->fdeltaZ[i][1]->Fill(deltaZ );
            }
          }
        }

      } // end loop over eflow branches

    }// end loop over entries

  }//end void

  //------------------------------------------------------------------------------

  void PrintHistograms(ExRootResult *result, MyPlots *plots)
  {
    result->Print("png");
  }

  //------------------------------------------------------------------------------

  void track_control_plots(const char *inputFile)
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
    gSystem->cd("PLOTS/10kevents/track_control_plots/");

    PrintHistograms(result, plots);

    result->Write("results.root");

    cout << "** Exiting..." << endl;

    delete plots;
    delete result;
    delete treeReader;
    delete chain;
  }
