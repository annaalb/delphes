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
  TH1* ftrackPT[2];
  TH1* ftrackEta[2];

  TH1* ftrackX[2];
  TH1* ftrackY[2];
  TH1* ftrackZ[2];
  TH1* ftrackT[2];
  TH1* ftrackTOuter[2];

  TH1* fdeltaT[2][3][2];
  TH1* fdeltaZ[2][3][2];

  TH1* fdeltaZ_PU_003[2];

};

//------------------------------------------------------------------------------

class ExRootResult;
class ExRootTreeReader;

//------------------------------------------------------------------------------

void BookHistogramsBasic(ExRootResult *result, MyPlots *plots)
{
  plots->fgenparticleT = result->AddHist1D(
    "genparticle_T", "genparticle T",
    "genparticle T [ns]", "number of genparticles", 100, -1, 1);

  plots->fvtxZ = result->AddHist1D(
    "vtx_z", "vertex Z",
    "vertex Z [m]", "number of vertices", 100, -0.25, 0.25);

  TString name[2] = {"no_smearing", "TimeSmearing"};

  for (size_t i = 0; i < 2; i++) {
  plots->ftrackPT[i] = result->AddHist1D(
    "track_pt_"+name[i], "track pt",
    "track pt [GeV]", "number of tracks", 100, 0.0, 10.0);
  plots->ftrackEta[i] = result->AddHist1D(
    "track_eta_"+name[i], "track eta",
    "track eta ", "number of tracks", 100, -5.0, 5.0);

  plots->ftrackT[i] = result->AddHist1D(
    "track_T_"+name[i], "track T",
    "track T [ns]", "number of tracks", 100, -1, 1);
    plots->ftrackTOuter[i] = result->AddHist1D(
      "track_TOuter_"+name[i], "track TOuter",
      "track TOuter [ns]", "number of tracks", 100, 0.0, 30);

        TString IsPU[2] = {"_signal", "_PU"};
        TString category[3] = {"_all", "_matched_PV", "_not_matched"};

      for (size_t j = 0; j < 3; j++) {
        for (size_t k = 0; k < 2; k++) {
      plots->fdeltaT[i][j][k] = result->AddHist1D(
        "track_timediff_"+name[i]+category[j]+IsPU[k], "track T - PV T",
        "track T - PV T [ns]", "number of tracks", 100, -1, 1);
      plots->fdeltaZ[i][j][k] = result->AddHist1D(
        "track_zdiff_"+name[i]+category[j]+IsPU[k], "track Z - vtx Z",
        "track Z - PV Z [m]", "number of tracks", 100, -0.5, 0.5);
      }
    }

  plots->fdeltaZ_PU_003[i] = result->AddHist1D(
    "track_PU_003_Z_minus_PV_Z_"+name[i], "track Z - vtx Z",
    "track Z - PV Z [m]", "number of tracks", 100, -0.25, 0.25);
}

  }

  //------------------------------------------------------------------------------

  void AnalyseEvents(ExRootTreeReader *treeReader, MyPlots *plots)
  {
    Long64_t allEntries = treeReader->GetEntries();

    cout << "** Chain contains " << allEntries << " events" << endl;

    TClonesArray *branchParticle = treeReader->UseBranch("mergerSignalParticle"); // particles are input array to tracks
    TClonesArray *branchTrack = treeReader->UseBranch("Track");
    TClonesArray *branchTrackSmeared = treeReader->UseBranch("TimeSmearing");
    TClonesArray *branchVtx = treeReader->UseBranch("Vertex");

    TClonesArray *branchTracks[2] = {branchTrack, branchTrackSmeared};

    Track *track;
    Candidate *candidate;
    GenParticle *genparticle;
    Vertex *vtx;
    Vertex *matched_vtx;

    Long64_t entry;

    Double_t trackT;
    Double_t vtxT;
    Double_t deltaT;

    Double_t trackZ;
    Double_t vtxZ;
    Double_t deltaZ;

    Double_t sum_eff=0;
    Double_t sum_mis=0;
    Double_t sum_purity=0;

    // Loop over all events
    for(entry = 0; entry < allEntries; ++entry)
    {
      // Load selected branches with data from specified event
      treeReader->ReadEntry(entry);

      // check particles
      for (size_t v = 0; v < branchParticle->GetEntriesFast(); v++) {
        genparticle = (GenParticle*) branchParticle->At(v);
        //cout << genparticle->IsPU << endl;
        // if (genparticle->IsPU == 0) {
        //   cout << genparticle->PID << endl;
        // }
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
      for (size_t i = 0; i < 2; i++) {
        Int_t n_signal = 0;
        Int_t n_PU = 0;
        Int_t nMatchedPV_signal=0;
        Int_t nMatchedPV_PU=0;
        Int_t nMatchedSV_signal=0;
        Int_t nMatchedSV_PU=0;
        Int_t nNotMatched_signal=0;
        Int_t nNotMatched_PU=0;
        TClonesArray* branch = branchTracks[i];
        for(Int_t k = 0; k < branch->GetEntriesFast(); ++k)
        {
          // Take track
          track = (Track*) branch->At(k);

          plots->ftrackPT[i]->Fill(track->PT);
          plots->ftrackEta[i]->Fill(track->Eta);
          trackT = track->T;
          trackZ = track->Z;
          plots->ftrackT[i]->Fill(trackT * 1000000000); // to have it in ns
          plots->ftrackTOuter[i]->Fill(track->TOuter * 1000000000); // to have it in ns
          double trackTOF = track->TOuter - track->T ;

           deltaT = (trackT - vtxT)* 1000000000;// to have it in ns
           deltaZ = (trackZ - vtxZ)/ 1000;// to have it in m
           // Get track and associated particle
           Track *t = static_cast<Track*>(branch->At(k));
           GenParticle *p = static_cast<GenParticle*>(t->Particle.GetObject());
           // cout << p->PID << "\n";
            //cout << p->IsPU << "\n";
           // signal timediff
           if (p->IsPU == 0) {
             ++ n_signal;
             plots->fdeltaT[i][0][0]->Fill(deltaT );
             plots->fdeltaZ[i][0][0]->Fill(deltaZ );
           }
           // PU timediff
           else if (p->IsPU == 1) {
             ++ n_PU;
             plots->fdeltaT[i][0][1]->Fill(deltaT );
             plots->fdeltaZ[i][0][1]->Fill(deltaZ );
            }
            if (p->IsPU == 1 && abs(deltaT) < 0.03){
              plots->fdeltaZ_PU_003[i]->Fill(deltaZ);
            }
            // get closest vertex (matching radius 0.3 cm = 3mm)
            matched_vtx = get_closest_vertex(track, branchVtx, 3);

            //  if (matched_vtx->Index == 0) { // a) matched to PV
            //if (abs(deltaZ) < 0.003) { // b) dz to PV < 0.003 m
            if (abs(deltaZ) < 0.003 && abs(deltaT) < 0.1) { // c) dz to PV < 0.003 m && dt < 0.1 ns
            //if (abs(deltaZ) < 0.002 && abs(deltaT) < 0.1) { // d) dz to PV < 0.001 m && dt < 0.1 ns
            //if (abs(deltaZ) < 0.001 && abs(deltaT) < 0.1) { // e) dz to PV < 0.001 m && dt < 0.1 ns
              if (p->IsPU == 0) {
                ++ nMatchedPV_signal;
                plots->fdeltaT[i][1][0]->Fill(deltaT );
                plots->fdeltaZ[i][1][0]->Fill(deltaZ );
              }
              else if (p->IsPU == 1) {
                ++nMatchedPV_PU;
                plots->fdeltaT[i][1][1]->Fill(deltaT );
                plots->fdeltaZ[i][1][1]->Fill(deltaZ );
               }
            }
            else if (matched_vtx->Index > 0) { // matched to PU vertex
              if (p->IsPU == 0) {
                ++ nMatchedSV_signal;
              }
              else if (p->IsPU == 1) {
                ++nMatchedSV_PU;
               }
            }
            else if (matched_vtx->Index == -1){ // not matched tracks
              if (p->IsPU == 0) {
                ++ nNotMatched_signal;
                plots->fdeltaT[i][2][0]->Fill(deltaT );
                plots->fdeltaZ[i][2][0]->Fill(deltaZ );

              }
              else if (p->IsPU == 1) {
                ++ nNotMatched_PU;
                plots->fdeltaT[i][2][1]->Fill(deltaT );
                plots->fdeltaZ[i][2][1]->Fill(deltaZ );
               }
            }
        } // end loop over tracks
        cout << "--------------------------------------------------" << endl;
        cout << "matched to primary vertex: Signal: "<<nMatchedPV_signal<<" PU: "<<nMatchedPV_PU << endl;
        cout << "matched to PU vertex: Signal: "<<nMatchedSV_signal<<" PU: "<<nMatchedSV_PU << endl;
        cout << "not matched track: Signal "<< nNotMatched_signal << " PU: "<<nNotMatched_PU << endl;

        cout << "Signal tracks " << n_signal << endl;
        cout << "PU tracks " << n_PU << endl;

        Double_t eff = (Double_t) nMatchedPV_signal/n_signal;
        Double_t mis = (Double_t) nMatchedSV_signal/n_signal;
        Double_t purity = (Double_t) nMatchedPV_signal/(nMatchedPV_signal+nMatchedPV_PU);

        cout << "efficiency " << eff << endl;
        cout << "misidentification rate " << mis << endl;
        cout << "Purity " << purity << endl;

        if (i == 1) { // only for time smearing
          sum_eff += eff;
          sum_mis += mis;
          sum_purity += purity;
        }
      } // end loop over time smeared

    }// end loop over entries
    cout << "------------------------------------" << endl;
    cout << "mean efficiency " << sum_eff/allEntries << endl;
    cout << "mean misidentification rate " << sum_mis/allEntries << endl;
    cout << "mean Purity "  << sum_purity/allEntries << endl;

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
    gSystem->cd("Plots/track_vertex_association/");

    PrintHistograms(result, plots);

    result->Write("results.root");

    cout << "** Exiting..." << endl;

    delete plots;
    delete result;
    delete treeReader;
    delete chain;
  }
