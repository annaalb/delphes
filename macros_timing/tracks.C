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
  TH1* ftrackPT[4];
  TH1* ftrackEta[4];

  TH1* ftrackZ[4];
  TH1* ftrackT[4];

  TH1* fdeltaT[4][4][2];
  TH1* fdeltaZ[4][4][2];

  TH1* eff;
  TH1* mis;
  TH1* purity;


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

  TString name[4] = {"smeared", "RecoPUTrack", "", "eflowTracks"};
  TString IsPU[2] = {"_signal", "_PU"};
  TString category[4] = {"_all", "_matched_PV", "_matched_PU", "_not_matched"};

  Int_t i = 0;
  //for (size_t i = 0; i < 4; i++) { // w/o and with timeSmearing applied
  plots->ftrackPT[i] = result->AddHist1D(
    "track_pt_"+name[i], "track pt",
    "track pt [GeV]", "number of tracks", 100, 0.0, 10.0);
  plots->ftrackEta[i] = result->AddHist1D(
    "track_eta_"+name[i], "track eta",
    "track eta ", "number of tracks", 100, -5.0, 5.0);

    plots->ftrackZ[i] = result->AddHist1D(
      "track_z_"+name[i], "track Z",
      "track Z [m]", "number of tracks", 1000, -0.1, 0.1);
  plots->ftrackT[i] = result->AddHist1D(
    "track_T_"+name[i], "track T",
    "track T [ns]", "number of tracks", 1000, -1, 1);

      for (size_t j = 0; j < 4; j++) { // matching categories
        for (size_t k = 0; k < 2; k++) { // signal or PU tracks
      plots->fdeltaT[i][j][k] = result->AddHist1D(
        "track_dt_"+name[i]+category[j]+IsPU[k], "track T - PV T",
        "track T - PV T [ns]", "number of tracks", 1000, -1, 1);
      plots->fdeltaZ[i][j][k] = result->AddHist1D(
        "track_dz_"+name[i]+category[j]+IsPU[k], "track Z - vtx Z",
        "track Z - PV Z [m]", "number of tracks", 1000, -0.1, 0.1);
      }
    }

    plots->eff = result->AddHist1D(
      "efficiency", "eff", "", "efficiency", 100, 0.0, 1);
      plots->mis = result->AddHist1D(
        "mis", "mis", "","mis", 100, 0.0, 1);
        plots->purity = result->AddHist1D(
          "purity", "purity", "","purity", 100, 0.0, 1);

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

    TClonesArray *branchVtx = treeReader->UseBranch("Vertex");

    TClonesArray *branchTracks[4] = {branchTrackSmeared, branchRecoPUTrack, branchTrack, branchEFlowTrack};

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

    Double_t dzsig;
    Double_t dtsig;

    Double_t sum_eff=0;
    Double_t sum_mis=0;
    Double_t sum_purity=0;

    Long64_t sum_entries;
    sum_entries = 100;
    // Loop over all events
    for(entry = 0; entry < sum_entries; ++entry)
    {
      cout << "Process event "<< entry << " of " << sum_entries << endl;
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
      Int_t i = 0;
      //for (size_t i = 0; i < 2; i++) {
        Int_t n_signal = 0;
        Int_t n_PU = 0;
        Int_t nMatchedPV_signal=0;
        Int_t nMatchedPV_PU=0;
        Int_t nMatchedSV_signal=0;
        Int_t nMatchedSV_PU=0;
        Int_t nNotMatched_signal=0;
        Int_t nNotMatched_PU=0;
        TClonesArray* branch = branchTracks[0];
        for(Int_t k = 0; k < branch->GetEntriesFast(); ++k)
        {
          // Take track
          track = (Track*) branch->At(k);

          plots->ftrackPT[i]->Fill(track->PT);
          plots->ftrackEta[i]->Fill(track->Eta);
          trackT = track->T;
          trackZ = track->Z;
          plots->ftrackT[i]->Fill(trackT * 1000000000); // to have it in ns

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

             //get closest vertex (matching radius 0.3 cm = 3mm)
             matched_vtx = get_closest_vertex(track, branchVtx, 3); // a) get closest vertex
              if (matched_vtx->Index == 0) { // a) matched to PV
          //  if (abs(deltaZ) < 0.003) { // b) dz to PV < 0.003 m
            //if (abs(deltaZ) < 0.003 && abs(deltaT) < 0.1) { // c) dz to PV < 0.003 m && dt < 0.1 ns
            //if (abs(deltaZ) < 0.002 && abs(deltaT) < 0.1) { // d) dz to PV < 0.002 m && dt < 0.1 ns
            //if (abs(deltaZ) < 0.001 && abs(deltaT) < 0.1) { // e) dz to PV < 0.001 m && dt < 0.1 ns
            //if (abs(deltaZ) < 10) { // f) dz to PV < 0.1 m

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
            // associate all other tracks to PU (no unassociated tracks left)
            //else{
              else if (matched_vtx->Index > 0) { // matched to PU vertex
              if (p->IsPU == 0) { // signal tracks
                ++ nMatchedSV_signal;
                plots->fdeltaT[i][2][0]->Fill(deltaT ); //plot dt for signal matched to PU
                plots->fdeltaZ[i][2][0]->Fill(deltaZ ); // dz
              }
              else if (p->IsPU == 1) { // PU tracks
                ++nMatchedSV_PU;
                plots->fdeltaT[i][2][1]->Fill(deltaT ); //plot dt for PU matched to PU
                plots->fdeltaZ[i][2][1]->Fill(deltaZ ); // dz
               }
            }
        } // end loop over tracks
         cout << "--------------------------------------------------" << endl;
        cout << "matched to primary vertex: Signal: "<<nMatchedPV_signal<<" PU: "<<nMatchedPV_PU << endl;
        cout << "matched to PU vertex: Signal: "<<nMatchedSV_signal<<" PU: "<<nMatchedSV_PU << endl;
        cout << "not matched track: Signal "<< nNotMatched_signal << " PU: "<<nNotMatched_PU << endl;
        //
        // cout << "Signal tracks " << n_signal << endl;
        // cout << "PU tracks " << n_PU << endl;

        Double_t eff = (Double_t) nMatchedPV_signal/n_signal;
        Double_t mis = (Double_t) nMatchedSV_signal/n_signal;
        Double_t purity = (Double_t) nMatchedPV_signal/(nMatchedPV_signal+nMatchedPV_PU);

          sum_eff += eff;
          sum_mis += mis;
          sum_purity += purity;
      //} // end loop over time smeared

    }// end loop over entries
    cout << "------------------------------------" << endl;
    cout << "mean efficiency " << sum_eff/sum_entries << endl;
    cout << "mean misidentification rate " << sum_mis/sum_entries << endl;
    cout << "mean Purity "  << sum_purity/sum_entries << endl;

    plots->eff->Fill(sum_eff/sum_entries);
    plots->mis->Fill(sum_mis/sum_entries);
    plots->purity->Fill(sum_purity/sum_entries);

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
    gSystem->cd("PLOTS/10kevents/track_vertex_association/scenario_a/");

    PrintHistograms(result, plots);

    result->Write("results.root");

    cout << "** Exiting..." << endl;

    delete plots;
    delete result;
    delete treeReader;
    delete chain;
  }
