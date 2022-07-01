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
  TH1* ftrackPT[7];
  TH1* ftrackEta[7];

  TH1* fdeltaT[6][4][3];
  TH1* fdeltaZ[6][4][3];

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

  TString IsPU[3] = {"_signal", "_PU", "_inclusive"};
  TString etabin[6] = {"", "_eta01p3", "_eta1p32", "_eta23", "_eta34", "_eta4"};
  TString ptbin[4] = {"_pt0to1","_pt1to2","_pt2to10","_pt10toInf"};

  for (size_t k = 0; k < 4; k++) { // different eta bins
  plots->ftrackPT[k] = result->AddHist1D(
    "No_signal_rej_pt"+ptbin[k], "track pt",
    "track pt [GeV]", "number of tracks", 200, 0.0, 20.0);
  plots->ftrackEta[k] = result->AddHist1D(
    "No_signal_rej_eta"+ptbin[k], "track eta",
    "track eta ", "number of tracks", 100, -5.0, 5.0);
    for (size_t e = 0; e < 6; e++) { // pt bins
      for (size_t i = 0; i < 3; i++) {  // signal or PU tracks or inclusive
      plots->fdeltaZ[e][k][i] = result->AddHist1D(
        "No_signal_rej_dz"+IsPU[i]+etabin[e]+ptbin[k], "track Z - vtx Z",
        "track Z - LV Z [m]", "number of tracks", 100, -0.002, 0.002);
      }
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
    TClonesArray *branchPUPPIeflow = treeReader->UseBranch("ParticleFlowCandidate"); // input to PUPPI jets

    TClonesArray *branchVtx = treeReader->UseBranch("Vertex");

    TClonesArray *branchEflow[1] = {branchPUPPIeflow};

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
    for(entry = 0; entry < 1000; ++entry)
    {
      cout << "Process event " << entry+1 << " of total " << allEntries << endl;
      // Load selected branches with data from specified event
      treeReader->ReadEntry(entry);

      vtxT = 0;
      // get PV
      for (size_t v = 0; v < branchVtx->GetEntriesFast(); v++) {
        vtx = (Vertex*) branchVtx->At(v);
        if (0 == vtx->Index) {
         vtxT = vtx->T;
         vtxZ = vtx->Z;
        }
      }

        TClonesArray* branch = branchEflow[0];
        for(Int_t k = 0; k < branch->GetEntriesFast(); ++k)
        {
        // Take track
        ParticleFlowCandidate *t = static_cast<ParticleFlowCandidate*>(branch->At(k));
        GenParticle *part = static_cast<GenParticle*>(t->Particles.At(0));

        if (t->Charge != 0) {

          Int_t e;
          Int_t p;
          if (abs(t->Eta) < 1.3) {e=1;};
          if (abs(t->Eta) > 1.3 && abs(t->Eta) < 2) {e=2;};
          if (2 < abs(t->Eta) && abs(t->Eta) < 3) {e=3;};
          if (3 < abs(t->Eta) && abs(t->Eta) < 4) {e=4;};
          if (abs(t->Eta) > 4) {e=5;};

          if (abs(t->PT) > 0 && abs(t->PT) < 1) {p=0;};
          if (abs(t->PT) > 1 && abs(t->PT) < 2) {p=1;};
          if (abs(t->PT) > 2 && abs(t->PT) < 10) {p=2;};
          if (abs(t->PT) > 10) {p=3;};

          // inclusive
          plots->ftrackPT[p]->Fill(t->PT);
          plots->ftrackEta[p]->Fill(t->Eta);
          trackZ = t->Z;

           deltaZ = (trackZ - vtxZ)/ 1000;// to have it in m

           plots->fdeltaZ[e][p][2]->Fill(deltaZ );
           plots->fdeltaZ[0][p][2]->Fill(deltaZ );
           if (part->IsPU == 0) {
             plots->fdeltaZ[e][p][0]->Fill(deltaZ );
             plots->fdeltaZ[0][p][0]->Fill(deltaZ );
           }
           else if (part->IsPU == 1) {
             plots->fdeltaZ[e][p][1]->Fill(deltaZ );
             plots->fdeltaZ[0][p][1]->Fill(deltaZ );
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

// root -l -b -q macros_timing/track_control_plots.C'("ROOTOUTPUT/VBF/EtaMax3/CMS_PhaseII_Snowmass_200PU_VBF_10k_ptmin10_EtaMax3.root", "PLOTS/study_timing_cut/VBF/track_control_plots_EtaMax3")'
  void Snowmass_origin_control_plots(const char *inputFile, const char *output)
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
    gSystem->cd(output);

    PrintHistograms(result, plots);

    result->Write("results.root");

    cout << "** Exiting..." << endl;

    delete plots;
    delete result;
    delete treeReader;
    delete chain;
  }
