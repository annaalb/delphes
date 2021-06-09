/*
Simple macro showing how to access branches from the delphes output root file,
loop over events, and plot simple quantities such as the jet pt and the di-electron invariant
mass.

root -l examples/Example1.C'("delphes_output.root")'
*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#else
class ExRootTreeReader;
class ExRootResult;
#endif

//------------------------------------------------------------------------------

void test_macro(const char *inputFile)
{
  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchEvent = treeReader->UseBranch("Event");

  TCanvas* canvas = new TCanvas("c", "c", 600, 600);
  // Book histograms
  TH1 *histJetPT = new TH1F("jet_pt", "jet P_{T}", 100, 0.0, 1000.0);
  TH1 *histJetEta = new TH1F("jet_eta", "jet #eta", 100, 0.0, 5.0);

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    // If event contains at least 1 jet
    if(branchJet->GetEntries() > 0)
    {
      for (size_t k = 0; k < branchJet->GetEntries(); k++) {// loop over jets
        Jet *jet = (Jet*) branchJet->At(k);
        // Plot jet transverse momentum
        histJetPT->Fill(jet->PT); //TODO draw in one plot PU and VBF jets
        histJetEta->Fill(jet->Eta);

      }
    }
  }

  // Show resulting histograms
  histJetEta->Draw();

canvas->SaveAs("Plots/test_VBF_jet_eta.eps");


}
