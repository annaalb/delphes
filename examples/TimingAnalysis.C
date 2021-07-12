/*
	Prints complete input particle arborescence for the first 100 events. Useful for debugging purposes.
	root -l examples/Example5.C'("delphes_output.root")'
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


	template<typename T>
	void CollectionFilter(const TClonesArray& inColl ,vector<T*>& outColl, Double_t ptMin=30, Double_t etaMax=2.5)
	{

	  const TObject *object;

	  for (Int_t i = 0; i < inColl.GetEntriesFast(); i++)
	  {

	   object = inColl.At(i);
	   const T *t = static_cast<const T*>(object);

	   if(t->P4().Pt() < ptMin) continue;
	   if(TMath::Abs(t->P4().Eta()) > etaMax) continue;
	   if(TMath::Abs(t->PID) != 11 && TMath::Abs(t->PID) != 13) continue;

	   outColl.push_back(t);

	  }
	}


	//------------------------------------------------------------------------------

	void TimingAnalysis(const char *inputFile)
	{
	  gSystem->Load("libDelphes");

	  // Create chain of root trees
	  TChain chain("Delphes");
	  chain.Add(inputFile);

	  // Create object of class ExRootTreeReader
	  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
	  Long64_t numberOfEntries = treeReader->GetEntries();

	  // Get pointers to branches used in this analysis
	  TClonesArray *branchTrack = treeReader->UseBranch("Track");

	  // Book histogram
	  TH1F *hDeltaT = new TH1F("hDeltaT", "hDeltaT", 100, -0.5e-9, 0.5e-9);

	  // define vector of tracks
	  vector<const Track*>     *tracks    = new vector<const Track*>();;

	  // define variables
	  Double_t recoVertexTime, recoTrackTime;
cout << " hi "<<numberOfEntries<< endl;
	  // Loop over all events
	  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
	  {
	    // Load selected branches with data from specified event
	    treeReader->ReadEntry(entry);

	    tracks    -> clear();

	    // Select 4 highest pT tracks with pT > 2 GeV
	    CollectionFilter(*branchTrack, *tracks , 1.0 , 4.5);
			cout << tracks->size() << endl;
	    if(tracks->size() < 4) continue;

	    // given that we have no vertexing, decide that time coordinate given by high pest pT track
	    recoVertexTime = tracks->at(0)->T;

	    // take the time at vertex of the 4th most energetic track
	    recoTrackTime = tracks->at(3)->T;

			cout << "recoVertexTime "<< recoVertexTime << endl;
			cout << "recoTrackTime" << recoTrackTime << endl;

	    // fill the difference
	    hDeltaT -> Fill(recoTrackTime - recoVertexTime);

	  } // end of event loop

	 // Plot the time difference observables

	 hDeltaT ->  Draw();
	 gPad -> SetLogy();

	}
