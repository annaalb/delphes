/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/** \class TrackPileUpSubtractor
 *
 *  Subtract pile-up contribution from tracks.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/TrackPileUpSubtractor.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootResult.h"

#include "TDatabasePDG.h"
#include "TFormula.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "TString.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

//------------------------------------------------------------------------------

TrackPileUpSubtractor::TrackPileUpSubtractor() :
  fFormula(0), fFormulaTime(0)
{
  fFormula = new DelphesFormula;
  fFormulaTime = new DelphesFormula;
}

//------------------------------------------------------------------------------

TrackPileUpSubtractor::~TrackPileUpSubtractor()
{
  if(fFormula) delete fFormula;
  if(fFormulaTime) delete fFormulaTime;

}

//------------------------------------------------------------------------------

void TrackPileUpSubtractor::Init()
{
  // import input array

  fVertexInputArray = ImportArray(GetString("VertexInputArray", "PileUpMerger/vertices"));
  fItVertexInputArray = fVertexInputArray->MakeIterator();

  // read resolution formula in m
  fFormula->Compile(GetString("ZVertexResolution", "0.001"));
  fFormulaTime->Compile(GetString("TVertexResolution", "0"));
  fEtaMax = GetDouble("EtaMax", 0.);

  fPTMin = GetDouble("PTMin", 0.);

  // import arrays with output from other modules

  ExRootConfParam param = GetParam("InputArray");
  Long_t i, size;
  const TObjArray *array;
  TIterator *iterator;

  size = param.GetSize();
  for(i = 0; i < size / 2; ++i)
  {
    array = ImportArray(param[i * 2].GetString());
    iterator = array->MakeIterator();

    fInputMap[iterator] = ExportArray(param[i * 2 + 1].GetString());
  }
}

//------------------------------------------------------------------------------

void TrackPileUpSubtractor::Finish()
{
  map<TIterator *, TObjArray *>::iterator itInputMap;
  TIterator *iterator;

  for(itInputMap = fInputMap.begin(); itInputMap != fInputMap.end(); ++itInputMap)
  {
    iterator = itInputMap->first;

    if(iterator) delete iterator;
  }

  if(fItVertexInputArray) delete fItVertexInputArray;
}

//------------------------------------------------------------------------------

void TrackPileUpSubtractor::Process()
{
  Candidate *candidate, *particle;
  map<TIterator *, TObjArray *>::iterator itInputMap;
  TIterator *iterator;
  TObjArray *array;
  Double_t z, zvtx = 0;
  Double_t t, tvtx = 0;
  Double_t pt, eta, phi, e;
  const Double_t c_light = 2.99792458E8;

  // find z position of primary vertex

  fItVertexInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItVertexInputArray->Next())))
  {
    if(!candidate->IsPU)
    {
      zvtx = candidate->Position.Z();
      tvtx = candidate->Position.T() * 1.0E-3 / c_light;
      // break;
    }
  }

  // loop over all input arrays
  for(itInputMap = fInputMap.begin(); itInputMap != fInputMap.end(); ++itInputMap)
  {
    iterator = itInputMap->first;
    array = itInputMap->second;

    // loop over all candidates
    iterator->Reset();
    while((candidate = static_cast<Candidate *>(iterator->Next())))
    {
      particle = static_cast<Candidate *>(candidate->GetCandidates()->At(0));
      const TLorentzVector &candidateMomentum = particle->Momentum;

      eta = candidateMomentum.Eta();
      pt = candidateMomentum.Pt();
      phi = candidateMomentum.Phi();
      e = candidateMomentum.E();


      z = particle->Position.Z();

      // apply pile-up subtraction
      // assume perfect pile-up subtraction for tracks outside fZVertexResolution

    //  if(candidate->Charge != 0 && candidate->IsPU && TMath::Abs(z - zvtx) > fFormula->Eval(pt, eta, phi, e) * 1.0e3)
      if(candidate->Charge != 0 && TMath::Abs(z - zvtx) > fFormula->Eval(pt, eta, phi, e) * 1.0e3)
      {
        candidate->IsRecoPU = 1;
      }
      else
      {
        candidate->IsRecoPU = 0;
        if(candidate->Momentum.Pt() > fPTMin) array->Add(candidate);
      }


      /// New part !
      Bool_t _debug=false;

      // z = candidate->InitialPosition.Z();
      // t = candidate->InitialPosition.T()* 1.0E-3 / c_light; // time in s
      //
      // // apply pile-up subtraction
      // // assume perfect pile-up subtraction for tracks outside fZVertexResolution
      // if (_debug && candidate->Charge != 0 && (abs(eta)<0||abs(eta)>4)) {
      //   std::cout << "TVertex Resolution "<< GetString("TVertexResolution", "0") << '\n';
      //   std::cout << "---- Module TrackPileUpSubtractor ----" << '\n';
      //   std::cout << "z = "<<z << '\n';
      //   std::cout << "zvtx = " << zvtx << '\n';
      //   std::cout << "abs(z - zvtx) = "<<TMath::Abs(z - zvtx) << '\n';
      //   std::cout << "fFormula->Eval(pt, eta, phi, e) * 1.0e3 = "<< fFormula->Eval(pt, eta, phi, e) * 1.0e3 << '\n';
      //   std::cout << "t = "<<t << '\n';
      //   std::cout << "tvtx = " << tvtx << '\n';
      //   std::cout << "abs(t - tvtx) = "<<TMath::Abs(t - tvtx) << '\n';
      //   std::cout << "fFormulaTime->Eval(pt, eta, phi, e) = "<< fFormulaTime->Eval(pt, eta, phi, e) << '\n';
      //   std::cout << "candidate->Momentum.Eta() = " << candidate->Momentum.Eta() << '\n';
      //   std::cout << "eta = candidate->particle->momentum.Eta() = "<< eta << '\n';
      //   std::cout << "EtaMax = " << fEtaMax << '\n';
      //   std::cout << "Eta diff = " << (TMath::Abs(eta)-fEtaMax) << '\n';
      // }



      // Bool_t charged = (candidate->Charge != 0);
      // Bool_t dz_smaller = (TMath::Abs(z - zvtx) < fFormula->Eval(pt, eta, phi, e) * 1.0e3);
      // Bool_t dt_smaller = (TMath::Abs(t - tvtx) < fFormulaTime->Eval(pt, eta, phi, e) );
      // Bool_t eta_smaller = (TMath::Abs(eta)<fEtaMax);
      //
      // if (TMath::Abs(eta)>=4) {eta_smaller=false;}
      //
      // // handle two eta regions:
      // // 1. (eta > etaMax) keep tracks if dz < x
      // // 2. (eta < etaMax) keep tracks if dz < x and dt < y
      // if(charged && ((!eta_smaller && dz_smaller) || (eta_smaller && dz_smaller && dt_smaller)) ) // include also signal track rejection // default Eta max = 0 (only dz cut)
      // {
      //   // if (_debug) {
      //   //      std::cout << "******************************** IsRecoPU = 0 ******************************" << '\n';
      //   //      std::cout << "charged "<< charged << '\n';
      //   //      std::cout << "eta_smaller "<< eta_smaller << '\n';
      //   //      std::cout << "dz_smaller "<< dz_smaller << '\n';
      //   //      std::cout << "dt_smaller "<< dt_smaller << '\n';
      //   //
      //   //    }
      //   candidate->IsRecoPU = 0;
      //   if(candidate->Momentum.Pt() > fPTMin) array->Add(candidate);
      // }
      // else{
      //   candidate->IsRecoPU = 1; // tracks that are rejected
      // }


    }
  }
}
