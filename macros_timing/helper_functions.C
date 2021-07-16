#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#endif

//-----------Helper functions-------------------------------------------------------------------

double get_distance(Jet *jet, GenParticle* genpart)
{
    double distance;
    double Delta_Eta;
    double Delta_Phi;
    Delta_Eta = jet->Eta - genpart->Eta;
    Delta_Phi = jet->Phi - genpart->Phi;
    distance = sqrt(pow(Delta_Eta,2) + pow(Delta_Phi,2));
    return distance;
}

double get_distance(Jet *jet1, Jet *jet2)
{
    double distance;
    double Delta_Eta;
    double Delta_Phi;
    Delta_Eta = jet1->Eta - jet2->Eta;
    Delta_Phi = jet1->Phi - jet2->Phi;
    distance = sqrt(pow(Delta_Eta,2) + pow(Delta_Phi,2));
    return distance;
}

double get_distance(Vertex *vertex, Track *track)
{
    double distance;
    distance = abs(vertex->Z - track->Z); // dz in mm
    return distance;
}

Jet* get_closest_jet(TClonesArray *branchJet, GenParticle* genpart, double matching_radius)
{
    Jet *matched_jet;
    Jet *jet;
    double distance;
    if(branchJet->GetEntriesFast() > 0)
    {
    for (size_t k = 0; k < branchJet->GetEntriesFast(); k++) { // loop over jets
        jet = (Jet*) branchJet->At(k);
        distance = get_distance(jet, genpart);
        if (distance < matching_radius) { // is the jet closer than the previous one?
            matched_jet = jet;
            matching_radius = distance; // new distance to matched jet
        }
    }// end loop over jets
    }
    return matched_jet;
}

Jet* get_closest_jet(TClonesArray *branchJet, Jet* jet_in, double matching_radius)
{
    Jet *matched_jet;
    Jet *jet;
    double distance;
    if(branchJet->GetEntriesFast() > 0)
    {
    for (size_t k = 0; k < branchJet->GetEntriesFast(); k++) { // loop over jets
        jet = (Jet*) branchJet->At(k);
        distance = get_distance(jet, jet_in);
        if (distance < matching_radius) { // is the jet closer than the previous one?
            matched_jet = jet;
            matching_radius = distance; // new distance to matched jet
        }
    }// end loop over jets
    }
    return matched_jet;
}

GenParticle* get_closest_particle(Jet* jet, std::vector<GenParticle*> genparts, double matching_radius)
{
    GenParticle *matched_particle;
    GenParticle *genpart;
    double distance;
    {
    for (size_t k = 0; k < genparts.size(); k++) { // loop over jets
        genpart = genparts[k];
        distance = get_distance(jet, genpart);
        if (distance < matching_radius) { // is the particle closer than the previous one?
            matched_particle = genpart;
            matching_radius = distance; // new distance to matched jet
        }
    }// end loop over jets
    }
    return matched_particle;
}

// for track vertex Association
Vertex* get_closest_vertex(Track* track, TClonesArray* branchVtx, double matching_radius)
{
    Vertex *matched_vertex;
    Vertex *vertex;
    double distance;
    matched_vertex->Index = -1;
    if(branchVtx->GetEntriesFast() > 0){
      for (size_t k = 0; k < branchVtx->GetEntriesFast(); k++) { // loop over vertices
          vertex = (Vertex*) branchVtx->At(k);
          distance = get_distance(vertex, track);
        if (distance < matching_radius) { // is the particle closer than the previous one?
            matched_vertex = vertex;
            matching_radius = distance; // new distance to matched vertex
        }
    }// end loop over vertices
    }
    return matched_vertex;
}

bool is_VBF(GenParticle *genparticle)
{
    if (abs(genparticle->PID) <7){
        bool from_beam_part=false;
        if((genparticle->M1 == 1 || genparticle->M1 == 2) && genparticle->M2 == 0){
            from_beam_part=true;
        }
        if((genparticle->M2 == 1 || genparticle->M2 == 2) && genparticle->M1 == 0){
            from_beam_part=true;
        }
        if (genparticle->Status==23 && from_beam_part) {
            return true;
        }
    }
    return false;
}

bool is_b(GenParticle *genparticle)
{
    if (abs(genparticle->PID) == 5 && genparticle->Status==23) {
        return true;}
    return false;
}

// id of jet constituents
// 0 - charged hadrons
// 1 - neutral hadrons
// 2 - gamma
// 3 - other
Int_t id(GenParticle *particle)
{
  switch (std::abs(particle->PID)) {
    case 211: // -+ pion
      return 0;
    case 321: // -+ kaon
      return 0;
    case 2212: // proton -+
      return 0;
    case 11: // elec
      return 3;
    case 13: // muon
      return 3;
    case 22: // gamma
      return 2;
    case 130: // K0
      return 1;
    case 111: // pion 0
      return 1;
    case 2112: // neutron
      return 1;
    default:
      return 3;
    }
}
