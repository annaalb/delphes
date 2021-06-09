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
