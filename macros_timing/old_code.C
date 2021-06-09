//------------match VBF genparts to jets---------------------------------
//   for (size_t l = 0; l < VBF_genparts.size(); l++) { //loop over VBF genparts
//     VBF_genpart = (GenParticle*) VBF_genparts.at(l);
//     //----match to AK4 jets
//     if(branchJet->GetEntriesFast() > 0)
//     {
//       VBF_matched_jet = get_closest_jet(branchJet, VBF_genpart, matching_radius_large);
//       plots->fVBFmatchedjetDeltaR[0]->Fill(get_distance(VBF_matched_jet, VBF_genpart));
//      if(VBF_matched_jet->PT>0){plots->fVBFmatchedjetPT[0]->Fill(VBF_matched_jet->PT);}
//       if(VBF_matched_jet->PT>20){plots->fVBFmatchedjetEta[0]->Fill(VBF_matched_jet->Eta);}
//     }
//
//     //----match to PUPPI jets
//     if(branchJetPUPPI->GetEntriesFast() > 0)
//     {
//         VBF_matched_jet = get_closest_jet(branchJetPUPPI,VBF_genpart, matching_radius_large);
//         plots->fVBFmatchedjetDeltaR[1]->Fill(get_distance(VBF_matched_jet, VBF_genpart));
//        if(VBF_matched_jet->PT>0){plots->fVBFmatchedjetPT[1]->Fill(VBF_matched_jet->PT);}
//         if(VBF_matched_jet->PT>20){plots->fVBFmatchedjetEta[1]->Fill(VBF_matched_jet->Eta);}
//     }//end if branch jet
//
//     //----match to Genjets
//     if(branchGenJet->GetEntriesFast() > 0)
//     {
//         VBF_matched_jet = get_closest_jet(branchGenJet, VBF_genpart, matching_radius_large);
//         plots->fVBFmatchedjetDeltaR[2]->Fill(get_distance(VBF_matched_jet, VBF_genpart));
//        if(VBF_matched_jet->PT>0){plots->fVBFmatchedjetPT[2]->Fill(VBF_matched_jet->PT);}
//        if(VBF_matched_jet->PT>20){plots->fVBFmatchedjetEta[2]->Fill(VBF_matched_jet->Eta);}
//     }//end if branch jet
//   } // end loop over VBF genparts
// //-----------------------------------------------------------------
// //------------match b genparts to jets---------------------------------
//   for (size_t l = 0; l < b_genparts.size(); l++) { //loop over VBF genparts
//     b_genpart = (GenParticle*) b_genparts.at(l);
//     //----match to AK4 jets
//     if(branchJet->GetEntriesFast() > 0)
//     {
//       b_matched_jet = get_closest_jet(branchJet, b_genpart, matching_radius_large);
//       plots->fbmatchedjetDeltaR[0]->Fill(get_distance(b_matched_jet, b_genpart));
//       if(b_matched_jet->PT>0){plots->fbmatchedjetPT[0]->Fill(b_matched_jet->PT);}
//       if(b_matched_jet->PT>20){plots->fbmatchedjetEta[0]->Fill(b_matched_jet->Eta);}
//     }
//
//     //----match to PUPPI jets
//     if(branchJetPUPPI->GetEntriesFast() > 0)
//     {
//         b_matched_jet = get_closest_jet(branchJetPUPPI,b_genpart, matching_radius_large);
//         plots->fbmatchedjetDeltaR[1]->Fill(get_distance(b_matched_jet, b_genpart));
//
//         if(b_matched_jet->PT>0){plots->fbmatchedjetPT[1]->Fill(b_matched_jet->PT);}
//         if(b_matched_jet->PT>20){plots->fbmatchedjetEta[1]->Fill(b_matched_jet->Eta);}
//     }//end if branch jet
//
//     //----match to Genjets
//     if(branchGenJet->GetEntriesFast() > 0)
//     {
//         b_matched_jet = get_closest_jet(branchGenJet, b_genpart, matching_radius_large);
//         plots->fbmatchedjetDeltaR[2]->Fill(get_distance(b_matched_jet, b_genpart));
//
//         if(b_matched_jet->PT>0){plots->fbmatchedjetPT[2]->Fill(b_matched_jet->PT);}
//         if(b_matched_jet->PT>20){plots->fbmatchedjetEta[2]->Fill(b_matched_jet->Eta);}
//     }//end if branch jet
//   } // end loop over b genparts
// //-----------------------------------------------------------------
// //------------match VBF genparts to jets---------------------------------
//   for (size_t l = 0; l < PU_genparts.size(); l++) { //loop over VBF genparts
//     PU_genpart = (GenParticle*) PU_genparts.at(l);
//     //----match to AK4 jets
//     if(branchJet->GetEntriesFast() > 0)
//     {
//       PU_matched_jet = get_closest_jet(branchJet, PU_genpart, matching_radius_large);
//       plots->fPUmatchedjetDeltaR[0]->Fill(get_distance(PU_matched_jet, PU_genpart));
//
//       if(PU_matched_jet->PT>0){plots->fPUmatchedjetPT[0]->Fill(PU_matched_jet->PT);}
//       if(PU_matched_jet->PT>20){plots->fPUmatchedjetEta[0]->Fill(PU_matched_jet->Eta);}
//     }
//
//     //----match to PUPPI jets
//     if(branchJetPUPPI->GetEntriesFast() > 0)
//     {
//         PU_matched_jet = get_closest_jet(branchJetPUPPI,PU_genpart, matching_radius_large);
//         plots->fPUmatchedjetDeltaR[1]->Fill(get_distance(PU_matched_jet, PU_genpart));
//
//         if(PU_matched_jet->PT>0){plots->fPUmatchedjetPT[1]->Fill(PU_matched_jet->PT);}
//         if(PU_matched_jet->PT>20){plots->fPUmatchedjetEta[1]->Fill(PU_matched_jet->Eta);}
//     }//end if branch jet
//
//     //----match to Genjets
//     if(branchGenJet->GetEntriesFast() > 0)
//     {
//         PU_matched_jet = get_closest_jet(branchGenJet, PU_genpart, matching_radius_large);
//         plots->fPUmatchedjetDeltaR[2]->Fill(get_distance(PU_matched_jet, PU_genpart));
//
//         if(PU_matched_jet->PT>0){plots->fPUmatchedjetPT[2]->Fill(PU_matched_jet->PT);}
//         if(PU_matched_jet->PT>20){plots->fPUmatchedjetEta[2]->Fill(PU_matched_jet->Eta);}
//     }//end if branch jet
//   } // end loop over PU genparts
//-----------------------------------------------------------------
