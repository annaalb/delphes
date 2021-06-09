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

class ExRootTreeReader;

void jet_response_plots(const char *inputFile)
{
  gSystem->Load("libDelphes");

  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile);

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);

  SetupGlobalStyle();
  //Simulation_label();

  TH1 *fmatchedresponse_AK4[6][3];
  TH1 *fmatchedresponse_PUPPI[6][3];

  THStack *hs[6][3];

  TString eta_min[6] = {"0", "1.5", "3", "4", "0", "0"};
  TString eta_max[6] = {"1.5", "3", "4", "5", "5", "5"};

  for (int l = 0; l < 3; l++) { // loop over jet collections bins
    for (int k = 0; k < 6; k++) { // loop over eta bins
      // jet response b matched jets
      TString hname_AK4 = TString::Format("%i_matched_response_AK4_bin_%i",l,k);
      TString hname_PUPPI = TString::Format("%i_matched_response_PUPPI_bin_%i",l,k);
      TString hname = TString::Format("%i_matched_response_bin_%i",l,k);

      // Book histograms
      fmatchedresponse_AK4[k][l] = new TH1F(hname_AK4, " matched jet pt / pt gen", 100, 0.0, 2.0);
      fmatchedresponse_PUPPI[k][l] = new TH1F(hname_PUPPI, " matched jet pt / pt gen", 100, 0.0, 2.0);

      fmatchedresponse_AK4[k][l]->SetLineColor(kRed);
      fmatchedresponse_AK4[k][l]->SetLineWidth(3);

      fmatchedresponse_PUPPI[k][l]->SetLineColor(kBlue);
      fmatchedresponse_PUPPI[k][l]->SetLineWidth(3);

      hs[k][l] = new THStack(hname,"");
      hs[k][l]->Add(fmatchedresponse_AK4[k][l], "sames");
      hs[k][l]->Add(fmatchedresponse_PUPPI[k][l],"sames");

    } // end loop over eta bins
  }

    TClonesArray *branchJet = treeReader->UseBranch("Jet");
    TClonesArray *branchJetPUPPI = treeReader->UseBranch("JetPUPPI");
    TClonesArray *branchJetAK8 = treeReader->UseBranch("JetAK8");
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");

    Long64_t allEntries = treeReader->GetEntries();

    cout << "** Chain contains " << allEntries << " events" << endl;

    Jet *jetAK8;

    GenParticle *genparticle;
    Jet *genjet;
    Jet *VBF_matched_jet;
    GenParticle *VBF_genpart;
    Jet *b_matched_jet;
    GenParticle *b_genpart;
    Jet *PU_matched_jet;
    GenParticle *PU_genpart;

    double matching_radius_large = 1;
    double matching_radius = 0.4;
    double gen_matching_radius = 0.2;

    Long64_t entry;

    Int_t nbins_eta = 6;
    Double_t n_eta_min[6] = {0, 1.5, 3, 4, 0, 0};
    Double_t n_eta_max[6] = {1.5, 3, 4, 5, 5, 5};

    Int_t i;
    // Loop over all events
    for(entry = 0; entry < allEntries; ++entry)
    {
      //cout << "new entry " << endl;
      std::vector<GenParticle*> VBF_genparts;
      std::vector<GenParticle*> b_genparts;
      std::vector<GenParticle*> PU_genparts;
      std::vector<GenParticle*> genparts;


      // Load selected branches with data from specified event
      treeReader->ReadEntry(entry);

      //------Analyse GenParticle---------
      if(branchParticle->GetEntriesFast() > 0)
      {
        for (size_t k = 0; k < branchParticle->GetEntriesFast(); k++) {
          genparticle = (GenParticle*) branchParticle->At(k);

          // select b quarks from higgs decay
          if (is_b(genparticle)){
              b_genparts.push_back(genparticle);
              genparts.push_back(genparticle);
               // plots->fgenparticlePT[0]->Fill(genparticle->PT);
               // plots->fgenparticleEta[0]->Fill(genparticle->Eta);
          }
          // select VBF jets
          if (is_VBF(genparticle)){
              VBF_genparts.push_back(genparticle);
              genparts.push_back(genparticle);
               // plots->fgenparticlePT[1]->Fill(genparticle->PT);
               // plots->fgenparticleEta[1]->Fill(genparticle->Eta);
          }
          // plots->fgenparticlePT[2]->Fill(genparticle->PT);
          // plots->fgenparticleEta[2]->Fill(genparticle->Eta);
          // plots->fgenparticleX[2]->Fill(genparticle->X);
          // plots->fgenparticleY[2]->Fill(genparticle->Y);
          // plots->fgenparticleZ[2]->Fill(genparticle->Z);

          //select genparticles from PU
          if (abs(genparticle->IsPU)==1) {
            PU_genparts.push_back(genparticle);
            // plots->fgenparticlePT[3]->Fill(genparticle->PT);
            // plots->fgenparticleEta[3]->Fill(genparticle->Eta);
          }
        } // end for genparticle branch
      }// end if genparticle branch not empty
      //---------Analyse GenJet-------
      if(branchGenJet->GetEntriesFast() > 0)
      {
        for (size_t k = 0; k < branchGenJet->GetEntriesFast(); k++) {
          genjet = (Jet*) branchGenJet->At(k);
           // plots->fgenjetPT->Fill(genjet->PT);
           // plots->fgenjetETA->Fill(genjet->Eta);
        }
      }

      // --------Analyse Ak4 jets --
      if(branchJet->GetEntriesFast() > 0)
      {
        Jet* jet;
        for (size_t k = 0; k < branchJet->GetEntriesFast(); k++) {
          jet = (Jet*) branchJet->At(k);
          if (jet->PT>30) {
          GenParticle* particle = get_closest_particle(jet, genparts, matching_radius);
          if (is_b(particle)) {
              //put into list of b jets
              genjet = get_closest_jet(branchGenJet, jet, gen_matching_radius);
              if(jet->PT>30 && genjet->PT>20 && genjet->PT<20000){
                for (size_t k = 0; k < nbins_eta; k++) {
                  if(n_eta_min[k] < abs(jet->Eta)&& abs(jet->Eta)<=n_eta_max[k]){fmatchedresponse_AK4[k][0]->Fill(jet->PT/genjet->PT);}
                }
              }
            }
          else if(is_VBF(particle)){
            // put into list of VBF jets
            genjet = get_closest_jet(branchGenJet, jet, gen_matching_radius);
            if(jet->PT>30 && genjet->PT>20){
              for (size_t k = 0; k < nbins_eta; k++) {
                if(n_eta_min[k] < abs(jet->Eta)&& abs(jet->Eta)<=n_eta_max[k]){fmatchedresponse_AK4[k][1]->Fill(jet->PT/genjet->PT);}
              }
            }

          }
          else{
            // rest of the jets
            genjet = get_closest_jet(branchGenJet, jet, gen_matching_radius);
            if(jet->PT>30 && genjet->PT>20){
              for (size_t k = 0; k < nbins_eta; k++) {
                if(n_eta_min[k] < abs(jet->Eta)&& abs(jet->Eta)<=n_eta_max[k]){fmatchedresponse_AK4[k][2]->Fill(jet->PT/genjet->PT);}
              }
            }
          }
        } // end if pt cut 30
        }
      }
      //-------- Analyse PUPPI jets -------------
      if(branchJetPUPPI->GetEntriesFast() > 0)
      {
        Jet* jet;
        Jet* genjet;
        for (size_t l = 0; l < branchJetPUPPI->GetEntriesFast(); l++) {
          int n=0;
          jet = (Jet*) branchJetPUPPI->At(l);
          if (jet->PT>30) {
          GenParticle* particle = get_closest_particle(jet, genparts, matching_radius);
          if (is_b(particle)) {
              //put into list of b jets
              genjet = get_closest_jet(branchGenJet, jet, gen_matching_radius);
              if(jet->PT>30 && genjet->PT>20 && genjet->PT<20000){
                for (size_t k = 0; k < nbins_eta; k++) {
                  if(n_eta_min[k] < abs(jet->Eta)&& abs(jet->Eta)<=n_eta_max[k]){
                    fmatchedresponse_PUPPI[k][0]->Fill(jet->PT/genjet->PT);
                    // if(fmatchedresponse_PUPPI[k][0]->GetBinContent(1) != 0 && n == 0){
                    //   cout << "jet: " << jet->PT << ", " << jet->Eta << endl;
                    //   cout << "genjet: " << genjet->PT << ", " << genjet->Eta << endl;
                    //   n=2;
                    // }
                  }
                }
              }
            }
          else if(is_VBF(particle)){
            // put into list of VBF jets
            // match to gen jet
            genjet = get_closest_jet(branchGenJet, jet, gen_matching_radius);
            if(jet->PT>30 && genjet->PT>20){
              for (size_t k = 0; k < nbins_eta; k++) {
                if(n_eta_min[k] < abs(jet->Eta)&& abs(jet->Eta)<=n_eta_max[k]){fmatchedresponse_PUPPI[k][1]->Fill(jet->PT/genjet->PT);}
              }
            }

            }
          else{
            // rest of the jets
            genjet = get_closest_jet(branchGenJet, jet, gen_matching_radius);
            if(jet->PT>30 && genjet->PT>20){
                for (size_t k = 0; k < nbins_eta; k++) {
                  if(n_eta_min[k] < abs(jet->Eta)&& abs(jet->Eta)<=n_eta_max[k]){fmatchedresponse_PUPPI[k][2]->Fill(jet->PT/genjet->PT);}
                }
            }
          }
        }//end if pt cut 30
        }
      }


    }// end loop over entries

    TLegend *legend;
    TCanvas* canvas = new TCanvas("c", "c", 800, 650);

    TPaveStats *stats;
    TPaveStats *stats2;

    SetupGlobalStyle();

    for (int l = 0; l < 3; l++) {
      for (int k = 0; k < 6; k++) { // loop over eta bins
        hs[k][l]->Draw("nostack");
        hs[k][l]->GetXaxis()->SetTitle("p_{T}^{reco} / p_{T}^{gen}");
        hs[k][l]->GetYaxis()->SetTitle("number of jets");
        //fmatchedresponse_AK4[k][l]->Draw();
        //first hist stats position
        canvas->Update();
        stats = (TPaveStats*)fmatchedresponse_AK4[k][l]->GetListOfFunctions()->FindObject("stats");
        stats->SetName("h1stats");
        stats->SetY1NDC(.75);
        stats->SetY2NDC(.85);
        stats->SetX1NDC(.65);
        stats->SetX2NDC(.85);
        stats->SetTextColor(2);
        stats->SetTextSize(0.03);
        stats->SetBorderSize(0);

        //fmatchedresponse_PUPPI[k][l]->Draw("SAMES");
        canvas->Update();
        stats2 = (TPaveStats*)fmatchedresponse_PUPPI[k][l]->GetListOfFunctions()->FindObject("stats");
         stats2->SetName("h1stats2");
         stats2->SetY1NDC(.65);
         stats2->SetY2NDC(.75);
         stats2->SetX1NDC(.65);
         stats2->SetX2NDC(.85);
         stats2->SetTextColor(4);
         stats2->SetTextSize(0.03);
         stats2->SetBorderSize(0);

        legend = new TLegend(0.2,0.7,0.4,0.85);
        legend->SetBorderSize(0);
        legend->SetTextSize(.04);
        legend->SetHeader(eta_min[k] + " < #eta <= "+eta_max[k]);
        legend->AddEntry(fmatchedresponse_AK4[k][l], "AK4 jet ", "l");
        legend->AddEntry(fmatchedresponse_PUPPI[k][l], "PUPPI jet ", "l");

        legend->Draw();


        gSystem->cd("Plots/jet_response");
        TString hname_b_stack = TString::Format("%i_matched_response_bin_%i.eps",l,k);
        canvas->SaveAs(hname_b_stack);
      }
  }

    cout << "** Exiting..." << endl;

    delete treeReader;
    delete chain;
  }
