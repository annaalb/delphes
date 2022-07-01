#! /bin/bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_99/x86_64-centos7-gcc10-opt/setup.sh
cd /afs/cern.ch/user/a/aalbrech/Delphes
source DelphesEnv.sh

root -l -b -q macros_timing/macro_QCD_jets.C'("ROOTOUTPUT/QCD_100kpileup/EtaMax3/CMS_PhaseII_Snowmass_200PU_QCD_10k_100kPU_EtaMax3.root", "PLOTS/study_timing_cut/QCD_100kPU/macro_QCD_jets_track_eta_EtaMax3")'
