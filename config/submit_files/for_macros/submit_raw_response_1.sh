#! /bin/bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_99/x86_64-centos7-gcc10-opt/setup.sh
cd /afs/cern.ch/user/a/aalbrech/Delphes
source DelphesEnv.sh

root -l -b -q macros_timing/macro_efficiency_purity.C'("ROOTOUTPUT/VBF/EtaMax3/CMS_PhaseII_Snowmass_200PU_VBF_10k_ptmin10_EtaMax3.root", "PLOTS/study_timing_cut/VBF/Raw_response/EtaMax3", "EtaMax3")'
