#! /bin/bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_99/x86_64-centos7-gcc10-opt/setup.sh
cd /afs/cern.ch/user/a/aalbrech/Delphes
source DelphesEnv.sh

root -l -b -q macros_timing/track_control_plots.C'("ROOTOUTPUT/QCD/EtaMax4/CMS_PhaseII_Snowmass_200PU_QCD_10k_EtaMax4.root", "PLOTS/study_timing_cut/QCD/track_control_plots_EtaMax4")'
