#! /bin/bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_99/x86_64-centos7-gcc10-opt/setup.sh
cd /afs/cern.ch/user/a/aalbrech/Delphes
source DelphesEnv.sh

./DelphesHepMC2 cards/CMS_PhaseII/CMS_PhaseII_Snowmass_200PU_dt_6.tcl ROOTOUTPUT/VBF/EtaMax4/CMS_PhaseII_Snowmass_200PU_VBF_1k_EtaMax4_6.root MCFILES/VBF_HH_bbbb_10k.hepmc
