#! /bin/bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_99/x86_64-centos7-gcc10-opt/setup.sh
cd /afs/cern.ch/user/a/aalbrech/Delphes
source DelphesEnv.sh

echo "Hello World Nr 1"

./DelphesHepMC2 cards/CMS_PhaseII/CMS_PhaseII_Snowmass_200PU.tcl ROOTOUTPUT/CMS_PhaseII_Snowmass_200PU_test.root MCFILES/VBF_HH_bbbb_10k.hepmc
