#!/bin/bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_99/x86_64-centos7-gcc10-opt/setup.sh
cd /afs/cern.ch/user/a/aalbrech/Delphes
source DelphesEnv.sh

echo $1
echo $2

./DelphesHepMC2 cards/CMS_PhaseII/CMS_PhaseII_Snowmass_200PU_dt_${1}.tcl ROOTOUTPUT/QCD_100kpileup/EtaMax${2}/CMS_PhaseII_Snowmass_200PU_QCD_100kpileup_EtaMax${2}_resolution_70ps_dummy_${1}.root MCFILES/QCD_HT100toInf_pythia8_events.hepmc
