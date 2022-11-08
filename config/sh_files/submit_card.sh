#! /bin/bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_99/x86_64-centos7-gcc10-opt/setup.sh
cd /afs/cern.ch/user/a/aalbrech/Delphes
source DelphesEnv.sh

argument1=${1}

./DelphesHepMC2 cards/CMS_PhaseII/CMS_PhaseII_Snowmass_200PU_dt_${2}_${1}.tcl ROOTOUTPUT/test_resolutions/QCD/${2}/CMS_PhaseII_Snowmass_200PU_QCD_100kpileup_${2}_${3}_${1}.root MCFILES/QCD/Event_$((argument1+1)).hepmc
