
        #! /bin/bash
        source /cvmfs/sft.cern.ch/lcg/views/LCG_99/x86_64-centos7-gcc10-opt/setup.sh
        cd /afs/cern.ch/user/a/aalbrech/Delphes
        source DelphesEnv.sh

        ./DelphesHepMC2 cards/CMS_PhaseII/CMS_PhaseII_Snowmass_200PU_dt_4.tcl ROOTOUTPUT/QCD_100kpileup/EtaMax4/CMS_PhaseII_Snowmass_200PU_QCD_100kpileup_EtaMax4_resolution_70ps_4.root MCFILES/QCD_HT100toInf_pythia8_events.hepmc
        