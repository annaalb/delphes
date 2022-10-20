import os
## define values: etamax, resolution, QCD or VBF, count = in how many parts should the sample be split
name = "dummy_script"
etamax = "4"
scenario = "resolution_70ps"
is_QCD = True
count = 10
if is_QCD:
    sample_name = "QCD_100kpileup"
    input = "MCFILES/QCD_HT100toInf_pythia8_events.hepmc"
else:
    sample_name = "VBF"
    input = "MCFILES/VBF_HH_bbbb_10k.hepmc"
##copy the delphes card with different seed and Nskip events TODO

## loop from 1 to 10 here to create sh files
for i in range(count):
    with open("config/sh_files/"+name+"_"+str(i)+".sh",'w+') as wrapper_script:
        wrapper_script.write("""
        #! /bin/bash
        source /cvmfs/sft.cern.ch/lcg/views/LCG_99/x86_64-centos7-gcc10-opt/setup.sh
        cd /afs/cern.ch/user/a/aalbrech/Delphes
        source DelphesEnv.sh

        ./DelphesHepMC2 cards/CMS_PhaseII/CMS_PhaseII_Snowmass_200PU_dt_"""+str(i)+""".tcl ROOTOUTPUT/"""+sample_name+"""/EtaMax"""+etamax+"""/CMS_PhaseII_Snowmass_200PU_"""+sample_name+"""_EtaMax"""+etamax+"""_"""+scenario+"""_"""+str(i)+""".root """+input+"""
        """)
with open("config/sub_files/"+name+".submit",'w+') as htc_config:
    htc_config.write("""
    universe          = vanilla
    notification      = Error
    notify_user       = anna.albrecht@desy.de
    executable = /afs/cern.ch/user/a/aalbrech/Delphes/config/sh_files/"""+name+"""_$(ProcId).sh
    output = /afs/cern.ch/user/a/aalbrech/Delphes/config/log/QCD/submit_card_$(ClusterId)_$(ProcId).out
    error = /afs/cern.ch/user/a/aalbrech/Delphes/config/log/QCD/submit_card_$(ClusterId)_$(ProcId).err
    log = /afs/cern.ch/user/a/aalbrech/Delphes/config/log/QCD/submit_card_$(ClusterId).log
    +JobFlavour = "microcentury"
    queue """+str(count)+"""
    """)
## submit everything to condor
string="condor_submit "+name+".submit"
print(string)
