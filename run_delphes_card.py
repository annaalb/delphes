import os
import sys
## define values: count = in how many parts should the sample be split
# the script loops over the different etamax possibilities and different resolutions. Copies the delphes cards and creates submit files for each configuration.
# run with python run_delphes_card.py

scenario = ["resolution_30ps", "resolution_50ps", "resolution_70ps"]
scenario_resolutions = ["30E-12", "50E-12", "70E-12"]

count_VBF = 2 # for full run this should be 4 * 2500 = 10 000 events
count = 2 # for full run this should be 400 * 2500 = 1 mio events

etamax_name = ["EtaMax0", "EtaMax3", "EtaMax4"]
card_prefix = "cards/CMS_PhaseII/CMS_PhaseII_Snowmass_200PU_dt_"
card_0 = card_prefix+scenario[0]+"_"+etamax_name[0]+"_0.tcl"

## copy card_0 with two different eta max configurations, the rest should stay the same
for e in range(2):
    print(e)
    card_eta_0 = card_prefix+scenario[0]+"_"+etamax_name[e+1]+"_0.tcl"
    with open(card_0, 'r') as firstfile, open(card_eta_0, 'w') as secondfile:
        for line in firstfile:
            if "set EtaMax { 0}" in line:
                secondfile.write("set EtaMax { "+ str(e+3)+"} \n")
            # write content to second file
            else:
                secondfile.write(line)
    print("created card: "+card_eta_0)
# copy card for each resolution (two copies "r" for all 3 eta max "e")
for e in range(3):
    for r in range(2):
        print(r)
        card_eta_0 = card_prefix+scenario[0]+"_"+etamax_name[e]+"_0.tcl"
        card_res_0 = card_prefix+scenario[r+1]+"_"+etamax_name[e]+"_0.tcl"
        with open(card_eta_0, 'r') as firstfile, open(card_res_0, 'w') as secondfile:
            for line in firstfile:
                if "set TimeResolutionForward " in line:
                    secondfile.write("set TimeResolutionForward "+ scenario_resolutions[r+1]+"} \n")
                    # write content to second file
                else:
                    secondfile.write(line)
        print("created card: "+card_res_0)

## then for each of these cards loop from 1 to count and copy the cards. (code below can stay the same)
for e in range(3):
    for r in range(3):
        ## loop from 1 to 'count' here
        for i in range(count):
            ##copy the delphes card with different seed and Nskip events
            if i!=0:
                card_i = card_prefix+scenario[r]+"_"+etamax_name[e]+"_"+str(i)+".tcl"
                # open both files
                with open(card_prefix+scenario[r]+"_"+etamax_name[e]+"_0.tcl",'r') as firstfile, open(card_i,'w') as secondfile:
                    # read content from first file
                    for line in firstfile:
                        if "set SkipEvents 0" in line:
                            secondfile.write("set SkipEvents "+str(2500*i)+" \n")
                        elif "set RandomSeed 1" in line:
                            secondfile.write("set RandomSeed "+str(i*10) +"\n")
                        # write content to second file
                        else:
                            secondfile.write(line)
        # submit file for VBF
        with open("config/sub_files/"+scenario[r]+"_"+etamax_name[e]+"_VBF.submit",'w+') as htc_config:
            htc_config.write("""
            notify_user       = anna.albrecht@desy.de
            executable = /afs/cern.ch/user/a/aalbrech/Delphes/config/sh_files/submit_card_VBF.sh
            output = /afs/cern.ch/user/a/aalbrech/Delphes/config/log/QCD/submit_card_$(ClusterId)_$(ProcId).out
            error = /afs/cern.ch/user/a/aalbrech/Delphes/config/log/QCD/submit_card_$(ClusterId)_$(ProcId).err
            log = /afs/cern.ch/user/a/aalbrech/Delphes/config/log/QCD/submit_card_$(ClusterId).log
            +JobFlavour = "microcentury"
            arguments         = $(ProcId) """+etamax_name[e]+""" """+scenario[r]+ """
            queue """+str(count_VBF)+"""
            """)
        # submit file for QCD
        with open("config/sub_files/"+scenario[r]+"_"+etamax_name[e]+".submit",'w+') as htc_config:
            htc_config.write("""
            notify_user       = anna.albrecht@desy.de
            executable = /afs/cern.ch/user/a/aalbrech/Delphes/config/sh_files/submit_card.sh
            output = /afs/cern.ch/user/a/aalbrech/Delphes/config/log/QCD/submit_card_$(ClusterId)_$(ProcId).out
            error = /afs/cern.ch/user/a/aalbrech/Delphes/config/log/QCD/submit_card_$(ClusterId)_$(ProcId).err
            log = /afs/cern.ch/user/a/aalbrech/Delphes/config/log/QCD/submit_card_$(ClusterId).log
            +JobFlavour = "microcentury"
            arguments         = $(ProcId) """+etamax_name[e]+""" """+scenario[r]+ """
            queue """+str(count)+"""
            """)

print("now run these commands:")
## submit everything to condor
for e in range(3):
    for r in range(3):
        string1="condor_submit config/sub_files/"+scenario[r]+"_"+etamax_name[e]+"_VBF.submit &"
        string="condor_submit config/sub_files/"+scenario[r]+"_"+etamax_name[e]+".submit &"
        print(string1)
        print(string)
