#!/bin/zsh

model=bmf # bmf or mcsm
target=lh2
tina=ttt30cm
gamma=cacao # tina, grape, dali, csi, ring, tilt, cacao
#
bmode=single # no, single, array
bevent=100000
barray=(18467 5657 371 934 467 82) # 50 MeV
#barray=(184670 56570  3710 9340 4670 820) # 50 MeV
#
#sloop="1st0 1st2 1st4 2nd0 2nd2 2nd4"
sloop="1st0 1st2 2nd0 2nd2"
eloop="50mev" #"45mev 50mev 55mev"
#
#loop="1st0 2nd0 1st2 2nd2"

############################################################
# interactive mode with a second argument
if [[ $# -gt 0 ]]
then
    echo "The second argument is detected. Run in the interactive mode."
    bmode=no
fi
############################################################

##############################
# filename
if [[ $bmode == "no" ]]
then
    fname=${gamma}_${tina}_${target}_${model}
else
    if [[ $bmode == "single" ]]
    then
	fname=${gamma}_${tina}_${target}_${model}
    else
	fname=${gamma}_${tina}_${target}_${model}_batch
    fi
fi
##############################

##########################################
# reaction model
cd reac
ln -sf 2nd0.channel.${model} 2nd0.channel
ln -sf 2nd2.channel.${model} 2nd2.channel
ln -sf 2nd4.channel.${model} 2nd4.channel
cd ..
##########################################

##########################################
# event number
if [[ $bmode == "no" ]]
then
    echo "interactive mode"
else # only for single
    echo "/tracking/verbose 0" > batch.mac
    echo "/run/beamOn ${bevent}" >> batch.mac
fi
##########################################


cat geometry/target.detector.${target} > det.detector # target selection
cat geometry/chamber.detector >> det.detector # add chamber

############################################################
# TiNA geometry
if [[ $tina == "ttt20cm" ]]
then
    cat geometry/TiNA.detector >> det.detector
elif [[ $tina == "ttt30cm" ]]
then
    cat geometry/TiNA.detector.30cm >> det.detector
elif [[ $tina == "tttyy1" ]]
then
    cat geometry/TiNA.detector.rev >> det.detector
fi
############################################################

############################################################
# gamma detector
if [[ $gamma == "grape" ]]
then
    cat geometry/GRAPE.detector >> det.detector # grape
elif [[ $gamma == "dali" ]]
then
    cat geometry/DALI2_sepa.detector >> det.detector # dali2
elif [[ $gamma == "csi" ]]
then
    cat geometry/csi.detector >> det.detector # csi
elif [[ $gamma == "ring" ]]
then
    cat geometry/csi.detector.ring >> det.detector # csi
elif [[ $gamma == "tilt" ]]
then
    cat geometry/csi.detector.tilt >> det.detector # csi
elif [[ $gamma == "cacao" ]]
then
    cat geometry/cacao.detector >> det.detector # CACAO
fi
############################################################

#########################################
# Test the Nan's detector
#cp geometry/nan.detector ./det.detector
#########################################

################################################################################
# interactive mode
if [[ $bmode == "no" ]]
then
    x="1st0"
    y="50mev"
    cat reac/80Sr.beam.${y} reac/${x}.channel > reac.reaction
    npsimulation -D det.detector -E reac.reaction
    exit
fi
################################################################################

i=1
for x in `echo $sloop`
do
    for y in `echo $eloop`
    do
	##############################
	# event number for each loop
	if [[ $bmode == "array" ]]
	then
	    echo "/tracking/verbose 0" > batch.mac
	    echo "/random/setSeeds 293949294 1929394929" >> batch.mac
	    echo "/run/beamOn ${barray[${i}]}" >> batch.mac
	fi
	##############################
	##############################
	# energy loop
	cat reac/80Sr.beam.${y} reac/${x}.channel > reac.reaction
	##############################
	npsimulation -D det.detector -E reac.reaction -B batch.mac -O Sr78_sim_${x}_${fname}_${y}.root
	npanalysis -T root/sim/Sr78_sim_${x}_${fname}_${y}.root SimulatedTree -O Sr78_ana_${x}_${fname}_${y}.root
    done
    i=$( expr $i + 1 )
done

########################################
# for batch -> making a merged file
if [[ $bmode == "array" ]]
then
    cd root/ana
    hadd -f Sr78_ana_${fname}_${eloop}.root Sr78_ana_????_${fname}_${eloop}.root
    cd ..
fi
########################################
