#!/bin/zsh
##### for simplicity, only 1000 events ####
#tee >batch.mac <<EOF
#/run/beamOn 10
#EOF
###########################################

#for x in alpha; do #proton deuton triton alpha
#    echo "simulation for $x"
export >export.txt
x=alpha
npsimulation -D detectors/stark_full.det -E sources/$x.source -O stark_$x.root -B batch.mac
npanalysis --last-sim -O stark_$x.root
#done

#cd ../../Outputs/Analysis
#hadd -f stark_pdta.root stark_proton.root stark_deuton.root stark_triton.root #stark_alpha.root
#cd ../../Projects/STARK

#root macro/stark_pdta.C
