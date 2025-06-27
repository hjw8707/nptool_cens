#!/bin/bash

##### for simplicity, only 1000 events ####
tee > batch.mac <<EOF
/run/beamOn 1000
EOF
###########################################

for x in proton deuton triton alpha
do
    echo "simulation for $x"
    npsimulation -D detectors/stark_full.det -E sources/$x.source -O stark_$x.root -B batch.mac
    npanalysis -T ../../Outputs/Simulation/stark_$x.root SimulatedTree -O stark_$x.root
done

cd ../../Outputs/Analysis
hadd -f stark_pdta.root stark_proton.root stark_deuton.root stark_triton.root stark_alpha.root
cd ../../Projects/STARK

root macro/stark_pdta.C 
