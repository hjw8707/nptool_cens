#!/bin/bash

##### for simplicity, only 10000 events ####
tee > batch.mac <<EOF
/run/beamOn 10000
EOF
###########################################

echo "simulation for $x"
npsimulation -D detectors/stark_full.det -E reactions/p208Pb_el.reac -O stark_elas.root -B batch.mac
npanalysis -T ../../Outputs/Simulation/stark_elas.root SimulatedTree -O stark_elas.root

root macro/stark_elas.C 
