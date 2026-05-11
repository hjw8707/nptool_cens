#!/bin/zsh

ens=(0289 0471 0701 0811 2251)
#cat geometry/target.detector geometry/TiNA.detector > det.detector
#cat geometry/target.detector geometry/TiNA.detector geometry/GRAPE.detector > det.detector
cat geometry/target.detector geometry/TiNA.detector geometry/DALI2_top_hole.detector > det.detector
#cat geometry/target.detector geometry/TiNA.detector geometry/DALI2_half.detector geometry/GRAPE_half.detector > det.detector
#npsimulation -D det.detector -E reac.reaction

for x in $ens
do
    sed -i 's/EnergyLow= .* keV/EnergyLow= ${x} keV/g' srcs/gamma.source
    sed -i 's/EnergyHigh= .* keV/EnergyHigh= ${x} keV/g' srcs/gamma.source
    npsimulation -D det.detector -E srcs/gamma.source -B batch.mac -O dali2_sim_${x}.root
    npanalysis --last-sim -O dali2_ana_${x}.root
done
