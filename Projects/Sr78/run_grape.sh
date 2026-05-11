#!/bin/zsh

y=0811
#cat geometry/target.detector geometry/TiNA.detector > det.detector
#cat geometry/target.detector geometry/TiNA.detector geometry/GRAPE.detector > det.detector
cat geometry/target.detector geometry/TiNA.detector geometry/DALI2.detector > det.detector
#cat geometry/target.detector geometry/TiNA.detector geometry/DALI2_half.detector geometry/GRAPE_half.detector > det.detector
#npsimulation -D det.detector -E reac.reaction

npsimulation -D det.detector -E gamma.source -B batch.mac -O dali2_sim_${y}.root
npanalysis --last-sim -O dali2_ana_${y}.root
