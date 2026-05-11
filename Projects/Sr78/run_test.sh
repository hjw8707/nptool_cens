#!/bin/zsh

#npsimulation -D det.detector -E srcs/gamma.source -B batch.mac -O cacao_sim_${x}.root
npanalysis --last-sim -O cacao_ana_${x}.root
