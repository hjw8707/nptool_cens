#!/bin/zsh
cat geometry/target.detector.ch2 geometry/TiNA.detector > det.detector
#npsimulation -D det.detector -E srcs/particle.source -B batch.mac -O Zr80_particle_tina.root
npanalysis --last-sim -O Zr80_particle_tina.root
