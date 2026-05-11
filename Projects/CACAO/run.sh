#!/bin/bash

python geometry/cacao.py > geometry/cacao.detector
cat chamber.detector > detector.det
cat geometry/cacao.detector >> detector.det
#cat cacao_single.detector >> detector.det

npsimulation -D detector.det -E gamma.source
#npsimulation -D detector.det -E gamma.source -B batch.mac -O cacao_single_sim.root
#npanalysis --last-sim -O cacao_single_ana.root