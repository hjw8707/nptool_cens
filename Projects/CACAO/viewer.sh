#!/bin/bash
set -e

cd "$(dirname "$0")"

mkdir -p root/sim root/ana
PYTHON_BIN=${PYTHON_BIN:-python3}
"$PYTHON_BIN" geometry/cacao.py > geometry/cacao.detector
cat chamber.detector > detector.det
cat geometry/cacao.detector >> detector.det
npsimulation -D detector.det -E gamma.source -O cacao_viewer.root -M geant4_vis.mac
