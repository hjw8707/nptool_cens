#!/bin/bash
set -e

cd "$(dirname "$0")"
mkdir -p root/sim root/ana
root -l -b -q make_natural_gamma_hist.C

npsimulation -D hpge.detector -E background_plus_co60_x-800.reaction -O background_plus_co60_viewer.root -M input/geant4_vis.mac
