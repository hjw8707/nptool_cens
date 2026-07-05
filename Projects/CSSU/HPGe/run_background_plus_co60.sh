#!/bin/bash
set -e

cd "$(dirname "$0")"

mkdir -p root/sim root/ana
root -l -b -q make_natural_gamma_hist.C
npsimulation -D hpge.detector -E background_plus_co60_x-800.reaction -B batch.mac -O background_plus_co60_x-800.root -N
root -l -b -q 'plot_hpge_energy.C("root/sim/background_plus_co60_x-800.root","root/ana/background_plus_co60_x-800")'
