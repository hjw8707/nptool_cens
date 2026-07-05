#!/bin/bash
set -e

cd "$(dirname "$0")"
mkdir -p root/sim root/ana

npsimulation -D detector.txt -E reaction_alpha_am241.txt -B batch.mac -O cssu_cf4_alpha.root -N
root -l -b -q 'plot_gas_energy.C("root/sim/cssu_cf4_alpha.root","root/ana/cssu_cf4_alpha")'
