#!/bin/bash
set -e

cd "$(dirname "$0")"

mkdir -p root/sim root/ana
root -l -b -q make_natural_gamma_hist.C
npsimulation -D hpge.detector -E natural_background.reaction -B batch.mac -O hpge_natural_background.root -N
