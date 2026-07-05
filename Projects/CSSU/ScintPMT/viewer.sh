#!/bin/bash
set -e

cd "$(dirname "$0")"

mkdir -p root/sim root/ana
npsimulation -D detector_viewer.txt -E reaction.txt -O scint_pmt_viewer.root -M geant4_vis.mac
