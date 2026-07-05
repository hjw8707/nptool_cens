#!/bin/bash
set -e

cd "$(dirname "$0")"
mkdir -p root/sim root/ana

npsimulation -D detector.txt -E reaction_alpha_am241.txt -O cssu_cf4_alpha_viewer.root -M geant4_vis.mac
