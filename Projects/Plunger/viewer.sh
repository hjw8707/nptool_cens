#!/bin/bash
set -e

cd "$(dirname "$0")"

mkdir -p root/sim root/ana
npsimulation -D Plunger.detector -E decay.reaction -O plunger_viewer.root -M geant4_vis.mac

