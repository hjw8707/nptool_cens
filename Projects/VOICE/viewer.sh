#!/bin/bash
set -e

cd "$(dirname "$0")"

mkdir -p root/sim root/ana
npsimulation -D VOICE_Kr86_vac.detector -E Kr86beamonly.reaction -O voice_viewer.root -M geant4_vis.mac

