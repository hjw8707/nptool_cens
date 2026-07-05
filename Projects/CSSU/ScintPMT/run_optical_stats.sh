#!/bin/bash
set -e

cd "$(dirname "$0")"

mkdir -p root/sim root/ana
NPS_OPTICAL_DEATH_STATS=1 npsimulation -D detector.txt -E reaction.txt -B batch.mac -O scint_pmt_optical_stats.root -N

