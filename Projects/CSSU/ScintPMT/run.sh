#!/bin/bash
set -e

cd "$(dirname "$0")"

mkdir -p root/sim root/ana
npsimulation -D detector.txt -E reaction.txt -B batch.mac -O scint_pmt_alpha.root -N
