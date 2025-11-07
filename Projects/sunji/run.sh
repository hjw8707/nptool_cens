#!/bin/bash
cat target.detector.lh2 > detector.det
cat si.detector >> detector.det

cat na21.beam > reaction.reac
cat dp_gs.channel >> reaction.reac

npsimulation -D detector.det -E reaction.reac -M startup.mac