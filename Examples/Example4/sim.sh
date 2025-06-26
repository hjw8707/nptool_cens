#!/bin/bash

for i in {1..10}
do
  npsimulation -D Example4.detector -E Example4.reaction --random-seed $i -O Example4_$i -B run.mac &
done
