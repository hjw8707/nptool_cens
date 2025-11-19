#!/bin/bash
make;
./run.sh 1;
./run.sh 2;
./run.sh 3;
./run.sh 4;

cd root/ana;
hadd -f sunji_ana.root dp_gs_ana.root dp_ex_ana.root dt_gs_ana.root d3He_gs_ana.root;
cd ../..;
