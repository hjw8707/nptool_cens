#!/bin/bash
make;
./run.sh 5;
./run.sh 6;

cd root/ana;
hadd -f sunji_aa.root aa_gs_ana.root aa_ex_ana.root;
cd ../..;
