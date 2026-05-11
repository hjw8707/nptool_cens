#!/bin/bash

ANALYSIS_ONLY=0

# 옵션 파싱
while getopts "a" opt; do
  case $opt in
    a)
      ANALYSIS_ONLY=1
      ;;
  esac
done

reaction=$1 # 1: dp_gs, 2: dp_ex, 3: dt_gs, 4: d3He_gs, 5: aa_gs, 6: aa_ex
nevent=10000

cat target.detector.lh2 > detector.det
cat si.detector.50mm_downonly >> detector.det

if [ -z "$reaction" ]; then
    # argument가 없으면 GUI 모드로 실행
    cat na21.beam > reaction.reac
    npsimulation -D detector.det -E reaction.reac -M startup.mac
    exit 0
fi

if [ "$reaction" -ge 1 ] && [ "$reaction" -le 4 ]; then
    cat na21.beam > reaction.reac
elif [ "$reaction" -ge 5 ] && [ "$reaction" -le 6 ]; then
    cat mg22.beam > reaction.reac
fi


tee > batch.mac <<EOF
/run/beamOn $nevent
EOF

if [ $reaction -eq 1 ]; then
cat dp_gs.channel >> reaction.reac
output=dp_gs
elif [ $reaction -eq 2 ]; then
cat dp_ex.channel >> reaction.reac
output=dp_ex
elif [ $reaction -eq 3 ]; then
cat dt_gs.channel >> reaction.reac
output=dt_gs
elif [ $reaction -eq 4 ]; then
cat d3He_gs.channel >> reaction.reac
output=d3He_gs
elif [ $reaction -eq 5 ]; then
cat aa_gs.channel >> reaction.reac
output=aa_gs
elif [ $reaction -eq 6 ]; then
cat aa_ex.channel >> reaction.reac
output=aa_ex
fi

if [ $ANALYSIS_ONLY -eq 0 ]; then
    npsimulation -D detector.det -E reaction.reac -M startup.mac -B batch.mac -O ${output}_sim.root
fi
npanalysis -T ./root/sim/${output}_sim.root SimulatedTree -O ${output}_ana.root