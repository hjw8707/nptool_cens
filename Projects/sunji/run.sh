#!/bin/bash
reaction=$1 # 1: dp_gs, 2: dp_ex, 3: dt_gs, 4: d3He_gs
nevent=10

if [ -z "$reaction" ]; then
    echo "Usage: $0 <reaction>"
    echo "  reaction: 1: dp_gs, 2: dp_ex, 3: dt_gs, 4: d3He_gs"
    exit 1
fi

tee > batch.mac <<EOF
/run/beamOn $nevent
EOF

cat target.detector.lh2 > detector.det
cat si.detector >> detector.det

cat na21.beam > reaction.reac

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
fi

npsimulation -D detector.det -E reaction.reac -M startup.mac -B batch.mac -O ${output}_sim.root