#!/bin/zsh

nevent=100000
ens=(471 2251)
dists=(115 135 155 175 195)
#cat geometry/target.detector geometry/TiNA.detector > det.detector
#cat geometry/chamber.detector geometry/target.detector.lh2 geometry/TiNA.detector.30cm geometry/cacao.detector > det.detector

if [[ $# -gt 0 ]]; then
    # 인자가 하나라도 있으면, npsimulation을 non-batch 모드로 실행하고 npanalysis는 실행하지 않음
    #for dist in $dists; do
        ./geometry/cacao/cacao.py > geometry/cacao.detector # cacao.py for the final configuration as a module
        cat geometry/chamber.detector geometry/target.detector.lh2 geometry/cacao.detector >det.detector
        sed -i "" "s|/run/beamOn .*|/run/beamOn ${nevent}|g" batch.mac
        sed -i "" "s/EnergyLow= .* keV/EnergyLow= ${x} keV/g" srcs/gamma.source
        sed -i "" "s/EnergyHigh= .* keV/EnergyHigh= ${x} keV/g" srcs/gamma.source
        npsimulation -D det.detector -E srcs/gamma.source -O cacao_sim_${x}_${dist}.root
    #done
else
    # 인자가 없으면 기존대로 batch 모드로 npsimulation과 npanalysis 실행
    #for dist in $dists; do
        ./geometry/cacao/cacao.py >geometry/cacao.detector
        cat geometry/chamber.detector geometry/target.detector.lh2 geometry/cacao.detector >det.detector

        for x in $ens; do
            # batch.mac 파일의 /run/beamOn 파라미터를 nevent 변수 값으로 변경
            sed -i "" "s|/run/beamOn .*|/run/beamOn ${nevent}|g" batch.mac
            sed -i "" "s/EnergyLow= .* keV/EnergyLow= ${x} keV/g" srcs/gamma.source
            sed -i "" "s/EnergyHigh= .* keV/EnergyHigh= ${x} keV/g" srcs/gamma.source
            npsimulation -D det.detector -E srcs/gamma.source -B batch.mac -O cacao_sim_${x}.root
            #npanalysis --last-sim -O cacao_ana_${x}.root
        done
    #done
fi
