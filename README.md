# CENS NPTool

## 프로젝트 설명

------
nptool은 저에너지 핵물리 실험을 위한 데이터 분석 및 몬테카를로 시뮬레이션 패키지입니다. 이 패키지는 복잡한 실험을 준비하고 분석하기 위한 통합 프레임워크를 제공하며, Geant4 및 ROOT 툴킷을 효율적으로 활용합니다. nptool은 사용자가 실험 데이터를 쉽게 처리하고 분석할 수 있도록 다양한 기능을 제공합니다. 이 프로젝트는 오픈 소스이며, 사용자와 개발자들이 자유롭게 기여할 수 있도록 설계되었습니다.

본 프로젝트는 CENS(센터명: Center for Exotic Nuclear Studies)에서 개발 및 유지보수하고 있는 NPTool의 공식 저장소입니다. CENS NPTool은 저에너지 핵물리 실험을 위한 데이터 분석과 몬테카를로 시뮬레이션을 지원하며, CENS의 연구 환경과 실험 요구에 최적화되어 있습니다. CENS 연구진 및 협력 연구자들이 실험 데이터를 효율적으로 처리하고 분석할 수 있도록 다양한 기능과 예제, 문서가 제공됩니다.

이 패키지는 다음과 같은 주요 기능을 포함합니다:

- **데이터 분석**: 실험 데이터를 처리하고 분석하는 데 필요한 다양한 도구를 제공합니다.
- **몬테카를로 시뮬레이션**: Geant4를 사용하여 복잡한 물리적 현상을 시뮬레이션할 수 있습니다.
- **유연한 구성**: 사용자가 필요에 따라 다양한 분석 및 시뮬레이션 설정을 조정할 수 있습니다.
- **커뮤니티 지원**: 사용자와 개발자들이 서로의 경험을 공유하고 문제를 해결할 수 있는 플랫폼을 제공합니다.

이 프로젝트는 저에너지 핵물리 분야의 연구자들에게 유용한 도구가 될 것입니다. 사용자는 이 패키지를 통해 실험 데이터를 보다 효율적으로 분석하고, 새로운 물리적 현상을 탐구할 수 있습니다.

nptool, which stands for Nuclear Physics Tool, is an open source and freely distributed data analysis and Monte Carlo simulation package for low-energy nuclear physics experiments. The NPTool package aims to offer an unified framework for preparing and analysing complex experiments, making an efficient use of Geant4 and ROOT toolkits.

- If you wish to contribute, contact Adrien Matta: <matta@lpccaen.in2p3.fr>
- For issues, use the git issue tracker or contact Adrien Matta: <matta@lpccaen.in2p3.fr>

## Useful links

- [nptool website](http://nptool.org): Find the latest information on the framework
- [nptool manual](http://nptool.org/manual): Find a more detailed manual on how to install and run the framework

## 목차

1. [코드 받기](#코드-받기)
2. [설치 및 환경설정](#설치-및-환경설정)
3. [NPLib 빌드](#nplib-빌드)
4. [NPSimulation 빌드](#npsimulation-빌드)
5. [벤치마크 및 예제](#벤치마크-및-예제)
6. [트릭 및 팁](#트릭-및-팁)
<!-- 7. [npreader](#npreader)
8. [npcalibration](#npcalibration) -->

## 코드 받기

### Git으로 받기

가장 추천하는 방법은 git을 사용하는 것입니다. 최신 버전의 코드를 쉽게 받을 수 있습니다. 먼저 git이 설치되어 있는지 확인하세요. 없다면 패키지 매니저를 통해 설치하세요. 그 다음, NPTool 패키지를 설치할 디렉토리로 이동하여 아래 명령어를 실행하세요.

```sh
git clone https://gitlab.in2p3.fr/np/nptool
```

이 명령어를 실행하면 nptool 폴더에 최신 버전의 nptool이 다운로드됩니다.

### GitHub에서 받기

또는 <https://gitlab.in2p3.fr/np/nptool> 페이지에서 Download ZIP 버튼을 클릭하여 압축 파일을 받을 수 있습니다. 받은 파일을 원하는 위치에 압축 해제하세요.

## 설치 및 환경설정

### 요구사항

NPTool 컴포넌트는 CMake 빌드 시스템을 사용하여 컴파일 및 설치됩니다. 시작하기 전에 CMake가 설치되어 있는지 확인하세요.

NPLib(코어 라이브러리)를 컴파일하려면 ROOT 5(5.34 테스트됨) 또는 6이 필요합니다. NPLib과 분석 프로젝트를 컴파일하는 데 충분합니다.

NPSimulation을 컴파일하려면 최신 Geant4(9.6, 10.1 테스트됨)가 필요합니다. GDML 지원이 필요하다면 Geant4를 GDML 옵션과 함께 설치해야 합니다.

### 환경 변수 설정

환경 변수(PATH, LD_LIBRARY_PATH 등)와 alias를 설정하려면 다음 스크립트를 실행하세요.

```sh
source <설치경로>/nptool/nptool.sh
```

`<설치경로>`는 NPTool 패키지를 압축 해제한 위치입니다. 이후 터미널을 재시작하세요.

이 명령어를 .profile, .bashrc, .zshrc 등에 추가하면 매번 입력할 필요가 없습니다.

## NPLib 빌드

NPLib은 NPTool 패키지의 핵심으로, 대부분의 실제 코드가 포함되어 있습니다. 독립적인 C++ 클래스 모음으로 구성되어 있으며, 프로그램과 매크로에서 사용할 수 있습니다.

1. NPLib 폴더로 이동:

    ```sh
    npl
    ```

2. CMake로 Makefile 생성:

    - 모든 검출기 컴파일:

    ```sh
    cmake ./
    ```

    - 일부 검출기만 컴파일:

    ```sh
    cmake ./ -DNPTOOL_DETLIST="DetFolder1 DetFolder2"
    ```

3. 컴파일 및 설치:

    ```sh
    make -jn install
    ```

    (j는 스레드 개수)

4. 검출기 추가 컴파일 시:

    ```sh
    nptool-cleaner
    cmake ./ -DNPTOOL_DETLIST="DetFolder1 DetFolder2 ..."
    make -jn install
    ```

5. Ninja 빌드 사용 시:

    ```sh
    cmake -GNinja ./
    ninja install
    ```

Ninja가 make보다 빠릅니다.

## NPSimulation 빌드

이 부분은 Geant4를 이용한 몬테카를로 시뮬레이션을 담당합니다. NPLib이 먼저 컴파일되어 있어야 하며, 그 후 NPSimulation을 컴파일할 수 있습니다.

1. NPSimulation 폴더로 이동:

    ```sh
    nps
    ```

2. CMake로 Makefile 생성:

    ```sh
    cmake ./
    ```

3. 컴파일 및 설치:

    ```sh
    make -jn install
    ```

실행 파일: `npsimulation`

사용 가능한 입력 플래그와 설명은 다음 명령어로 확인할 수 있습니다:

```sh
npsimulation -h
```

## 벤치마크 및 예제

벤치마크와 예제를 실행하려면 추가 파일이 필요합니다. $NPTOOL 디렉토리에서 다음 명령어를 실행하세요.

```sh
git clone https://github.com/adrien-matta/NPData
```

### 벤치마크

벤치마크는 설치 또는 업그레이드의 무결성 확인, CPU 성능 비교 등에 유용합니다. 두 가지 주요 벤치마크가 제공됩니다.

1. cats (빔 트래커 데이터 분석)
2. gaspard (실리콘 어레이 시뮬레이션)

각 벤치마크는 결과를 그림으로 출력하며, 참조 결과와 비교할 수 있습니다.

#### cats 벤치마크 실행

```sh
cd $NPTOOL/Benchmarks/cats
npanalysis -D benchmark_cats.detector -C calibration.txt -R RunToTreat.txt -O benchmark_cats
```

#### gaspard 벤치마크 실행

```sh
cd $NPTOOL/Benchmarks/gaspard
npsimulation -D benchmark_gaspard.detector -E 132Sndp_benchmark.reaction -O benchmark_gaspard -B batch.mac
```

결과 확인:

```sh
root -l ShowResult.C
```

### 예제

예제는 여러 검출기를 조합한 복잡한 분석 사례를 다룹니다. Example1 실행 예시는 다음과 같습니다.

```sh
npsimulation -D Example1.detector -E Example1.reaction -O Example1
```

이후 GUI 또는 프롬프트에서 이벤트를 생성:

```text
> run/beamOn/ 10000
> exit
```

분석:

```sh
npp Example1
cmake ./
make -jn
npanalysis -R RunToTreat.txt -O Example1
```

결과 확인:

```sh
root -l ShowResult.C
```

## 트릭 및 팁

- _npsimulation_과 _npanalysis_는 어느 디렉토리에서나 실행할 수 있습니다.
- _npanalysis_는 현재 디렉토리에서 분석 라이브러리(_libNPAnalysis_)를 자동으로 로드합니다. 없으면 PhysicsTree만 생성합니다.
- 최근 시뮬레이션을 빠르게 분석하려면:

```sh
npanalysis --last-sim
```

- _npsimulation_은 배치 모드(-B 플래그, UI 없이)로도 실행할 수 있습니다.

```sh
npsimulation -D Example1.detector -E Example1.reaction -B path/to/macro.mac -O FileName
```
