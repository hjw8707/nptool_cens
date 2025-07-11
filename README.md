# CENS NPTool

## Administrators

- Jongwon Hwang (CENS, <jwhwang@ibs.re.kr>)

nptool은 저에너지 핵물리 실험을 위한 데이터 분석 및 몬테카를로 시뮬레이션 패키지입니다. 이 패키지는 복잡한 실험을 준비하고 분석하기 위한 통합 프레임워크를 제공하며, Geant4 및 ROOT 툴킷을 효율적으로 활용합니다. nptool은 사용자가 실험 데이터를 쉽게 처리하고 분석할 수 있도록 다양한 기능을 제공합니다. 이 프로젝트는 오픈 소스이며, 사용자와 개발자들이 자유롭게 기여할 수 있도록 설계되었습니다.

nptool is a data analysis and Monte Carlo simulation package for low-energy nuclear physics experiments. This package provides an integrated framework for preparing and analyzing complex experiments, efficiently utilizing the Geant4 and ROOT toolkits. nptool offers various features to help users easily process and analyze experimental data. This project is open source and designed to allow users and developers to freely contribute.

본 프로젝트는 CENS(센터명: Center for Exotic Nuclear Studies)에서 개발 및 유지보수하고 있는 NPTool의 공식 저장소입니다. CENS NPTool은 저에너지 핵물리 실험을 위한 데이터 분석과 몬테카를로 시뮬레이션을 지원하며, CENS의 연구 환경과 실험 요구에 최적화되어 있습니다. CENS 연구진 및 협력 연구자들이 실험 데이터를 효율적으로 처리하고 분석할 수 있도록 다양한 기능과 예제, 문서가 제공됩니다. 특히, CENS에서 개발된 검출기에 대한 시뮬레이션 코드들이 추가되어 있습니다. CENS의 모든 유저들이 직접 검출기를 추가하거나 코드를 수정할 수 있습니다.

This project is the official repository of NPTool, developed and maintained by CENS (Center for Exotic Nuclear Studies). CENS NPTool supports data analysis and Monte Carlo simulation for low-energy nuclear physics experiments, optimized for the research environment and experimental requirements of CENS. Various features, examples, and documentation are provided to help CENS researchers and collaborators efficiently process and analyze experimental data. In particular, simulation codes for detectors developed at CENS are included. All CENS users can directly add detectors or modify the code.

현재 이용할 수 있는 CENS 검출기의 목록은 다음과 같습니다.

The following CENS detectors are currently available:

- ASGARD (Clover HPGe Array)
- Khala (LaBr3 Array, the part of IDATEN)
- Fatima (LaBr3 Array, the part of IDATEN)
- STARK (Si Array)
- STARKjr (Si Array)
- VOICE (Active Target TPC)
- CACAO (CsI:Tl Array)
- Plunger (Plunger for lifetime measurement)

## 검증된 실행 환경 (Verified Environment)

- **OS**: macOS 15.5
- **Compiler**: Clang 17.0.0
- **CMake**: 3.31.6
- **ROOT**: 6.35.01
- **Geant4**: 11.3.1

## 유용한 링크 (Useful Links)

- [nptool 웹사이트](http://nptool.org): 프레임워크에 대한 최신 정보를 확인할 수 있습니다. (Official website of nptool)
- [nptool 매뉴얼](http://nptool.org/manual): 설치 및 실행 방법에 대한 자세한 매뉴얼을 확인할 수 있습니다. (Detailed manual for installation and execution)

---

## 목차 (Table of Contents)

- [설치 방법 (Installation Guide)](#설치-방법-installation-guide)
  - [의존성 설치 (Install Dependencies)](#1-의존성-설치-install-dependencies)
  - [소스 코드 받기 (Get-the-Source-Code)](#2-소스-코드-받기-get-the-source-code)
  - [설치 및 환경설정 (Installation-and-Setup)](#3-설치-및-환경설정-installation-and-setup)
  - [NPLib 빌드 (Building-NPLib)](#4-nplib-빌드-building-nplib)
  - [NPSimulation 빌드 (Building-NPSimulation)](#5-npsimulation-빌드-building-npsimulation)
  - [벤치마크 (Benchmarks)](#6-벤치마크-benchmarks)
- [이용 방법 (Usage Guide)](#이용-방법-usage-guide)
  - [프로그램 컨셉 (Program Concept)](#1-프로그램-컨셉-program-concept)
  - [각 프로그램의 사용법 (Usage of each program)](#2-각-프로그램의-사용법-usage-of-each-program)
  - [프로젝트 디렉토리 (Project Directory)](#3-프로젝트-디렉토리-project-directory)
  - [예제 (Examples)](#4-예제-examples)
  - [검출기별 메뉴얼 (Detector-specific manuals)](#5-검출기별-메뉴얼-detector-specific-manuals)

## 설치 방법 (Installation Guide)

### 1. 의존성 설치 (Install Dependencies)

nptool을 설치하기 전에 다음 필수 소프트웨어가 시스템에 설치되어 있어야 합니다.

Before installing nptool, make sure the following required software is installed on your system.

- **CMake** (최소 3.10 이상 권장 / recommended: 3.10 or higher)
- **ROOT** (6.x 버전 권장 / recommended: version 6.x)
- **Geant4** (11.x 버전 권장 / recommended: version 11.x)
- **C++ 컴파일러** (예: gcc, clang / e.g., gcc, clang)

각 소프트웨어의 설치 방법은 공식 홈페이지 또는 패키지 매니저(예: Homebrew, apt, yum 등)를 참고하세요.

For installation instructions, please refer to the official website of each software or use a package manager (e.g., Homebrew, apt, yum, etc.).

### 2. 소스 코드 받기 (Get the Source Code)

아래 명령어를 사용하여 nptool 소스 코드를 클론하세요.

Clone the nptool source code using the command below.

```sh
git clone https://github.com/hjw8707/nptool_cens
```

이 명령어를 실행하면 nptool 폴더에 최신 버전의 nptool이 다운로드됩니다.

By running this command, the latest version of nptool will be downloaded into the nptool folder.

### 3. 설치 및 환경설정 (Installation and Setup)

#### 환경 변수 설정 (Setting Environment Variables)

환경 변수(PATH, LD_LIBRARY_PATH 등)와 alias를 설정하려면 다음 스크립트를 실행하세요.

To set environment variables (such as PATH, LD_LIBRARY_PATH) and aliases, run the following script:

```sh
source <설치경로/Installation Path>/nptool/nptool.sh
```

`<설치경로>`는 NPTool 패키지를 압축 해제한 위치입니다. 이후 터미널을 재시작하세요.

`<Installation Path>` is the directory where you extracted the NPTool package. After running the above command, please restart your terminal.

이 명령어를 .profile, .bashrc, .zshrc 등에 추가하면 매번 입력할 필요가 없습니다.

If you add this command to your `.profile`, `.bashrc`, or `.zshrc` file, you won't need to enter it every time.

### 4. NPLib 빌드 (Building NPLib)

NPLib은 NPTool 패키지의 핵심 라이브러리로, 대부분의 실제 코드가 포함되어 있습니다. NPLib은 독립적인 C++ 클래스 모음으로 구성되어 있으며, 프로그램과 매크로에서 사용할 수 있습니다.

NPLib is the core library of the NPTool package, containing most of the actual code. It consists of a collection of independent C++ classes and can be used in programs and macros.

#### 1. NPLib 폴더로 이동 (Move to the NPLib folder)

```sh
npl
```

#### 2. CMake로 Makefile 생성 (Create Makefile with CMake)

```sh
cmake ./
```

#### 3. 컴파일 및 설치 (Compile and Install)

```sh
make -jn install
```

(j는 스레드 개수 (j is the number of threads))

#### 4. 검출기 추가 컴파일 시 (Compile additional detectors)

```sh
nptool-cleaner
cmake ./ -DNPTOOL_DETLIST="DetFolder1 DetFolder2 ..."
make -jn install
```

#### 5. Ninja 빌드 사용 시 (Use Ninja build)

```sh
cmake -GNinja ./
ninja install
```

Ninja가 make보다 빠릅니다. (Ninja is faster than make)

### 5. NPSimulation 빌드 (Building NPSimulation)

이 부분은 Geant4를 이용한 몬테카를로 시뮬레이션을 담당합니다. NPLib이 먼저 컴파일되어 있어야 하며, 그 후 NPSimulation을 컴파일할 수 있습니다.

This part is responsible for Monte Carlo simulations using Geant4. NPLib must be compiled first, and then NPSimulation can be compiled.

#### 1. NPSimulation 폴더로 이동 (Move to the NPSimulation folder)

```sh
nps
```

#### 2. CMake로 Makefile 생성 (Create Makefile with CMake)

```sh
cmake ./
```

#### 3. 컴파일 및 설치 (Compile and Install)

```sh
make -jn install
```

실행 파일: `npsimulation` (Executable: `npsimulation`)

사용 가능한 입력 플래그와 설명은 다음 명령어로 확인할 수 있습니다:

Available input flags and descriptions can be checked with the following command:

```sh
npsimulation -h
```

### 6. 벤치마크 (Benchmarks)

벤치마크는 설치 또는 업그레이드의 무결성 확인, CPU 성능 비교 등에 유용합니다. 두 가지 주요 벤치마크가 제공됩니다.

Benchmarks are useful for checking the integrity of installation or upgrade, comparing CPU performance, etc. Two main benchmarks are provided.

- cats (빔 트래커 데이터 분석)
- gaspard (실리콘 어레이 시뮬레이션)

각 벤치마크는 결과를 그림으로 출력하며, 참조 결과와 비교할 수 있습니다.

Each benchmark outputs results as images and can be compared with reference results.

#### 1. cats 벤치마크 실행 (Run the cats benchmark)

```sh
cd $NPTOOL/Benchmarks/cats
npanalysis -D benchmark_cats.detector -C calibration.txt -R RunToTreat.txt -O benchmark_cats
```

#### 2. gaspard 벤치마크 실행 (Run the gaspard benchmark)

```sh
cd $NPTOOL/Benchmarks/gaspard
npsimulation -D benchmark_gaspard.detector -E 132Sndp_benchmark.reaction -O benchmark_gaspard -B batch.mac
```

#### 3. 결과 확인 (Check the results)

```sh
root -l ShowResult.C
```

---

## 이용 방법 (Usage Guide)

### 1. 프로그램 컨셉 (Program Concept)

nptool은 저에너지 핵물리 실험 데이터를 효율적으로 분석하고 실험 환경을 시뮬레이션할 수 있는 통합 프레임워크입니다. 이 프로그램의 핵심 컨셉은 "실험 데이터의 수집 → 시뮬레이션 → 분석 → 결과 도출"의 전체 과정을 하나의 환경에서 일관되게 처리할 수 있도록 지원하는 것입니다.

nptool is an integrated framework designed to efficiently analyze low-energy nuclear physics experimental data and simulate experimental environments. The main concept of the program is to enable the entire cycle—"data acquisition → simulation → analysis → result output"—to be handled seamlessly within a single environment.

#### 주요 흐름 (Main Flow)

1. **설정 파일 준비 (Prepare Configuration Files)**
   사용자는 detector, reaction, calibration 등 실험에 필요한 설정 파일을 작성합니다. 이 파일들은 실험 환경(검출기 배열, 반응 조건 등)과 분석 조건을 정의합니다.

   Users prepare configuration files (such as detector, reaction, and calibration) that define the experimental environment (detector array, reaction conditions, etc.) and analysis conditions.

2. **시뮬레이션 실행 (Run Simulation)**
   `npsimulation` 프로그램을 통해 Geant4 기반의 몬테카를로 시뮬레이션을 수행합니다. 이 과정에서 실제 실험과 유사한 이벤트 데이터를 생성하며, 결과는 ROOT 파일로 저장됩니다.

   Run a Monte Carlo simulation based on Geant4 through the `npsimulation` program. During this process, events similar to actual experiments are generated, and the results are saved as ROOT files.

3. **데이터 분석 (Analyze Data)**
   `npanalysis` 프로그램을 사용하여 시뮬레이션 또는 실제 실험에서 얻은 데이터를 분석합니다. 분석 과정에서는 사용자가 직접 작성한 분석 클래스(Analysis class)를 로드하여, 원하는 물리량을 추출하거나 추가적인 데이터 처리를 수행할 수 있습니다.

   Analyze the data obtained from simulation or actual experiments using the `npanalysis` program. During the analysis process, the user can load the analysis class (Analysis class) they have written to extract desired physical quantities or perform additional data processing.

4. **결과 확인 및 시각화 (Check the results and visualize)**
   분석 결과는 ROOT 파일로 저장되며, 사용자는 ROOT 환경에서 히스토그램, 그래프 등 다양한 방식으로 결과를 시각화하고 해석할 수 있습니다.

   The analysis results are saved as ROOT files, and users can visualize and interpret the results in various ways using the ROOT environment, such as histograms and graphs.

### 2. 각 프로그램의 사용법 (Usage of each program)

#### **npsimulation**

**npsimulation**은 Geant4를 이용한 시뮬레이션 실행 프로그램으로 일반적으로 다음과 같이 사용됩니다.

**npsimulation** is a program for running Geant4-based Monte Carlo simulations. It is typically used as follows.

```sh
npsimulation -D <detector_file> -E <reaction_file> (-B <batch_file> -O <output_file>)
```

<detector_file>은 검출기 설정 파일, <reaction_file>은 반응 조건 설정 파일, <batch_file>은 배치 모드 설정 파일, <output_file>은 출력 파일 이름입니다.

<detector_file> is the detector configuration file, <reaction_file> is the reaction condition configuration file, <batch_file> is the batch mode configuration file, and <output_file> is the output file name.

#### **npanalysis**

**npanalysis**은 실험 데이터 또는 시뮬레이션 데이터를 분석하는 프로그램으로 일반적으로 다음과 같이 사용됩니다.

**npanalysis** is a program for analyzing experimental or simulation data. It is typically used as follows.

```sh
npanalysis -T <tree_name> <file_name> -O <output_file>
```

또는 간단하게 마지막으로 실행한 시뮬레이션 파일을 분석할 수 있습니다.

Or, you can analyze the last simulation file simply.

```sh
npanalysis --last-sim -O <output_file>
```

### 3. 프로젝트 디렉토리 (Project Directory)

nptool의 프로그램들은 어느 디렉토리에서도 실행할 수 있습니다. 그러나 일반적으로 프로젝트 디렉토리 (Projects) 내에 관련된 하위 디렉토리를 생성하고 그 밑에서 실행하는 것을 권장합니다.

nptool's programs can be run from any directory. However, it is recommended to create a related subdirectory under the project directory (Projects) and run it there.

```sh
cd $NPTOOL/Projects/<project_name>
```

### 4. 예제 (Examples)

예제는 여러 검출기를 조합한 복잡한 분석 사례를 다룹니다. Example1 실행 예시는 다음과 같습니다.

Examples deal with complex analysis cases combining multiple detectors. The following is an example of running Example1.

```sh
npsimulation -D Example1.detector -E Example1.reaction -O Example1
```

이후 GUI 또는 프롬프트에서 이벤트를 생성.

Then, create events in the GUI or prompt.

```text
> run/beamOn/ 10000
> exit
```

#### 분석 (Analysis)

```sh
npp Example1
cmake ./
make -jn
npanalysis -R RunToTreat.txt -O Example1
```

#### 결과 확인 (Check the results)

```sh
root -l ShowResult.C
```

### 5. 검출기별 메뉴얼 (Detector-specific manuals)

검출기별 메뉴얼은 각 검출기의 설정 및 사용법을 설명합니다.

Detector-specific manuals explain the setup and usage of each detector.

- [Plunger](Documentation/detectors/Plunger.md)
