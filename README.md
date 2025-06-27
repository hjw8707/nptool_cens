CENS NPTool
======

## 프로젝트 설명

nptool은 저에너지 핵물리 실험을 위한 데이터 분석 및 몬테카를로 시뮬레이션 패키지입니다. 이 패키지는 복잡한 실험을 준비하고 분석하기 위한 통합 프레임워크를 제공하며, Geant4 및 ROOT 툴킷을 효율적으로 활용합니다. nptool은 사용자가 실험 데이터를 쉽게 처리하고 분석할 수 있도록 다양한 기능을 제공합니다. 이 프로젝트는 오픈 소스이며, 사용자와 개발자들이 자유롭게 기여할 수 있도록 설계되었습니다.

본 프로젝트는 CENS(센터명: Center for Exotic Nuclear Studies)에서 개발 및 유지보수하고 있는 NPTool의 공식 저장소입니다. CENS NPTool은 저에너지 핵물리 실험을 위한 데이터 분석과 몬테카를로 시뮬레이션을 지원하며, CENS의 연구 환경과 실험 요구에 최적화되어 있습니다. CENS 연구진 및 협력 연구자들이 실험 데이터를 효율적으로 처리하고 분석할 수 있도록 다양한 기능과 예제, 문서가 제공됩니다.

이 패키지는 다음과 같은 주요 기능을 포함합니다:

- **데이터 분석**: 실험 데이터를 처리하고 분석하는 데 필요한 다양한 도구를 제공합니다.
- **몬테카를로 시뮬레이션**: Geant4를 사용하여 복잡한 물리적 현상을 시뮬레이션할 수 있습니다.
- **유연한 구성**: 사용자가 필요에 따라 다양한 분석 및 시뮬레이션 설정을 조정할 수 있습니다.
- **커뮤니티 지원**: 사용자와 개발자들이 서로의 경험을 공유하고 문제를 해결할 수 있는 플랫폼을 제공합니다.

이 프로젝트는 저에너지 핵물리 분야의 연구자들에게 유용한 도구가 될 것입니다. 사용자는 이 패키지를 통해 실험 데이터를 보다 효율적으로 분석하고, 새로운 물리적 현상을 탐구할 수 있습니다.

nptool, which stands for Nuclear Physics Tool, is an open source and freely
distributed data analysis and Monte Carlo simulation package for low-energy
nuclear physics experiments. The NPTool package aims to offer an unified
framework for preparing and analysing complex experiments, making an efficient
use of Geant4 and ROOT toolkits.

- If you wish to contribute, contact Adrien Matta: <matta@lpccaen.in2p3.fr>
- For issues, use the git issue tracker or contact  Adrien Matta : <matta@lpccaen.in2p3.fr>

# Usefull links

- [nptool website](http://nptool.org) : Find the latest informations on the framework
- [nptool manual](http://nptool.org/manual) : Find a more detail manual on how to install and run the framework

# Contents

1. [Getting the code](#getting-the-code)
1. [Using git](#using-git)
2. [Downloading from Git Hub](#downloading-from-git-hub)
2. [Setup](#setup)
1. [Requirements](#requirements)
2. [Preparing the build](#preparing-the-build)
3. [Building NPLib](#building-nplib)
4. [Building NPSimulation](#building-npsimulation)
3. [Benchmarks and Examples](#benchmarks-and-examples)
1. [Benchmarks](#benchmarks)
2. [Examples](#examples)
4. [Tricks](#tricks)

## Getting the code

### Using git

The recommended method to obtain the source code is to use git. This is an
easy way to access the last version of the code. First make sure you have git
installed. If not, use your package manager to get it. Then go to the directory
where you want to install the NPTool package and do:

```
git clone https://gitlab.in2p3.fr/np/nptool
```

This will create the nptool folder with the latest version of nptool.

### Downloading from Git Hub

Alternatively, you can browse the following page <https://gitlab.in2p3.fr/np/nptool>,
and click the Download ZIP button on the right side of the page. Then, unzip
the archive at the desire location

## Setup

### Requirements

NPTool components are compiled and installed using the CMake build system,
so be sure to have a working CMake installation before starting.

In order to compile NPLib, the NPTool core libraries, ROOT 5 (tested with 5.34)
or 6 should be installed. This is sufficient to compile NPLib and any analysis project.

In order to compile NPSimulation, a recent installation of Geant4 (tested
with version 9.6 and 10.1) is needed. If you want to use GDML support in
NPTool, Geant4 should be installed with GDML support.

### Preparing the build

To set the needed environment variables, PATH and LD\_LIBRARY\_PATH, and
aliases, source the following script doing:

```
source <pathname>/nptool/nptool.sh
```

where `<pathname>` is the location where you unpacked the NPTool package.
Then, restart your terminal.

You should typically add the previous command line to your .profile,
.bashrc or .tcshrc file.

### Building NPLib

NPLib is the core of the NPTool package, holding most of the actual code. It is
made of a collection of stand alone C++ classes that can be used in programs
and macros.

First, go to the NPLib folder by using the command:

```
npl
```

In order to prepare the compilation CMake must be run to generate the Makefile.
If no arguments are given to CMake, all detectors will be compiled. If you wish
to limit the number of detectors to be compiled, specify the detector folder
name (respecting the case). Note that more than one detector can be specified.

All detectors compiled:

```
cmake ./
```

**OR** some detectors compiled:

```
cmake ./ -DNPTOOL_DETLIST="DetFolder1 DetFolder2"
```

Then, the whole NPLib can be compiled with _n_ threads using:

```
make -jn install
```

If you wish to recompile in order to get support for more detectors, do:

```
nptool-cleaner
cmake ./ -DNPTOOL_ETLIST="DetFolder1 DetFolder2 ..."
make -jn install
```

If you have google ninja build installed you can alternatively ask CMake to
generate the ninja.build file:

```
cmake -GNinja ./
ninja install
```

Compilation using Ninja is faster than using make.

### Building NPSimulation

This part of the package relies on Geant4 to perform Monte Carlo simulation.
NPLib needs first to be compiled and configured correctly before NPSimulation
can be compiled. The compilation is done as follow:

```
nps
cmake ./
make -jn install
```

This will produce the _npsimulation_ executable. For a detailed list of the
available input flags and their meaning, run the following command:

```
npsimulation -h
```

## Benchmarks and Examples

You need to download additional files to be able to run the benchmarks and
the examples. In the $NPTOOL directory, do the following:

```
git clone https://github.com/adrien-matta/NPData

```

### Benchmarks

Benchmarks play an important role to check the installation or upgrade integrity
of the NPTool package. They are also useful for comparing CPU performances on
different computer platforms. So far two benchmarks are included in the NPTool
package. The first one (_cats_) analyses experimental data from a beam tracker
detector using the _npanalysis_ facility, while the second one (_gaspard_) runs
a silicon array simulation using the _npsimulation_ facility and display some
basics control spectra. Each benchmark produces figures that can be compared to
the reference figures provided in NPTool. These two benchmarks cover all the core
functionalities of NPTool's framework.

The first benchmark can be run with the following commands:

```
cd $NPTOOL/Benchmarks/cats
npanalysis -D benchmark_cats.detector -C calibration.txt -R RunToTreat.txt -O benchmark_cats
```

For the second benchmark do:

```
cd $NPTOOL/Benchmarks/gaspard
npsimulation -D benchmark_gaspard.detector -E 132Sndp_benchmark.reaction -O benchmark_gaspard -B batch.mac
```

In both cases, the results can be displayed and compared to reference results
using the following command:

```
root -l ShowResult.C
```

### Examples

With respect to benchmarks, examples deal with more complex analysis cases
where several detectors are present. In this context, physical information
resulting from the combination of information from several detectors can be
calculated.

As a standardized test case, you can run Example1 by the following command:

```
npsimulation -D Example1.detector -E Example1.reaction -O Example1
```

This will open the _npsimulation_ GUI (if you are using Qt) or the prompt
terminal. In either case events can be generated using the following
command:

```
> run/beamOn/ 10000
> exit
```

This will simulate the <sup>11</sup>Li(d,<sup>3</sup>He)<sup>10</sup>He-><sup>8</sup>He+n+n
reaction and produce an output ROOT file located in $NPTOOL/Outputs/Simulation/Example1.root.

The `Example1.detector` file located in $NPTOOL/Inputs/DetectorConfiguration
and the `Example1.reaction` file located in $NPTOOL/Inputs/EventGenerator are
self explanatory thanks to the use of understandable tokens.

The simulated file can be analysed using the following commands:

```
npp Example1
cmake ./
make -jn
npanalysis -R RunToTreat.txt -O Example1
```

This will produce an analysed tree located in the $NPTOOL/Outputs/Analysis/Example1.root file.
Note that since the input files needed by _npsimulation_ are stored in the
ROOT file, _npanalysis_ use these automatically as inputs.

The results can be displayed using the command:

```
root -l ShowResult.C
```

You should be able to see the light particle identification, the light particle kinematical line and the associated excitation energy spectrum fitted by a gaussian.

The Example1 input files and NPAnalysis project are simple basis that can be used to start doing your own simulations.

## Tricks

- _npsimulation_ and _npanalysis_ facilities can be run from any directory.
- _npanalysis_ looks in the current directory for an analysis library
  (_libNPAnalysis_) to load and use. If not present, the analysis is limited
  to build the PhysicsTree.
- To perform a quick analysis of the last MOnte Carlo simulation, do:

```
npanalysis --last-sim
```

Any additional flag can be used

- _npsimulation_ can be used in batch mode (with no UI) if the user provides
a Geant4 macro file specified with the -B flag

```
npsimulation -D Example1.detector -E Example1.reaction -B path/to/macro.mac -O FileName
```

## npreader

- _npreader_ is a new executable, working similarly to _npanalysis_. It also uses _Analysis.cxx_/_Analysis.cxx_ to convert raw data to physics tree.
- _npreader_ uses the TTreeReader.next() method to retrieve data event by event. A software trigger condition can be set to filter the data.
  Depending on the trigger condition, a small part of the data can be extracted by the _UnallocateBeforeBuild()_ method to check the condition before building detectors.
  The method then returns _true_ or _false_ depending on this trigger, and will continue the analysis in case of _true_, but will stop and go to the next event in case of _false_.
  Conditions can also be set after building detectors but before treating the analysis using the _UnallocateBeforeTreat()_ method
  or after treating the analysis but before writing the tree using the _FillOutputCondition()_ method.
  These method are useful to treat and write only relevant information which accelerates significantly the execution of the analysis.
- An example of such analysis using the new npreader utility can be found in project E805.
- Currently, the _npreader_ utility is only implemented for MUST2, Exogam, CATS and TAC detectors.

The syntax is the same as npanalysis:

```
npreader -D /path/to/detector.detector -E /path/to/reaction.reaction -R RunToTreat.txt -O OutputFile -C /path/to/calibration.txt
```

## npcalibration

- _npcalibration_ is an executable intended for detector calibrations. It does not rely on _Analysis.cxx_,
  and only uses standard calibration functions in several detector physics libraries.
- _npcalibration_ requires configuration files. Examples of such configuration files can be found
  in project E805 DoCalibration folder.
- The motivation behind _npcalibration_ is to get standard calibration files that can be directly used
  by _npanalysis_ or _npreader_. A root file is also created after the calibration has been performed, it is recommended to write
  the relevant histograms and fit functions used during the calibration to this file to check that everything
  worked as expected.
- Currently _npcalibration_ supports MUST2 (DSSD energy and CsI energy) and Exogam (Crystals energy) calibration.
- Output calibration files can be found in the Calibration folder. This folder can be sourced in your
  project.config.
- Some calibrations require histogram CUTS (for instance particle identification in CsI detectors of MUST2).
  A path to these CUTS can also be sourced in your project.config. (Example can be found in Project E805: project.config)

The syntax is the following:

```
npcalibration -DC /path/to/DoCalibration.do -R RunToTreat.txt -O OutputFile
```

If necessary, a preliminary calibration file can be added (for instance, DSSD energy calibration is necessary for MUST2
CsI calibration)

```
npcalibration -DC /path/to/DoCalibration.do -R RunToTreat.txt -O OutputFile -C /path/to/calibration.txt
```

Exogam calibration requires the CUBIX software <https://cubix.in2p3.fr/>. If CUBIX is installed, it needs to be sourced
in your environment. To compile nptool with CUBIX, add to your cmake command:

```
-DCUBIX=1
```
