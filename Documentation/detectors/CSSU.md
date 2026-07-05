# CSSU

`CSSU` is a simulation-oriented detector package for CSSU detector studies. It
currently provides two detector types:

- `ScintillatorBar`: a BC400-like plastic scintillator bar with one cylindrical
  PMT counter at each end.
- `GasBox`: a CF4 gas volume with optional separated active slices.

## Requirements

Optical photon generation for `ScintillatorBar` requires optical physics:

```txt
OpticalPhysics 1
```

This can be enabled in `NPSimulation/PhysicsListOption.txt` or in the
project-specific physics-list option file.

## ScintillatorBar Input

```txt
CSSU
 Type= ScintillatorBar
 POS= 0 0 50 mm
 Length= 200 mm
 Width= 20 mm
 Thickness= 40 mm
 PMTLength= 100 mm
 PMTDiameter= 40 mm
 ScintillationYield= 10000
 AttenuationLength= 2000 mm
 ScintillatorRIndex= 1.58
```

`ScintillationYield` is photons per MeV and is written without a unit.
`ScintillatorRIndex` is unitless and defaults to `1.58`. `PMTThickness` and
`PMTFace` are accepted as legacy aliases for `PMTLength` and `PMTDiameter`.
The bar axis is the global X axis:

```txt
PMT 1 -- scintillator bar -- PMT 2
  -X                         +X
```

## GasBox Input

```txt
CSSU
 Type= GasBox
 POS= 0 0 100 mm
 Length= 120 mm
 Width= 120 mm
 Thickness= 200 mm
 PressureTorr= 200
 StepLimit= 10 mm
 SegmentsZ= 4
 SegmentGapZ= 5 mm
 EnergyThreshold= 0.001 MeV
 ResoEnergy= 1 keV
```

`GasBox` creates a non-sensitive CF4 gas mother volume at `POS`. The active gas
slices are placed inside that mother volume along Z. `SegmentGapZ` leaves
non-sensitive CF4 gas between active slices, so particles continue losing energy
in the gaps but only active slices are scored.

The active slice deposits are summed event-by-event. The `CSSU` output branch
stores one total gas energy per `GasBox`, not one entry per slice.

Gas parameters:

- `PressureTorr`: CF4 pressure in Torr. Density and pressure scale from the
  `200 Torr` reference used in `ejungwoo/cssu_geant4_simulation`.
- `StepLimit`: Geant4 maximum step in the gas volumes.
- `SegmentsZ`: number of active slices along Z. `ActiveSegmentsZ` and
  `NbSlices` are accepted as aliases.
- `SegmentGapZ`: gap between active slices. `GapZ` is accepted as an alias.
- `EnergyThreshold`: threshold applied after summing active-slice energy.
- `ResoEnergy`: Gaussian sigma applied to the summed active energy.

## Output

The simulation writes a `CSSU` branch containing `TCSSUData`.

Energy entries contain:

- detector number
- measured deposited energy
- time

For `ScintillatorBar`, the energy is the bar energy deposit. For `GasBox`, it is
the event-by-event total of all active slices after threshold and resolution.

PMT entries, used by `ScintillatorBar`, contain:

- detector number
- PMT number (`1` for left, `2` for right)
- optical photon count
- first photon time

## Examples

- `Projects/CSSU/ScintPMT`: plastic scintillator and PMT optical-photon example.
- `Projects/CSSU/GasCF4`: CF4 gas detector based on
  `ejungwoo/cssu_geant4_simulation`.

Run the gas example:

```sh
cd Projects/CSSU/GasCF4
./run.sh
```

Open the gas viewer:

```sh
./viewer.sh
```

## Current Scope

Implemented:

- `ScintillatorBar` geometry
- `GasBox` CF4 mother gas volume
- separated active gas slices with non-sensitive gas gaps
- event-level summed gas energy response
- BC400-like optical scintillator material
- left/right cylindrical PMT sensitive volumes
- scintillator energy scorer
- PMT optical photon count scorer

Not implemented yet:

- arbitrary detector rotation
- gas ionization charge, drift, diffusion, avalanche gain, or electronics
- wire, strip, pad, or Micromegas readout geometry
- reflective wrapping geometry
- PMT quantum efficiency, gain, or transit-time spread
