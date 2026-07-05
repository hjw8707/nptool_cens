# CSSU CF4 Gas Detector Example

This example ports the simple gas volume from `ejungwoo/cssu_geant4_simulation`
into the NPTool `CSSU` detector package.

The detector is a `120 x 120 x 200 mm` CF4 gas box at `200 Torr` with a
`10 mm` step limit. It is centered at `z=100 mm`, so an alpha source at the
origin travels along `+Z` into the entrance plane. `SegmentsZ` splits the active
gas volume into separated sensitive slices along `Z`; `SegmentGapZ` leaves CF4
gas between active slices, so particles continue losing energy in the gaps but
those gaps are not recorded as active hits. The active slice deposits are summed
event-by-event and stored as one total gas detector energy in the `CSSU` branch.

Run:

```sh
source ../../../nptool.sh
./run.sh
```

Viewer:

```sh
./viewer.sh
```

Main files:

- `detector.txt`: `CSSU Type= GasBox` detector block.
- `reaction_alpha_am241.txt`: 5.443 MeV alpha source along `+Z`.
- `plot_gas_energy.C`: creates deposited-energy histograms from the `CSSU` branch.

Outputs:

- `root/sim/cssu_cf4_alpha.root`
- `root/ana/cssu_cf4_alpha.root`
- `root/ana/cssu_cf4_alpha.png`
