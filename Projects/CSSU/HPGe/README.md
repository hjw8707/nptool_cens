# CSSU HPGe Natural Background Test

This directory contains a minimal coaxial HPGe simulation for ambient gamma-ray background tests.

Run from the repository environment:

```sh
source ../../../nptool.sh
./run.sh
```

Main files:

- `hpge.detector`: one `Coaxial_Germanium` detector at 150 mm on the +X axis. `ResoEnergy` is Gaussian sigma, so `0.85 keV` is about `2 keV` FWHM.
- `natural_background.reaction`: isotropic gamma event-generator input at x=-400 mm using a ROOT energy histogram.
- `background_plus_co60_x-800.reaction`: natural background at x=-400 mm plus Co-60 1173/1332 keV gamma lines from x=-800 mm.
- `make_natural_gamma_hist.C`: creates `natural_gamma_lines.root` with representative K-40, U-chain, and Th-chain gamma lines.
- `run_background_plus_co60.sh`: runs the two-source example and writes an HPGe deposited-energy histogram.
- `run_co60_position_scan.sh`: scans Co-60 source positions outside the detector. Override defaults with `POSITIONS="-800 -1200 -1600" EVENTS=50000 CONE_DEG=3 ./run_co60_position_scan.sh`.
- `viewer_hpge_co60.sh`: opens the Geant4 viewer with the natural source at x=-400 mm and Co-60 at x=-800 mm.
- `plot_hpge_energy.C`: converts `SimulatedTree` HPGe energy deposits into ROOT and PNG histograms.
- `plot_hpge_scan_overlay.C`: overlays the scan histograms without normalizing counts.
- `batch.mac`: runs 10000 events.

Output:

- `root/sim/hpge_natural_background.root`: simulation output containing the `Coaxial_Germanium` branch.
- `root/sim/background_plus_co60_*.root`: scan outputs.
- `root/ana/background_plus_co60_*.root` and `.png`: deposited-energy histograms.
- `root/ana/background_plus_co60_overlay.png`: count-preserving position comparison.
