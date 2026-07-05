# CSSU ScintPMT

This example shoots an approximate 241Am alpha source into a long plastic
scintillator bar and counts optical photons at two PMT faces.

## Files

- `detector.txt`: CSSU scintillator bar, 200 mm x 20 mm x 40 mm, with one
  cylindrical PMT photon counter on each long-axis end.
- `detector_viewer.txt`: same geometry for the OpenGL viewer.
- `reaction.txt`: one 5.486 MeV alpha per event, emitted toward the bar.
- `reaction_with_am241_gamma.txt`: same alpha plus one 59.54 keV gamma per
  event for response tests.
- `PhysicsListOption.txt`: enables Geant4 optical physics.
- `batch.mac`: runs 1 event by default because optical photon tracking is
  expensive.
- `geant4_vis.mac`: viewer macro for Geant4 OpenGL visualization.
- `viewer.sh`: runs `detector_viewer.txt` and `reaction.txt` with
  `geant4_vis.mac`.
- `run_optical_stats.sh`: runs one event and writes optical photon death
  statistics to `optical_photon_death_stats.txt`.
- `vis.mac`: compatibility wrapper that executes `geant4_vis.mac`.

## Expected photon numbers

Real 241Am produces one alpha per decay. The dominant alpha line is about
5.486 MeV, with lower-intensity alpha lines near 5.443 MeV and 5.388 MeV. It
also emits a 59.54 keV gamma in about 36% of decays.

The default `reaction.txt` does not generate the nuclear gamma; it generates
only the alpha. The PMT signals in this example are optical scintillation
photons from the plastic, not gamma rays.

In this CSSU implementation, optical photons are generated approximately as:

```text
number of scintillation photons = ScintillationYield * deposited energy
```

The default `detector.txt` uses `ScintillationYield=1000` photons/MeV so the
example can be run quickly. With full 5.486 MeV alpha energy deposit, this
produces an upper limit of about 5,500 optical photons per event.

For a plastic-scintillator-like light yield, change `ScintillationYield` to
`10000`. Then the no-quenching upper limit is about 55,000 optical photons per
event, but the Geant4 optical tracking will be much slower.

The PMT branches count how many of those optical photons reach each PMT face.
This first implementation does not yet model alpha quenching, PMT quantum
efficiency, photocathode gain, or realistic electronics.

## Run

```sh
source ../../../nptool.sh
./run.sh
```

The output ROOT file is written to:

```text
Projects/CSSU/ScintPMT/root/sim/scint_pmt_alpha.root
```

To open the Geant4 viewer:

```sh
./viewer.sh
```

To check where optical photons are killed:

```sh
./run_optical_stats.sh
cat optical_photon_death_stats.txt
```

To include the 59.54 keV gamma test source, run:

```sh
npsimulation -D detector.txt -E reaction_with_am241_gamma.txt -B batch.mac -O scint_pmt_alpha_gamma.root -N
```

Remember that `reaction_with_am241_gamma.txt` emits the gamma every event. For
241Am-like gamma rates, scale gamma-sensitive observables by about 0.36.
