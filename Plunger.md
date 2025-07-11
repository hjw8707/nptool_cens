# Plunger

## 1. Introduction

Plunger is the measurement system mainly composed of a target and a stopper used to measure the lifetime of a nucleus in an excited state.

## 2. Detector file (.detector)

The Plunger section of the `.detector` file defines the geometry and materials for the main components of the plunger system. Each component is specified with a block starting with its name, followed by its parameters:

- **Plunger Target**
  - `R`: Radius of the target (e.g., `20 mm`)
  - `Thickness`: Thickness of the target foil (e.g., `0.05 mm`)
  - `Z`: Position along the beam axis (e.g., `-10 mm`)
  - `Material`: Material of the target (e.g., `Al` for aluminum)

- **Plunger Stopper**
  - `R`: Radius of the stopper (e.g., `20 mm`)
  - `Thickness`: Thickness of the stopper (e.g., `1 mm`)
  - `Z`: Position along the beam axis (e.g., `0.65 mm`)
  - `Material`: Material of the stopper (e.g., `Al`)

- **Plunger Chamber**
  - `R`: Inner radius of the chamber (e.g., `100 mm`)
  - `Thickness`: Thickness of the chamber wall (e.g., `1 mm`)
  - `Material`: Material of the chamber (e.g., `Al`)
  - `PipeR`: Radius of the beam pipe (e.g., `30 mm`)
  - `PipeZ0`: Start position of the pipe along the beam axis (e.g., `-120 mm`)
  - `PipeZ1`: End position of the pipe along the beam axis (e.g., `120 mm`)

Example:

```txt
Plunger Target
 R = 20 mm
 Thickness = 0.05 mm
 Z = -10 mm
 Material = Al

Plunger Stopper
 R = 20 mm
 Thickness = 1 mm
 Z = 0.65 mm
 Material = Al

Plunger Chamber
 R = 100 mm
 Thickness = 1 mm
 Material = Al
 PipeR = 30 mm
 PipeZ0 = -120 mm
 PipeZ1 = 120 mm
```
