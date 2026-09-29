# Brief M4 — Calorimetry and thermal design

Branch: `claude/lucid-davinci-gel4vu-m4`. Read `_common.md` first.

## Questions
1. What calorimeter design gives the best power resolution σ_P (mW) and systematic accuracy (% of input) for a **closed** electrolytic cell? The cell has 30–100 mL electrolyte, an internal recombiner and 0.5–10 W input (config C1 and the electrolyte side of C3).
2. How position-dependent is the calibration constant? Heat released at the cathode vs the recombiner vs Joule heating in the electrolyte. Which geometry makes it position-independent to < 0.1 %?
3. Quantify the Shanahan "calibration constant shift" (CCS) critique for each design and state the design features that bound it.
4. For gas-phase heated reactors (C4/C5: 200–300 °C, 10–100 W heater), what accuracy is realistic?
5. Can a heat channel usefully be added to the DFM cell (C3) without compromising the charged-particle detectors?

## Scope and methods
- Compare isoperibolic (Fleischmann–Pons Dewar), mass-flow, and Seebeck envelope (heat-flux, thermopile or Peltier-module wall) calorimeters, and a **twin differential** arrangement (D₂O cell vs H₂O control cell in the same thermal environment).
- Build a lumped plus spatially resolved thermal model: axisymmetric conduction with an effective convection coefficient in the electrolyte, and the envelope/thermopile. Compute:
  - the sensitivity (µV/W or K/W)
  - the time constants
  - the noise floor from sensor noise and environmental drift (specify the required enclosure stability, e.g. ±0.01 K)
  - the calibration method (resistive heater at the cathode position and at the recombiner position; pulse and step calibrations)
- Recombiner: heat of recombination at a given current, and the effect of recombination inefficiency on apparent excess power (open vs closed cells). Include cited thermoneutral potentials for D₂O and H₂O (1.527 V for D₂O vs 1.481 V for H₂O).
- Error budget table for each design, including electrical-power measurement with AC ripple from pulsed or modulated currents. This matters because we plan current modulation.

## Required outputs
- The recommended calorimeter geometry: dimensions, materials, sensor types and placement, calibration heaters.
- σ_P and systematic accuracy for each option. The design choice with its reasoning.
- An error budget and the pre-registered calibration protocol.
- Guidance for a heat–⁴He co-measurement: the all-metal sealed cell volume, headspace sampling, and why the cell must avoid glass (He permeation).
