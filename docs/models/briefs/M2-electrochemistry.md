# Brief M2 — Electrochemical cell geometry, current distribution and loading

Branch: `claude/lucid-davinci-gel4vu-m2`. Read `_common.md` first.

## Questions
1. Which electrode geometry maximises the **fraction of cathode area at D/Pd ≥ 0.90 (ideally ≥ 0.95)** while avoiding current hot-spots, excessive ohmic heating and bubble screening?
2. For the **detector-facing membrane (DFM, config C3)**, what anode geometry and masking give a current distribution within ±5 % over a 10–25 mm diameter Pd disk that forms one wall of the cell? The other face of the disk is in vacuum.
3. What operating window (current density, electrolyte, temperature) maps to the loading required?

## Scope
- Geometries:
  - (a) coaxial wire cathode (diameter 0.25–2 mm, length 1–5 cm) inside a cylindrical Pt mesh or helix anode (radius 5–20 mm)
  - (b) planar foil between two parallel anode meshes
  - (c) DFM disk cathode as a cell wall, with a disk, mesh or ring anode at gap g = 2–30 mm, with or without an insulating PTFE collar/shield over the disk edge
  - (d) sphere-in-sphere reference
- Primary current distribution (Laplace). Secondary distribution with hydrogen-evolution kinetics on Pd in 0.1 M LiOD/D₂O: Butler–Volmer/Tafel. Take the exchange current density and Tafel slope from literature with citations; compute the Wagner number. Use axisymmetric finite differences or finite volumes in `scipy.sparse`.
- Loading model: map local current density → overpotential → effective D fugacity (Volmer–Tafel and Volmer–Heyrovsky mechanisms, both cited) → loading x through the high-pressure Pd–D isotherm (Baranowski-type data). Use literature correlations of x vs current density (Akita, Kunimatsu, McKubre/SRI, Crouch-Baker) as validation. Quantify the widely reported need for ≳ 10³–10⁴ atm equivalent fugacity to reach x ≥ 0.9.
- Bubble void fraction and its effect on conductivity (Bruggeman). Ohmic heating and electrolyte temperature rise. iR drop and cell voltage. Required supply voltage.
- Surface additives reported to raise loading (Al, Si, thiourea, Pt deposits): quantify only if a literature model exists; otherwise list them.
- Internal recombiner (Pt or Pd on alumina): recommended location relative to electrodes and detectors, gas-flow geometry, heat load. Provide the stoichiometric recombination heat at a given current: 1.48 V thermoneutral equivalent for H₂O and the D₂O value, with citations.

## Required outputs
- Maps of current density and loading for each geometry. A table of area fraction with x ≥ 0.90 and x ≥ 0.95 versus geometry and total current.
- The recommended geometry with dimensions and tolerances, for both the DFM (C3) and the coaxial calorimetric cell (C1).
- The operating window: current-density ramp schedule (for example 10 → 100 → 500 mA/cm²), electrolyte concentration, temperature.
- Interfaces: supply the back-face fugacity/loading boundary condition to M3 (permeation) as a table of x_in vs current density.
