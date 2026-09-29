# Brief M3 — Deuterium transport, permeation, phase transition and mechanics

Branch: `claude/lucid-davinci-gel4vu-m3`. Read `_common.md` first.

## Questions
1. Loading time versus characteristic dimension for wires, foils, films and nanoparticles, including the α/β miscibility gap and the moving phase boundary.
2. For the **detector-facing membrane (C3)**, loaded at the back face at x_in = 0.85–0.97 with the front face in vacuum: steady-state and transient profiles, permeation flux J, and **loading at the exit face** as functions of thickness (5–200 µm) and exit-face recombination rate constant k_r. Compare clean Pd, PdO, Au overlayer, CaO or Pd/CaO multilayer (Iwamura) and Ni/Cu multilayer. Which exit-face barrier keeps the entire membrane at x ≥ 0.9 while maintaining a large flux?
3. Mechanics: when does loading, deloading or gradient cycling crack or buckle the membrane or film? Use Vegard strain (Δa/a ≈ 3.5 % α→β; partial molar volume of H/D in Pd ≈ 1.7 cm³/mol, cited), yield strength of annealed vs cold-worked Pd, and a free-standing vs supported membrane with edge clamping. Consider the pressure differential across a vacuum-backed membrane (electrolyte at 1 atm on one side).
4. Two regimes, with a protocol for each:
   - (a) **high-loading steady state**: avoid cracks
   - (b) **controlled non-equilibrium**: repeated α/β cycling to create fresh crack faces and vacancies on purpose, for fracto-emission and defect-site generation. Estimate new crack area per cycle, and vacancy concentration (Fukai superabundant vacancies: conditions, and whether room-temperature electrodeposition or electrolysis produces them, with cited estimates).
5. Temperature dependence (20–90 °C) of diffusion, solubility and the miscibility gap. Isotope effect D vs H, because the H₂O control must match.

## Scope and methods
- Composition-dependent chemical diffusion coefficient (Darken factor from the isotherm), cited D₀ and Eₐ for D in Pd.
- 1D finite-volume solver with a nonlinear isotherm. 2D axisymmetric solver for the membrane with a clamped edge if needed.
- Surface kinetics: absorption/desorption via Pick's model or equivalent, with rate constants cited for each surface condition. Where k_r is unknown, scan it over orders of magnitude and state the required value.
- Stress: thin-film/membrane elasticity (Stoney for supported films, plate theory for clamped membranes), hydrogen-induced gradient stress, and critical conditions for dislocation generation and cracking.

## Required outputs
- A loading-time table (t for 90 % of equilibrium) vs dimension.
- A DFM design chart: thickness × exit-face k_r → exit-face loading and flux, with the recommended region marked.
- Recommended membrane thickness, diameter, clamping and support geometry (for example a perforated support grid on the vacuum side, with open fraction versus stress), annealing state and grain size.
- Two cycling protocols (current-amplitude waveform, period, number of cycles), with predicted crack density and loading history.
- Interfaces: take the back-face x_in vs current density from M2. If M2 is not yet available, use x_in ∈ {0.85, 0.90, 0.95, 0.97} as the parameter.
