# Brief M5 — Detection geometry, backgrounds and statistics

Branch: `claude/lucid-davinci-gel4vu-m5`. Read `_common.md` first.

## Questions
1. For reaction products born at depth z in the active layer, what fraction escapes and reaches each detector, and with what energy spectrum? Products: 3.02 MeV p, 1.01 MeV T, 0.82 MeV ³He, and α particles over a range of energies. Materials: Pd, PdD₀.₉, PdO, Au, Ni, Cu, CaO, D₂O, H₂O, Mylar, Al, D₂ gas, Si.
2. For the **detector-facing membrane (C3)**, what is the optimal Si spectrometer geometry, including a ΔE–E telescope for particle ID? Specify area, distance, collimation, thicknesses and vacuum level. Give the expected background in the 2.6–3.1 MeV proton window, and the minimum detectable reaction rate.
3. What CR-39 stack design (Mylar/Al filter thicknesses) discriminates p, T and α? What background and artifact controls does it need?
4. What neutron detection system gives the best 2.45 MeV sensitivity for a tabletop cell? Options: moderated ³He or BF₃ bank, EJ-309/stilbene with PSD for energy information, and bubble detectors. What background is achievable at sea level with feasible shielding and a muon veto?
5. Statistics: pre-registered analysis, an on/off current-modulation (lock-in) test, D vs H control cells, blind analysis, and look-elsewhere correction across channels. What run time reaches 5σ for source rates 10⁻⁴–10⁰ s⁻¹?

## Scope and methods
- Stopping powers: implement Bethe–Bloch with shell/Barkas corrections, or use NIST PSTAR/ASTAR tabulations and Bragg's rule for compounds, and cite the mean excitation energies. Validate against NIST and SRIM values in the literature within ~5 %. Compute ranges and energy loss.
- Monte Carlo (numpy) of isotropic emission from depth distributions (uniform over the top 1, 10 or 100 µm, or δ at the surface) through overlayers to a detector of given area and distance. Output spectra, including a detector resolution of ~20–30 keV FWHM for Si.
- Background model for Si in vacuum: cosmic muons (MIP energy deposit per thickness), alpha emitters in construction materials (U/Th chains, ²¹⁰Po), and noise. Show the benefit of ΔE–E coincidence, pulse-shape cuts and anticoincidence veto. Cite measured backgrounds from low-background charged-particle work where possible (Lipson, Kasagi, the Google/Nature 2019 team).
- Neutrons: a simple MC for moderated ³He efficiency versus HDPE thickness and geometry, validated against published efficiency curves. Scintillator PSD options. Cosmic neutron background with its barometric and solar modulation.
- Statistics: Currie L_C and L_D, Feldman–Cousins, and Bayesian alternatives. The power of a modulation test at duty cycle 50 % and period T_mod. Worked examples.
- Controls: identical H-loaded twin, Pt dummy cathode, detector blank, and source calibrations (²⁴¹Am α, and a ²⁵²Cf neutron source if available, or AmBe; flag the licensing needs).

## Required outputs
- Escape-fraction and spectrum tables/plots versus overlayer thickness and source depth. The **maximum useful active-layer thickness** for charged-particle detection.
- The recommended detector layout for C3 and for C1 (distances, areas, thicknesses, shielding, veto), with MDA per channel.
- A pre-registration template: observables, windows, cuts, stopping rules and the control schedule.
