# Brief M1 — Site physics, non-thermal energy sources, and per-configuration yield ladder

Branch: `claude/lucid-davinci-gel4vu-m1`. Read `_common.md` first. Read these literature digests in full:
- `docs/research/R3-low-energy-nuclear-physics.md`: screening data, threshold resonance, Dubey 2025, ARPA-E, NASA LCF
- `docs/research/R4-theories-geometric-predictions.md`

M0 established the governing facts:
- The rate is dominated by the rarest high-enhancement sites.
- The exponent gap from molecular D₂ to neutron-detectable is ~100 e-folds.
- Only two levers act on the exponent directly: effective screening, and the relative kinetic energy E of the pair.

Your job is to make both levers concrete for real site classes and real "cold" non-equilibrium processes. Then turn them into a **yield ladder** for each candidate configuration C1–C6 in `docs/01-master-plan.md`.

## Part A — Site catalogue (Pd–D; also Ti–D, Zr–D, Ni–H/D where data exist)
For each site class, compile from DFT/experimental literature (with citations):
- minimum D–D separation
- local electron density (a screening proxy)
- vibrational energy
- binding energy
- achievable number density (cm⁻³ or cm⁻²) and how geometry/processing controls it

Site classes:
- bulk octahedral and tetrahedral sites
- vacancy–D_n clusters (superabundant vacancies)
- divacancies/voids holding D₂ at high pressure
- dislocation cores
- grain boundaries
- free surface and subsurface sites
- crack faces/nanogaps (Storms' "nuclear active environment")
- oxide/metal interfaces (PdO/Pd, CaO/Pd, ZrO₂/Pd)
- nanoparticles of 2–10 nm

Estimate U_e,eff for each with at least two screening models:
- Thomas–Fermi/Lindhard with the local electron density
- an empirical scaling calibrated to the accelerator-measured U_e values in R3, mapped to thermal energies with the Yukawa-WKB method of `sim/m0_rate_budget.py` (import and reuse its functions)

Report the rate per pair and per cm³ of active material. State plainly where the numbers are guesses.

## Part B — Non-thermal ("cold-compatible") energy sources for the pair
Quantify each with equations and numbers. For each: energy distribution of deuterons produced, number per event and per cm² per cycle, and the in-flight thick-target D–D yield. Use the R3 screening-enhanced cross-sections and PdD stopping power. The processes:
1. **Fracto-emission / crack-tip charge separation** in PdD, TiD₂ and LiD during α/β cycling. Fields across nascent cracks (literature: 10⁷–10⁹ V/m), gap widths, D₂ pressure in the gap, mean free path, and the resulting deuteron energy spectrum. Compare predicted neutrons per cycle with published fracto-fusion claims (Derjaguin/Klyuev, Menlove/LANL Ti chips, Lipson) and with nulls.
2. **Desorption/phase-transition transients** (Lipson's Pd/PdO heterostructure proton claims).
3. **Electrochemical double-layer** fields: show the magnitude, which is expected to be negligible.
4. **Optical/plasmonic near fields at nanotips** acting on deuterons and electrons: show the magnitude.
5. **Cosmic-ray-driven baselines**, which are real and irreducible:
   - cosmic-neutron knock-on deuterons fusing in flight (NASA-LCF-like)
   - muon-catalysed fusion from stopped cosmic μ⁻ in PdD and D₂O

   Compute the expected rates in a 1 cm³ PdD cathode plus 50 mL D₂O at sea level. This is the "cold" floor any experiment sits on, and it tells us what the D-vs-H control must cancel.
6. **Near-threshold resonance** (Czerski et al.): if real, by what factor does it raise the rates above? What products and signatures follow: e⁺e⁻ pairs (511 keV), and the internal-pair / E0 channel?

## Part C — Yield ladder per configuration
For each configuration C1–C6, combine:
- the number of each site class that the configuration creates (use M3's crack/vacancy estimates if available on branch `claude/lucid-davinci-gel4vu-m3`, otherwise estimate)
- the non-equilibrium processes the configuration drives, and at what rate
- the fraction of products that reach detectors: use M5 if available on branch `claude/lucid-davinci-gel4vu-m5`, otherwise use M0's assumptions

Produce expected detected events per day under three scenarios:
- (i) standard physics
- (ii) "accelerator screening applies at the site", with R3 values
- (iii) a hypothetical anomaly: an extra ×10ⁿ enhancement localised at site class k, for each k

The key output is a matrix [configuration × site class] of **sensitivity to a localised anomaly**: the minimal enhancement factor at site class k that configuration C would detect at 5σ in 30 days. This matrix is the main quantitative input to the lead's decision.

## Outputs
Tables, figures, the matrix above, and design recommendations: which site classes to maximise, by what processing and geometry, and which drive protocol maximises non-thermal energies without violating "cold".
