# M5 — Detection geometry, backgrounds and statistics

Code: [`sim/m5_stopping.py`](../../sim/m5_stopping.py), [`sim/m5_escape.py`](../../sim/m5_escape.py), [`sim/m5_si_telescope.py`](../../sim/m5_si_telescope.py), [`sim/m5_cr39.py`](../../sim/m5_cr39.py), [`sim/m5_neutron.py`](../../sim/m5_neutron.py), [`sim/m5_stats.py`](../../sim/m5_stats.py) · Raw outputs: `figs/m5_*.txt` · Regenerate everything with `pip install numpy scipy matplotlib pycatima`, then run each script (`python3 sim/m5_stats.py` also re-runs the telescope model). The scripts use fixed seeds. The Monte Carlo (MC) numbers carry about ±3 % statistical noise for charged particles and about ±10 % for neutrons.

## Summary for lead

1. **Silicon ΔE–E telescope on C3 is the instrument. Everything else is a cross-check.** Baseline: Ø20 mm membrane; 25 µm/600 mm² ΔE plus 500 µm E, 4 mm away, ≤10⁻⁴ mbar. It detects **0.22 of emitted 3.02 MeV protons** (0.11 per fusion). Background in the proton window is **0.01–0.03/day**, set by cosmic-neutron (n,p) reactions inside the silicon (uncertain ×3). A single Si detector without particle ID sees 3.7/day.
2. **30-day 5σ reach:** **7.6×10⁻⁵ fusions/s** against an equal-time H control, **3×10⁻⁵** with background characterised over ≥10× the live time. That is 20–90× better than M0 assumed.
3. **³He bank** (18 tubes, 4–5 cm HDPE in front, Cd plus 20 cm borated shield): ε = 0.23, background ≈0.35 cps (M0 assumed 0.05). Reach is 2×10⁻² fusions/s, **≈300× worse than Si**. EJ-309 is 3× worse again; bubble detectors are useless for physics. Barometric drift is 6σ/hPa over 30 days. **Pre-register that protons below ~10⁻² fusions/s come without detectable neutrons.**
4. **Useful depth:** 3 MeV protons reach the 2.6–3.1 MeV peak only from the **top 5 µm**, and PID-window counting only from the top **16–23 µm**. Tritons come from ≤3.7 µm, ³He from ≤0.9 µm. Front overlayers must be ≤100 nm. 1 µm Ni/Cu kills ³He.
5. **Statistics trap:** at near-zero background, an equal-time control needs **≥18–22 signal counts** for 5σ (binomial limit 0.5ⁿ). Budget background time at ≥10× the live time, which cuts that to ~8. Current modulation needs T_mod ≥ 10 τ_response and is ~5× slower, so it is the confirmatory test.
6. **CR-39** channel C (6 µm Mylar + 55 ± 2 µm Al) is α-free by construction. Reach is ~1.5×10⁻³ fusions/s in C1, but a 30 µm electrolyte film kills it. It is auxiliary only.

## 1. Questions answered

| # | Question (brief) | Where |
|---|---|---|
| Q1 | Escape fraction and spectra vs birth depth, overlayers and materials | §5.1–5.3, figs `m5_escape_*` |
| Q2 | C3 Si spectrometer and ΔE–E telescope geometry, background in 2.6–3.1 MeV, minimum detectable rate | §5.4–5.6 |
| Q3 | CR-39 stack (Mylar/Al) for p/t/α discrimination; controls | §5.7 |
| Q4 | Best 2.45 MeV neutron system; achievable sea-level background | §5.8 |
| Q5 | Statistics: Currie, Feldman–Cousins, Bayes, modulation, blind analysis, look-elsewhere, run time to 5σ | §5.9, §9 |

## 2. Model and equations

**Stopping powers.** The electronic stopping cross-section per atom, S_e(Z₁,Z₂,T/A) [eV/(10¹⁵ at/cm²)], comes from the Ziegler (SRIM-type) elemental parameterisation. The implementation used is catima/pycatima 1.98 (<https://pypi.org/project/pycatima/>). Compounds use Bragg's rule, S = Σ nᵢ Sᵢ. Isotopes (d, t, ³He) use velocity scaling. Nuclear stopping uses the ZBL universal formula. Mass stopping is converted with ρ. The CSDA range is R(E) = ∫ dE / S(E), and residual energy is E_out = R⁻¹(R(E₀) − L).

As an independent check, I implemented Bethe–Bloch from first principles:

S = K z² (Z/A) β⁻² [ln(2mc²β²γ²/I) − β² − C/Z + L₂]

with K = 0.307075 MeV cm²/mol, the Barkas–Berger shell correction C, and the Bloch term L₂ = −y² Σ 1/[n(n²+y²)], where y = zα/β.

Energy straggling uses Bohr's formula, σ² = 4π z² e⁴ n_e L. This is an upper bound at low velocity.

**Escape and transport.** Emission is isotropic from depth z, and particles follow straight CSDA paths. A particle leaving at direction cosine μ traverses z/μ and then tᵢ/μ in each overlayer. It then crosses the vacuum gap, the detector dead layer (50 nm Si) and the ΔE and E layers, with Gaussian resolution added. The analytic limit for a δ source with no threshold is f(z) = ½(1 − z/R). For a uniform source over [0, t] it is f = ½(1 − t/2R) when t ≤ R, and R/(4t) when t > R.

For a uniform reaction density, the signal per unit area is ∝ ∫₀ᵗ f dz. This saturates at t = R_eff = R(E₀) − R(E_min) and reaches 90 % at 0.684 R_eff.

**Particle identification (PID).** An event is a proton if both ΔE and E exceed 100 keV and 0.75 ≤ ΔE/ΔE_p(E_tot, normal) ≤ 2.0. The window is E_tot ∈ [1.41, 3.10] MeV (the PID window) or [2.6, 3.1] MeV (the peak window).

**Neutron transport.** A vectorised MC tracks neutrons through elastic scattering, isotropic in the CM frame, on H (Gammel σ_np), C and O (coarse ENDF digitisation, ±15 %). Below 0.2 eV a one-group thermal model applies. It uses an isotropic transport cross-section per H of 28.6 b, calibrated to the water thermal diffusion coefficient D = 0.16 cm, and 1/v absorption averaged over a Maxwellian (×0.886). ³He absorption uses σ(0.0253 eV) = 5333 b × (1/v).

The geometry is a cylindrical HDPE annulus around a void cavity containing the cell, with 1-inch tubes on a ring. A Cd sheet is optional, as is an outer borated shield that absorbs thermal neutrons.

The external cosmic background is simulated as an isotropic inward flux: current φ/4 per cm² over the outer surface, with ×1.5 added for the unmodelled end faces.

**Statistics.**

- Currie: L_C = 1.645√(2B), L_D = 2.71 + 4.65√B (paired blank).
- Known-background discovery: n_crit is the smallest n with P(N ≥ n | B) ≤ 2.87×10⁻⁷, and s solves P(N ≥ n_crit | B + s) = 0.5.
- On/off with t_off = τ t_on: the Li & Ma (1983) eq. 17 Asimov median.
- Feldman–Cousins: the likelihood-ratio-ordered Poisson belt.
- Bayes: flat prior on s ≥ 0, with 1 − CL = P(N ≤ n | b + s)/P(N ≤ n | b).
- Modulation, for a first-order response lag τ to a 50 % square wave of period T: c = ⟨on⟩ − ⟨off⟩ = 1 − (4τ/T) tanh(T/4τ). Derived analytically and checked against the ODE: 0.3966 vs 0.3965.
- Look-elsewhere: Šidák, p_local = 1 − (1 − p_global)^(1/N).

## 3. Parameters

| Parameter | Value | Units | Source | Uncertainty |
|---|---|---|---|---|
| Electronic stopping (validation reference) | NIST PSTAR (p, ≤2 MeV), ASTAR (α, ≤10 MeV) | MeV cm²/g | NIST tables as transcribed in Geant4 `G4PSTARStopping.cc`/`G4ASTARStopping.cc` (<https://github.com/Geant4/geant4>, originals at <https://physics.nist.gov/PhysRefData/Star/Text/>) | 2–5 % (NIST) |
| Electronic stopping (engine) | Ziegler/SRIM-type coefficients | eV/(10¹⁵ at/cm²) | catima (<https://pypi.org/project/pycatima/>) | ±5 % vs NIST above 0.5 MeV/u (§4) |
| Mean excitation energies I | Pd 470, Ag 470, Ni 311, Cu 322, Au 790, Al 166, Si 173, O 95, Ca 191, H 19.2 | eV | ICRU 37/49 via <https://physics.nist.gov/PhysRefData/Star/Text/method.html> | ±2–5 % |
| ρ(PdD₀.₉) | 10.90 | g/cm³ | 4 f.u./(N_A a³), with a = 4.04 Å (β-phase a = 4.02–4.025 Å at onset: <https://inis.iaea.org/records/9hbv3-2bz85>) | ±2 % |
| ρ Pd, PdO, Au, Ni, Cu, CaO, Al, Si, Mylar, D₂O | 12.02, 8.3, 19.32, 8.908, 8.96, 3.34, 2.699, 2.329, 1.397, 1.104 | g/cm³ | CRC Handbook; NIST | <1 % |
| CR-39 | C₁₂H₁₈O₇, 1.31 | g/cm³ | <https://www.tasl.co.uk/tastrak-padc.php> | — |
| Si detector resolution | 25 keV (E), 50 keV (ΔE) | FWHM | ORTEC ULTRA class; ΔE from preamp noise slope for a ~0.6 nF quadrant (assumption) | ×1.5 |
| Si entrance window | 50 nm Si-equivalent | — | ORTEC ULTRA spec (<https://www.ortec-online.com/products/radiation-detectors/silicon-charged-particle-radiation-detectors/si-charged-particle-radiation-detectors-for-alpha-spectroscopy/ultra>) | ±50 % (irrelevant for p) |
| Detector-intrinsic α background | ≤24 counts/day, 3–8 MeV, 450 mm² | — | ULTRA-AS warranty (brochure <https://www.ortec-online.com/-/media/ametekortec/brochures/u/ultra_ultra-as-a4.pdf>, via search snippet) | upper limit |
| Thin ΔE availability | 20 µm, ±0.5 % rms uniformity | — | Micron Semiconductor, <https://www.aanda.org/articles/aa/full_html/2021/06/aa39754-20.html> | large-area 25 µm must be quoted by vendor |
| α → p misID floor (edge and partial-charge events) | 10⁻⁴ | — | **assumption**; the MC Gaussian tails give 0 in 10⁵ | ×10 |
| Muon flux (horizontal) | 1 | cm⁻² min⁻¹ | PDG Cosmic-ray review <https://pdg.lbl.gov/2022/reviews/rpp2022-rev-cosmic-rays.pdf> | ±10 % |
| Cosmic neutron flux, >10 MeV (NYC) | 3.6×10⁻³ | cm⁻² s⁻¹ | Gordon et al., IEEE TNS 51, 3427 (2004), quoted in <https://arxiv.org/pdf/2605.16534> | ±10 % |
| Cosmic neutron flux, 1–10 MeV / 0.4 eV–1 MeV / thermal | 9 / 10 / 4 | cm⁻² h⁻¹ | Gordon 2004 spectrum shape (my coarse integration) | ±50 % |
| Neutron barometric coefficient | 0.66–0.82 (use 0.72) | %/hPa | <https://doi.org/10.3390/s26030925>; NMDB | ±10 % |
| Solar modulation of neutron flux | 10–20 % over a cycle; Forbush decreases up to ~10–30 % over days | — | same | — |
| σ(³He(n,p)) at 0.0253 eV | 5333 | b | ENDF/B-VIII | 0.2 % |
| σ_np elastic | Gammel formula | b | e.g. Marion & Young, *Nuclear Reaction Analysis* (1968) | ~2 % below 20 MeV |
| ²⁸Si(n,p) averaged over >5 MeV | 0.25 | b | ENDF/B-VIII shape (assumption for the average) | ±50 % |
| Water thermal D, L; Fermi age | 0.16 cm, 2.85 cm; 27 cm² | — | Lamarsh, *Introduction to Nuclear Reactor Theory*, Tables 5.2/5.3 | — |
| Optimum HDPE in front of ³He, bare source | 4–5 cm | — | Rees & Czirr, NIM A 691 (2012) (<https://www.osti.gov/pages/biblio/1178530>) | — |
| EJ-309 composition | H 5.43×10²², C 4.35×10²² | cm⁻³ | Eljen data sheet <https://eljentechnology.com/products/liquid-scintillators/ej-301-ej-309> | 1 % |
| EJ-309 proton light | L = 0.817E − 2.63(1 − e^(−0.297E^0.9)) | MeVee | Enqvist et al., NIM A 715, 79 (2013) | ±10 % |
| EJ-309 γ singles / PSD misID at 100 keVee | 150 cps / 10⁻³ | — | **assumption** (typical) | ×3 |
| BD-PND sensitivity | 1.9–3.7 (use 3) | bubbles/µSv | <https://www.osti.gov/etdeweb/servlets/purl/20409858> (search snippet) | — |
| h*(10) at 2.5 MeV | 416 | pSv cm² | ICRP 74 | 5 % |
| Sea-level neutron H*(10) | ~9 | nSv/h | UNSCEAR 2000 (annual neutron component ~0.08 mSv) | ±50 % |
| CR-39 registration threshold | LET ≥ 12 keV/µm (p ≲ 4 MeV at 6.25 N NaOH, 70 °C, 6 h); critical angle 50° (p), 70° (α) | — | **assumption**, consistent with <https://www.osti.gov/servlets/purl/1076447> and Sci. Rep. 7, 2152 (2017) (critical energy 21–22.5 MeV at 16–24 h etch) | ±30 % |
| CR-39 background | ~30 tracks/cm² per 28 d (radon-dominated) | — | <https://arxiv.org/pdf/1806.06567> | site-dependent |
| 6 µm Mylar cut-offs | p 0.45, t 0.55, ³He 1.40, α 1.45 | MeV | Mosier-Boss et al., EPJ AP 46, 30901 (2009) | — |
| Lipson Au/Pd/PdO:D proton emission | (4.0 ± 1.0)×10⁻³ | p/s (4π) | <https://www.osti.gov/etdeweb/biblio/20845769> | as reported |
| Kasagi TiDₓ ΔE–E products | p to ~17 MeV, α to ~6.5 MeV, t and ³He at 4.75 MeV | — | <https://www.osti.gov/etdeweb/biblio/109126>; JCF3 abstracts | as reported |

## 4. Verification

**Stopping versus NIST** (`m5_stopping.txt`, `figs/m5_stopping_validation.png`):

- Protons at 0.5–2 MeV: within **±4.5 %** for Al, Si, Cu, Ag, Au, Pt, Mylar, water and air.
- Protons at 0.1–0.3 MeV: Si +9 to +12 %, Mylar −10 %, air −17 %. This is a known SRIM–PSTAR difference near the stopping maximum. It affects only the last ~1 µm of range.
- α particles at 2–10 MeV: within **±5 %** for all materials, except Mylar and water at 8–10 MeV (+4.4 to +5.1 %; Bragg's rule with gas-phase H/O).
- **Pd** has no NIST table. Against Ag × (Z/A) scaling (same I = 470 eV), Pd comes out +4 to +7 %. I therefore assign **±6 %** to all Pd/PdD ranges.
- The independent **Bethe–Bloch** check agrees with the engine to **≤4 %** at ≥10 MeV/u for Si, Al, Cu, Pd and Au. Below ~5 MeV/u the Barkas–Berger shell formula fails (η < 0.1), which is why low-energy values rely on the empirical tables.
- Literature cross-checks: 6 µm Mylar stops p < 0.46, t < 0.54, ³He < 1.42 and α < 1.44 MeV (published: 0.45/0.55/1.40/1.45). The CSDA range of 3.02 MeV p in Si is 91.5 µm (PSTAR at 3 MeV: 92.7 µm). The 5.486 MeV α range is 27.3 µm in Si (SRIM ~28) and 4.25 cm in air (4.1–4.2).

**Escape MC versus analytic** (`m5_escape.txt`): δ sources at 0, 0.25R, 0.5R and 0.9R reproduce ½(1 − z/R) to ≤0.002 absolute. The only exception is ³He at 0.9R, where the 1 keV floor gives 0.037 vs 0.050. Uniform sources over 0.5R, R and 3R match to ≤0.001. Disk-geometry MC matches the on-axis solid-angle formula (0.3077 vs 0.3077). The straight-line approximation overestimates depths by ≤5 % (the CSDA/projected-range detour factor).

**Neutron MC** (`m5_neutron.txt`):

- Water thermal diffusion length: 2.86–2.88 cm (Lamarsh 2.85; D was the calibration input, so this checks the absorption).
- HDPE thermal diffusion length: 2.33 cm, a prediction (literature ~2.1–2.3).
- Water Fermi age to 1.46 eV: 25.3 cm² for 2.0 MeV and 29.9 cm² for 2.45 MeV (fission-spectrum reference 27 cm²).
- The optimum front moderator comes out at 5 cm (ε = 0.237), with a plateau from 4 to 6 cm. Rees & Czirr's published optimum for bare sources is 4–5 cm.

**Statistics:**

- FC intervals match Feldman & Cousins Table IV (e.g. b = 0.5: [0, 1.94], [0, 3.86], [0.61, 6.92]; b = 1, n = 0: 1.61).
- The exception is n = 0 at b ≥ 2, where my belt gives 1.08 (b = 2) and 0.95 (b = 3) against published 1.26 and 1.08. The difference comes from how the table resolves ties in the ordering. **For quoted limits, use the published FC table values.**
- The Poisson n_crit agrees with direct summation.

## 5. Results

### 5.1 Ranges (CSDA, µm; `m5_stopping.txt`)

| product | Pd | PdD₀.₉ | PdO | Au | Ni | Cu | CaO | D₂O | Mylar | Al | Si | CR-39 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| p 3.02 | 31.3 | 33.3 | 39.8 | 28.8 | 32.8 | 34.9 | 63.2 | 149 | 114 | 81.9 | 91.5 | 118 |
| t 1.01 | 4.72 | 4.82 | 5.67 | 5.45 | 5.20 | 5.90 | 7.86 | 15.9 | 12.0 | 10.2 | 10.2 | 12.3 |
| ³He 0.82 | 1.46 | 1.47 | 1.74 | 1.80 | 1.63 | 1.88 | 2.40 | 4.73 | 3.49 | 2.97 | 2.98 | 3.57 |
| α 5.30 (²¹⁰Po) | 9.8 | 10.2 | 12.1 | 10.0 | 10.3 | 11.3 | 18.3 | 40.6 | 31.3 | 23.7 | 25.9 | 32.0 |
| α 8.78 (²¹²Po) | 19.7 | 20.9 | 24.8 | 19.0 | 20.7 | 22.3 | 38.7 | 89.1 | 68.5 | 50.0 | 55.6 | 70.4 |
| α 12 | 31.4 | 33.4 | 39.9 | 29.1 | 33.0 | 35.2 | 63.3 | 149 | 114 | 82.0 | 91.4 | 117 |
| t 4.75 / ³He 4.75 | 35.2 / 9.3 | 37.1 / 9.8 | — | — | — | — | — | — | — | — | 96.9 / 25.3 | — |
| p 14.7 | 412 | 444 | — | — | — | — | — | 2387 | 1833 | 1241 | 1353 | 1893 |

In D₂ gas at 1 bar the ranges are 726 mm (p), 60 mm (t) and 17 mm (³He). Stopping power at birth in PdD₀.₉ is 57.5 keV/µm (p), 189 keV/µm (t) and 719 keV/µm (³He).

### 5.2 Escape fraction vs depth, overlayer and material (fraction of 4π; `m5_escape.txt`)

Source in PdD₀.₉ under 20 nm PdO:

| depth (µm) | p 3.02 (E > 50 keV) | t 1.01 | ³He 0.82 | p in 2.6–3.1 MeV |
|---|---|---|---|---|
| 0 | 0.500 | 0.498 | 0.493 | 0.499 |
| 1 | 0.484 | 0.383 | 0.069 | 0.427 |
| 2 | 0.470 | 0.269 | 0 | 0.355 |
| 5 | 0.423 | 0 | 0 | 0.139 |
| 10 | 0.348 | 0 | 0 | 0 |
| 20 | 0.196 | 0 | 0 | 0 |
| 30 | 0.045 | 0 | 0 | 0 |

Effect of overlayers on a surface source (p / t / ³He / p-window):

- **PdO 100 nm:** 0.499 / 0.490 / 0.465 / 0.494
- **Au 50 nm:** 0.499 / 0.495 / 0.481 / 0.495
- **CaO 100 nm:** 0.499 / 0.493 / 0.474 / 0.497
- **Ni 1 µm:** 0.485 / 0.391 / **0.110** / 0.427
- **Cu 1 µm:** 0.486 / 0.404 / **0.157** / 0.430
- **10 µm D₂O film:** 0.466 / 0.155 / **0** / 0.350
- **6 µm Mylar:** 0.473 / 0.232 / **0** / 0.382

Figures: `figs/m5_escape_vs_overlayer.png` (all 8 overlayer materials, 4 products) and `figs/m5_escape_vs_depth.png`.

### 5.3 Maximum useful active-layer thickness (PdD₀.₉, uniform reaction density)

| observable | R_eff (µm) | 90 % of the saturated signal at |
|---|---|---|
| p in the 2.6–3.1 MeV peak | — | **5.0 µm** |
| p in the telescope PID window (1.41–3.1 MeV) | 22.9 | **16 µm** |
| p, any energy > 0.3 MeV (single Si) | 31.9 | 22 µm |
| t > 0.2 MeV | 3.7 | 2.5 µm |
| ³He > 0.2 MeV | 0.9 | 0.6 µm |
| α 5.3 / α 12 (> 0.5 MeV) | 9.1 / 32.3 | 6.2 / 22 µm |
| p 14.7 (D–³He secondary) | 438 | 300 µm |

For bulk-distributed sites, PID-window counting collects **3.2×** more protons than peak-window counting.

Detected spectra through the baseline geometry are in `figs/m5_escape_spectra.png`. For a surface δ source the proton peak sits at 3.017 MeV with the full ~25 keV resolution. For a U(0–10 µm) source the mean is 2.51 MeV, with 48 % in 2.6–3.1 MeV. For U(0–100 µm) the mean is 1.84 MeV, with 21 % in the peak window. Tritons and ³He from U(0–1 µm) have mean energies of 0.84 and 0.39 MeV.

### 5.4 C3 telescope geometry (`m5_si.txt` §1, `figs/m5_si_geometry.png`)

The table gives the detected fraction of emitted p (4π, PID-accepted) for a U(0–5 µm) source:

| membrane a (mm) | ΔE area (mm²) | d = 2 | 3 | 4 | 6 | 8 | 12 | 20 mm |
|---|---|---|---|---|---|---|---|---|
| 10 | 300 | 0.168 | 0.161 | 0.152 | 0.133 | 0.111 | 0.076 | 0.037 |
| 10 | 450 | 0.219 | 0.208 | 0.194 | 0.172 | 0.152 | 0.109 | 0.055 |
| 10 | **600** | 0.227 | 0.225 | **0.216** | 0.200 | 0.177 | 0.133 | 0.071 |
| 10 | 900 | 0.228 | 0.228 | 0.226 | 0.222 | 0.208 | 0.172 | 0.099 |
| 5 | 450 | 0.226 | 0.227 | 0.228 | 0.216 | 0.187 | 0.125 | 0.059 |
| 15 | 600 | 0.165 | 0.162 | 0.157 | 0.145 | 0.132 | 0.107 | 0.063 |
| 15 | 1200 | 0.226 | 0.226 | 0.226 | 0.214 | 0.199 | 0.171 | 0.114 |

These values use the final PID band (0.75–2.0). The saturation value, ≈0.227, is the escape-limited fraction times the PID acceptance.

**Rule:** the detector radius should be at least the membrane radius plus the gap. Beyond that, extra area buys nothing, but background scales with area.

Scaling the membrane up while keeping the detector fixed loses efficiency in direct proportion. Use **several telescopes (a 2×2 array) for membranes larger than ~Ø20 mm**.

ΔE thickness trade-off (from `m5_si.txt` §2):

| ΔE thickness | lowest identifiable p | largest α stopped in ΔE | ΔE deposited by a 3.02 MeV p |
|---|---|---|---|
| 10 µm | 0.83 MeV | 2.62 MeV | 206 keV |
| 15 µm | 1.05 MeV | — | 312 keV |
| 20 µm | 1.24 MeV | — | 422 keV |
| **25 µm** | **1.41 MeV** | 5.17 MeV | **535 keV** |
| 40 µm | 1.87 MeV | — | 898 keV |
| 65 µm | 2.50 MeV | 9.70 MeV (all natural α) | — |

I chose 25 µm. It gives a 3 MeV proton ΔE signal of ≥10× the ΔE noise while keeping the proton window 1.41–3.1 MeV wide. Tritons (1.01 MeV) and ³He (0.82 MeV) stop in the ΔE, and are measured there alone ("ΔE-only" events: ΔE fired, E below threshold).

### 5.5 Particle identification (`figs/m5_si_pid.png`)

The MC used 40 000–200 000 events per species. Fraction passing the proton cut in the PID window:

- 3.02 MeV p: 0.84
- ²¹⁰Po α (surface), 8.78 MeV α (bulk), 12 MeV α, 1.01 MeV t, 0.82 MeV ³He, 4.75 MeV t and 4.75 MeV ³He: **0**
- 3 MeV d: 0.15

With the band cap at 3.0, 4.75 MeV tritons leak at 8 % and deuterons at 32 %. **That is why the upper edge was lowered to 2.0**, at a cost of 20 % in efficiency (0.267 → 0.217).

No angular collimator is needed. The α/p separation holds at every angle. A physical angle-limiting collimator at 50° would cost 20 % of the signal (0.174 vs 0.217).

**Signal efficiency** (detected/emitted p) by source distribution, peak / PID window:

| source | peak window | PID window |
|---|---|---|
| surface | 0.220 | 0.220 |
| U(0–5 µm) | 0.206 | 0.217 |
| U(0–10 µm) | 0.117 | 0.210 |
| U(0–20 µm) | 0.059 | 0.168 |
| U(0–33 µm) | 0.035 | 0.103 |
| U(0–100 µm) | 0.012 | 0.034 |

The triton and ³He ΔE-only windows each detect about 0.31 of those emitted from the surface.

### 5.6 Background budget, C3 telescope (counts/day, sea level; `m5_si.txt` §5)

| source | single Si, 2.6–3.1 MeV, no PID | telescope, 2.6–3.1 MeV | telescope, PID window |
|---|---|---|---|
| cosmic muons, no veto (Landau with δ-ray escape) | 0.36 | 0 | 0 |
| detector-intrinsic α (spec maximum; misID 10⁻⁴) | 3.2 | 0.001 | 0.003 |
| membrane bulk U + Th at 1 ppb each | 0.007 | ~0 | ~0 |
| membrane surface ²¹⁰Po at 1 /cm²/day | — | 3×10⁻⁵ | 1×10⁻⁴ |
| n–p recoils from residual H (H/D = 1 %) | 3×10⁻⁵ | 2×10⁻⁵ | 6×10⁻⁵ |
| n–d recoils misidentified as p | 0.002 | 0.0005 | 0.0012 |
| **²⁸Si(n,p) inside ΔE and E** | 0.086 | **0.009** | **0.023** |
| **Total** | **3.66** | **0.010** | **0.028** |
| with a 99 % muon veto and 10 cm HDPE (fast n × 0.6) | 3.26 | 0.007 | 0.018 |

**D-specific versus H-specific real protons.** In a PdH control, n–p recoils give **0.006/day of real protons** in the PID window. The D cell has no equivalent. This makes the H control slightly *more* background-rich, which is conservative, but it must be modelled in the D−H subtraction.

**Beam–target D–D induced by recoil deuterons** is negligible: of order 10⁻¹¹ s⁻¹.

**The residual background is identical in D, H and off states** (Si(n,p) is intrinsic to the silicon). It costs statistics only; it cannot fake a D-specific signal.

**The dominant artifact risk is electronic.** Current modulation couples into the detector electronics as EMI, and the membrane heats with current, which changes Si leakage current and noise. Schenkel et al. (<https://arxiv.org/abs/1905.03400>) had Si diodes swamped by EMI during plasma pulses. These are handled by ΔE–E time coincidence, digitised pulse shapes and a **covered blank telescope** (§6, recommendation 6), not by the background model.

**Vacuum** (`m5_si.txt` §6): residual gas loss over a 10 mm path at 1 mbar is 0.12 keV (p), 0.53 keV (t) and 2.1 keV (³He); at 10⁻² mbar it is negligible. The D₂ gas load from permeation through a 3.14 cm² membrane is 6×10⁻⁵ to 6×10⁻³ mbar L/s for fluxes of 10¹⁵–10¹⁷ D/cm²/s. With 50 L/s of D₂ pumping speed that gives 10⁻⁶ to 10⁻⁴ mbar. The binding requirement is not energy loss but safe biasing: never bias the detectors through the Paschen region (~10⁻²–10 mbar).

### 5.7 CR-39 stack (`m5_cr39.txt`, `figs/m5_cr39_filters.png`)

All four channels sit behind 6 µm Mylar in electrolyte (C1/C2); in vacuum (C3) there is no barrier. Residual energy at the CR-39 in MeV, normal incidence / 45°, with `*` marking a registered track:

| species | A: barrier only | B: +12 µm Al | C: +55 µm Al | D: +100 µm Al |
|---|---|---|---|---|
| p 3.02 | 2.92/2.88* | 2.64/2.47* | 1.35/0* | 0 |
| t 1.01 | 0.54/0.30* | 0 | 0 | 0 |
| ³He 0.82 | 0 (stopped by the Mylar; visible only without a barrier) | 0 | 0 | 0 |
| ²²²Rn–²¹²Po α (5.5–8.8) | yes* | yes* | **0** | 0 |
| α 12 | 11.6* | 10.5* | 4.3* (0 at 45°) | 0 |
| α 16 | yes* | yes* | yes* | 5.8* |
| p 14.7 | not registered at a 6 h etch (LET 4 keV/µm); use a 16–24 h etch or Si | | | |

**Identification logic:**

- 3.02 MeV p gives the signature **ABC–**.
- Tritons give **A–––** (³He only without a barrier). t and ³He are then separated by pit size: ³He LET ≈ 4× t.
- All natural α give **AB––**, so **channel C is α-free by construction**.
- Energetic α above 14 MeV give ABCD.
- **Channel D is the blind control for C:** a real 3 MeV proton line must appear in C and vanish in D.

**Tolerance on the C filter.** R(8.785 MeV α, Al) is 50.1 µm, but after the 6 µm Mylar the α needs only ~46 µm of Al, so 55 µm has 9 µm of margin. At the other end, the proton must still reach the CR-39 with E ≥ 0.5 MeV; at 70 µm of Al it arrives at only 0.64 MeV. Hence **55 ± 2 µm Al**, or 84 µm of all-Mylar, which stops α ≤ 9.96 MeV.

**Efficiency** (channel C per emitted p, U(0–5 µm) source): 0.097 with no electrolyte film, 0.065 behind a 10 µm film, **0.010 behind 30 µm**, and 0 behind 100 µm. Multiply by the fraction of cathode area the CR-39 covers.

### 5.8 Neutron detection (`m5_neutron.txt`, `figs/m5_neutron_he3.png`)

**³He bank: efficiency vs HDPE thickness in front of the tubes** (cavity r = 8 cm, 18 × 1-inch × 40 cm tubes at 4 atm):

| front HDPE | 1 | 2 | 3 | 4 | 5 | 6 | 8 | 10 cm |
|---|---|---|---|---|---|---|---|---|
| ε | 0.14 | 0.15 | 0.20 | 0.22 | 0.24 | 0.23 | 0.21 | 0.18 |

- Back moderator 2 / 4 / 6 / 10 cm: ε = 0.11 / 0.18 / 0.22 / 0.23.
- 10 atm tubes: 0.24. Tube pitch 3 cm: 0.23.
- Cavity radius 5 / 8 / 12 / 16 / 20 cm: ε = 0.26 / 0.24 / 0.21 / 0.17 / 0.15.

**Background.** Unshielded, the bank sees ~1.8 cps. With 1 mm Cd plus a borated-HDPE shield of 10 / 20 / 30 cm it sees 0.60 / **0.34** / 0.18 cps of cosmic neutrons. Adding ~10 % muon-induced and 0.005 cps intrinsic gives **≈0.35–0.39 cps** at 20 cm.

This is 7× the 0.05 cps M0 assumed. It is also uncertain ±50 %, because of the cosmic spectrum and because >20 MeV spallation is not modelled. Treat 0.35 cps as a model estimate. **Measure the actual background in the lab before the design is frozen.**

**Barometric systematic.** Over 30 days the bank accumulates ~8×10⁵ background counts. A 1 hPa change moves the rate by 0.72 %, which is **6σ**. Correct for pressure with a logged β fitted in situ, and rely on on/off and D/H differencing.

**EJ-309** (two 5″×5″ cells at 10 cm): geometric acceptance 0.078 per cell, intrinsic efficiency above 0.1 MeVee 0.71, so ε = 0.11 per neutron. Background in the 0.1–0.8 MeVee neutron band is 0.19 (cosmic n) + 0.15 (γ leakage) = 0.34 cps per cell. Its 30-day reach is 3× worse than the ³He bank. What it offers is energy information (the 2.45 MeV recoil edge at 0.73 MeVee) and fast timing. It needs a muon veto.

**Bubble detectors** (BD-PND at 3 cm): 1.1×10⁻⁵ bubbles per emitted neutron, background ≈0.65 bubbles/day. A 30-day 5σ result needs **~2 fusions/s**. Use them as γ-blind dosimetry and safety cross-checks only.

### 5.9 Minimum detectable rates and run time to 5σ (`m5_stats.txt`, `figs/m5_stats_runtime.png`)

All values are median 5σ reach in D–D fusions/s. A proton is emitted in 50 % of fusions and a neutron in 50 %.

| channel (per fusion ε, background) | 1 day | 14 days | 30 days |
|---|---|---|---|
| **C3 Si telescope, PID window** (0.109, 0.028/day), equal-time control | 1.9×10⁻³ | 1.5×10⁻⁴ | **7.6×10⁻⁵** |
| same, background known from 10× live time | 6×10⁻⁴ | 6×10⁻⁵ | **3.3×10⁻⁵** |
| C3 single Si, no PID (0.103, 3.7/day) | 3.1×10⁻³ | 5.1×10⁻⁴ | 3.3×10⁻⁴ |
| CR-39 channel C, C1 geometry (0.0097 with 30 % coverage and a 10 µm film; 3 tracks/cm² per run) | 2.4×10⁻² | 2.6×10⁻³ | 1.5×10⁻³ |
| ³He bank (0.11, 0.31–0.39 cps) | 0.12 | 3.2×10⁻² | **2.2×10⁻²** |
| EJ-309 pair (0.056, 0.67 cps) | 0.36 | 9.5×10⁻² | 6.5×10⁻² |
| BD-PND pair | — | — | ~2 |

**Run time to 5σ for the C3 telescope** against an equal-time control:

| rate (fusions/s) | 10⁻⁴ | 10⁻³ | 10⁻² | 10⁻¹ | 1 |
|---|---|---|---|---|---|
| days | 22 | 2.0 | 0.19 | 0.019 | 0.002 |

For comparison, the ³He bank needs 1.4×10⁶ days at 10⁻⁴, 140 days at 10⁻², 1.4 days at 10⁻¹ and 0.02 days at 1 fusion/s.

**Signal counts needed for 5σ** at 50 % power:

| B (counts) | 0 | 0.1 | 1 | 10 | 100 |
|---|---|---|---|---|---|
| known background | 1.7 | 4.6 | 8.7 | 19.7 | 54.7 |
| equal-time on/off | — | 18.7 | 22.0 | 35.9 | 83.6 |

**Current modulation.** The contrast is c = 0.32, 0.61 and 0.80 at T/τ = 5, 10 and 20.

With T = 1 h and τ = 5 min, and a signal present only while the drive is on, the C3 telescope needs **11 days at 10⁻³ fusions/s** (versus 2 days for D vs H). The modulation test costs about ×5 in time. It is nonetheless the strongest defence against slow drifts and against cell-specific artifacts.

## 6. Design recommendations for the lead

1. **Primary H1 channel: a Si ΔE–E telescope in vacuum facing the C3 membrane.**
   - ΔE: **25 ± 2 µm**, **600 mm²** (Ø27.6 mm), segmented into **4 quadrants** to keep each at ≈0.6 nF and ~50 keV FWHM.
   - E: **500 µm**, ≥700 mm² (Ø30 mm), ion-implanted with ≤50 nm windows, ULTRA-AS class (certified low-α).
   - Gap between ΔE and E: **1.5 mm**.
   - Membrane-to-ΔE distance: **4 ± 1 mm**. Efficiency falls about 4 % per extra mm (0.227 at 2 mm, 0.200 at 6 mm), so go as close as membrane bowing and thermal load allow.
   - PID band: **0.75–2.0 × ΔE_p(E_tot)**.
   - Proton windows: **1.41–3.10 MeV (PID)** and **2.60–3.10 MeV (peak)**.
   - Expected performance: ε = 0.22 per emitted p, background 0.01–0.03/day, 30-day reach 3–8×10⁻⁵ fusions/s.
2. **Membrane active area Ø20 mm per telescope.** For more area, tile telescopes (2×2 for Ø40 mm). Do not use a larger single detector at the same distance: efficiency drops ∝ 1/area once the membrane is larger than the detector.
3. **Collimation.** Use no angle-limiting collimator. Use an aperture mask of **low-α electroformed Cu or Si, 1 mm inside the ΔE active edge**, to keep events off the edges (partial charge collection causes misID). Wetted and in-view seals must be metal or PTFE, not elastomers, because Viton and other elastomers contain U/Th.
4. **Vacuum:** ≤10⁻⁴ mbar, with a turbo pump of ≥50 L/s for D₂ and a bias interlock that allows bias only below 10⁻⁴ mbar. Log pressure as a covariate, because it tracks the D permeation flux.
5. **Active layer (hand-off to M1/M3/M6):**
   - Candidate sites must lie **within 5 µm of the front face** to appear in the 3 MeV peak, or **within 16–23 µm** for PID-window counting.
   - A PdD layer thicker than **~30 µm** adds nothing to the charged-particle signal.
   - Keep the front overlayer at **≤100 nm of oxide or Au** (≤1 % loss for p, ≤7 % for ³He); keep it ≤50 nm if the ³He line is an observable.
   - Avoid ≥1 µm Ni/Cu caps on the detector side (C5-style multilayers lose 78–84 % of ³He and ~20 % of t).
6. **Controls built into the geometry:**
   - **Blank telescope:** a second, identical telescope covered with **100 µm Al** (stops 3 MeV p; R = 82 µm). It sees all internal, cosmic and EMI backgrounds and no membrane products. It runs throughout.
   - **H twin:** an identical PdH membrane cell with an identical telescope.
   - **Pt (non-absorbing) dummy membrane** with identical current modulation. It tests EMI and thermal pickup.
   - **Detector swaps:** swap the telescopes between the D and H cells at the midpoint of the run, as an ABBA sequence.
7. **Calibration.**
   - ²⁴¹Am check source at the exempt quantity (US 10 CFR 30.71 Schedule B: 0.01 µCi) plus a pulser, for energy scale and dead layer.
   - The PID band must be calibrated with **real 1.5–3 MeV protons and deuterons** before the physics runs. Options are a d(d,p) or proton beam at a partner accelerator lab (off-site, so the charter is not violated), or ⁶LiF + thermal neutrons for 2.73 MeV tritons.
   - Neutron efficiency: use ²⁵²Cf or a D–D neutron generator at a licensed partner facility. **Do not acquire AmBe or ²⁵²Cf; they need a specific licence** (R6).
8. **Neutron bank around the C3 vacuum chamber and the C1 cell** (cavity r ≤ 8 cm):
   - 18 × 1-inch × 40 cm ³He tubes at 4 atm, at 4 cm pitch.
   - **HDPE 4.5 ± 1 cm** in front of the tube axes and **6 cm behind**, then a 1 mm Cd sheet, then **20 cm borated HDPE** (30 cm if the floor load allows: 0.18 cps).
   - A plastic-scintillator muon veto on top.
   - Log pressure and temperature every 10 s.
   - Expected: ε = 0.23, background ≈0.35 cps, 30-day reach 2×10⁻² fusions/s.
   - This bank is also the R6 safety trip. **Keep lead away from it:** lead produces muon-induced neutrons.
9. **C1 layout.** Charged products cannot leave the electrolyte (range in D₂O is 149 µm for p and 16 µm for t). C1's H1 channels are therefore:
   - the **neutron bank** above (keep the cell plus calorimeter within an r ≤ 8 cm cavity; at r = 16 cm ε falls to 0.17), and
   - **CR-39 4-channel strips** (A/B/C/D above) held against ≥30 % of the cathode behind 6 µm Mylar, with the electrolyte film minimised to **≤10 µm** by a clamped geometry.
   Expected C1 reach: n at 2×10⁻², CR-39 at ~1.5×10⁻³ fusions/s. That is 20–300× worse than C3.
10. **EJ-309:** two 5″×5″ cells as a secondary, energy-resolving neutron channel with PSD and veto. Bubble detectors for safety only.
11. **CR-39 protocol** (whenever CR-39 is used):
    - Pre-etch 1 h and scan every chip before exposure, and register the existing pits.
    - Store in N₂.
    - Use a 6 µm Mylar chemical barrier in electrolyte.
    - Run equal-area blank chips in the same electrolyte with no current, and in the H cell.
    - Etch in the same bath; read blind (chip IDs coded by a third party); analyse both faces and sequential etches.
    - Channel D must be empty of the "3 MeV p" population.
12. **Statistics:** pre-register one primary test: C3 PID window, D minus H, current on. Accumulate background exposure (blank telescope, off periods, H-control, pre-run) of **≥10× the live time**. Then 5σ needs ~8–10 signal counts rather than ~20.

## 7. Sensitivities and uncertainties

- **Si(n,p) internal background (±×3)** sets the telescope background. Even at ×3 (0.08/day), the 30-day reach moves only from 7.6 to ~9×10⁻⁵ with an equal-time control. The binomial count limit dominates, so **the recommendation is robust**.
- **α → p misidentification floor (assumed 10⁻⁴).** At 10⁻², the intrinsic α alone would add 0.3/day, and the telescope would degrade toward ~2×10⁻⁴ in 30 days. **Calibrating the misID with ²⁴¹Am and ²²⁸Th through the full PID chain is mandatory.** It could flip the "no angle collimator" decision only if misID is dominated by oblique tracks.
- **ΔE noise.** At 100 keV instead of 50 keV FWHM, p/α separation is unaffected but p/d separation degrades. This matters only if recoil deuterons are relevant (≤0.003/day).
- **Stopping in Pd (±6 %)** moves all depths and R_eff by ±6 %. No recommendation changes.
- **Cosmic neutron spectrum (±50 %)** scales the ³He bank background and the telescope (n,p) term linearly. The Si-vs-n conclusion (≥100× advantage for Si) holds for any plausible value.
- **Active-depth distribution.** This is the dominant physics unknown. If sites are bulk-uniform over >30 µm, the telescope detects 0.03–0.10 per emitted p instead of 0.22. It would then collect ~5× more in the PID window than in the peak window. **Do not pre-register the peak window alone.**
- **Response lag τ of any signal to current (M3).** If τ exceeds T_mod/5, the modulation test loses more than 70 % of its contrast. Choose T_mod from M3's τ.
- **CR-39 registration model:** the 12 keV/µm threshold and critical angles are assumptions, ±30 % on efficiency. CR-39 does not drive any primary decision.
- **Straight-line transport** neglects multiple scattering. That is fine for p and t (detour factor ≤1.05). For ³He at its end of range the error is larger, but ³He escapes only from ≤1 µm anyway.

## 8. Open questions and hand-offs

- **M3:** membrane thickness under a 1 atm ΔP with a support grid. The grid shadows the telescope, so a grid opening fraction ≥0.9 is needed. Also: bowing, which sets d_min; D permeation flux, for the gas load (≤10¹⁷ D/cm²/s is fine); and **τ_response to current steps**, for T_mod.
- **M2:** Joule and electrochemical heating of the membrane, and how much it varies with current modulation. Si leakage current doubles for every ~7 °C, so keep the telescope ΔT below 1 °C during modulation, or cool the detectors to −20 °C.
- **M1/M6:** the predicted depth distribution of candidate sites relative to the front face. It decides the peak-vs-PID emphasis and the required layer thickness.
- **M4/C1:** the size of the calorimeter jacket. It limits the neutron-bank cavity radius (ε 0.24 at r = 8 cm, 0.15 at 20 cm).
- **R5:** confirm the ULTRA-AS background spec from the datasheet. The brochure could not be fetched from this environment; the value here comes from a search snippet.
- **Lead:** an **on-site measurement** of neutron background and fast-n flux is worth more than further modelling. A single ³He tube in HDPE over one week settles the ±50 %.
- **Not modelled:** the ⁴He channel (H2), γ spectroscopy, and e⁺e⁻ (511 keV) coincidence for a hypothetical ⁴He* resonance. A 2″ NaI adds sensitivity to 23.8 MeV γ for R3/M1 if wanted.

## 9. Pre-registration template (fill in, sign, hash and timestamp before the first D run)

```
PRE-REGISTRATION  v__  date ____  SHA-256 of this file posted to ____ before run 1
1. Hypotheses
   H1a (primary): D-loaded C3 membrane emits 3.02 MeV protons at rate r_p > 0 while current is ON.
   Null: r_p(D,on) = r_p(H,on) = r_p(D,off) (within the background model).
2. Primary observable (N_trials = 1)
   Counts in telescope T1 (facing D membrane): coincidence dE>100 keV AND E>100 keV, |t_dE - t_E| < 100 ns,
   pulse-shape chi2 < cut (fixed from 241Am/pulser data), ratio dE/dE_p(E_tot) in [0.75, 2.00],
   E_tot in [1.41, 3.10] MeV, event not within 1 ms of a veto hit, detector temperature within +/-1 C.
3. Background estimate (fixed before unblinding)
   B = pooled rate in the same window from: blank telescope (100 um Al), H-cell telescope,
   D-cell OFF periods, pre-run (no loading); each weighted by live time; total background live time >= 10x.
4. Test statistic and threshold
   Li-Ma Z (on/off with tau = t_bkg/t_on) or Poisson known-B if tau >= 10. Discovery: Z >= 5.0 (global;
   primary only). Secondary channels (t window 0.80-1.05 MeV dE-only; 3He 0.60-0.85 MeV dE-only;
   n bank; EJ-309 recoil edge; alpha 9-20 MeV; p 3.5-16 MeV incl. punch-through): N_trials = 10,
   local threshold 5.43 sigma each; reported as evidence only.
5. Confirmatory test (required for a claim)
   Current modulation: square wave, 50 % duty, period T_mod = max(1 h, 10 tau_M3), phase randomised
   per cycle by the DAQ from a sealed seed; demodulated excess must be >= 3 sigma with correct sign,
   and D-specific (H cell shows < 1 sigma).
6. Consistency checks (pre-declared, must all pass)
   (a) energy spectrum consistent with the depth-convolved 3.02 MeV response (KS p > 0.01);
   (b) excess in T1 and T2 (if two telescopes face the D membrane) consistent within 2 sigma;
   (c) blank telescope shows no ON/OFF modulation (< 2 sigma);
   (d) n bank: observed n/p ratio consistent with 1 within its (large) error, OR n upper limit consistent
       with the p-inferred rate given the ~300x lower n sensitivity (NO claim is weakened by absent n
       below 1e-2 fusions/s);
   (e) no correlation of excess with detector temperature, leakage current, vacuum pressure spikes,
       or electrochemical cell voltage transients (declared regression; |coef| < 2 sigma).
7. Blinding
   Primary-window events hidden (box) until cuts frozen on: calibration data, blank telescope, H cell,
   and 10 % random "open" sample of D data. Cell labels D/H coded by a third party (both cells run
   identically; D2O vs H2O filling logged in sealed record).
8. Stopping rules
   Fixed live time: 30 days D-ON per campaign (ADR-001). Early stop only for safety trip (n bank > 10x bkg)
   or hardware failure (declared, run excluded in full). No extension based on observed significance.
9. Control schedule (per 30-day campaign)
   Day -7..0: pre-run (no loading), both telescopes + blank; 241Am/pulser calibration.
   Day 1..30: D and H cells loaded in parallel; modulation per item 5; telescopes swapped D<->H at day 15;
              blank telescope continuous; Pt dummy membrane run for 3 days before and after.
   Day 31..37: post-run background (unloaded); calibration repeated; CR-39 chips etched blind.
10. Reporting
   Report Z, FC 90 % interval on r_p (published FC table values), Bayesian 90 % UL (flat prior),
   efficiency with uncertainties, and the full background table - regardless of outcome.
```
