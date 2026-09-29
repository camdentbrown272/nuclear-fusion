# M2 — Electrochemical cell geometry, current distribution and loading

Code: [`sim/m2_common.py`](../../sim/m2_common.py) (isotherm and kinetics), [`sim/m2_fv.py`](../../sim/m2_fv.py) (finite-volume solver), [`sim/m2_geom.py`](../../sim/m2_geom.py) (geometries), [`sim/m2_current.py`](../../sim/m2_current.py), [`sim/m2_loading.py`](../../sim/m2_loading.py), [`sim/m2_cell.py`](../../sim/m2_cell.py).
Raw outputs: [`figs/m2_current.txt`](figs/m2_current.txt), [`figs/m2_loading.txt`](figs/m2_loading.txt), [`figs/m2_cell.txt`](figs/m2_cell.txt), [`figs/m2_m3_interface.txt`](figs/m2_m3_interface.txt).
To reproduce every number: `pip install numpy scipy matplotlib`, then `python3 sim/m2_current.py && python3 sim/m2_loading.py && python3 sim/m2_cell.py` (about 80 s total).

## Summary for lead

1. **Geometry does not create high loading; the surface does.** x rises only +0.03–0.05 per decade of current. x ≥ 0.90 needs a vacuum-annealed + etched ("good") surface at **≥ 86 mA/cm²**. Typical surfaces need about 2 A/cm². x ≥ 0.95 needs recombination poisoned about 10–30× and ≥ 57 mA/cm²; no geometry reaches it otherwise.
2. **Geometry's job is to leave no area under- or over-polarised.** At loading currents Wa = 0.001–0.1, so the current distribution is primary. Flush disks and free wire tips give 3–7× edge peaks. Fixes: insulating walls meeting the electrode at 90° (tube cell, end plates), or a PTFE collar deeper than the disk radius.
3. **C1:** Ø1.0 mm Pd wire, 30 mm long, between PTFE end plates. Pt helix anode at R = 10 mm spanning the same plates; eccentricity ≤ 2 mm. Uniformity is exact. 6.0 V and 1.7 W at 300 mA/cm² in 0.1 M LiOD.
4. **C3 (DFM):** Ø20 mm disk as the floor of a tube whose bore equals the wetted diameter. Pt mesh anode at g = 4.0 ± 0.2 mm, **1.0 M LiOD**: ±2.2 %, 3.2 V and 2.0 W at 200 mA/cm². 0.1 M LiOD with a 10–15 mm gap needs 16–46 V and is not viable. Seal edge flush to ±0.1 mm; mask r > 9.5 mm. Alternative: Ø10 mm disk with a ring anode at g ≥ 1.1 × radius.
5. **Biggest risk (M3/M1):** a clean Pd vacuum face drains the membrane. x_out ≥ 0.85 needs exit desorption k_des ≤ 2×10⁻¹⁰ mol cm⁻² s⁻¹ atm⁻¹, about 10⁸× below clean Pd.
6. x = 0.90 needs a fugacity of 3×10³–3×10⁴ atm (0.12 eV). Only about 0.06–0.15 V of the 0.3–0.55 V overpotential loads the lattice.

---

## 1. Questions answered

1. Which electrode geometry maximises the area fraction at x ≥ 0.90 / 0.95 without hot spots, excessive ohmic heating or bubble screening? (§5.1–5.3, §6)
2. For the DFM (C3), which anode and masking give ±5 % over a 10–25 mm disk? (§5.2, §6.2)
3. Which operating window (current density, electrolyte, temperature) gives the required loading? (§5.4–5.6, §6.3)
4. Interface: entry-side x_in(i, permeation flux) for M3. (§5.7, [`figs/m2_m3_interface.txt`](figs/m2_m3_interface.txt))

## 2. Model and equations

### 2.1 Field (current distribution)
- In the electrolyte: ∇·(κ∇Φ) = 0, with Φ in V and κ in S/cm. Axisymmetric (r, z) for the wire, DFM and sphere; planar (x, y) per unit depth for the foil. Finite volumes on non-uniform tensor grids refined at edges (`m2_fv.py`, `scipy.sparse`).
- Materials: electrolyte, Pd cathode (metal at 0 V), Pt anode (Dirichlet Φ = V_a), insulator (zero flux). The cell walls and free surface are insulating.
- **Primary**: Φ = 0 on the cathode.
- **Secondary**: on each cathode face, h(Φ_P − Φ_s) = i(η) with η = −Φ_s and h = κ/(Δ/2). i(η) comes from the kinetic model (§2.2). Newton iteration on Φ, with an outer root-find on V_a so that the total current is imposed.
- Wagner number: Wa = κ (dη/di) / L.
- **Bubbles**: Bruggeman κ_eff = κ(1−ε)^1.5. Void fraction by drift flux: ε = j_g/(j_g + u_b), with superficial gas velocity j_g = i V_m/(2F).
  - Wire: the plume width grows as d(z) = 0.5 mm + 0.1 (z − z_tip).
  - DFM: a uniform D₂ column in the gap.
  - Adherent-bubble coverage: Θ ≈ 0.023 (j/[A m⁻²])^0.3.

### 2.2 Loading chain i → η → f → x
Volmer–Tafel–Heyrovsky kinetics with Langmuir adsorption (y = θ/(1−θ)), φ = Fη/RT and η measured against the reversible D₂ electrode:

- Volmer: v_V = k_V/(1+y) · [e^{−βφ} − (y/K) e^{(1−β)φ}]
- Heyrovsky: v_H = k_H/(1+y) · [y e^{−β_Hφ} − K p e^{(1−β_H)φ}]
- Tafel: v_T = k_T/(1+y)² · [y² − K² p]

Detailed balance holds exactly: at φ = −½ ln p, y = K√p.

Steady state: v_V = v_H + 2v_T + J_abs, with i = F(v_V + v_H). J_abs is the permeation drain.

Surface/subsurface equilibrium gives the **effective fugacity f = y²/K² (atm)**. Two limits follow:
- Volmer–Tafel: f ≈ i/c_T
- Volmer–Heyrovsky: f → f_max = (k_V/(k_H K))²

The model is parameterised by observables:

| Parameter | Controls |
|---|---|
| i₀V | η only (the Tafel line); has no effect on f |
| c_T | recombination activity |
| f_max | the surface ceiling on loading |
| K | coverage scale |

**Isotherm (β-PdD).** A mean-field lattice gas:

½RT ln f = Δh + W(x − 0.7) − TΔs + RT ln[x/(1−x)]

It is anchored at x = 0.67 at 1 atm and x = 0.90 at F90 = 10⁴ atm (band 3×10³–3×10⁴). Δh = −17.5 kJ/mol D; W = 33 kJ/mol follows from the anchors. It is valid only for the β phase (x ≳ 0.6).

### 2.3 Cell budget
- V_cell = E_rev + |η_c| + η_a + I R_Ω(ε), with E_rev(D₂O) = 1.2615 V.
- Open-cell heat: I(V − E_tn). Closed cell: all of I·V becomes heat, of which I·E_tn is released at the recombiner.
- DFM temperature bound: 1-D conduction only, with q‴ = i²/κ_eff in the gap, an adiabatic (vacuum-backed) membrane and the anode plane at bulk temperature. Convection makes the real rise smaller.
- Permeation illustration: J = D c_Pd (x_in − x_out)/L = k_des f_out.

## 3. Parameters

Note on sources: publisher and repository sites (ScienceDirect, OSTI, lenr-canr, arXiv, JCMNS, Wikipedia) were blocked from this session. Values below come from search-result abstracts (marked **abs**) or from the cited paper as recalled (marked **rec**). Anything **rec** is bracketed in §7.

| Parameter | Value | Units | Source | Uncertainty |
|---|---|---|---|---|
| κ(LiOD/D₂O) | 60.56c − 14.25c² + 2.514cT − 0.5459c²T (0.1 M, 25 °C: 12.1; 1 M: 95.5) | mS/cm | https://jcmns.org/api/v1/articles/72600-conductivity-and-molar-conductivity-of-liod-heavy-water-solution.pdf (abs) | ±5 %; validity at 1 M not verified |
| E_tn(D₂O), ΔH_f | 1.5267 V, −294.60 kJ/mol | | https://jcmns.org/article/72550-the-thermoneutral-potential-in-electrochemical-calorimetry-for-the-pd-d2o-system.pdf (abs) | ±0.001 V |
| E_tn(H₂O) | 1.4812 V (ΔH_f −285.83 kJ/mol) | | https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185 (rec) | ±0.001 V |
| E_rev(D₂O) | 1.2615 V (ΔG_f −243.44 kJ/mol) | | https://webbook.nist.gov/cgi/cbook.cgi?ID=C7789200 (rec) | ±0.003 V |
| Isotherm anchor x(1 atm) | 0.67 | | PdH ≈ 0.70 at 1 atm: https://en.wikipedia.org/wiki/Palladium_hydride (abs); PdD lower (rec) | ±0.02 |
| Isotherm anchor f(x = 0.90) | 10⁴ (band 3×10³–3×10⁴) | atm | Baranowski high-pressure data as used in the LENR literature (rec); https://www.sciencedirect.com/science/article/abs/pii/002250889090070Z | factor 3 |
| Check: electrolytic loading equivalent | 80–150 bar | bar | Baranowski, Filipek, Raczyński: https://www.osti.gov/etdeweb/biblio/76885 (abs) | — |
| Check: "1 V ≈ 800 atm" | 800 | atm | Berlinguette et al., Nature 644, 640 (2025): https://www.nature.com/articles/s41586-025-09042-7 ; https://science.ubc.ca/news/2025-08/researchers-use-electrochemistry-boost-nuclear-fusion-rates (abs) | press statement |
| Check: x = 1 formation pressure | 2.7 GPa (PdD), 1.9 GPa (PdH) | | search abstract, https://arxiv.org/pdf/1708.09316 (abs) | — |
| Plateau Δh (PdD) | −17.5 | kJ/mol D | ≈ −35 kJ/mol D₂ (rec) | ±3 |
| Mechanism VT → VH; maximum loading depends on symmetry factors | — | | Zhang et al.: https://www.lenr-canr.org/acrobat/ZhangWSthemaximum.pdf ; https://www.researchgate.net/publication/245146846 (abs) | — |
| H vs D kinetics, alkaline Pd | D slower | | https://www.sciencedirect.com/science/article/pii/0022072896045895 (abs) | — |
| Tafel slope | 0.114–0.118 (model) | V/dec | 102–130 mV/dec on Pd in base: https://www.sciencedirect.com/science/article/abs/pii/036031999400112D (abs); 120 mV theory: https://pmc.ncbi.nlm.nih.gov/articles/PMC11019648/ (abs) | ±20 % |
| i₀V | 3×10⁻⁵ | A/cm² | 2.2×10⁻⁴ for Pd in acid (H): https://pmc.ncbi.nlm.nih.gov/articles/PMC10979421/ (abs); lower for D and in base (assumption) | ×10 (changes V only) |
| c_T (poor / typical / good / exceptional) | 10⁻³ / 10⁻⁴ / 10⁻⁵ / 10⁻⁶ | A cm⁻² atm⁻¹ | fitted to validation bands (§4.3) | class definition |
| f_max (same classes) | 3×10³ / 3×10⁴ / 3×10⁵ / 3×10⁶ | atm | fitted | class definition |
| K, β, β_H | 3×10⁻³, 0.5, 0.5 | | assumption | §7 |
| Thiourea retards Tafel | qualitative | | https://www.sciencedirect.com/science/article/abs/pii/002207289303067Y (abs) | — |
| D in β-PdD, 300 K | 2.5×10⁻⁷ | cm²/s | H in Pd 4.2×10⁻⁷ at 20 °C, https://www.sciencedirect.com/science/article/abs/pii/0360319985900667 (abs); isotope factor about 0.6 (rec) | ×2 |
| c_Pd | 0.1129 | mol/cm³ | ρ = 12.02 g/cm³ | — |
| Bubble coverage Θ ∝ j^0.3 | prefactor 0.023 | | Vogt & Balzer, Electrochim. Acta 2005: https://www.sciencedirect.com/science/article/abs/pii/S001346860400948X (exponent abs, prefactor rec) | ×2 |
| Bubble-swarm rise velocity u_b | 1 (0.3–3) | cm/s | Stokes 0.14–5 cm/s for 50–300 µm (computed) plus plume entrainment | ×3 |
| OER on Pt: b, i₀ | 0.12 V/dec, 10⁻⁶ A/cm² | | assumption (η_a ≈ 0.44–0.6 V) | ±0.15 V |
| k_th, c_p of D₂O | 0.595 W/m/K, 4.21 J/g/K | | rec | ±3 % |

## 4. Verification

1. **Coax between plates** (R = ln(R_a/a)/(2πκL)). Error is +0.77 %, +0.65 %, +0.59 % at grid factors 1, 0.5, 0.25; the distribution is exactly uniform.
2. **Newman disk** (i/i_avg = ½(1−r²)^−½). At r = 0, 0.5, 0.8, 0.9 the numerical values are 0.502/0.575/0.814/1.131 against analytic 0.500/0.573/0.811/1.125. Error is ≤ 0.6 % on the fine grid; the resistance is within 1 % of 1/(4κR_d).
3. **1-D secondary check** (tube, g = 10 mm): V_a matches iR + |η| within 0.4–0.5 % at 0.01–0.5 A/cm².
4. **Eccentric cylinders**: exact bipolar solution against the small-e limit 1 + 2ae/R_a², agreeing to 10⁻³.
5. **Kinetic limits.** The Tafel slope comes out at 116 and 114 mV/dec. The Tafel share of recombination is 98 % at 1 mA/cm² and 56 % at 1 A/cm² (the VT→VH transition). Detailed balance is exact at η = 0 (f = 1 atm, x = 0.67).
6. **Validation bands** (shaded boxes in `m2_loading_vs_i.png`):
   - untreated Pd: x = 0.75–0.82 → matches the typical/poor classes
   - vacuum-annealed + etched: 0.91–0.93 → good class at 0.1–0.3 A/cm²
   - SRI cathodes: 0.85–0.95 at 0.1–0.6 A/cm²
   - Baranowski 80–150 bar → x = 0.79–0.81
   - "800 atm" → x = 0.844

   These bands define the surface classes. They are **calibration, not an independent test.**

## 5. Results

### 5.1 Loading vs current density (entry side, 25 °C, no drain)

| i (mA/cm²) | poor | typical | good | exceptional |
|---|---|---|---|---|
| 10 | 0.733 | 0.792 | 0.851 | 0.910 |
| 100 | 0.786 | 0.845 | 0.903 | 0.961 |
| 300 | 0.810 | 0.868 | 0.924 | 0.974 |
| 1000 | 0.832 | 0.889 | 0.942 | 0.979 |
| i for x ≥ 0.90 | never | 2.2 A/cm² | **86 mA/cm²** | 6.6 mA/cm² |
| i for x ≥ 0.95 | never | never | 2.6 A/cm² | **57 mA/cm²** |

- **Fugacity for x = 0.90:** 3×10³–3×10⁴ atm, i.e. Δμ_D = 0.12 eV. For x = 0.95: 4×10⁴–7×10⁵ atm.
- **Nernstian comparison.** Treating the whole overpotential as loading, f = exp(−2Fη/RT), would give 10¹⁰–10¹⁶ atm; the model gives 10²–10³.
- **Mechanism.** The overpotential that actually loads the lattice is η_H = −0.06 to −0.15 V. The rest is Volmer activation (see `m2_loading_vs_i.png`).
- **Recombination poisons** (factor P on k_T, k_H), for the typical surface at 100 mA/cm²:

  | P | 1 | 3 | 10 | 30 |
  |---|---|---|---|---|
  | x | 0.845 | 0.877 | 0.910 | 0.942 |

### 5.2 Current distribution by geometry
Primary distribution at 100 mA/cm²; loading uses the good-surface class (full tables in `m2_current.txt`). Figures: `m2_coax_maps.png`, `m2_dfm_profiles.png`, `m2_dfm_sweeps.png`, `m2_foil.png`, `m2_area_fraction.png`.

| Geometry | max/mean | min/mean | area within ±5 % | V_a at 100 mA/cm² (0.1 M) |
|---|---|---|---|---|
| A1 coax, free tip (a = 0.5, R_a = 10 mm) | 2.87 (tip) | 0.95 | 90 % | 1.6 V |
| A2 **coax between PTFE end plates** | 1.000 | 1.000 | 100 % | 1.7 V |
| A2b same + bubble plume, 500 mA/cm² | 1.18 | 0.98 | 92 % | 7.5 V |
| A3 coax, tip in PTFE sleeve | 1.34 | 0.97 | 88 % | 1.6 V |
| A5 coax, a = 1 mm, R_a = 20 mm, free tip | 2.67 | 0.92 | 22 % | 2.7 V |
| A7 coax, anode half the wire length | 2.20 | 0.87 | 32 % | 1.8 V |
| B1 foil, open edges | 7.51 | 0.72 | 6 % | — |
| B2 **foil framed in insulating slot** | 1.000 | 1.000 | 100 % | — |
| C1 DFM flush in wide floor | 7.30 | 0.56 | 4 % | 4.5 V |
| C2 DFM, collar 5 mm | 1.04 | 0.93 | 88 % | 7.6 V |
| C3 DFM, collar 10 mm | 1.005 | 0.988 | 100 % | 11.7 V |
| C4/C5 **DFM tube + mesh anode** | 1.000 | 1.000 | 100 % | 8.8 V (g = 10) / 2.9 V (g = 3) |
| C6 DFM tube + ring anode, g = 10 mm | 1.022 | 0.949 | 99 % | 9.7 V |
| C11 **DFM R_d = 5 mm + ring anode, g = 7 mm** | 1.005 | 0.988 | 100 % | 6.3 V |
| C12 DFM mesh g = 4 mm, 0.2 mm sag | 1.022 | 0.988 | 100 % | — |
| C9 seal overhang 0.5 mm (0.2 mm crevice) | 1.11 | 0.03 | 3 % | — |
| C10 seal recess 0.5 mm (flush annulus) | 3.33 | 0.91 | 15 % | — |

**Area fraction with x ≥ 0.90.** Because x depends on log i, this fraction is simply the fraction of area with i ≥ i₉₀ (86 mA/cm² for a good surface). The same holds for x ≥ 0.95 with i₉₅.

- **Uniform geometries** (A2, B2, C3–C5, C11): the fraction jumps from 0 to 100 % at i_avg = i₉₀.
- **Nonuniform geometries** (C1, B1): the fraction is 36 % at 100 mA/cm² (good surface) and reaches 100 % only at 300 mA/cm². Even then, 90 % of that area sits at 0.56× the mean current, so the flush disk wastes about 1.8× current and heat for the same loading.
- **x ≥ 0.95**: reached nowhere below 0.5 A/cm² except with the exceptional surface class. In the tables, the good surface reaches only 1 % (flush disk at 0.5 A/cm²).

### 5.3 Wagner number (typical kinetics, 0.1 M)

| i (mA/cm²) | κR_ct | Wa (L = 0.5 mm wire) | Wa (L = 10 mm disk) |
|---|---|---|---|
| 1 | 6.1 mm | 12 | 0.61 |
| 10 | 0.61 mm | 1.2 | 0.061 |
| 100 | 0.060 mm | 0.12 | 0.006 |
| 500 | 0.012 mm | 0.024 | 0.001 |

At loading currents the distribution is primary, so geometry must be uniform by construction.

### 5.4 Tolerances (primary distribution)

- **DFM ring anode in a tube.** The smallest gap giving ±5 % is g_min = 1.0–1.1 × R_d (5, 8, 11 and 13 mm for R_d = 5, 7.5, 10, 12.5 mm).
- **Collar in a wide cell.** A well depth ≥ 0.75–1.0 × R_d is needed. For R_d = 10 mm: at h = 7.5 mm the range is [0.969, 1.015], at 10 mm [0.988, 1.005].
- **Mesh anode sag/tilt.** At g = 4 mm: sag 0.2 mm → ±2.2 %, 0.5 mm → [0.970, 1.057]. At g = 3 mm: 0.3 mm → [0.970, 1.050].
- **Mesh or helix pitch.** Ripple is 2e^(−2πd/p): 0.37 % at p = d, negligible for p ≤ d/2.
- **Coax eccentricity.** For ±5 %: e ≤ 4 mm (a = 0.5 mm, R_a = 10 mm), 2.3 mm (a = 1 mm, R_a = 10 mm), 1.1 mm (a = 0.5 mm, R_a = 5 mm).
- **Seal step** (R_d = 10 mm tube):
  - −0.1 mm (a flush insulating annulus 0.1 mm wide): 95 % of the area within ±5 %, with a rim peak of 1.8×.
  - +0.1 mm (a crevice): 93.5 % within ±5 %, and the crevice bottom falls to 0.58×.
  - ±0.5 mm: only 3–15 % within ±5 %.

### 5.5 Cell voltage, heat and temperature (from `m2_cell.txt`)

**C1**, 1 mm wire × 30 mm, R_a = 10 mm:

| i (mA/cm²) | 0.1 M: V_cell / P_in / open-cell heat | 1.0 M: V_cell / P_in |
|---|---|---|
| 100 | 3.4 V / 0.32 W / 0.17 W | 2.3 V / 0.21 W |
| 300 | 6.0 V / 1.7 W / 1.3 W | 2.7 V / 0.8 W |
| 500 | 8.7 V / 4.1 W / 3.4 W | 3.1 V / 1.5 W |

**C3**, Ø20 mm (I = 0.314 A per 100 mA/cm²). ΔT is the conduction-only upper bound.

| Configuration | i (mA/cm²) | V_cell | P_in | ΔT bound |
|---|---|---|---|---|
| ring g = 15 mm, 0.1 M | 100 | 16.3 V | 5.1 W | 164 K ✗ |
| ring g = 15 mm, 0.1 M | 300 | 45.9 V | 43 W | ✗ |
| mesh g = 3 mm, 0.1 M | 100 | 4.8 V | 1.5 W | 7 K |
| mesh g = 3 mm, 0.1 M | 200 | 7.5 V | 4.7 W | 28 K |
| **mesh g = 4 mm, 1.0 M** | 200 | 3.2 V | 2.0 W | 8 K |
| mesh g = 3 mm, 1.0 M | 300 | 3.4 V | 3.2 W | 11 K |
| ring g = 15 mm, 1.0 M | 200 | 6.0 V | 3.8 W | 91 K |


**Temperature.** Loading falls about 0.006 per +10 K at fixed fugacity (for example, good surface at 100 mA/cm²: 0.915 at 5 °C, 0.903 at 25 °C, 0.883 at 60 °C). This is an upper bound on x at high T. Conductivity falls about 40 % from 25 to 5 °C.

**Bubbles.** The D₂ column void is 1.3 % at 100 mA/cm² and 6 % at 500 mA/cm² (u_b = 1 cm/s), or 4 % / 17 % at u_b = 0.3 cm/s. For the wire, the plume at 500 mA/cm² reduces the within-±5 % area from 100 % to 92 %. Adherent coverage Θ = 0.18–0.30 raises the local current on the free area by 1.2–1.4×, which changes x by only +0.004–0.005.

**Recombiner.** Heat is 1.527 W/A for D₂O and 1.481 W/A for H₂O (a 45 mW/A difference in the light-water control). Gas flow is 7.6 mL/min D₂ plus 3.8 mL/min O₂ per ampere. A stoichiometric headspace stores 136 J per 20 mL.

**Loading time.** A 1 mm wire needs 77 C/cm to reach x = 0.9, and 90 % diffusional equilibration takes about 1.4 h. A 50 µm membrane needs 49 C/cm², which is 41 min at 20 mA/cm² at 100 % efficiency; its diffusion time is 100 s.

### 5.6 Sensitivity of x(i) (good surface, 100 mA/cm²; baseline 0.903)

| Assumption | Effect on x(100 mA/cm²) |
|---|---|
| Isotherm band (F90 = 3×10³ / 3×10⁴) | 0.929 / 0.882 |
| c_T ×3 / ÷3 | 0.879 / 0.925 |
| f_max ×10 / ÷10 | 0.910 / 0.889 |
| β_H = 0.60 | loading peaks at 0.887 near 0.3–0.5 A/cm², then falls (the maximum reported by Zhang et al.) |
| β_H = 0.40 | monotonic; 0.974 at 1 A/cm² |
| K ×10 | 0.952 |
| K ÷10 | 0.894 |
| i₀V | no effect on x |

### 5.7 DFM with permeation (hand-off to M3)
Entry-side x_in(i, i_p/i) for four surface classes is in `m2_m3_interface.txt` (and `f_in_with_drain()` in the code). A drain of 10 % of the current lowers x_in by 0.003; a 50 % drain lowers it by 0.02–0.03.

In the coupled illustration (good surface, 50 µm), both the flux J and x_out are set by the exit-face desorption coefficient k_des:

| k_des (mol cm⁻² s⁻¹ atm⁻¹) | 10⁻¹¹ | 10⁻⁹ | 10⁻⁷ | 10⁻³ |
|---|---|---|---|---|
| x_out at 300 mA/cm² | 0.903 | 0.832 | 0.733 | 0.497 |
| J (D cm⁻² s⁻¹) | 7×10¹⁶ | 3×10¹⁷ | 6×10¹⁷ | 1.3×10¹⁸ |

The clean-Pd estimate of k_des (2sZ) is 0.03–2.6, so clean Pd sits at or beyond the right-hand end of this table. For x_out ≥ 0.85, k_des must be ≤ 1–6×10⁻¹⁰ across L = 25–100 µm and i = 100–300 mA/cm².

## 6. Design recommendations for the lead

### 6.1 C1 coaxial calorimetric cell
1. **Cathode:** Pd wire Ø1.00 ± 0.02 mm, wetted length 30.0 ± 0.5 mm.
   - The ends pass through PTFE or PCTFE end plates with an interference fit (hole Ø0.95 mm). This avoids crevices; a 0.1 mm crevice leaves an annulus at 0.58× current.
   - Reason: end plates give exact uniformity (§5.2 A2), whereas a free tip has a 2.9× hot spot.
   - Alternative Ø0.5 mm: loads 4× faster (0.35 h) at half the power.
2. **Anode:** Pt wire helix or cylindrical mesh at R_a = 10 ± 1 mm, spanning exactly the same plates. Pitch ≤ 5 mm (ripple < 0.01 %). Concentricity ≤ 2 mm (±2 %). An anode shorter than the wire drops the within-±5 % area to 32 %, so do not shorten it.
3. **Supply:** 0.1 M LiOD for lineage comparability, giving V_cell 6.0 V and 1.7 W at 300 mA/cm². Use a supply rated ≥ 20 V / 0.5 A compliance for bubbles and aging. Use 0.5 M if M4 wants less than 1 W of input at 300 mA/cm².
4. **Recombiner:** hydrophobic (PTFE-bonded) Pt or Pd on alumina, in the headspace ≥ 30 mm above the electrolyte and directly above the anode, so gases rise through it by buoyancy. Condensate should drip back to the electrolyte away from the cathode. Headspace ≤ 20 mL (≤ 136 J stored). Heat load is 1.527 W/A, inside the calorimetric boundary (M4).

### 6.2 C3 detector-facing membrane
5. **Cell:** the disk is horizontal, wetted face up (cell floor), with the vacuum and detectors below. D₂ bubbles leave vertically; a ceiling-mounted disk would trap gas.
6. **Disk and walls:** wetted Ø20.0 mm (R_d = 10 mm). The cell is a PCTFE or PTFE tube with bore 20.00 ± 0.05 mm, sealed so the gasket edge is flush with the bore within ±0.1 mm. The detector aperture or collimator should exclude r > 9.5 mm, where the ±0.1 mm seal tolerance leaves a residual 0.58–1.8× rim.
7. **Anode:** planar Pt mesh (≥ 50 % open, pitch ≤ 1 mm) spanning the bore at g = 4.0 ± 0.2 mm, flat or parallel within 0.2 mm (±2.2 %).
   - **Electrolyte: 1.0 M LiOD** (V_cell 3.2 V, 2.0 W and ΔT bound 8 K at 200 mA/cm²).
   - With 0.1 M LiOD, use g = 3 mm and i ≤ 100 mA/cm² (4.8 V, 1.5 W).
   - Downside: D₂ bubbles pass through the anode (parasitic oxidation; Pt redeposition on the membrane).
   - Alternative that separates the gas streams: Ø10 mm disk, Pt ring anode (≥ 1 mm wire or band) at the wall, g = 7 ± 0.5 mm (±1.2 %). Area and detector solid angle drop 4× (M5 to weigh this).
8. **Do not use** a flush disk in a wide floor (7× edge peak) or a collar shallower than R_d.
9. **Exit (vacuum) face:** it must carry a desorption barrier with k_des ≤ 2×10⁻¹⁰ mol cm⁻² s⁻¹ atm⁻¹, or x_out stays below 0.85 whatever the electrochemistry does. This is the dominant C3 design variable (M3/M1).

### 6.3 Operating window (both cells)
10. **Surface preparation** is mandatory (vacuum anneal plus acid etch). Without it x ≥ 0.9 is unreachable.
11. **Ramp** (hold each step until the resistance ratio or x plateaus):

    | Step | Current density | Hold (wire / membrane) | Purpose |
    |---|---|---|---|
    | 1 | 5–10 mA/cm² | ≥ 4 h / ≥ 1.5 h | through α→β, M3 cracking limits |
    | 2 | 50 mA/cm² | ≥ 2 h / ≥ 0.5 h | |
    | 3 | 100 mA/cm² | | x ≈ 0.90 for a good surface |
    | 4 | 200–300 mA/cm² | | C1 up to 500 only if x has not saturated; C3 ≤ 300 |

    Above 300 mA/cm² loading gains are ≤ 0.01–0.02 while heat grows as i².
12. **Temperature:** 10–25 °C electrolyte during loading, with 5 °C as a hard floor (D₂O freezes at 3.8 °C). Each 10 K lower gains about +0.006 in x.
13. **Additives:** thiourea, and the As/S/CN⁻ class, act as recombination poisons. In this model a 10× poisoning factor gives +0.065 in x at 100 mA/cm². Al and Si (from glass or added) and Pt deposits have no quantitative literature model, so they are listed only. Pt deposits are expected to lower loading because they catalyse recombination. Any additive must also go into the H₂O control.

## 7. Sensitivities and uncertainties (what could flip a recommendation)
- **Surface class and isotherm anchor.** The "86 mA/cm² for x = 0.9" figure moves between about 30 and 300 mA/cm² across the F90 band and c_T ×3. The geometric recommendations do not change, because they depend only on uniformity.
- **Sign of β_H − β.** If β_H > β, loading peaks at 0.3–0.5 A/cm² and falls at higher current. Measure in situ; never run above the peak.
- **1 M LiOD.** The conductivity fit may not be valid at 1 M, and a possible chemical effect of Li concentration on loading is not modelled. If the lead insists on 0.1 M, C3 must use g ≤ 3 mm and i ≤ 100 mA/cm².
- **ΔT bounds** are conduction-only. Real bubble-driven convection lowers them, so M4 must compute actual temperatures.
- **Bubble parameters** (u_b, Θ prefactor) are uncertain to ×3. They affect voltage by ≤ 10 % and loading negligibly.
- **Permeation** (§5.7) uses constant D, a linear gradient, and a β-only isotherm. x_out < 0.6 is outside the model's validity.

## 8. Open questions / hand-offs
- **M3:** use `m2_m3_interface.txt` and `f_in_with_drain()` as the entry boundary condition. Model concentration-dependent D, α/β cracking during step 1, and the exit-face k_des. M3 owns this number, and it decides whether C3 can reach x_out ≥ 0.85.
- **M4:** heat loads (C1 0.2–4 W; C3 1–3 W plus recombiner 1.53 W/A), real ΔT with convection, and recombiner placement inside the calorimeter.
- **M5:** C3 aperture r ≤ 9.5 mm; Ø20 vs Ø10 mm disk trade (mesh vs ring anode); D₂ gas load on the vacuum side, up to 10¹⁷–10¹⁸ D cm⁻² s⁻¹ if the exit is unpassivated.
- **M1/M6:** the exit-face barrier (PdO/CaO/multilayer) must be judged both for site physics and for k_des.
- **Experiment:** measure c_T and f_max for the actual cathode lot, using the loading-vs-current curve at 10/50/100 mA/cm² with the resistance ratio, before committing to a current schedule.
