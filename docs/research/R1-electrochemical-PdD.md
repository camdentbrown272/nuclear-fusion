# R1 — Electrochemical Pd–D (Fleischmann–Pons-type) Systems: Geometry, Materials, and What Correlates with Reported Positive Results

*Research brief for the LENR cell-geometry design effort. Scope: electrolytic loading of deuterium into palladium (bulk cathodes and co-deposited films). Plasma, fusor, and beam-driven systems are out of scope except where they bear directly on Pd–D loading physics.*

---

## 0. How reliable this document is (read first)

- **Network limits during this session:** direct fetching of primary PDFs was blocked by the environment's egress policy for every host tried (lenr-canr.org, nature.com, jcmns.org, arxiv.org, osti.gov, currentscience.ac.in, rle.mit.edu, wikipedia.org). Literature search worked for a limited number of queries before the session's search budget ran out. So every claim below comes from (a) search-result abstracts and snippets that confirm specific facts, or (b) the author's knowledge of the primary literature.
- **The dagger mark (†)** flags numbers that could **not** be checked against the primary source in this session. Treat them as best recollection. They must be checked against the cited PDF before anyone uses them in a quantitative model. Numbers without † were confirmed by at least one retrieved snippet or are simple derived calculations.
- **Evidence grades used throughout:**
  - **A** = independently replicated in peer-reviewed mainstream venues.
  - **B** = replicated by several LENR groups but not accepted by the mainstream.
  - **C** = single group (or one group plus close collaborators).
  - **D** = anecdotal or disputed.
  - **N** = a credible null result against the claim.

---

## 1. Executive summary

1. **Reports of excess heat in electrolytic Pd–D cells come from several independent LENR laboratories**: Fleischmann–Pons, SRI, ENEA, China Lake/NRL, Energetics, Storms, and IMRA Japan. Typical sizes are 5–30 % of input power, at mW to W scale, with rare episodes claimed above 100 %. None has been reproduced on demand, and the mainstream does not accept them. **Grade B for "reported by multiple groups"; N for any strong claim**, given the 1989 Caltech, Harwell and MIT nulls, the 1992–97 Japanese NHE program, and the 2019 Google program.
2. **A high average loading, D/Pd ≳ 0.85–0.90 measured by resistance, is the most consistently reported necessary condition.** No well-characterised positive result below about 0.83 is known. Loading is necessary but clearly not sufficient: NHE and others reached about 0.9 with no heat. **B (necessary); N (sufficient).**
3. **SRI's empirical rule:** excess power appears only above a loading threshold x₀ ≈ 0.83–0.875 and a current-density threshold i₀ ≈ 0.1–0.4 A cm⁻². It scales roughly as (x − x₀)²(i − i₀), and later versions multiply by the deuterium flux |dx/dt|. Loading must be held for ≥ 10 diffusion time constants, which in practice meant weeks. **C** (SRI; similar thresholds reported by IMRA, ENEA and Kunimatsu).
4. **Heat–helium-4 correlation.** China Lake found He-4 correlated with excess heat in 18 of 21 runs, at 10¹⁰–10¹² He s⁻¹ W⁻¹. SRI's M4 cell recovered ≈ 104 ± 10 %† of the 23.85 MeV/He expectation after flushing the cathode. ENEA and Rome report similar values. This is the most physically meaningful signature claimed. However, most of the measured He sat *below* the atmospheric level (5.24 ppm), so a leak cannot be ruled out. **B− / D (disputed; Jones & Hansen 1995).**
5. **Metallurgy dominates reproducibility, but only one group claims control of it.** ENEA reports that annealed, etched 50 µm foils with grain size < 100 µm and a predominant ⟨100⟩ texture reach D/Pd > 0.95 in > 90 % of samples. ENEA also reports that excess heat correlates with the surface-roughness power spectral density (PSD) at micron to sub-micron scales. **C.**
6. **"Triggers"** (laser stimulation of Au-plated cathodes, dual-laser beat frequencies of 8.3, 15.3 and 20.4 THz, SuperWave current modulation, magnetic fields) are single-group claims with at most partial replication by close collaborators. **C/D.**
7. **Nuclear radiation is at least ~10¹⁰ below what heat of fusion origin would imply.** The 1989 Fleischmann–Pons neutron and gamma claim was withdrawn in practice after Petrasso et al. **N.** The SPAWAR CR-39 "triple tracks" from co-deposition are disputed because the detector sat in chemical contact with the electrolyte. **D.**
8. **Google-funded program** (Berlinguette et al., *Nature* 2019). It could not reproducibly make or hold the claimed D/Pd regime. Best loading was H/Pd = 0.96 ± 0.02 with *light* hydrogen. It found no anomalous heat or nuclear products. **N, but the program did not test the x > 0.9 regime squarely.** Its 2025 follow-up (*Nature*) showed electrochemical loading raising *beam-target* D–D fusion rates by 15 %. That result is solid but is not LENR.
9. **The artifacts are real and documented in mainstream journals** (**A**): recombination and Faradaic efficiency below 100 % in open cells, calibration-constant shifts, stored chemical energy, He contamination, tritium contamination and electrolytic enrichment, and neutron-detector noise. A credible design has to remove them *by construction*.
10. **What this means for design:** the highest-value choices are:
    - a closed, He-tight cell with an internal recombiner and flow or Seebeck calorimetry;
    - a thin, pre-screened, annealed Pd cathode in a geometry that gives uniform current density and uniform loading, with resistance used to verify x ≥ 0.95;
    - a small gas headspace, so that about 1 W of 24 MeV/He heat drives He-4 *above* atmospheric concentration within about a day.

---

## 2. Table of experiments

Legend: n.m. = not measured; JM = Johnson Matthey; † = not re-verified in this session.

| Group | Year | Cathode material | Shape | Dimensions (mm) | Surface area | Anode geometry | Electrolyte | Current density (mA cm⁻²) | Max D/Pd | Duration | Claimed signal & magnitude | Grade | Reference |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Fleischmann, Pons, Hawkins (Utah/Southampton) | 1989 | Pd (JM†) | rods; sheet; cube | rods Ø1, Ø2, Ø4 × 100 (measured on 12.5 mm segments and rescaled†); cube 10×10×10 | Ø4×100 rod ≈ 12.6 cm²; cube 6 cm² | Pt wire helix wound on a glass "cage", concentric with the cathode | 0.1 M LiOD in 99.5 % D₂O / 0.5 % H₂O | 8, 64, 512 (cube ≈ 125†) | n.m. | ~10²–10³ h† | Excess enthalpy 10–111 % of input†, up to ≈ 26 W cm⁻³ (scaled†). Neutrons ≈ 4×10⁴ s⁻¹†. Tritium. Cube "meltdown" anecdote | Heat C; neutrons **N** | JEAC 261 (1989) 301 |
| Fleischmann, Pons, Anderson, Li, Hawkins | 1990 | Pd | rods | Ø1, Ø2, Ø4 × 12.5† | 0.4–1.6 cm² | Pt helix, concentric | 0.1 M LiOD | 64–1024† | n.m. | weeks–months | Excess heat, including bursts. General Electric (GE) reanalysis (Wilson 1992) found smaller but non-zero values | C | JEAC 287 (1990) 293 |
| Pons & Fleischmann | 1993 | Pd | rod | Ø2 × 12.5† | ≈ 0.85 cm² | Pt helix, 4 Dewar cells | 0.1 M LiOD | stepped up to boiling† | n.m. | weeks, then boil-off | ≈ 3.7 kW cm⁻³ during boil-off (≈ 145 W from 0.039 cm³†). "Heat after death" ≈ 3 h | C/D (Morrison critique) | Phys. Lett. A 176 (1993) 118 |
| CEA Grenoble (Lonchampt et al.) | 1996–98 | Pd built to F-P spec | rod | ≈ Ø2 × 12.5† | ≈ 0.85 cm² | Pt helix | 0.1 M LiOD | F-P protocol | n.m. | weeks | Partial reproduction of boil-off excess heat | C | Lonchampt ICCF-6 |
| SRI (McKubre, Crouch-Baker, Tanzella; EPRI-funded) | 1989–94 | Pd wire (JM/Engelhard†) | wire | Ø1 × ~30 (some Ø3)† | ≈ 0.94 cm² | Concentric Pt† | 1.0 M LiOD, some with Al/Si additives† | ~100–1000† | 0.95–1.0 (4-wire resistance) | weeks–months | Typically 5–30 % of input. "~MJ over ~5 days" in best cells | B/C | McKubre ICCF-3; JEAC 368 (1994) 55 |
| SRI "M4" (heat–He) | 1994† | Pd | wire† | † | † | Concentric Pt† | LiOD, closed cell | stepped | > 0.9 | weeks | He-4 correlated with excess heat; ≈ 104 ± 10 %† of 23.85 MeV/He after flushing | C | Hagelstein et al. 2004 |
| IMRA Japan (Kunimatsu et al.) | 1992 | Pd | rod† | † | † | Partially immersed fuel-cell (gas-diffusion) anode that consumes D₂, so pressure change gives loading | LiOD, closed | † | **0.89** | days–weeks | Excess heat only above x ≈ **0.83** | C | Kunimatsu ICCF-3 |
| NHE (Japan, MITI) | 1992–97 | Pd, several suppliers† | rods/plates | various | various | Pt | LiOD, closed | up to ~1000† | ≥ 0.9 in some cells† | months | **No reproducible excess heat**; program terminated 1998 | **N** | Kasagi, country history |
| ENEA Frascati (Violante, Apicella, Sarto, Castagna, Lecci) | 2003–12 | Pd, ENEA-processed from commercial lots | foil | ~0.05 thick; cm-scale strips† | ~1–5 cm²† | Pt, parallel† | 0.1 M LiOD | ~10–500† | ≥ **0.95** (> 90 % of samples) | days–weeks | Excess power in a subset of cathodes, up to ~100 %+ of low input†. He-4 claims | C (B− with SRI/NRL) | Violante ICCF-12; Curr. Sci. 2015 |
| ENEA (Preparata, De Ninno, Del Giudice) | 1997–2004 | Pd | thin wire | Ø0.05 × ~600† | ~1 cm² | Coaxial, plus longitudinal electromigration voltage ("Coehn effect")† | LiOD† | † | ≈ 1 claimed (resistance) | days | Heat plus He-4 | C | ENEA reports† |
| INFN-LNF (Celani et al.) | 1996 | Pd | thin wire | Ø0.05–0.25† | < 1 cm² | Pt | Low-conductivity electrolyte, millisecond high-power pulses | very high peak | ≈ 1 (R/R₀) claimed | days | Excess heat claims | C | OSTI 434378 |
| NAWC China Lake (Miles, Bush et al.) | 1990–95 | Pd (JM), Pd–B (NRL), Pd–Ce | rods/wires | ~Ø1–6†, cm long | ~1–4 cm²† | Pt spiral | 0.2 M LiOD† | ~100–300† | n.m. | days–weeks | Excess 0.05–0.5 W†. He-4 at 10¹⁰–10¹² s⁻¹ W⁻¹. **18/21** runs correlated; **12/12** no-heat runs had no He | C (JEAC); B− with SRI/ENEA He | JEAC 346 (1993) 99; DTIC ADA315020 |
| NRL (Imam) plus Miles | 1997–2001 | Pd–B (~0.1–0.75 wt % B)† | rod | † | † | Pt | LiOD | † | n.m. | weeks | Excess heat in most Pd–B cathodes (≈ 7/8 reported†) | C | NRL/MR/6320-01-8526 |
| Energetics Technologies (Dardik et al.) | 2003–08 | Pd foil | foil | † | † | Pt | 0.1 M LiOD† | "SuperWave" nested sinusoidal modulation | ≥ 0.9† | days–weeks | Up to ~1 MJ in one cell ("Cell 64"†); peak gain ~25×† | C (B− via SRI/ENEA) | OSTI etdeweb 20813957 |
| SPAWAR (Szpak, Mosier-Boss) | 1991–2009 | Pd co-deposited on Au, Ag, Pt, Ni | wire or screen substrate | substrate Ø ~0.25–1† | ~0.1–1 cm² | Pt wire | ~0.03–0.05 M PdCl₂ + 0.3 M LiCl in D₂O† | stepped from ~1 mA to ~100 mA total† | "high" (claimed) | hours–days | IR hot spots, heat, tritium, X-rays, CR-39 pits and "triple tracks" | C/D | JEAC 302 (1991) 255; Naturwiss. 96 (2009) 135 |
| Letts & Cravens | 2003 | Cold-rolled Pd, Au-plated | foil | ~0.15 g Pd | ~0.5 cm²† | Pt | LiOD† | † | high† | days | Excess ~0.1–1 W† switched on by a 30 mW, 670 nm† laser | C | Cravens & Letts ICCF-10 |
| Letts, Cravens, Hagelstein | 2008–12 | Same as above | foil | ~0.15 g | † | Pt | LiOD | † | † | days | Excess-power response at beat frequencies **8.3, 15.3, 20.4 THz**† | C | ACS Sourcebook 2 (2009); arXiv cond-mat/0603213 |
| Storms (LANL/ENECO) | 1993–2007 | Pd sheet and wire, many lots | sheet/wire | † | † | Pt | LiOD, closed with recombiner | ~100–500† | ≥ 0.9 (screened samples) | weeks | Tens of % of input†. Only a minority of samples ever worked† | C | Fusion Technol. 23 (1993) 230† |
| Osaka (Takahashi) | 1992 | Pd | plate | 25 × 25 × 1† | ~12.5 cm² | Pt | LiOD† | high/low pulsed | n.m. | months | 1–10 W claimed† | C | ICCF-3† |
| Caltech (Lewis), Harwell (Williams), MIT (Albagli) | 1989–90 | Pd | rods/wires | various | various | Pt | LiOD | up to ~500 | ≈ 0.75–0.8 where measured† | days–weeks | No heat, neutrons or tritium | **N** (weak: low x, short runs) | Nature 340 (1989) 525; Nature 342 (1989) 375 |
| Utah Physics (Salamon et al.) | 1989–90 | Fleischmann–Pons's own cells | – | – | – | – | – | – | – | ~weeks† | No neutrons or gammas above limits | **N** | Nature 344 (1990) 401 |
| Google program (UBC, MIT, LBNL) | 2015–19 | Pd foils, wires, films, powders | various | various | various | Pt | LiOD and others; polymer and ceramic proton conductors | various | **H/Pd 0.96 ± 0.02** (light H); claimed D/Pd regime not sustained | months | No excess heat, no nuclear products | **N** (regime not reached) | Nature 570 (2019) 45 |
| UBC "Thunderbird" (Chen, Berlinguette et al.) | 2025 | Pd foil membrane target | foil | † | † | D⁺ plasma-immersion ion implantation on one face; electrolysis on the other | D₂O electrolyte† | † | raised by electrolysis | hours | **+15 % D–D neutron rate** when electrochemically loaded; energy out ≪ energy in | A− (one group, *Nature*; beam-target, not LENR) | Nature (2025) s41586-025-09042-7 |

---

## 3. Quantitative relationships

### 3.1 SRI (McKubre) empirical excess-power law

**Original form** (ICCF-3 1992; EPRI TR-104195 1994):

```
P_xs = M · (x − x₀)² · (i − i₀)        for x > x₀ and i > i₀; otherwise 0
```

- P_xs is excess power in W.
- x is the *average* bulk D/Pd, from 4-terminal resistance.
- i is the average cathodic current density in A cm⁻².
- M is a cell- and cathode-specific constant in W cm² A⁻¹. It is not predicted by the model; it absorbs "hidden variables" such as surface state and metallurgy.
- **Fitted constants:** one secondary source quoting McKubre gives **x₀ = 0.832 and i₀ = 0.4 A cm⁻²**. Later SRI summaries put the loading threshold at **x₀ ≈ 0.875–0.90**, and current thresholds in other cells were lower, roughly 0.1–0.25 A cm⁻²†. Conclusion: **x₀ is 0.83–0.90 and i₀ is 0.1–0.4 A cm⁻², and both are cathode-dependent.**

**Flux-extended form** (SRI, ICCF-10 onward; JCMNS "one perspective … SRI"):

```
P_xs = M′ · (x − x₀)² · (i − i₀) · |dx/dt|      (equivalently · |i_D|)
```

- i_D is the net interfacial deuterium flux expressed as a current density.
- SRI's figure "Correlation of excess power and deuterium flux for 1 mm dia. Pd cathode" is the empirical basis. Heat bursts coincided with periods of *changing* loading, such as current steps, rather than with static high loading.
- This is the justification for current modulation: Energetics' SuperWave, Takahashi's high/low pulsing, and SRI's current ramps.

**Initiation criterion** (SRI; restated by Cravens & Letts as their "enabling criteria"). All four conditions must hold at once:

1. high loading;
2. t_init ≥ 10 τ_D, where τ_D ∝ a² (a = radius or half-thickness);
3. i > i₀, where the current threshold is uncorrelated with bulk loading;
4. a non-zero D flux.

SRI noted that in practice loading "needed to be maintained for several weeks to a month."

### 3.2 Diffusion time constants (derived; order of magnitude)

These use D_D ≈ 5×10⁻⁷ cm² s⁻¹ at about 300 K. Dilute-phase tracer values are 4–6×10⁻⁷, from D₀ ≈ 1.7×10⁻³ cm² s⁻¹ and E_a ≈ 0.21 eV†. The time constant is τ ≈ a²/D.

| Cathode | a (cm) | τ_D | 10 τ_D |
|---|---|---|---|
| 50 µm foil (half-thickness 25 µm) | 2.5×10⁻³ | ~12 s | ~2 min |
| Ø1 mm wire (SRI) | 0.05 | ~1.4 h | ~14 h |
| Ø2 mm rod (F-P 1993) | 0.10 | ~5.6 h | ~2.3 d |
| Ø4 mm rod (F-P 1989) | 0.20 | ~22 h | ~9 d |
| 1 cm cube | 0.50 | ~6 d | ~58 d |

**Implication:** the observed incubation periods (weeks) are much longer than 10 τ_D for thin cathodes. So bulk diffusion is **not** what limits initiation. Something slower controls it: evolution of the surface film (Li, Pt or Al deposits), impurity accumulation, defect or vacancy formation, or H/D isotopic exchange. This is a geometry-independent time cost the designer must budget for.

### 3.3 Loading thermodynamics and kinetics

- **Equivalent D₂ fugacity vs absorption overpotential** (Nernst):
  f_D₂ = exp(−2F·η_abs / RT). At 298 K, fugacity rises **10× for every −29.6 mV** of η_abs.
  - x ≈ 0.9 needs equivalent fugacities of order 10³–10⁴ atm†, i.e. η_abs ≈ −90 to −120 mV.
  - x ≥ 0.95 needs more than 10⁴ atm†.
  - η_abs is only the part of the measured overpotential that drives absorption. It is set by the ratio of Volmer (adsorption) to Tafel/Heyrovsky (desorption) rates, not by the total overpotential.
- **Loading vs current density (empirical):** x rises roughly linearly with log i, then saturates at a **cathode-specific plateau of 0.85–0.97** above about 50–200 mA cm⁻²†.
  - The plateau is set by surface recombination kinetics and by leakage through cracks, not by i. Pushing i higher mostly adds D₂ evolution and Joule heat.
  - "Bad" cathodes saturate at 0.7–0.85 regardless of current (Storms; NHE; Google program).
  - W.-S. Zhang's analysis of the maximum loading ratio treats this as a Volmer/Tafel/Heyrovsky balance.
- **Charge needed to load (derived):**
  Q = x · F · n_Pd, with n_Pd = ρ/M = 12.02/106.42 = 0.113 mol cm⁻³.
  For x = 0.9 that is **9.8 kC cm⁻³ ≈ 2.7 A·h cm⁻³**. Cravens & Letts recommend passing **10⁷ C per mole of Pd**, about **115× the stoichiometric charge**. Coulombic loading efficiency near saturation is therefore about 1 % or less.
- **Resistance as a loading gauge:** R/R₀ rises to a maximum of about 1.8–2.0 near x ≈ 0.7–0.75, then *falls* to about 1.5–1.7 as x approaches 0.95–1.0†. The curve is double-valued, so the branch has to be tracked continuously from the start of loading. It is also temperature-sensitive and is raised artificially by cracking.

### 3.4 Calorimetric energy balance

- **Open cell:**
  P_xs = P_out − (E_cell − E_tn)·I − P_heater, with E_tn(D₂O, 25 °C) ≈ **1.527 V** (1.481 V for H₂O).
- **Closed cell with complete recombination:**
  P_xs = P_out − E_cell·I − P_heater.
- **Recombination / Faradaic-efficiency artifact (open cell):**
  ΔP_artifact = (1 − η_F)·E_tn·I.
  Example: I = 0.5 A and η_F = 0.9 gives **76 mW**. At low current density, η_F can fall well below 0.9 (Jones et al. 1995).

### 3.5 Energy scales that decide credibility (derived)

- **Chemical ceiling of the cathode:** D stored in PdD₁.₀ burned to D₂O gives ≈ 0.056 mol D₂ cm⁻³ × 294 kJ mol⁻¹ ≈ **17 kJ cm⁻³**. The hydride formation enthalpy (ΔH ≈ −35 kJ mol⁻¹ D₂†) adds about 2 kJ cm⁻³. Any claim much larger than **~20 kJ per cm³ of Pd** (and much larger than any plausible Li or electrolyte reaction) is non-chemical *if* the calorimetry is right.
- **He-4 yield at 23.85 MeV:** 1 W ↔ **2.62×10¹¹ He s⁻¹**, i.e. 2.62×10¹¹ He J⁻¹. China Lake's 10¹⁰–10¹² s⁻¹ W⁻¹ brackets this value.
- **He-4 accumulation in a closed headspace:**
  C_He ≈ 843 · f_gas · P_xs[W] / V_gas[cm³] ppmv per day, where f_gas is the fraction released from the lattice (SRI M4 suggests about 0.6 before flushing†).
  Example: 1 W into 100 cm³ gives ≈ **8 ppmv per day**, well above air (5.24 ppm) within a day.
- **Separating He-4 from D₂ in a mass spectrometer:** the masses are 4.00260 u and 4.02820 u, so a resolving power m/Δm ≥ **~160** is needed, or D₂ must first be removed with a getter.
- **Expected neutrons if the heat came from conventional D–D fusion:** ≈ **8.5×10¹¹ n s⁻¹ per W** (50 % branch; mean Q ≈ 3.65 MeV). Reported LENR neutron rates are ~10⁻¹–10¹ s⁻¹ above background, a deficit of ≥ 10¹⁰. Hagelstein (2010) used the absence of secondary signatures to argue that any charged products carry less than about 20 keV† each.

---

## 4. Geometry-relevant findings

Each parameter lists the direction of the reported effect and a confidence grade. Grade letters as in §0; "Eng." means established electrochemical engineering or metallurgy, independent of the LENR question.

| # | Parameter | Direction of reported effect | Evidence / notes | Confidence |
|---|---|---|---|---|
| G1 | **Cathode thickness / diameter** | Thinner (50 µm foils, 1 mm wires) reaches high, uniform x faster and more reliably than 4 mm rods or cubes | τ_D ∝ a² (§3.2). Thick cathodes have large core-to-surface loading gradients and more stress cracking. SRI and ENEA moved to thin cathodes for this reason. F-P's 1989 volumetric scaling (larger rods, more W cm⁻³) was never confirmed | Loading: Eng./A. Heat: C |
| G2 | **Current-density uniformity** (anode–cathode symmetry) | Uniform i gives uniform x, and SRI's law depends on *average* x | F-P: concentric Pt helix around a rod. SRI: concentric anode around a wire. Foils should face anodes on both sides. Edges and corners concentrate current and crack first | Eng. (A) for loading; C for heat |
| G3 | **Edges, ends, corners** | Detrimental: current crowding causes local over-evolution, delamination and cracking | The F-P cube is an extreme case. Masking wire ends and foil edges is standard practice† | Eng.; C |
| G4 | **Surface roughness PSD** (ENEA) | Heat-producing cathodes show enhanced PSD at micron to sub-micron lateral scales, roughly 10⁵–10⁷ m⁻¹ spatial frequency; exact band **not verified**† | ENEA ties this to surface-plasmon-polariton (SPP) coupling. For 633 nm light at a Pd/D₂O interface, the SPP wavelength is about 0.45 µm, i.e. k ≈ 1.4×10⁷ m⁻¹, so roughness must have Fourier components at that scale. ENEA's AFM correlation study has not been replicated independently | C |
| G5 | **Grain size** | < 100 µm reported favourable | ENEA SEM/AFM/EBSD work. Very large grains from over-annealing were reported inactive† | C |
| G6 | **Crystallographic texture** | Predominant ⟨100⟩ orientation reported favourable | ENEA | C |
| G7 | **Cracks and voids** | Strongly detrimental to loading: crack surfaces leak D₂. Some theories nonetheless place the "nuclear active environment" in cracks (Storms) | The α→β transition expands the lattice by about 10 % in volume†. Repeated load cycles degrade maximum x | Loading: A. NAE-in-cracks: D |
| G8 | **Thermo-mechanical history** (anneal vs cold-work) | Evidence conflicts. ENEA: cold-roll, anneal at ~850–900 °C† for ~1 h in vacuum, then etch. Others report some cold-worked lots working | Annealing lowers dislocation density and raises grain size, so it trades off against G5 | C/D |
| G9 | **Surface area-to-volume ratio** | Unknown; key open question (surface vs bulk) | ENEA's thin foils, SPAWAR co-deposits and Letts' plated films point to a near-surface effect. F-P 1989 claimed a volume effect | D |
| G10 | **Co-deposited dendritic Pd** (fractal, high area) | Claimed to produce hot spots and nuclear signatures within hours to days, with no incubation | SPAWAR; Letts' modified protocol (lower PdCl₂, hundreds of mA cm⁻²) claimed heat "in all" runs. Substrate choice matters: Au and Ag favoured over Pt and Ni† | C/D |
| G11 | **Anode–cathode spacing and anode area** | Close spacing lowers cell voltage and Joule heat but increases O₂ crossover, i.e. unmeasured recombination at the cathode. Small, high-current-density Pt anodes dissolve and plate Pt onto the cathode | Pt is an excellent Tafel recombination catalyst, so Pt on the cathode lowers x. This effect is well established electrochemically | Eng. (A); heat link C |
| G12 | **Cathode orientation and bubble release** | Vertical cathode with free bubble detachment gives steadier, more uniform surface activity | Bubble shielding and local current redistribution at > 100 mA cm⁻² | Eng. |
| G13 | **Axial voltage drop in long thin wires** | Non-uniform current along the wire unless fed from both ends. ENEA/Preparata turned this into electromigration on purpose | Example: Ø50 µm × 600 mm wire has R of order 30–60 Ω (loaded). Electromigration (Coehn effect) was claimed to push x to about 1 | C |
| G14 | **Headspace volume and cell materials** (for He) | Small, all-metal or He-impermeable volume makes the He signal decisive | Borosilicate glass is He-permeable. SRI's Case replication used a metal-sealed, He-leak-tight vessel. See §3.5 for the ppm/day calculation | A (physics); C for results |
| G15 | **Calorimeter geometry** | Closed cell with recombiner inside a flow or Seebeck calorimeter removes the main open-cell artifacts | Evolution: F-P partially silvered Dewar (isoperibolic), then SRI isothermal mass-flow, then Seebeck envelopes (Storms and others) | A (artifact physics) |
| G16 | **Optical access plus Au overlayer** (laser trigger) | Claimed positive with p-polarised light at specific beat frequencies | Letts/Cravens/Hagelstein; ENEA single-laser work. Requires a window and a thin Au electrodeposit on the cathode | C/D |
| G17 | **Static magnetic field** | Claimed needed for the laser effect (~0.05 T order†) and to affect co-deposition (SPAWAR) | Single groups | D |
| G18 | **Alloying** (B, Ce, Ag, Rh) | Pd–B reported more reproducible than pure Pd (NRL/Miles). Ag shrinks the miscibility gap, so less cracking but a lower maximum x | C | C |
| G19 | **Isotopic purity** | H contamination above about 1 % reported to suppress excess heat. H₂O cells serve as a (claimed) null control | F-P, McKubre, Storms | B |
| G20 | **Temperature** | Higher T raises claimed excess power (F-P positive feedback; Storms) but lowers equilibrium x unless D₂ over-pressure is applied (SRI pressurised cells†) | Direction for heat: C. Direction for loading: A | C |

---

## 5. Known artifacts and failure modes specific to electrochemical cells

1. **Recombination and Faradaic efficiency below 100 % (open cells).**
   - Open-cell accounting assumes all electrolysis gas leaves the cell. Recombination at the cathode (including O₂ crossover reduced there), on the cell walls, or in the headspace makes the cell look like a heat source (§3.4).
   - Jones et al. (J. Phys. Chem. 99, 6966, 1995) measured Faradaic efficiency below 100 %, enough to explain many small claims. Kreysa et al. (JEAC 266, 437, 1989) showed the heat from D₂ + O₂ recombining on Pd.
   - **Mitigation:** closed cell with an internal recombiner, sized and verified for 100 % recombination, with gas-volume monitoring.
2. **Calibration-constant shift** (Shanahan, Thermochim. Acta 387, 95, 2002; later exchanges with Marwan et al., J. Environ. Monit. 2010).
   - A shift of about 1–3 % in the calorimeter constant, for example when the heat-release location moves because of a change in recombination or electrolyte level, can mimic 1–3 % "excess".
   - **Mitigation:** calibrate with a heater located where the cathode is, using heat distributions like the operating ones; recalibrate during runs; use a flow calorimeter whose constant depends on flow rate and c_p, not on heat distribution.
3. **Electrolyte level, top-up and evaporation.** These change the heat-transfer coefficient in isoperibolic Dewar cells. F-P's analysis relies on accurate level modelling.
4. **Boil-off accounting.** Liquid entrainment (droplets, foam) during vigorous boiling means less vaporisation enthalpy than assumed, which overstates excess (Morrison critique of Pons & Fleischmann 1993).
5. **Power measurement with modulated or pulsed current** (SuperWave, pulsing).
   - The mean of the product ⟨V·I⟩ must be sampled at adequate bandwidth, with true-RMS and 4-wire sensing at the cell.
   - Ripple and lead resistance must be accounted for.
   - Any stimulus power (lasers, RF) must be counted.
6. **Stored chemical energy and "heat after death".** Deloading D oxidising in air, or D₂ + O₂ recombining on the drying cathode, can release up to about 20 kJ per cm³ of Pd (§3.5). Claims must exceed this by a wide margin.
7. **Mixing and thermal gradients.** Positional temperature errors in unstirred cells (GE critique, Wilson et al. 1992). Mitigation: stirring, or flow calorimetry.
8. **Loading-measurement artifacts.**
   - R/R₀ is double-valued and temperature-dependent (§3.3).
   - Cracks raise R and imitate changes in loading.
   - Loading is non-uniform between surface and core and along the wire.
   - Gas-pressure loading gauges (Kunimatsu) are sensitive to leaks and to recombination.
9. **Helium contamination.** Air contains 5.24 ppm He. Glass is He-permeable; O-ring seals leak. China Lake's He sat mostly below ambient, so Jones & Hansen (J. Phys. Chem. 99, 6973, 1995) could attribute it to leaks. He is also retained in Pd, which lowers apparent yields (SRI M4 needed flushing).
10. **Tritium.**
    - T is electrolytically enriched from the D₂O inventory (the T/D separation factor is of order 2†).
    - Some Pd stock contained tritium from processing (the Texas A&M-era controversies of 1989–90).
    - Mitigation: measure T in fresh D₂O, in the electrolyte, in recombined gas and in unused cathode material.
11. **Neutron and gamma detectors.**
    - Microphonics and electrical-noise pickup in BF₃ and He-3 tubes; temperature drift; cosmic-ray showers that look like "bursts"; radon.
    - Gamma-line misassignment sank the F-P 1989 neutron claim: the reported 2.22 MeV line had the wrong energy and lacked a Compton edge (Petrasso et al., Nature 339, 183, 1989).
12. **CR-39 in electrolyte contact.** Kowalski (EPJ Appl. Phys. 44, 287, 2008) argued that pits can come from chemical or mechanical attack (D₂, O₂, Cl₂ bubbles; Pd particulates) when the detector touches the cathode and electrolyte.
    - Fair to both sides: "triple tracks" from ¹²C(n,n′)3α need neutrons above about 8 MeV (Q = −7.27 MeV), so they cannot come from 2.45 MeV D–D neutrons. SPAWAR invoked secondary D–T 14.1 MeV neutrons.
    - Mitigation: physical barriers such as 6–25 µm Mylar or PE, H₂O and blank controls, blind etching and counting.
13. **"Transmutation" products on cathodes.** Electrodeposition concentrates ppb-level impurities (Cu, Zn, Fe, Pt, Si, Ca, Mg) from the electrolyte, glass and anode onto the cathode surface. Isotopic anomalies are required before any transmutation claim is taken seriously.
14. **Pt anode dissolution.** Pt re-deposits on the cathode, lowers x over time, and changes the surface catalysis mid-run. This is one likely reason for drifting M in SRI's law and for loss of reproducibility.
15. **Selection effects.** Groups report successful cathodes and runs. Denominators (all cathodes tried) are rarely published. Cravens' "Factors affecting the success rate…" is one of the few attempts to tabulate them.

---

## 6. Design implications, ranked by expected value

Expected value = (effect on the probability of a positive, *credible* result) × (feasibility). The top items raise credibility regardless of whether the effect is real. That is intentional: a null result from a well-built cell is also a valuable, publishable outcome.

1. **Build a closed, He-tight cell with an internal recombiner inside a flow or Seebeck calorimeter.**
   - Use all-metal or glass-free wetted seals where possible (PTFE-lined metal body, metal gaskets) and a small headspace (≤ ~50–100 cm³).
   - Target calorimetric uncertainty ≤ 1 % of input or ≤ 20 mW, whichever is smaller, with in-situ heater calibration at the cathode location.
   - *Why:* this removes artifacts 1, 2, 4, 6 and 9 at once. It is also the only way to test the one signature with a quantitative nuclear prediction (He-4 at about 2.6×10¹¹ He J⁻¹). At 1 W into 100 cm³, He-4 exceeds ambient within about a day (§3.5), which defeats the leak argument.
2. **Use thin, pre-screened, metallurgically controlled Pd cathodes, run several in parallel.**
   - Geometry: 50–100 µm foil strips (about 1 × 2 cm) or Ø0.5–1 mm wires (about 3 cm long). These give τ_D from seconds to about 1 h, so loading is uniform.
   - Process: anneal around 850–900 °C† under vacuum, aiming for grains < 100 µm; etch in a controlled way; characterise before and after by AFM (PSD), EBSD (texture) and ICP-MS (impurities).
   - **Pre-screen** every cathode with an H₂O loading test. Keep only samples that reach H/Pd ≥ 0.95 with little dimensional swelling (the Storms and ENEA criteria).
   - *Why:* the NHE and Google nulls, and ENEA's positive claims, all point to material state as the dominant hidden variable (M in SRI's law). Parallel cells (for example 4–8) give statistics and denominators.
3. **Make current density and loading uniform by geometry.**
   - Wire on the axis of a cylindrical Pt mesh anode, or a foil centred between two parallel Pt meshes.
   - Mask ends and edges with PTFE or heat-shrink.
   - Anode–cathode gap of about 5–15 mm, a compromise between Joule heat and O₂ crossover.
   - Anode area ≥ 3× cathode area, to keep anode potential and Pt dissolution low.
   - Mount the cathode vertically, with free bubble release; feed long wires from both ends.
   - *Why:* SRI's law uses *average* x. Non-uniform x wastes the loading budget and starts cracks at edges.
4. **Verify loading continuously with 4-wire AC resistance** (temperature-compensated; branch tracked from x = 0).
   - Load slowly through the α→β transition at about 5–20 mA cm⁻² to limit cracking, then ramp to 100–500 mA cm⁻².
   - Pass at least 10⁷ C per mol Pd, and hold x ≥ 0.95 for more than 10 τ_D and, in practice, weeks.
   - Record at ≥ 1 Hz so |dx/dt| can be computed.
   - *Why:* this is the best-supported necessary condition (**B**) and lets the team test SRI's equation directly.
5. **Control the electrolyte chemistry.**
   - 0.1–1.0 M LiOD in ≥ 99.9 % D₂O; keep H/(H+D) < 0.5 at % and monitor it.
   - Test ~100–200 ppm Al (and/or Si) as a loading promoter per SRI practice†. Avoid Cl⁻ in bulk-cathode cells.
   - Prefer PTFE or quartz over borosilicate if Si is not deliberately added.
   - *Why:* surface films set the loading plateau (§3.3).
6. **Program deuterium flux.** Once x ≥ 0.95, superimpose current steps, ramps or nested modulation to generate |dx/dt| without dropping mean x below threshold, and measure input power with high-bandwidth ⟨V·I⟩.
   - *Why:* the flux term in SRI's law, plus Energetics and Takahashi claims (**C**). It costs little given item 4.
7. **Operating temperature and pressure.** After loading, test elevated temperature (roughly 50–90 °C) under D₂ over-pressure (a few bar up to tens of bar, if the vessel allows) so loading is not lost.
   - *Why:* claimed positive temperature dependence (**C**), balanced against loading physics (**A**).
8. **Nuclear diagnostics, in priority order:**
   - (a) On-line He-4 with m/Δm ≥ 160 or D₂ getter removal; periodic cathode flushing by reverse polarity.
   - (b) Tritium in D₂O, electrolyte and recombined gas, plus a pre-run cathode assay.
   - (c) Neutron counting (He-3 tubes in moderator plus a liquid scintillator for pulse-shape discrimination), with a matched-background cell.
   - (d) CR-39 *behind* a 6–25 µm Mylar/PE barrier, with H₂O and blank controls, etched and counted blind.
   - *Why:* He-4 is the only product with a quantitative, heat-linked expectation. Radiation searches are cheap add-ons but historically low-yield (§3.5).
9. **Optional, low-cost trigger ports:** an optical window with a thin Au electrodeposit option and tunable dual lasers able to reach beat frequencies of 8–21 THz; and a Helmholtz coil (up to about 0.1 T). Treat these as factorial variables, not the core design.
   - *Why:* low prior (**C/D**) but cheap if designed in from the start.
10. **Co-deposition side-line** (SPAWAR / Letts modified protocol): Pd from dilute PdCl₂ onto Au or Ag wire at high current density.
    - Use it as a fast, cheap screening arm (hours to days, no incubation) inside the same calorimeter and He-tight design, never with CR-39 in contact with the electrolyte.
    - *Why:* it may shortcut the weeks-long initiation, but its credibility is lower (**C/D**).

**What not to prioritise:** large cathodes such as cubes or thick rods (slow, cracked, non-uniform loading); open Dewar cells (artifact-prone); glass-only He sampling; unshielded CR-39 in contact with the electrolyte; neutron detection as the primary diagnostic.

---

## 7. Open questions

1. **What is M?** Which measurable property (surface PSD band, vacancy concentration, impurity or Li/Pt film, grain size or texture, dislocation density) decides whether a cathode above x₀ produces heat? No group has published a blind, predictive pre-run test.
2. **Surface or bulk?** Does excess power scale with area, volume, or something else? This decides whether to optimise for S/V (films, foils) or for Pd mass.
3. **Is D/Pd ≥ 0.95 reproducibly achievable** outside ENEA's process? The Google program managed 0.96 only with light H; NHE reported about 0.9 but no heat.
4. **Superabundant vacancies:** do vacancy-ordered phases (Fukai-type Pd₃VacD₄), which form at high D chemical potential or during co-deposition, play a role in the "initiation time"? This is a hypothesis only.
5. **Heat–He at decisive concentrations:** can any group produce He-4 well above 5.24 ppm in a sealed metal cell, correlated with heat, under blind protocols? This is the single most important experiment.
6. **ENEA's PSD band:** what is the exact spatial-frequency range, and does it survive independent AFM analysis on independently made cathodes? (The band was not verified in this session.)
7. **Laser and THz triggering:** do the 8.3, 15.3 and 20.4 THz beat-frequency responses replicate outside the Letts lab, and do they line up with measured PdD phonon dispersion at x ≈ 0.95?
8. **The flux term:** is |dx/dt| causal, or a marker of surface restructuring? A controlled test would compare constant x with modulated i against constant i with modulated x.
9. **Denominators:** what fraction of *all* cathodes tried by SRI, ENEA, Energetics and China Lake ever produced heat? Published figures are selective.
10. **ARPA-E LENR program (2023–):** results relevant to electrochemical Pd–D had not surfaced in accessible sources at the time of writing. The eight teams were mostly gas-loading, nanoparticle and diagnostics projects.

---

## 8. References

URLs marked [S] were confirmed to exist by the search engine in this session; their contents were not fully read. Items without a URL are standard citations recalled from the primary literature.

**Fleischmann–Pons and critiques**
1. M. Fleischmann, S. Pons, M. Hawkins, "Electrochemically induced nuclear fusion of deuterium," *J. Electroanal. Chem.* 261 (1989) 301; erratum 263 (1989) 187.
2. M. Fleischmann, S. Pons, M.W. Anderson, L.J. Li, M. Hawkins, "Calorimetry of the palladium–deuterium–heavy water system," *J. Electroanal. Chem.* 287 (1990) 293–348.
3. S. Pons, M. Fleischmann, "Calorimetry of the Pd–D₂O system: from simplicity via complications to simplicity," *Phys. Lett. A* 176 (1993) 118.
4. J. Rothwell, "Review of the calorimetry of Fleischmann and Pons," LENR-CANR. https://lenr-canr.org/acrobat/RothwellJreviewofth.pdf [S]
5. G. Lonchampt, L. Bonnetain, P. Hieter, "Reproduction of Fleischmann and Pons experiments," ICCF-6 (1996). https://www.lenr-canr.org/acrobat/LonchamptGreproducti.pdf [S]
6. R.H. Wilson et al. (GE), "Analysis of experiments on the calorimetry of LiOD–D₂O electrochemical cells," *J. Electroanal. Chem.* 332 (1992) 1.
7. D.R.O. Morrison, "Comments on claims of excess enthalpy by Fleischmann and Pons using simple cells made to boil," *Phys. Lett. A* 185 (1994) 498.
8. R.D. Petrasso et al., "Problems with the γ-ray spectrum in the Fleischmann et al. experiments," *Nature* 339 (1989) 183.
9. N.S. Lewis et al. (Caltech), *Nature* 340 (1989) 525; D.E. Williams et al. (Harwell), *Nature* 342 (1989) 375; D. Albagli et al. (MIT), *J. Fusion Energy* 9 (1990) 133.
10. M.H. Salamon et al., "Limits on the emission of neutrons, γ-rays, electrons and protons from Pons/Fleischmann electrolytic cells," *Nature* 344 (1990) 401.
11. G. Kreysa, G. Marx, W. Plieth, "A critical analysis of electrochemical nuclear fusion experiments," *J. Electroanal. Chem.* 266 (1989) 437.
12. J.E. Jones, L.D. Hansen, S.E. Jones, D.S. Shelton, J.M. Thorne, "Faradaic efficiencies less than 100 % during electrolysis of water can account for reports of excess heat in 'cold fusion' cells," *J. Phys. Chem.* 99 (1995) 6966.
13. S.E. Jones, L.D. Hansen, "Examination of claims of Miles et al. in Pons–Fleischmann-type cold fusion experiments," *J. Phys. Chem.* 99 (1995) 6973.
14. K.L. Shanahan, "A systematic error in mass flow calorimetry demonstrated," *Thermochim. Acta* 387 (2002) 95.
15. J. Marwan, M.C.H. McKubre, F.L. Tanzella, P.L. Hagelstein, M.H. Miles, M.R. Swartz, E. Storms, Y. Iwamura, P.A. Mosier-Boss, L.P.G. Forsley, "A new look at low-energy nuclear reaction (LENR) research: a response to Shanahan," *J. Environ. Monit.* 12 (2010) 1765. https://www.lenr-canr.org/acrobat/MarwanJanewlookat.pdf [S]

**SRI / McKubre**
16. M.C.H. McKubre, S. Crouch-Baker, A.M. Riley, S.I. Smedley, F.L. Tanzella, "Excess power observations in electrochemical studies of the D/Pd system; the influence of loading," *Frontiers of Cold Fusion* (ICCF-3, 1992), p. 5. https://lenr-canr.org/acrobat/McKubreMCHexcesspowe.pdf [S]
17. M.C.H. McKubre et al., "Isothermal flow calorimetric investigations of the D/Pd and H/Pd systems," *J. Electroanal. Chem.* 368 (1994) 55.
18. M.C.H. McKubre et al., "Development of advanced concepts for nuclear processes in deuterated metals," EPRI TR-104195 (1994).
19. M.C.H. McKubre et al., "The need for triggering in cold fusion reactions," ICCF-10 (2003). https://www.researchgate.net/publication/241489694_The_Need_for_Triggering_in_Cold_Fusion_Reactions [S]
20. M.C.H. McKubre, "Cold fusion (LENR): one perspective on the state of the science," ICCF-15 (2009). https://lenr-canr.org/acrobat/McKubreMCHcoldfusionb.pdf [S]
21. M.C.H. McKubre, "Cold Fusion, LENR, CMNS, FPE: one perspective on the state of the science based on measurements made at SRI," *J. Condensed Matter Nucl. Sci.* 4 (2011) 17. https://jcmns.org/article/72124-cold-fusion-lenr-cmns-fpe-one-perspective-on-the-state-of-the-science-based-on-measurements-made-at-sri/attachment/150110.pdf [S]
22. M.C.H. McKubre, "Cold fusion: comments on the state of scientific proof," *Curr. Sci.* 108 (2015) 495. https://www.currentscience.ac.in/Volumes/108/04/0495.pdf [S]
23. P.L. Hagelstein, M.C.H. McKubre, D.J. Nagel, T.A. Chubb, R.J. Hekman, "New physical effects in metal deuterides" (submission to the 2004 DOE review). https://www.researchgate.net/publication/255610539_New_physical_effects_in_metal_deuterides [S]
24. U.S. DOE, "Report of the Review of Low Energy Nuclear Reactions" (Dec. 2004). https://www.lenr-canr.org/acrobat/DOEreportofth.pdf [S]

**IMRA / NHE / Japan**
25. K. Kunimatsu et al., "Deuterium loading ratio and excess heat generation during electrolysis of heavy water by a palladium cathode in a closed cell using a partially immersed fuel cell anode," ICCF-3 (1992). https://lenr-canr.org/acrobat/KunimatsuKdeuteriuml.pdf [S]
26. J. Kasagi et al., "Country history of Japanese work on cold fusion." https://lenr-canr.org/acrobat/KasagiJcountryhis.pdf [S]
27. W.-S. Zhang, "The maximum hydrogen (deuterium) loading ratio in the Pd|H₂O (D₂O) electrochemical system." https://www.lenr-canr.org/acrobat/ZhangWSthemaximum.pdf [S]
28. "Effects of D/Pd ratio and cathode pretreatments on excess heat in closed Pd|D₂O+D₂SO₄ electrolytic cells," *J. Condensed Matter Nucl. Sci.* https://jcmns.org/article/72434-effects-of-d-pd-ratio-and-cathode-pretreatments-on-excess-heat-in-closed-pd-d2o-d2so4-electrolytic-cells.pdf [S]

**ENEA / INFN**
29. V. Violante et al., "Progress in excess of power experiments with electrochemical loading of deuterium in palladium," ICCF-12 (2006) 55. https://ui.adsabs.harvard.edu/abs/2006cmns...12...55V/abstract [S]
30. V. Violante et al., "Reproducibility of excess of power and evidence of ⁴He in palladium foils loaded with deuterium." https://www.researchgate.net/publication/253408661 [S]
31. V. Violante et al., "Review of materials science for studying the Fleischmann and Pons effect," *Curr. Sci.* 108 (2015) 540†.
32. F. Celani et al., "Deuterium overloading of palladium wires by means of high power ms pulsed electrolysis and electromigration: suggestions of a 'phase transition' and related excess heat" (1996). https://www.osti.gov/etdeweb/biblio/434378 [S]

**China Lake / NRL / helium**
33. B.F. Bush, J.J. Lagowski, M.H. Miles, G.S. Ostrom, "Helium production during the electrolysis of D₂O in cold fusion experiments," *J. Electroanal. Chem.* 304 (1991) 271.
34. M.H. Miles, R.A. Hollins, B.F. Bush, J.J. Lagowski, R.E. Miles, "Correlation of excess power and helium production during D₂O and H₂O electrolysis using palladium cathodes," *J. Electroanal. Chem.* 346 (1993) 99.
35. M.H. Miles, "Anomalous effects in deuterated systems," NAWCWPNS TP 8302 (1996). https://apps.dtic.mil/sti/pdfs/ADA315020.pdf [S]
36. M.H. Miles, "Correlation of excess enthalpy and helium-4 production: a review," ICCF-10. https://lenr-canr.org/acrobat/MilesMcorrelatioa.pdf [S]
37. M.H. Miles et al., "Consistency of helium production with the excess power…" (2024). https://www.sciencedirect.com/science/article/am/pii/S1572665724007641 [S]; https://lenr-canr.org/acrobat/MilesMconsistenc.pdf [S]
38. M.H. Miles, M. Fleischmann, M.A. Imam, "Calorimetric analysis of a heavy water electrolysis experiment using a Pd–B alloy cathode," NRL/MR/6320-01-8526 (2001).
39. D.J. Nagel, P.L. Hagelstein et al. on SRI Case heat–He; see refs. 20 and 23.

**Energetics / SuperWave**
40. I. Dardik et al., "Excess heat in electrolysis experiments at Energetics Technologies," ICCF-11 (2004). https://www.osti.gov/etdeweb/biblio/20813957 [S]
41. I. Dardik et al., "Ultrasonically-excited electrolysis experiments at Energetics Technologies." https://www.lenr-canr.org/acrobat/DardikIultrasonic.pdf [S]
42. M.C.H. McKubre, F.L. Tanzella et al., SRI replications of Energetics SuperWave cells, ICCF-13/14 (2007–08)†.

**SPAWAR co-deposition**
43. S. Szpak, P.A. Mosier-Boss, J.J. Smith, "On the behavior of Pd deposited in the presence of evolving deuterium," *J. Electroanal. Chem.* 302 (1991) 255.
44. S. Szpak, P.A. Mosier-Boss, M.H. Miles, M. Fleischmann, "Thermal behavior of polarized Pd/D electrodes prepared by co-deposition," *Thermochim. Acta* 410 (2004) 101.
45. P.A. Mosier-Boss et al., "Use of CR-39 in Pd/D co-deposition experiments," *Eur. Phys. J. Appl. Phys.* 40 (2007) 293.
46. P.A. Mosier-Boss et al., "Triple tracks in CR-39 as the result of Pd–D co-deposition: evidence of energetic neutrons," *Naturwissenschaften* 96 (2009) 135.
47. L. Kowalski, "Comment on 'Use of CR-39 in Pd/D co-deposition experiments'," *Eur. Phys. J. Appl. Phys.* 44 (2008) 287; reply by Mosier-Boss et al., ibid. 291.

**Letts / Cravens / Hagelstein / Storms**
48. D. Cravens, D. Letts, "The enabling criteria of electrochemical heat: beyond reasonable doubt," ICCF-10 (2003). https://www.lenr-canr.org/acrobat/CravensDtheenablin.pdf [S]
49. D. Cravens, "Factors affecting the success rate of heat generation in CF cells." https://www.lenr-canr.org/acrobat/CravensDfactorsaff.pdf [S]
50. D. Letts, D. Cravens, P.L. Hagelstein, "Dual laser stimulation and optical phonons in palladium deuteride," in *Low-Energy Nuclear Reactions and New Energy Technologies Sourcebook* Vol. 2, ACS Symp. Ser. 1029 (2009)†.
51. P.L. Hagelstein et al., "On the laser stimulation of low-energy nuclear reactions in deuterated palladium," arXiv:cond-mat/0603213. https://arxiv.org/pdf/cond-mat/0603213 [S]
52. D. Letts, P.L. Hagelstein, "Modified Szpak protocol for excess heat," *J. Condensed Matter Nucl. Sci.* 6 (2012) 44†.
53. P.L. Hagelstein, "Constraints on energetic particles in the Fleischmann–Pons experiment," *Naturwissenschaften* 97 (2010) 345.
54. E. Storms, "Status of cold fusion (2010)," *Naturwissenschaften* 97 (2010) 861.
55. E. Storms, *The Science of Low Energy Nuclear Reaction* (World Scientific, 2007).
56. D.E. Hubler, "Anomalous effects in hydrogen-charged palladium — a review," *Surf. Coat. Technol.* 201 (2007) 8568†.

**Google-funded program and follow-ups**
57. C.P. Berlinguette, Y.-M. Chiang, J.N. Munday, T. Schenkel, D.K. Fork, R. Koningstein, M.D. Trevithick, "Revisiting the cold case of cold fusion," *Nature* 570 (2019) 45. https://www.nature.com/articles/s41586-019-1256-6 [S]
58. Nature news, "Google revives controversial cold-fusion experiments" (2019). https://www.nature.com/articles/d41586-019-01683-9 [S]
59. "Coming in from the cold," *Nature Materials* editorial (2019). https://www.nature.com/articles/s41563-019-0530-1 [S]
60. "Materials advances result from study of cold fusion," *MRS Bulletin*. https://www.cambridge.org/core/journals/mrs-bulletin/article/materials-advances-result-from-study-of-cold-fusion/6DF79E2F5317AE9DCA6AE9284132F288 [S]
61. J.D. Benck et al., "Producing high concentrations of hydrogen in palladium via electrochemical insertion from aqueous and solid electrolytes," *Chem. Mater.* 31 (2019) 4234†.
62. K. Chen, C.P. Berlinguette et al., "Electrochemical loading enhances deuterium fusion rates in a metal target," *Nature* (2025). https://www.nature.com/articles/s41586-025-09042-7 [S]; commentary: https://physicsworld.com/a/electrochemical-loading-boosts-deuterium-fusion-in-a-palladium-target/ [S]; https://www.chemistryworld.com/news/electrochemistry-offers-modest-boost-to-deuterium-fusion-reaction/4022070.article [S]
63. ARPA-E LENR Exploratory Topic (2023, $10 M, 8 teams). https://www.greencarcongress.com/2023/02/20230218-lenr.html [S]

**Reviews and context**
64. "A review of experiments reporting non-conventional phenomena in nuclear matter aiming at identifying common features in view of possible interpretation," arXiv:2308.13533. https://arxiv.org/pdf/2308.13533 [S]
65. "Gatekeeping: a partial history of cold fusion," arXiv:2601.09996 (2026). https://arxiv.org/pdf/2601.09996 [S]
66. V.A. Chechin et al., "Critical review of theoretical models for anomalous effects in deuterated metals," arXiv:nucl-th/0303057. https://arxiv.org/pdf/nucl-th/0303057 [S]
67. Y. Fukai, N. Okuma, "Formation of superabundant vacancies in Pd hydride under high hydrogen pressures," *Phys. Rev. Lett.* 73 (1994) 1640.
