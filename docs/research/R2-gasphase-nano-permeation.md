# R2: Gas-phase loading, nanostructured materials, and deuterium permeation/flux systems

*Research note for the LENR geometric-configuration design project. Compiled 2026-09-29.*

**How this was sourced, and its limits.** Web search worked. Direct retrieval of primary PDFs did not: the network egress policy blocked lenr-canr.org, arXiv, IOPscience, Frontiers, J-STAGE, OSTI and jcfrs.org. So the numbers here come from two places: (a) abstracts and snippets returned by search, which usually quote the paper's own abstract, and (b) the author's prior knowledge of the literature. Values marked **(v)** come from (b) or from secondary summaries. Check them against the primary PDF before using them in a design calculation. Unmarked values were seen in search-returned text from the primary source or its abstract.

**Evidence grades.** A = independently replicated in mainstream peer-reviewed venues. B = replicated by several LENR groups but not accepted by the mainstream. C = a single group or collaboration. D = anecdotal or disputed. N = a credible null result against the claim.

---

## 1. Executive summary

1. **Nanoscale hydride thermodynamics is well established, and it cuts against "super-loading" claims. [A]** In Pd particles smaller than about 3–5 nm, the α/β miscibility gap narrows or disappears at room temperature and the plateaus slope. The bulk-like H capacity *falls* as size shrinks (roughly 0.3–0.5 H/Pd at 1 bar for particles of about 2–3 nm). Only surface and subsurface sites add capacity, in proportion to the surface fraction. So reported D/Pd ≈ 1.1–3 in ZrO₂–nano-Pd (Arata; Kitamura 2009) is best explained by support, oxide or Zr side reactions, not by real lattice super-loading. **[claims: C/D]**
2. **Superabundant vacancies (Fukai) are real. Any link to nuclear effects is conjecture. [A for SAV; D for relevance]** Vacancy fractions of about 10–30 at.% form in Pd and Ni under GPa-level H pressure above ~600 °C, and near surfaces during electrodeposition or cathodic charging. Each vacancy traps up to 6 H in its neighbouring octahedral sites. The H–H spacing in such a cluster stays around 2.5 Å, far larger than the separation fusion would need.
3. **Pd/CaO multilayer D₂-permeation transmutation (Iwamura, MHI) is the best-replicated gas-phase nuclear-ash claim, and it is contested. [B, contested]** The claimed reactions are Cs→Pr, Sr→Mo and Ba→Sm. Pr reached about 10¹⁴ cm⁻² (v) at MHI and rose from ≤2×10¹¹ to 1.6×10¹² cm⁻² in Toyota's ICP-MS replication (JJAP 2013). H₂ controls stayed null, and Pr scaled with D flux. NRL's replication failed, and NRL alleged Pr contamination in the MHI lab, but that critique was never published as a paper.
4. **Ni–Cu nanolayer films on Ni (Tohoku/Clean Planet) are the most carefully instrumented recent heat claim, peer-reviewed in JJAP 2024. [C]** The structure is 6×[Cu 2 nm / Ni 14 nm] on Ni plates 25×25×0.1 mm. The samples are H₂-loaded at about 270 °C, then evacuated and heated. Claimed excess is 1–6 W, at >10 keV per absorbed H (JJAP) and ≥410 ± 108 keV/H (radiant calorimetry), with no neutrons or gammas. **No lab outside the Tohoku–Clean Planet–Kobe network has published a replication.**
5. **Ni-based nanocomposite powders (PdNi/ZrO₂, CuNi/ZrO₂) give heat at 200–300 °C in both H₂ and D₂. [B]** Claimed: 5–25 W per 100–200 g, about 1.3–1.5 W/g-Ni, 15 eV to 2.1 keV per H/D. This was reproduced qualitatively between Kobe and Tohoku under NEDO (2015–2017). Similar results come from Biberian (CleanHME, Ni alloys on alumina) and Beiting (Aerospace Corp., ZrO₂–Ni–Pd), but it is not accepted outside the field. Because H and D perform alike, D+D fusion is not what these results point to.
6. **A pattern recurs across groups: the signal appears in the desorption or flux phase, not during static loading. [C, trending B]** Examples: Iwamura's films heat up after evacuation plus a temperature ramp. Biberian sees more heat while pumping H out. Nissan's DSC shows sustained heat at desorption temperatures (300–500 °C) but not at absorption temperatures. Pr yield scales with D flux. NTT (1990) reported out-diffusion bursts. **Flux, not static loading ratio, is the best control variable for gas-phase designs.**
7. **Arata–Zhang (DS-cathode; ZrO₂–nano-Pd, He-4): [C/D].** Kitamura partly reproduced the absorption anomaly. The small "post-saturation" heat is plausibly chemical: PdO reduction, H/D exchange with support hydroxyls, or Zr hydriding. He-4 was never measured with a leak-tight, high-resolution, independently run protocol.
8. **Mizuno's Pd-rubbed Ni mesh: [D].** Claimed 250 W to 3 kW excess, measured by air-flow calorimetry. Replication attempts in 2019–2020 gave small, transient or no excess heat.
9. **Celani constantan wires and Brillouin Q-pulse tubes: [C/D]. Parkhomov and Rossi's E-Cat: [D/N].** Celani's claimed excess is 5–10 W. Brillouin's COP is 1.2–1.6 in SRI's contracted tests. Rossi's 2014 Lugano COP of about 3.5 depended on using an alumina emissivity value wrong for the camera's band; corrected analyses give COP ≈ 1 (v). Biberian ran more than 20 Parkhomov-type tests and saw no excess within ±2 W. **Rossi's claims are not credible and should carry zero design weight.**
10. **Credible or neutral null context.** The Google-funded programme (Berlinguette et al., *Nature* 2019) found no anomalous heat in the Pd–H configurations it tested **[N, for those configurations (v)]**. We found no public positive heat claim from the EU HERMES project (Pd–D, synchrotron and EPR diagnostics) **[weak N; absence of evidence]**. CleanHME (2020–2025) claims multi-lab excess heat and particle emission. Its main output in mainstream journals, though, is accelerator electron-screening physics (Czerski et al., screening energy ≈340 eV in ZrD₂), which is a beam experiment and does not show a heat source.

---

## 2. Table of experiments

Abbreviations: RT = room temperature; RC = reaction chamber; sccm = standard cm³/min; xs = excess.

| # | Group, year | Material / structure | Geometry / dimensions | Gas & pressure | Temp. | Flux | Loading | Duration | Claimed signal & magnitude | Grade | Ref |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | Arata & Zhang (Osaka), 1997–2003 | "DS-cathode": Pd-black sealed inside a Pd shell; D₂O electrolysis outside, D permeates inward | Closed Pd cylinder, O.D. ~20 mm, wall ~2 mm (v); Pd-black ~1–4 g, tens-of-nm grains (v) | Internal D₂ builds from permeation; high internal pressure claimed (v) | ~RT–80 °C | Set by electrolysis current; not quantified | "Spill-over" high D/Pd claimed | 10²–10³ h | W-scale xs over months; ⁴He and ³He in Pd-black by high-resolution QMS (v) | C/D | [R1][R2] |
| 2 | SRI (McKubre et al.) analysis of Arata DS cathodes, ~2000–2003 | Arata-supplied DS cathodes | As #1 | As #1 | As #1 | n/a | n/a | n/a | ³He/⁴He in internal gas; ³He partly attributed to tritium decay (v) | C | [R3] |
| 3 | Arata & Zhang, 2005–2008 (public demonstration 22 May 2008) | ZrO₂–nano-Pd from oxidized amorphous Zr₆₅Pd₃₅ ribbon | Pd particles ~5 nm in ZrO₂ matrix (v); ~7 g powder (v) | D₂ vs H₂ admitted to evacuated cell; atm to tens of atm (v) | Rose to ~70 °C during uptake (v) | Static | D/Pd ≈ 3 claimed (v) | ~50 h+ | Inner cell stays ~1–2 °C above wall for ≥50 h with D₂, not H₂ (v); ⁴He claimed | C/D | [R1][R4] |
| 4 | Kitamura, Takahashi et al. (Kobe), 2009 | Pd·ZrO₂ nanopowder (commercial) | Pd ~10 nm (v); a later summary says ~100 nm | D₂ or H₂ up to ~1 MPa (v); twin "A1/A2" flow calorimeters (v) | RT | Static | **D/Pd = 1.1, H/Pd = 1.1 ± 0.3** | Hours to days | Absorption energy **2.4 ± 0.2 eV/D, 1.8 ± 0.4 eV/H**; "significantly positive" heat after saturation with D₂ only | C (reproduces Arata's absorption anomaly; heat marginal) | [R5] |
| 5 | Kidwell et al. (NRL), 2009–2015 | Pd nanoparticles on oxide supports | nm-scale Pd on Al₂O₃/ZrO₂ (v) | D₂ vs H₂ pressurization | RT–moderate | Static | n/a | Hours | Larger heat with D₂ than H₂; patented (US 9,182,365). H/D exchange with support hydroxyls flagged as a chemical explanation (v) | C/D | [R6] |
| 6 | Kitamura, Takahashi A. & K., Technova, Nissan (Kobe "C1" system), 2012–2018 | PNZ (Pd₁Ni₇₋₁₀/ZrO₂), CNZ (Cu₁Ni₇/ZrO₂), CNS (CuNi/SiO₂) | Ni-alloy nanoparticles ~several–20 nm (v) in ZrO₂; powder grains <500 µm with ZrO₂-rich surface; **100–200 g**; **RC 500 cm³** | H₂ or D₂, ~0.1–1 MPa (v) | **RC 200–300 °C** (up to ~450 °C) | Static plus temperature cycling | Low (Ni-based) | Weeks per run | **5–25 W** (PNZ6: 120 g, 25 W); IJHE 2018: ΔT 15–16 °C ≈ 11–12 W ≈ **1.3–1.5 W/g-Ni**; oil-flow calorimetry | B | [R7][R8][R9] |
| 7 | NEDO "MHE" project (Tohoku, Kobe, Nagoya, Kyushu, Nissan, Technova), 2015–2017 | PNZ, CNZ from one shared preparation route | As #6; RC 500 cm³ at both Tohoku and Kobe | H₂/D₂ | 200–400 °C | Static | n/a | Weeks | Qualitative Kobe–Tohoku reproducibility "good"; **15 eV to 2.1 keV per H/D**; one secondary source says 11 of 11 runs gave net heat | B | [R9][R10] |
| 8 | Nissan (with Tohoku), 2019–2023 | PNZ alloy | mg-scale DSC samples | H₂ | **Sustained heat at 300–500 °C** | Desorption regime | n/a | Hours to days | Heat appears at desorption temperatures, **not** at absorption temperatures | C | [R11] |
| 9 | Tohoku et al., 2024 (arXiv) | CuNi/ZrO₂ nanocomposites | As #6 | H₂ | Elevated | n/a | n/a | n/a | ³He detected after runs; high mass resolution required (see §5) | C | [R12] |
| 10 | Biberian, Valat, Sjöberg et al. (CleanHME), 2022–2024 | Ni and Ni–Cu nanoparticles from hydrotalcite precursors on amorphous Al₂O₃; also CNZ | Powder; Calvet (Seebeck) calorimeter calibrated with alumina in H₂ up to 850 °C | H₂ | Up to ~850 °C | **More heat on pumping out H₂** | n/a | Days | Excess "100–1000×" the expected chemical heat; magnitude of order W (v) | C (B together with #6/#7) | [R13][R14] |
| 11 | Beiting (Aerospace Corp.), 2017–2018 | ZrO₂–Ni–Pd particles (Zr 50.5, Ni 24.9, Pd 2.03, O 22.0 wt%) | **Identical test and reference cells**, 10 cm³ powder each, 2 thermocouples per cell | H₂ (v) | ~300 °C (v) | Static | n/a | **950 h** | "Heat not chemical in origin"; with 10 g added Sm₂Co₁₇, **3.5 MJ ± 6% in 950 h (≈1 W mean)** | C | [R15] |
| 12 | Iwamura, Sakano, Itoh (MHI), 2002–2014 | Pd/[CaO/Pd]×5 multilayer on bulk Pd; Cs, Sr or Ba deposited on the D₂-facing surface | Pd plate 25×25×0.1 mm (v); top Pd 40 nm, CaO 2 nm, Pd 18 nm ×5 (v) | **D₂ ~1 atm upstream, vacuum downstream** | ~70 °C (v) | **2–3 sccm**, ≈1–5×10¹⁷ D cm⁻² s⁻¹ (derived, §3) | Low (α/β Pd) | ~1 week | Cs→Pr up to ~10¹⁴ cm⁻² (v); Sr→Mo with ⁹⁶Mo enrichment; Ba→Sm. **Pr ∝ D flux.** No Pr without CaO or with H₂ (v) | B (contested) | [R16][R17][R18][R19] |
| 13 | Hioki et al. (Toyota Central R&D), 2013 | Cs-ion-implanted Pd/CaO multilayer | Iwamura-type | D₂ vs H₂ permeation | As #12 | Permeation | n/a | 3 D₂ runs | **Pr ≤2.0×10¹¹ → 1.6×10¹² cm⁻²** (ICP-MS, detection limit ~10¹⁰ cm⁻²); **no rise with H₂** | Supports #12 | [R20] |
| 14 | Grabowski, Kidwell et al. (NRL), 2008–2009 | MHI structures, joint effort | As #12 | D₂ | As #12 | Permeation | n/a | n/a | **Failed to replicate.** Alleged Pr contamination in the MHI lab. Unpublished (oral only, ICCF-15) | N (weak, unpublished) | [R21][R22] |
| 15 | Iwamura, Itoh, Kasagi, Clean Planet (Tohoku CLEAR), 2017–2024 | Ni-based multilayer: Cu/Ni, sometimes with CaO or Y₂O₃ layers (v) | **Bulk Ni 25×25×0.1 mm; 6 × [Cu 2 nm / Ni 14 nm]** (≈96 nm stack); samples sit on a ceramic heater in a vacuum chamber (v) | H₂ loading at sub-atmospheric pressure (v), then evacuation | Load at ~250–270 °C for ~16 h; then rapid heating to ~500–900 °C (v) | Out-flux during heating (est. 10¹³–10¹⁵ H cm⁻² s⁻¹, §3) | Small (Ni) | Days to weeks | **>10 keV per absorbed H (JJAP 2024); no γ or n** | C | [R23][R24] |
| 16 | Kasagi, Itoh, Iwamura et al. (Tohoku/Clean Planet), 2023–2025 | NiCu multilayer vs Ni-only film | As #15 | As #15 | As #15 | Desorption | n/a | 80 h; ≥215 h | Photon calorimetry (0.3–5.5 µm, 3 detectors): **4–6 W; 460 ± 120 kJ in 80 h; ≥410 ± 108 keV/H**. 2025: **~1.1 W** NiCu, Ni-only slightly less; decay time >2000 h | C | [R25][R26] |
| 17 | ICCF-26 report (Iwamura group), 2025 | Cu/Ni multilayer | IR imaging | As #15 | — | — | — | — | **Micron-scale hot spots; local temperature claimed above Ni's melting point** | C/D | [R27] |
| 18 | Mizuno (Hydrogen Engineering Application & Development Co.), 2017 | Ni mesh cathode activated by D₂ glow discharge; Pd rod anode wrapped in Pd wire | Mesh wound inside a cylindrical steel reactor (v) | **D₂ 100–300 Pa** | >200 °C | Static (plus discharge) | n/a | 35 days | **80 W in → 78 W xs; 108 MJ** (air-flow calorimetry) | D | [R28] |
| 19 | Mizuno & Rothwell, 2019 (ICCF-22) | Ni mesh with Pd deposited by rubbing ("burnishing") | Mesh area ~10² cm² (v); mg-scale Pd (v) | D₂ ~10²–10³ Pa (v) | 200–450 °C (v) | Static | n/a | Days | **50 W in → ~250 W xs; 300 W in → 2–3 kW** | D | [R29] |
| 20 | Mizuno-type replications (Hokkaido Univ. of Science, Zhang, others), 2019–2020 | As #19 | As #19 | As #19 | As #19 | — | — | Hours | "Not going well"; small xs that fades after 2–3 h; one positive report from a university close to Mizuno | D/N (weak) | [R30][R31] |
| 21 | Celani et al. (INFN-LNF), 2011–2020 | Constantan Cu₅₅Ni₄₄Mn₁; surface "flash-oxidised" (up to 20 kVA/g pulses) into sub-µm Cu/Ni oxide texture; later coatings and Capuchin knots | **Wire Ø 200 µm, length ~100–160 cm** | H₂/D₂, often mixed with Ar/Xe, ~0.1–1 MPa (v) | Effect rises above ~400 °C | DC/pulsed longitudinal current (electromigration) | Low | Days to weeks | 5–10+ W xs (v); anomalous resistance behaviour | C/D | [R32][R33] |
| 22 | MFMP replication of Celani, 2012–2013 | Celani-supplied wire | As #21 | H₂ | As #21 | — | — | 48 h | Initially **6 W (12.5% of input)**, later **5.3%**; not established as robust | D | [R34] |
| 23 | Brillouin Energy / SRI (Tanzella), 2016–2018 | "Hydrogen Hot Tube" with Ni-based core, driven by high-voltage "Q-pulses" | Tubular core (v) | H₂, several bar (v) | ~300–600 °C (v) | Pulsed | n/a | Weeks | **COP 1.2–1.6**, one core 1.91–2.08; several W | C/D | [R35][R36] |
| 24 | Parkhomov, 2014–2016; replications | Ni powder + LiAlH₄ in alumina tube | ~1 g charge (v) | H₂ from LiAlH₄ | 1100–1300 °C | — | — | Hours to days | COP ~2.5 claimed (water-evaporation calorimetry, v). **Biberian: >20 runs, no xs within ±2 W** | D/N | [R37] |
| 25 | Rossi E-Cat; Levi et al. "Lugano", 2013–2014 | Ni + LiAlH₄ in alumina reactor | ~1 g fuel (v) | — | ~1200–1400 °C claimed | — | — | 32 days | COP ~3.2–3.6 from IR thermography; Ni-62/Li isotopic shifts in samples Rossi supplied (v) | D/N | [R38][R39] |
| 26 | Li X.Z. et al. (Tsinghua), 2003 | Pd tube, D₂ inside, pumped vacuum outside | Thin-walled Pd tube (v) | D₂ ~1 atm | ~RT–moderate | Controlled permeation | α/β | Days | Heat flow correlated with "abnormal" D flux; sub-W to W (v) | C/D | [R40] |
| 27 | Yamaguchi & Nishioka (NTT), 1990 | Pd plate, one face Au-coated (diffusion barrier), D-loaded, then out-diffused in vacuum | mm-thick Pd plate (v) | D₂ loading, then vacuum | RT | **Driven out-diffusion** | High (v) | Hours | Heat bursts, neutron bursts; later ⁴He claims (v) | D | [R41] |
| 28 | Berlinguette et al. (Google-funded; UBC, MIT, LBNL), 2015–2019 | Electrochemical Pd; Pd in D plasma; metal–H systems at elevated temperature (v) | Various | D₂/H₂ | Up to ~10³ °C (v) | Various | Pd loading hard to push above ~0.9 | 420+ experiments (v) | **No anomalous heat** found | N (for its configurations) | [R42] |
| 29 | Czerski et al. (Szczecin, CleanHME), 2024–2025 | d+d on ZrD₂ and Zr targets; UHV accelerator (eLBRUS, 10⁻¹⁰ mbar, up to 1 mA) | Beam on target | — | — | Beam | — | — | Screening energy **340 eV** (vs ~100 eV defect-free); signatures of a new near-threshold d+d channel (PRX 2025) | Mainstream-published, single group (beam-driven; mechanism only) | [R43][R44] |

### 2.1 Context notes not in the table

- **Japanese industry–academia status, 2020–2026.**
  - Tohoku's Condensed Matter Nuclear Reaction Division (CLEAR, founded 2015 with Clean Planet funding) and Clean Planet are the centre of activity. Clean Planet has announced a boiler partnership (Miura Co., 2019) and investors including Mitsubishi Estate (v). **No independent test data from these partners are public.**
  - Toyota Central R&D contributed the 2013 transmutation replication and NEDO-era materials work. Technova (K. Takahashi) and Kobe (Kitamura) continue the powder work. Nissan has published DSC studies.
  - ICCF-26 (Morioka, 2025; Clean Planet was a sponsor) reported "confirmed reproducibility" of Cu–Ni/ZrO₂ heat. Every cited confirmation comes from inside the same network.
  - Iwate (Narita) has run smaller Pd/Ni desorption and permeation heat studies [R51].
- **Flux hypothesis lineage.**
  - Preparata (coherence-QED) argued that driving D along Pd with a longitudinal current (the Coehn effect) raises local loading and flux. ENEA and Celani used 50 µm Pd wires with longitudinal current. They reported high loading and small heat claims (C, v).
  - NASA Glenn (Fralick et al.) noted an unexpected temperature rise during D₂ outgassing of Pd in 1989. They later reported surface "transmutations" in pressure-cycled Pd–Ag (NASA TM 2015; IJHE 2020, v). That work is single-group and EDS-based: D.
  - Biberian and Armanet (ICCF-8, 2000) reported heat during D diffusion through Pd tubes (C/D, v).

---

## 3. Quantitative relationships (with units)

### 3.1 Energy bookkeeping (what a "credible nuclear" claim must add up to)
- 1 W = 6.24×10¹⁸ eV s⁻¹.
- D+D→⁴He (Q = 23.8 MeV) gives **2.6×10¹¹ ⁴He s⁻¹ W⁻¹**, or about 2.3×10¹⁶ He per W·day (≈0.85 µL STP).
- The conventional branches would give **1.5×10¹² T s⁻¹ W⁻¹** (D+D→T+p, 4.03 MeV) or **1.9×10¹² n s⁻¹ W⁻¹** (D+D→³He+n, 3.27 MeV). Watt-level heat with no neutrons or tritium, as Tohoku reports, rules out conventional d+d branches by a factor of >10⁸–10¹⁰ (v, depends on detector sensitivity).
- Chemical ceiling per H: H₂ absorption into hydrides gives ≤0.2–0.9 eV/H (Pd ≈ 0.2 eV/H; ZrH₂ ≈ 0.9 eV/H). H + ½O → water gives ~1.2–2.5 eV/H. So **claims above ~10 eV/H exceed chemistry**, provided the H count and the time integral are right.
- Claimed ratios: 2.4 ± 0.2 eV/D (Kitamura 2009, first phase, chemical range); 15 eV to 2.1 keV/H (PNZ/CNZ); >10 keV/H (JJAP 2024); ≥410 ± 108 keV/H (Kasagi, 460 kJ). Implied number of H: 460 kJ / 410 keV ≈ **7×10¹⁸ H atoms (≈12 µmol)**. Because that denominator is so small, the per-H figure is only as good as the absorption measurement.

### 3.2 Pd–H(D) thermodynamics (bulk vs. nano)
- Van 't Hoff: ln p_pl[bar] ≈ −ΔH/(RT) + ΔS/R, with ΔH ≈ −39 kJ/mol H₂ and ΔS ≈ −93 J mol⁻¹ K⁻¹ (bulk PdH, absorption, v). This gives p_pl ≈ 0.01–0.02 bar at 25 °C and ≈20 bar at 300 °C.
- Bulk critical point: T_c ≈ 566 K (293 °C), p_c ≈ 20 bar, x_c ≈ 0.25–0.27. **Kitamura-type reactors at 200–300 °C run Pd near or above T_c, where there is no distinct β phase.**
- Pd–D plateau pressure is about 3–5× higher than Pd–H at RT (inverse isotope effect, v). β-phase x ≈ 0.6–0.7 at 1 bar RT. x → 0.95+ needs ~10⁴ bar equivalent fugacity.
- **Size effect [A]:** T_c and the gap width fall with size. For d ≲ 3–5 nm there is no two-phase region at RT, plateaus slope, and capacity drops (≈0.3–0.5 H/Pd at 1 bar for ~2.6 nm vs ≈0.6–0.7 for bulk; Yamauchi 2008 (v)). Surface H adds roughly the surface fraction f_s (~0.6 at 2 nm, ~0.25 at 5 nm). Single-particle studies (Baldi 2014, Syrenova 2015, Griessen 2016) show the absorption-plateau hysteresis comes from coherency-strain nucleation barriers. The desorption plateau is nearly size-independent.
- **Implication:** an honest D/Pd > 1 at ≤1 MPa in 5–10 nm Pd is not expected. Bookkeeping artifacts that can produce it include Zr/ZrO₂ side reactions. With Zr/Pd = 1.86 in Zr₆₅Pd₃₅-derived material, about 30% residual metallic Zr forming ZrD₂ would add ~1.1 D/Pd. PdO reduction consumes D₂ into D₂O. Dead-volume and temperature errors in Sieverts measurements also contribute.

### 3.3 Ni, Cu, constantan
- Ni hydride needs ≳0.6 GPa at RT. Solution is endothermic (ΔH_s ≈ +0.15–0.17 eV/H, v). Sieverts' law: x = K_s√p, with x ~10⁻⁵–10⁻³ H/Ni at 1 bar between RT and 800 °C (v). At 10²–10³ Pa, x is 3–10× smaller.
- Cu has lower solubility still (ΔH_s ≈ +0.4–0.5 eV/H, v).
- **Consequence:** in Ni–Cu films and Ni powders, bulk H is tiny. Trapping at interfaces, vacancies, dislocations and oxide boundaries likely dominates. This keeps the chemical background small, which helps credibility, but the H count is error-prone.
- Diffusion: D_H(Ni) ≈ 6.9×10⁻³ exp(−0.41 eV/kT) cm² s⁻¹ (v). That is ~1×10⁻⁶ cm² s⁻¹ at 270 °C and ~8×10⁻⁵ cm² s⁻¹ at 800 °C. √(Dt) after 16 h at 270 °C is ~2.5 mm, far more than the 0.1 mm substrate. **The Ni substrate saturates during loading and acts as an H reservoir that drains outward through the nanolayers on heating.** Reservoir estimate: 0.56 g Ni per plate (5.7×10²¹ Ni atoms) at x ~10⁻⁴–10⁻³ holds 6×10¹⁷–6×10¹⁸ H per plate, consistent with the ~7×10¹⁸ H implied in §3.1.
- Estimated out-flux on heating: ~10¹⁸ H leaving 12.5 cm² over 10²–10⁴ s gives **~10¹³–10¹⁵ H cm⁻² s⁻¹** (derived estimate).

### 3.4 Permeation flux (Pd membranes)
- J = (Φ/L)(√p₁ − √p₂). Pd permeability Φ at ~350 K is of order 10⁻⁹–10⁻⁸ mol H₂ m⁻¹ s⁻¹ Pa⁻⁰·⁵ (v). For L = 0.1 mm and p₁ = 1 atm this gives J ~10¹⁷–10¹⁸ H cm⁻² s⁻¹.
- Iwamura's 2–3 sccm equals 0.9–1.3×10¹⁸ D₂ s⁻¹ (1 sccm = 4.48×10¹⁷ molecules s⁻¹), or 1.8–2.7×10¹⁸ D s⁻¹. Over an effective area of about 4–6 cm² (v) that is **J_D ≈ 3×10¹⁷–7×10¹⁷ D cm⁻² s⁻¹**. Over one week (6×10⁵ s), the total fluence is ≈2–4×10²³ D cm⁻².
- **Transmutation yield per permeated D:** MHI ~10¹⁴ Pr cm⁻² / ~3×10²³ D cm⁻² ≈ **3×10⁻¹⁰** (v). Toyota 1.4×10¹² / (similar fluence, v) ≈ **5×10⁻¹²**. Toyota's yield was about 10²× lower.
- Diffusion time-lag through 0.1 mm Pd at 70 °C: L²/(6D) ≈ 10–20 s (D_H(Pd) ≈ 1×10⁻⁶ cm² s⁻¹ at 343 K). **Permeation is limited by surface kinetics.** Iwamura documented "self-poisoning" that lowered permeation rate over time [R19].
- For comparison, electrochemical cathodes: 100 mA cm⁻² corresponds to 6.2×10¹⁷ D cm⁻² s⁻¹ arriving at the surface, but most recombines to gas. Net permeation fluxes are typically 10¹⁴–10¹⁶ D cm⁻² s⁻¹ (v).
- McKubre's electrochemical excess-power correlation includes an explicit flux term: P_xs = M(x − 0.875)²(i − i₀)|i_D| (v). This is the only quantitative "flux matters" law in the field, and it comes from electrochemistry (see R1).

### 3.5 Calorimetric sensitivities for gas-phase/vacuum set-ups
- Radiative loss: P_rad = εσA(T⁴ − T₀⁴). For a 25×25 mm plate (two faces, A = 12.5 cm²) at 800 °C, the blackbody limit is 93 W, or 18.7 W at ε = 0.2.
  - **dP/dε ≈ 0.93 W per 0.01 change in ε.**
  - dP/dT ≈ 0.07 W/K at ε = 0.2, so **1 W ≈ 14 K** of sample temperature.
  - In a temperature-based calorimeter, a 5% relative emissivity drop, for example when H₂ reduces a surface oxide (oxide ε ~0.5–0.9 vs clean metal ε ~0.1–0.2), mimics ~1 W of "excess".
- Spectral coverage: a 0.3–5.5 µm photon calorimeter captures only F(0→λT) ≈ 72% of blackbody power at 1073 K (λT = 5900 µm·K) and ≈55% at 800 K. The remaining 28–45% has to be modelled, which assumes grey-body behaviour.
- Gas thermal conductivity at ~300 K (W m⁻¹ K⁻¹): H₂ 0.187, He 0.157, D₂ ~0.14, N₂ 0.026, Ar 0.018. **Switching H₂ ↔ D₂ ↔ vacuum changes conductive loss and thermocouple readings.** Evacuation, the step that precedes Iwamura's heat phase, cuts conduction and raises sample temperature at constant heater power. Calibration must therefore be done in the same gas state.
- Mass-spectrometric separation: ⁴He (4.0026 u) vs D₂ (4.0282 u) needs m/Δm ≳ 160. ³He (3.0160 u) vs HD (3.0219 u) or H₃⁺ (3.0235 u) needs **m/Δm ≳ 510**. Most RGAs cannot do either without special high-resolution operation.

---

## 4. Geometry-relevant findings

| Parameter | Association claimed/observed | Direction | Confidence | Basis |
|---|---|---|---|---|
| **Metal feature size** (particle diameter, layer thickness) | Signals are reported only for nm-scale structures: 2–20 nm particles in ZrO₂/SiO₂/Al₂O₃; 2 nm/14 nm Cu/Ni layers; 2 nm CaO/18 nm Pd layers | Smaller and more interfaces appear better; no dose–response curve published | Low–moderate (C) | #6, #7, #12, #15 |
| **Cu:Ni layer ratio** in multilayers | Cu 2 nm / Ni 14 nm (1:7) is the standard; NiCu beat Ni-only, but Ni-only still gave "comparable but slightly smaller" heat (2025) | Cu helps modestly; not essential | Low (C) | #15, #16 |
| **Number of layer periods** | 6 periods (~96 nm) standard; 5 periods in Pd/CaO | No systematic scan published | Unknown | #12, #15 |
| **Oxide spacer layers** (CaO, Y₂O₃, ZrO₂) | Pd/CaO: no Pr without CaO (v). Ni systems: ZrO₂ matrix prevents sintering at 200–400 °C | Needed for transmutation (claim); needed for nanostructure stability (well founded) | Moderate for stability (A-level materials science); C for nuclear role | #6, #12 |
| **Reservoir + active skin** (thick substrate plus nano-skin) | 0.1 mm Ni or Pd substrate feeding flux through a ~100 nm active stack | Asymmetric geometry is common to all "flux" successes | Moderate (C/B pattern) | #12, #15 |
| **Flux direction/phase** | Heat during evacuation plus heating (out-flux); Pr at the D₂-entry surface under in-flux | Flux phase matters more than static loading | Moderate (C → B) | #8, #10, #12, #15, #27 |
| **Temperature** | Ni systems need 200–900 °C: powders 200–300 °C, films 500–900 °C, constantan >400 °C, Nissan DSC 300–500 °C. Pd/CaO at ~70 °C | Higher T helps Ni systems (Arrhenius-like diffusion/desorption) | Moderate as an empirical window | #6, #8, #15, #21 |
| **Powder mass / active area** | 100–200 g powders give 5–25 W (~1.3 W/g-Ni); 12.5 cm² films give 1–6 W (~0.1–0.5 W/cm²) | Power roughly scales with active mass/area; saturation unknown | Low | #6, #15, #16 |
| **Pd distribution on Ni** (Mizuno) | Smaller, more uniform rubbed-on Pd particles give more heat (claim) | Finer is better | Very low (D) | #19 |
| **Wire diameter, knots, surface texture** (Celani) | 200 µm, ~1–1.6 m wires; sub-µm oxidised texture; knots as local current/thermal-gradient concentrators | Rough and knotted better (claim) | Very low (D) | #21 |
| **Micron-scale hot spots** | IR imaging shows µm-scale local heating above Ni's melting point (claim) | Suggests localised active sites (NAE-type), not a uniform bulk effect | Low (C/D); needs independent imaging | #17 |
| **Sealed shell with enclosed nanopowder** (Arata DS) | Permeation into enclosed Pd-black builds internal pressure | Unclear | Very low (C/D) | #1 |
| **Isotope (H vs D)** | Ni systems: H ≈ D. Pd/CaO: D only | Isotope-agnostic in Ni argues against d+d; transmutation claim is D-specific | Moderate for Ni observation | #6, #12, #13 |
| **Particle size < ~5 nm in Pd** | Capacity falls, gap vanishes | Does *not* raise loading | High (A) | [R45–R48] |

---

## 5. Known artifacts and failure modes (ranked by relevance to gas-phase/nano designs)

1. **Emissivity drift in vacuum/radiative set-ups.** H₂ reduces air-formed NiO/CuO, Cu inter-diffuses into Ni (seen above ~400–600 °C in synchrotron XRD [R24]), and surfaces roughen. Each changes ε. At 800 °C an emissivity change of 0.01 is worth ~0.9 W, the same size as the claimed signals (§3.5). Photon calorimetry is more robust but has to extrapolate 28–45% of the spectrum beyond 5.5 µm.
2. **Gas-state-dependent heat transfer.** Conduction changes when H₂ is swapped for D₂ (0.19 vs 0.14 W m⁻¹ K⁻¹) and when gas is evacuated. Excess inferred from temperature is only valid against a calibration in the identical gas and pressure state. "Heat on evacuation" is precisely the configuration most exposed to this artifact.
3. **Thermocouple placement and drift.** Type-K sensors drift several K per 10³ h in reducing atmospheres at 800 °C (v). At 0.07 W/K that is 0.1–0.5 W. Sensors that contact the heater rather than the sample report heater temperature. A few inner and outer TCs cannot map gradients in 100 g powder beds.
4. **Heater calibration drift and input-power metrology.** Heater resistance ages, and V·I sampled slowly misreads pulsed or AC drive. Brillouin's Q-pulses and Mizuno's glow discharge need ≥100 MHz-bandwidth, true-RMS power measurement.
5. **Chemical heat mis-assigned in nanopowders.** PdO + H₂ → Pd + H₂O releases ~1.6 eV per molecule. CuO + H₂ releases ~0.9 eV per molecule. Metallic Zr hydriding releases ~0.9 eV/H. H/D exchange with ZrO₂ hydroxyls creates spurious D₂ > H₂ differences via zero-point-energy effects. Water release and re-adsorption in the support add more. Low-end claims (15–100 eV/H) sit within 1–2 orders of magnitude of these and are vulnerable to integration over weeks of baseline.
6. **Loading-ratio bookkeeping.** Normalising uptake "per Pd" when the matrix holds 1.9 Zr per Pd; Sieverts dead-volume and temperature errors; spill-over to the support. **"D/Pd > 1" in nano-Pd/ZrO₂ should be treated as an artifact until phase-resolved evidence exists (in-situ XRD/neutron).**
7. **Transmutation analytics.**
   - XPS: Pr 3d₅/₂ (~933 eV) overlaps Cu 2p₃/₂ (~932.7 eV).
   - ICP-MS in a Pd matrix: ¹⁰⁵Pd³⁶Ar⁺ sits at m/z 141, where ¹⁴¹Pr is monoisotopic. ¹⁰⁴Pd³⁶Ar⁺ = 140 (Ce). ⁵⁸Ni⁴⁰Ar⁺ = 98 and ⁶⁰Ni⁴⁰Ar⁺ = 100 (Mo, Ru). ⁵⁶Fe⁴⁰Ar⁺ = 96 (Mo).
   - SIMS in CaO-containing stacks: **(⁴⁰Ca)₂¹⁶O⁺ = 96, which would mimic "⁹⁶Mo enrichment."** We could not verify whether the original Sr→Mo analysis corrected for this.
   - ⁹⁶Zr and ⁹⁶Ru isobars in Zr-containing systems.
   - Environmental Pr/Mo contamination (NRL's allegation).
   - Cs→Pr involves two monoisotopic species, so isotope ratios cannot separate contamination from transmutation.
8. **Helium analysis.** Air contains 5.24 ppm He. He permeates glass and elastomer seals. Standard QMS cannot separate ⁴He from D₂ or ³He from HD/H₃⁺ (m/Δm ≈ 160 and 510 needed). Natural D/H (1.56×10⁻⁴) produces HD at m/z 3 even in "pure H₂" runs.
9. **IR thermography misuse.** Rossi's Lugano test used a broadband total emissivity instead of the 7.5–13 µm band emissivity (alumina ε_band ≈ 0.9+). That overestimated temperature and output (v). Any camera-based power claim needs band-specific ε and a dummy-reactor calibration.
10. **Air-flow calorimetry.** Anemometer placement, non-uniform velocity profiles, and calibration at low power with a different temperature distribution. This was Mizuno's method [R29][R50].
11. **Nanostructure evolution.** Sintering, Cu–Ni inter-diffusion, and oxide-layer breakup at 300–900 °C. These make runs history-dependent and can produce apparent "decay" or "activation" of heat.
12. **Radiation detector artifacts.** EMI from discharges and pulsed heaters, radon variations, cosmic-ray bursts, and CR-39 handling. Relevant to CleanHME's "particle emission" claims.
13. **Non-independence.** Shared material sources and shared analysis teams within one network (NEDO/Tohoku/Kobe/Clean Planet) inflate apparent replication. SRI's Brillouin work was a paid contract.

---

## 6. Design implications, ranked by expected value

EV is judged as (probability of a signal) × (credibility if seen) ÷ cost/complexity.

**EV-1. Replicate the Tohoku/Clean Planet Ni–Cu multilayer protocol with emissivity-proof calorimetry. (Highest EV.)**
- *Geometry:* bulk Ni 25×25×0.1 mm substrates. Sputter 6×[Cu 2.0 nm / Ni 14.0 nm] on both faces with in-situ thickness control of ±0.2 nm. Build twin active and dummy samples from the same sputter batch. The dummy is either never H-exposed or loaded with Ar.
- *Protocol:* bake out, then H₂ loading at 250–270 °C for ~16 h, then evacuate and ramp to 600–900 °C. Repeat load/desorb cycles.
- *Calorimetry:* surround the heater-sample assembly with a water-cooled shroud (flow calorimetry, P = ṁc_pΔT), so the measured output does not depend on ε. Add in-situ two-colour pyrometry to track ε. Calibrate in vacuum after identical thermal and H history. Target ≤0.1 W uncertainty at ~20 W input. The claimed 1–6 W would then be a 10–60σ signal.
- *Add-ons:* calibrated RGA/QMS measurement of H out-flux (to test heat vs flux); µm-resolution IR microscopy (to test the hot-spot claim); a Ni-only control; neutron/γ monitors.
- *Why first:* the signal is claimed in a peer-reviewed mainstream journal, the dimensions are well specified, the set-up is small and cheap, and the chemical background is intrinsically tiny because so little H is present. A clean null is also valuable.

**EV-2. Pd/CaO permeation transmutation with isotopically tagged targets and interference-aware analysis.**
- *Geometry:* Pd 25×25×0.1 mm with Pd 40 nm / [CaO 2 nm / Pd 18 nm]×5. Deposit the target on the D₂-facing surface. D₂ at 1 atm upstream, vacuum downstream, ~70 °C. Hold J_D ≈ 3–7×10¹⁷ D cm⁻² s⁻¹ with a mass-flow meter and log it.
- *Controls (same batch):* H₂ permeation, a no-CaO stack, and a no-flux (static D₂) run.
- *Target choice:* pick a target whose claimed +4d/+6d products have low natural abundance. For example, ⁸⁶Sr-enriched Sr predicts ⁹⁴Mo under the claimed ΔA = +8 rule. That defeats the contamination critique; Cs→Pr cannot.
- *Analysis:* blinded HR-ICP-MS (resolving Pd-Ar polyatomics) at two external labs, plus TOF-SIMS with a Ca₂O⁺ check at m/z 96, plus a pre-run survey of lab contamination.
- *Expected signal* if the claim holds: 10¹²–10¹⁴ product atoms cm⁻² against a ~10¹⁰ cm⁻² detection limit. This has the highest credibility payoff of any gas-phase option because it gives an isotopic signature. The risk is that the ~10²× yield spread between MHI and Toyota suggests poorly controlled hidden variables such as surface poisoning and CaO stoichiometry.

**EV-3. Ni-alloy nanocomposite powders (CNZ/PNZ) at 100 g scale with helium and water accounting.**
- *Set-up:* RC ~500 cm³, 100–200 g of Cu₁Ni₇/ZrO₂ or Pd₁Ni₁₀/ZrO₂, RC at 200–300 °C, flow calorimetry at ≤1% of heater power (heater ~80–130 W; claimed excess 5–25 W).
- *Controls:* ZrO₂-only and Ar-fill blanks. Quantify oxygen (water produced) to bound the chemical term.
- *Helium:* for D₂ runs, use a sealed all-metal system with getter and high-resolution QMS. 10 W for one week via ⁴He would give ~1.6×10¹⁸ He (~60 µL STP) in 500 cm³, which is easy to detect. For H₂ runs, test the ³He claim at m/Δm ≥ 510.
- *Why third:* this is the most-reproduced heat claim in-network (grade B), but the chemical confounders are larger and the H/D-agnostic behaviour gives no obvious ash to look for.

**EV-4. Apply flux cycling across all designs.** Make the primary independent variable controlled out-flux or through-flux, not maximum loading. Use repeated load (low T) / desorb (high T) cycles, or permeation with a controlled Δp. Log J(t) and P_xs(t) at the same time. A reproducible correlation between P_xs and J that survives the gas-state calibration controls would be the most informative single result.

**EV-5. Instrumentation and credibility layer (applies to all).**
- Twin active and dummy cells from one batch.
- Blinded sample coding and third-party post-mortem analysis.
- Pre-registered analysis of excess.
- Power metrology at ≥100 MHz for any pulsed drive.
- Continuous neutron and γ spectroscopy, since expected nulls are also informative.
- No elastomer seals in helium experiments.

**Low or negative EV (do not pursue as primary):**
- Mizuno Pd-on-Ni mesh with air-flow calorimetry.
- Celani wires.
- Brillouin-type pulsed tubes (hard input-power metrology).
- Parkhomov/Rossi LiAlH₄–Ni above 1000 °C (thermocouple and heater failures, IR-emissivity pitfalls, nulls).
- Any design whose central claim is "D/Pd > 1 in nano-Pd".
- Arata-style He-4 claims without leak-tight high-resolution mass spectrometry.

---

## 7. Open questions

1. Does excess heat in Ni–Cu films scale with interface area (number of periods, Cu thickness) or with substrate reservoir size (thickness)? No systematic dimensional scan has been published.
2. Is the signal controlled by out-flux J, by temperature, or by their product? The desorption-phase association is consistent across groups but confounded with the evacuation calibration artifact.
3. What is the ash in H–Ni systems? ³He (CNZ, 2024) is the only candidate so far. Is it above HD/H₃⁺ interference at m/Δm ≥ 510?
4. Why do Pd/CaO Pr yields differ by ~10² between MHI and Toyota? What are the surface poisoning, CaO stoichiometry and flux-uniformity variables?
5. Did the Sr→Mo SIMS analysis exclude (⁴⁰Ca)₂¹⁶O⁺ at m/z 96?
6. Will any lab outside the Tohoku–Clean Planet–Kobe network publish a replication of the multilayer heat? This matters more than any further in-network result.
7. Can micron-scale hot spots above Ni's melting point be confirmed post-mortem, for example by melt craters in SEM, and with independent IR calibration?
8. What did HERMES and CleanHME actually measure? Their final deliverables and data sets (Zenodo/CORDIS) need to be read in full. We could not retrieve them.
9. Are SAV-rich structures made by high-pressure anneal or co-deposition more active than as-sputtered films? No direct test exists.
10. Why do H and D behave similarly in Ni powders but not in Pd/CaO transmutation? Are these the same phenomenon at all?

---

## 8. References (URLs)

Primary PDFs were not retrievable in this session (egress policy). URLs are as returned by search.

- [R1] Rothwell, J., "Report on Arata's paper and lecture about his 'Solid Fusion' reactor" (2008). https://www.lenr-canr.org/acrobat/RothwellJreportonar.pdf
- [R2] Chubb, T.A., "In honor of Yoshiaki Arata." https://lenr-canr.org/acrobat/ChubbTAinhonorofy.pdf
- [R3] McKubre, M. et al., SRI analyses of Arata DS-cathodes, ICCF-8/ICCF-10 proceedings (v). Index: https://lenr-canr.org/
- [R4] Arata & Zhang 2008 demonstration coverage. https://www.researchgate.net/publication/252590253_The_Establishment_of_Solid_Nuclear_Fusion_Reactor
- [R5] Kitamura, A. et al., "Anomalous effects in charging of Pd powders with high density hydrogen isotopes," *Phys. Lett. A* 373 (2009) 3109, doi:10.1016/j.physleta.2009.06.061. https://www.osti.gov/etdeweb/biblio/21377082
- [R6] Kidwell, D. et al. (NRL), US Patent 9,182,365, "Excess enthalpy upon pressurization of nanosized metals with deuterium." https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/9182365
- [R7] Kitamura, A. et al., "Excess heat evolution from nanocomposite samples under exposure to hydrogen isotope gases," *Int. J. Hydrogen Energy* 43 (2018) 16187. https://www.sciencedirect.com/science/article/abs/pii/S0360319918320925 ; https://www.researchgate.net/publication/326597398
- [R8] Takahashi, A. et al., JCMNS 33 (ICCF-22) paper; "Enhancement of excess thermal power in interaction of nano-metal and H(D)-gas." https://www.researchgate.net/publication/343893948 ; https://www.researchgate.net/publication/344754960
- [R9] Kitamura/Takahashi et al., "Anomalous heat effects induced by metal nano-composites and hydrogen gas," JCMNS. https://jcmns.org/article/72496-anomalous-heat-effects-induced-by-metal-nano-composites-and-hydrogen-gas.pdf
- [R10] Asia Times, "Cold fusion 2: Japan wins with systematic method" (2019, secondary). https://asiatimes.com/2019/11/cold-fusion-japan-takes-the-lead/
- [R11] "A study of sustained heat generation from Pd-Ni-Zr alloys in hydrogen from the viewpoint of hydrogen absorption," JCMNS 37 (2023) 111–123. https://jcmns.org/article/124623-a-study-of-sustained-heat-generation-from-pd-ni-zr-alloys-in-hydrogen-from-the-viewpoint-of-hydrogen-absorption.pdf
- [R12] "Detections of He-3 in Ni-based binary metal nanocomposites with Cu in zirconia exposed to hydrogen gas at elevated temperatures" (2024). https://arxiv.org/pdf/2409.05382
- [R13] Biberian, J.-P. et al., "Excess heat in nanoparticles of nickel alloys in hydrogen," JCMNS 38 (2024) 186–195. https://jcmns.org/article/124956-excess-heat-in-nanoparticles-of-nickel-alloys-in-hydrogen
- [R14] CleanHME project site, final event (Jan 2025), publications; CORDIS 951974. https://cleanhme.eu/ ; https://www.cleanhme.eu/?p=710 ; https://cleanhme.eu/?page_id=27 ; https://cordis.europa.eu/project/id/951974/results
- [R15] Beiting, E.J., "Investigation of the nickel-hydrogen anomalous heat effect," Aerospace Corp. ATR-2017-01760. https://www.lenr-canr.org/acrobat/BeitingEinvestigat.pdf
- [R16] Iwamura, Y., Sakano, M., Itoh, T., "Elemental analysis of Pd complexes: effects of D₂ gas permeation," *Jpn. J. Appl. Phys.* 41 (2002) 4642. Bibliography: https://newenergytimes.com/v2/reports/Iwamura-LENR-Bibliography.pdf
- [R17] Iwamura et al., "Low energy nuclear transmutation in condensed matter induced by D₂ gas permeation through Pd complexes: correlation between deuterium flux and nuclear products" (ICCF-10). https://lenr-canr.org/acrobat/IwamuraYlowenergyn.pdf ; https://lenr-canr.org/acrobat/IwamuraYobservatiob.pdf
- [R18] Iwamura et al., "Recent advances in deuterium permeation transmutation experiments," JCMNS 10 (2013) 63; "Increase of reaction products…," JCMNS 13 (2014) 242; 2012 ANS paper. https://jcmns.org/article/72216-recent-advances-in-deuterium-permeation-transmutation-experiments/attachment/150206.pdf ; https://jcmns.org/article/72253.pdf ; https://newenergytimes.com/v2/conferences/2012/ANS2012W/2012Iwamura-ANS-LENR-Paper.pdf
- [R19] "Effects of self-poisoning of Pd on the deuterium permeation rate and surface elemental analysis for nuclear transmutation," JCMNS. https://jcmns.org/article/72161-effects-of-self-poisoning-of-pd-on-the-deuterium-permeation-rate-and-surface-elemental-analysis-for-nuclear-transmutation/attachment/150147.pdf
- [R20] Hioki, T. et al., "ICP-MS study on the increase in the amount of Pr atoms for Cs-ion-implanted Pd/CaO multilayer complex with deuterium permeation," *Jpn. J. Appl. Phys.* 52 (2013) 107301. News: https://news.newenergytimes.net/2013/10/22/journal-publishes-toyotas-independent-replication-of-mitsubishi-lenr-transmutation/
- [R21] Grabowski, K. et al. (NRL), "Evaluation of the claim of transmutation of cesium…" https://www.lenr-canr.org/acrobat/GrabowskiKevaluation.pdf
- [R22] New Energy Times, "Latest NRL salvo attacks validity of Mitsubishi LENR research" (2016). https://news.newenergytimes.net/2016/05/12/latest-nrl-salvo-attacks-validity-of-mitsubishi-lenr-research/
- [R23] Iwamura, Y. et al., "Anomalous heat generation that cannot be explained by known chemical reactions produced by nano-structured multilayer metal composites and hydrogen gas," *Jpn. J. Appl. Phys.* 63 (2024) 037001. https://iopscience.iop.org/article/10.35848/1347-4065/ad2622
- [R24] "Changes in the structure and composition of nano-sized Ni-Cu multilayer films with increasing temperature in an atmosphere with and without hydrogen," *Front. Mater.* (2024). https://www.frontiersin.org/journals/materials/articles/10.3389/fmats.2024.1407810/full
- [R25] Kasagi, J. et al., "Photon radiation calorimetry for anomalous heat generation in NiCu multilayer thin film during hydrogen gas desorption," arXiv:2311.18347; JCMNS. https://arxiv.org/abs/2311.18347 ; https://jcmns.org/article/134004-photon-radiation-calorimetry-for-anomalous-heat-generation-in-nicu-multilayer-thin-film-during-hydrogen-gas-desorption ; slides: http://ikkem.com/iccf23/orppt/ICCF23-IA-17%20Kasagi.pdf
- [R26] Kasagi, J. et al., "Measurement of radiant spectrum for excess heat generation in NiCu and Ni thin film during hydrogen gas desorption," arXiv:2509.13847 (2025). https://arxiv.org/abs/2509.13847
- [R27] ICCF-26 (Morioka, 26–30 May 2025) summaries and abstracts. https://www.sciengine.com/JMCC/doi/10.16084/j.issn1001-3555.2025.06.009 ; https://iccf26.org/img_common/abstract/iccf26_abstract.pdf
- [R28] Mizuno, T., "Observation of excess heat by activated metal and deuterium gas," JCMNS 25 (2017). https://www.lenr-canr.org/acrobat/MizunoTpreprintob.pdf
- [R29] Mizuno, T. & Rothwell, J., "Increased excess heat from palladium deposited on nickel" (ICCF-22, 2019); "Excess heat from palladium deposited on nickel," JCMNS. https://www.lenr-canr.org/acrobat/MizunoTincreasede.pdf ; https://jcmns.org/article/72489-excess-heat-from-palladium-deposited-on-nickel.pdf
- [R30] E-Cat World, "Jed Rothwell reports on status of Mizuno replication efforts" (Dec 2019, secondary). https://e-catworld.com/2019/12/04/jed-rothwell-reports-on-status-of-mizuno-replication-efforts/
- [R31] E-Cat World, "Excess heat production reported in new Mizuno replication (Hokkaido University of Science)" (secondary). https://e-catworld.com/2019/12/10/excess-heat-production-reported-in-new-mizuno-replication-hokkaido-university-of-science/
- [R32] Celani, F. et al., JCMNS 27 (2018) 9–21; JCMNS 33 (2020) 46–73. https://jcmns.org/article/72475.pdf ; https://jcmns.org/article/72549-progress-toward-an-understanding-of-lenr-ahe-effects-in-coated-constantan-wires-in-d2-atmosphere-dc-ac-voltage-stimulation.pdf
- [R33] Celani et al., "First evaluation of coated constantan wires incorporating Capuchin knots…" https://www.researchgate.net/publication/338701975
- [R34] MFMP, "Celani's wire excess heat effect replication," JCMNS; "Final report on calorimetry-based excess heat trials using Celani treated NiCuMn wires," JCMNS. https://jcmns.org/article/72341.pdf ; https://jcmns.org/article/72376
- [R35] SRI International (Tanzella), "March 2018 isoperibolic hydrogen hot tube reactor studies." https://brillouinenergy.com/newwebsite/wp-content/uploads/2018/12/SRI_Technical_Report.pdf
- [R36] "Critique of the SRI Brillouin HHT report." https://coldfusioncommunity.net/critique-of-the-sri-brillouin-hht-report/
- [R37] Parkhomov-type replication nulls (Biberian, >20 runs, ±2 W), as reported in search-returned summaries. Primary not retrieved; see JCMNS index https://lenr-canr.org/wordpress/?page_id=1495
- [R38] Levi, G. et al., "Indication of anomalous heat energy production in a reactor device," arXiv:1305.3913. https://arxiv.org/pdf/1305.3913
- [R39] Ericsson, G. & Pomp, S., "Comments on the report 'Indications of anomalous heat energy production…'," arXiv:1306.6364. https://arxiv.org/pdf/1306.6364. Lugano (2014) emissivity critiques and Rossi v. Darden (S.D. Fla., 2016–2017, settled) are cited from memory (v).
- [R40] Li, X.Z. et al., "Correlation between abnormal deuterium flux and heat flow in a D/Pd system," *J. Phys. D: Appl. Phys.* 36 (2003) 3095 (v).
- [R41] Yamaguchi, E. & Nishioka, T., "Cold nuclear fusion induced by controlled out-diffusion of deuterons in palladium," *Jpn. J. Appl. Phys.* 29 (1990) L666 (v).
- [R42] Berlinguette, C.P. et al., "Revisiting the cold case of cold fusion," *Nature* 570 (2019) 45. https://www.nature.com/articles/s41586-019-1256-6
- [R43] Czerski, K. et al., "Observation of thermal deuteron-deuteron fusion in ion tracks," arXiv:2409.02112; Dubey, R., Czerski, K. et al., *Materials* 18 (2025) 1331; *Phys. Rev. X* 15 (2025), "Experimental signatures of a new channel of the deuteron-deuteron reaction at very low energy." https://arxiv.org/abs/2409.02112 ; https://doi.org/10.3390/ma18061331
- [R44] HERMES project (H2020 FET 952184) site, CORDIS, Zenodo; Winzeler, H.B., JCMNS 41 (2026) 153–163. https://hermesproject.eu/ ; https://cordis.europa.eu/project/id/952184/results ; https://zenodo.org/communities/hermes-h2020/about ; https://jcmns.org/article/163383.pdf
- [R45] Pundt, A. & Kirchheim, R., "Hydrogen in metals: microstructural aspects," *Annu. Rev. Mater. Res.* 36 (2006) 555; Sachs, C. et al., *Phys. Rev. B* 64 (2001) 075408 (v).
- [R46] Yamauchi, M. et al., "Nanosize effects on hydrogen storage in palladium," *J. Phys. Chem. C* 112 (2008) 3294 (v).
- [R47] Baldi, A. et al., *Nat. Mater.* 13 (2014) 1143; Syrenova, S. et al., *Nat. Mater.* 14 (2015) 1236; Griessen, R. et al., *Nat. Mater.* 15 (2016) 311 (v).
- [R48] Fukai, Y. & Ōkuma, N., "Formation of superabundant vacancies in Pd hydride under high hydrogen pressures," *Phys. Rev. Lett.* 73 (1994) 1640; Fukai, Y., *J. Alloys Compd.* 356–357 (2003) 263; *The Metal–Hydrogen System* (Springer, 2005) (v).
- [R49] Yamaura, S. et al., "Hydrogen absorption of nanoscale Pd particles embedded in ZrO₂ matrix prepared from Zr–Pd amorphous alloys," *J. Mater. Res.*; Zr–Pd–Ni oxide composites, *Mater. Trans.* 44 (2003) 696. https://www.cambridge.org/core/journals/journal-of-materials-research/article/abs/hydrogen-absorption-of-nanoscale-pd-particles-embedded-in-zro2-matrix-prepared-from-zrpd-amorphous-alloys/083A5AFCB431307DC33C2DEF8B8739D6 ; https://www.jstage.jst.go.jp/article/matertrans/44/4/44_4_696/_article
- [R50] "Basics of air-flow calorimetry," JCMNS. https://jcmns.org/article/72559-basics-of-air-flow-calorimetry.pdf
- [R51] Narita, S., "Summary of LENR research in Japan," ARPA-E LENR Workshop (Oct 2021). https://arpa-e.energy.gov/sites/default/files/migrated/2021LENR_workshop_Narita.pdf
- [R52] Japan CF-Research Society abstracts (JCF20/21/23). https://jcfrs.org/JCF21/jcf21-abstracts.pdf ; https://jcfrs.org/wordpress/wp-content/uploads/2023/03/jcf23-abstracts-updated.pdf
- [R53] Clean Planet Inc., technology page (corporate claims). https://www.cleanplanet.co.jp/en/technology/
- [R54] ICCF-25 Book of Abstracts (Szczecin, 2023). https://newenergytimes.com/v2/conferences/2023/ICCF25/ICCF-25-Book-of-Abstracts-2023.07.04.pdf
- [R55] Review: "A review of experiments reporting non-conventional phenomena in nuclear matter…," arXiv:2308.13533. https://arxiv.org/pdf/2308.13533
