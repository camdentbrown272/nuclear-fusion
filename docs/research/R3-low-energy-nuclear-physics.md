# R3 — Low-energy D–D (and p–D, p–Li) nuclear physics in condensed matter

*Scope: the established and frontier nuclear physics of very-low-energy fusion reactions in metals: accelerator screening data, the proposed ⁴He threshold resonance and new e⁺e⁻ channel, ARPA-E and NASA programs, rate theory, and structural dependences. The aim is to find what a cold, lattice-based device could credibly measure.*

**How the sources were checked (read first).** Live web access in this session was badly limited. Publisher, arXiv and OSTI pages were blocked by the network proxy, and the shared web-search budget ran out after a handful of queries. Every factual claim therefore carries one of three tags:

- **[v]**: confirmed this session from search-engine abstracts or summaries of the primary source.
- **[m]**: from my memory of the primary literature and **not re-verified** this session. Treat the numbers as ±1 unit in the last digit at best, and check them against the paper before any quantitative design use.
- **[calc]**: computed here. The code is in the Appendix. It is consistent with `sim/m0_rate_budget.py` and `docs/models/M0-rate-budget.md`.

Evidence grades: **A** = independently replicated, mainstream; **B** = replicated by several groups but debated; **C** = single group; **D** = anecdotal or disputed; **N** = credible null.

---

## 1. Executive summary

1. **Screening enhancement at keV energies is real, but its size is not understood. Existence A−, values B.** Four or more independent groups (Bochum, Tohoku, TU Berlin/Szczecin, JINR) see d(d,p)t yields in metals at E_cm ≲ 10 keV above bare-nucleus expectations. The fitted screening energies are U_e ≈ 100–800 eV. Theory gives 20–30 eV (gas, Thomas–Fermi) to about 110 eV (the best dielectric-function model). For the same metal the value differs 2–3× between groups: Pd is 310 eV at Tohoku and 800 eV at Bochum.
2. **The accelerator U_e cannot be a static potential acting on thermal deuterons. Grade A for the logic, N for the nulls. [calc]** If Pd's U_e = 800 eV applied to room-temperature PdD, 1 cm³ would emit about 7×10¹⁷ n/s, roughly 0.8 MW/cm³. With a Yukawa-softened potential (U_eff = 331 eV) it is still about 2×10⁹ n/s/cm³. Bulk null searches (≲10⁻²³–10⁻²⁵ fusions/pair/s) bound U_eff(thermal) ≲ 150–165 eV.
3. **Measurable D–D in metals requires keV relative energies. [calc]** At E_cm = 1 keV the screening enhancement is ×4 for U_e = 100 eV, ×48 for 300 eV and ×3×10³ for 800 eV. A 1 mA D⁺ beam on PdD at 5 keV gives about 11–240 n/s, depending on U_e. At 2 keV it gives 10⁻⁴–0.8 n/s. **The screening lever is only visible at 1–3 keV.**
4. **Czerski's threshold resonance and new e⁺e⁻ channel. Grade C. [v]** The group proposes a 0⁺ single-particle d+d state in ⁴He at E_x ≈ 23.85 MeV, right at the d+d threshold (23.847 MeV), with Γ_p ≈ 40 meV. Dubey et al. (PRX 15, 041004, 2025) report e⁺e⁻ pairs (continuum up to 22.8 MeV), 511 keV annihilation and bremsstrahlung from D+D on a ZrD₂ plate down to E_d = 5 keV. They infer Γ_e⁺e⁻ ≥ 10 Γ_p, and that below about 5 keV D+D proceeds *mainly* through e⁺e⁻ emission. This comes from a single group with no independent replication.
5. **"Thermal fusion in ion tracks". Grade C. [v]** The Szczecin group reports a yield plateau below E_d ≈ 2.5 keV in ZrD₂ (2024), extended to Ti and Pd (arXiv, May 2026). They interpret it as thermal-spike fusion. Beam-purity and background artifacts have not been excluded by an independent group.
6. **Structure matters, but the evidence is weak.** Reported U_e dependences, with grades:
   - Host electronic character (metal vs insulator): **B**.
   - D concentration (low U_e at high loading): **C/B**.
   - Surface contamination and oxides: **B**. Every group treats this as a major systematic.
   - Vacancy-type defects: **C**.
   - Temperature (T^−½): **C, disputed**.
   - There are **no** credible accelerator data on cracks or nanostructure. This is a gap.
7. **NASA "lattice confinement fusion" (2020). Grade C, unreplicated.** Fast neutrons were reported from bremsstrahlung-irradiated ErD₃ and TiD₂. The reacting deuterons are recoils from MeV photoneutrons, making this *hot* fusion inside a cold lattice.
8. **ARPA-E LENR program (2023–): 8 teams, $10 M. Provisionally N.** I know of no peer-reviewed confirmation of anomalous nuclear products by mid-2026. This is *unverified*, because the search budget ran out. Related results:
   - UBC, Nature 2025: electrochemical D loading raised a beam-target D–D rate by about 15%. **C**, and the mechanism is conventional. [m]
   - Google program, Nature 2019: no anomalies. **N**. [m]
9. **Muon-catalysed fusion (A) shows the one lever that works.** Shrinking the D–D separation 207× raises the rate by about 73 orders of magnitude (3×10⁻⁶⁴ → ~10⁹–10¹² s⁻¹). Electronic and lattice effects change U_eff by a few-fold at most. The rate scales as U_eff^(18–85).
10. **Design bottom line: the highest-expected-value configuration** has these parts:
    - A thin, well-characterised deuterided film (Pd, Zr, Ti), with defects and oxide treated as controlled variables.
    - UHV-clean surfaces.
    - A sub-10 keV deuteron probe (glow discharge or ion source).
    - A multi-channel detector set: Si for p/t/³He, a 511 keV coincidence pair plus an MeV-electron telescope, and pulse-shape-discriminating neutron counters.
    - "Beam-off", purely cold runs with the same detectors, plus H-loaded controls.

---

## 2. Screening measurements

### 2.1 Table of reported screening energies

For a thick-target yield the screening energy U_e is defined through the enhancement f(E) = σ_s(E)/σ_b(E) ≈ exp(πη U_e/E) (Assenbaum, Langanke & Rolfs 1987). E is the centre-of-mass energy. E_d is the laboratory deuteron energy, with E_cm = E_d/2.

| Host / target | U_e (eV) | Energy range | Group (lead) | Year | Notes | Tag |
|---|---|---|---|---|---|---|
| D₂ gas | 25 ± 5 | E_cm ≈ 1.6–130 keV | Bochum (Greife et al.) | 1995 | Reference value. Adiabatic limit ≈ 20 eV | [m] |
| Ta | 309 ± 12 | E_d ≈ 5–60 keV | Bochum (Raiola et al., EPJA 13, 377) | 2002 | First large metallic value | [m] |
| **Pd** | **800 ± 90** | E_d ≈ 5–60 keV | Bochum (Raiola et al., EPJA 19, 283) | 2004 | Largest in the 58-sample survey | [m] |
| Pt | 670 ± 50 | same | Bochum | 2004 | Reported T^−½ decrease 20→340 °C (J. Phys. G 31, 1141, 2005) | [m] |
| Co | 640 ± 70 | same | Bochum | 2004 | Also used in the temperature study | [m] |
| Al | 520 ± 50 | same | Bochum | 2004 | Conflicts with TU Berlin (190 eV) | [m] |
| Cu / Fe / Ni | 470±50 / 460±60 / 380±40 | same | Bochum | 2004 | | [m] |
| Ag / Au / W / Ta | 330±40 / 280±50 / 250±30 / 270±30 | same | Bochum | 2004 | Ta revised down from 309 eV | [m] |
| Be | 180 ± 40 | same | Bochum | 2004 | | [m] |
| Group 3–4 metals and lanthanides (Sc, Ti, Y, Zr, Hf, La, …) at 20 °C | ≲ 30 (gas-like) | same | Bochum | 2004 | High D solubility (hydride formation). Became large at ~200 °C once x fell to a few % (EPJA 27 s01, 79, 2006) | [v] qualitative, [m] numbers |
| Insulators and semiconductors (C, Si, Ge, oxides) | ≲ 30–60 | same | Bochum | 2004 | "Small effects in insulators, semiconductors, lanthanides" | [v] |
| Pd | 310 ± 30 | E_d ≈ 2.5–10 keV | Tohoku (Kasagi, Yuki, … with Lipson) | 2002 | JPSJ 71, 2881 | [m] |
| PdO | ≈ 600 | same | Tohoku | 2002 | Oxide surface roughly doubles U_e vs Pd | [m] |
| Fe / Au | ≈ 200 / ≈ 70 | same | Tohoku | 2002 | | [m] |
| Au/Pd/PdO heterostructure | enhancement ~×10–50 at 2.5 keV | E_d ≈ 2.5–10 keV | Yuki, Kasagi, Lipson et al. (JETP Lett. 68, 823) | 1998 | Reported as an enhancement, not a clean U_e | [m] |
| Yb | "large" (value not verified) | keV | Tohoku | ~2004–08 | I could not retrieve the number | — |
| Al | 190 ± 15 | E_d ≈ 5–60 keV, UHV | TU Berlin (Czerski et al., EPL 54, 449) | 2001 | UHV, controlled surfaces | [m] |
| Zr | 297 ± 8 | same | TU Berlin | 2001 | | [m] |
| Ta | 322 ± 15 | same | TU Berlin | 2001 | | [m] |
| Al, Zr, Ta (systematics) | as above | 5–60 keV | Huke, Czerski et al. (PRC 78, 015803) | 2008 | Surface C/O layers and D-profile evolution strongly bias U_e. UHV and in-situ cleaning required | [m] |
| ZrD₂ (UHV) | 105 ± 15 (fit including the threshold resonance) | down to E_d ≈ 5 keV | Szczecin (Czerski et al.) | 2020–22 | Theory 112 eV. The resonance absorbs the excess | [v] value; [m] which paper and fit |
| ZrD₂ | ≈ 340 (thick-target fit) | E_cm ≥ 0.675 keV | Szczecin (arXiv:2409.02112) | 2024 | Yield falls 7 orders of magnitude. Plateau below E_d ≈ 2.5 keV | [v] |
| Zr with O/C contamination | varies; linked to vacancy defects (PAS, XRD) | keV | Szczecin (Kowalska, Targosz-Ślęczka et al., Materials 2023, 2025) | 2023–25 | Defects are claimed to increase U_e | [v] titles and abstracts only |
| Ti, Pd | measured; values not retrieved | down to ~1 keV | Szczecin (arXiv:2605.27438) | 2026 | "Different thermal characteristics and screening energies" | [v] |
| ZrD₂, TiD₂ (pulsed Hall plasma accelerator) | ~100–300, large errors | E_d ≈ 2–7 keV | JINR Dubna / Tomsk (Bystritsky et al.) | 2008–13 | Low confidence in my recall | [m] |
| Li metal; PdLi alloys (⁶,⁷Li(p,α)) | ~1,000–4,000 | E_p ≈ tens of keV | Bochum/Lisbon (Cruz et al., PLB 624, 181; J. Phys. G 2008) | 2005–08 | Largest values reported for any system | [m] |

**Reading the table honestly.**

- The *existence* of enhanced yields at E_cm ≲ 5 keV in metals is replicated, and no group reports gas-like yields for clean Pd or Ta at 2–5 keV.
- The *numbers* are not reproducible between laboratories. Al is 520 vs 190 eV, Pd 800 vs 310, and ZrD₂ ≈ 30 vs 297–340 vs 105 (the last with a resonance in the fit).
- The spread tracks the experimental regime: UHV or not, implanted vs pre-loaded deuterium, lowest energy reached, and whether stopping powers and D profiles were measured.

### 2.2 Theory and the "screening puzzle"

**Expected values.**

- *Adiabatic (atomic/molecular) limit* for d+d: U_e ≈ 20–28 eV. The D₂ gas measurement (25 ± 5 eV) agrees. [m]
- *Thomas–Fermi conduction-electron screening* [calc]: for one conduction electron per Pd atom (n = 6.8×10²² cm⁻³), k_F = 1.26 Å⁻¹ and λ_TF = 0.57 Å, so U_e ≈ e²/λ_TF ≈ **25 eV**. Three electrons per atom gives 30 eV. Metals should add only tens of eV to the atomic value.
- *Dielectric-function models* (Czerski, Huke et al., EPJA 27 s01, 83, 2006) add static and dynamic contributions. They reach about **100–112 eV** for Zr. [v/m]
- *Classical Debye model* (Raiola et al. 2004/2006): U_e = 2.09×10⁻¹¹ (Z₁Z₂)(n_eff ρ_a/T)^½ eV, with ρ_a in m⁻³ and T in K. It reproduces the *magnitude* of the Bochum numbers with n_eff taken from Hall coefficients, and it predicts the T^−½ law. It is physically untenable, for two reasons [calc]:
  - For Pd, reproducing 800 eV needs R_D ≈ 0.018 Å and n_eff ≈ 6.4 electrons per atom.
  - A Debye sphere then contains about 10⁻⁵ electrons. Conduction electrons are degenerate (E_F ≫ kT), so kT is not the right screening energy.
  The model is a mnemonic, not an explanation.

**The puzzle.** Measured U_e in metals is 3–30× above credible theory. Leggett & Baym (PRL 63, 191, 1989) gave a rigorous bound: equilibrium many-body screening in a metal cannot raise D–D tunnelling anywhere near the rates claimed in 1989. [m, standard result] That bound applies to *static, ground-state* configurations, and it is the right frame for §7.3. Candidate explanations, none established:

1. **Low-energy stopping-power errors.** Hydrogen stopping in transition metals below the Bragg peak is known only to roughly 10–20%. [calc] At E_cm = 5 keV, a 20% yield normalisation error mimics ΔU_e ≈ 130 eV. At 2.5 keV the same error mimics only about 46 eV. **Lower-energy data are intrinsically more robust.**
2. **Deuteron depth profiles.** Implanted D diffuses, accumulates at surfaces or oxide interfaces, and is desorbed by the beam. The thick-target analysis assumes a known uniform profile.
3. **Surface contamination.** C and O layers shift the beam energy and trap D. Huke et al. (2008) and the later Szczecin papers show this changes fitted U_e substantially.
4. **Nuclear-physics mimicry.** A near-threshold resonance (Czerski, §3) raises low-energy yields with a steeper energy dependence than screening.
5. **Dynamic or target-specific effects**: thermal spikes (§3.3), channelling, and non-equilibrium electron densities along the ion track.
6. **Bare S-factor normalisation.** Gas-target reactions such as ³He(d,p) and ⁶,⁷Li(p,α) also show U_e above the adiabatic limit (the "lithium puzzle"). Trojan-Horse-Method bare S-factors partly reconcile these. [m]

The screening enhancement is an established *phenomenon* with a *disputed magnitude and mechanism*. That combination is what makes it both interesting and treacherous for a device.

### 2.3 Reported dependences (details in §8)

- Metal vs insulator: large vs small. **B.**
- Temperature: Bochum reports U_e(Pt, Co) ∝ T^−½. Others attribute this to temperature-dependent D concentration and profiles. **C, disputed.**
- D concentration: at room temperature, hydride-forming metals (Ti, Zr, …) show small U_e at high x (Bochum) but about 300 eV in Zr under UHV at TU Berlin. **Contradictory.**
- Oxide surface: PdO ≈ 2× Pd (Tohoku). **C.**
- Lattice defects: U_e rises with vacancy-type defect density measured by positron annihilation spectroscopy (PAS) (Szczecin). **C.**

---

## 3. Threshold resonance and new-channel results

### 3.1 The proposed 0⁺ threshold resonance in ⁴He (Czerski et al., 2016–2024)

- **Origin.** EPL 113, 22001 (2016) [m]: thick-target d(d,p)t yields on Zr under UHV rose faster at the lowest energies than any constant U_e could explain. The authors proposed a **0⁺ single-particle (d+d) resonance at the d+d threshold**. Later cross-section analyses kept the proposal: Acta Phys. Pol. B 51 (2020); PRC 106, L011601 (2022), "Deuteron–deuteron nuclear reactions at extremely low energies". [v]
- **Parameters** [v]:
  - Excitation energy E_x ≈ **23.85 MeV** in ⁴He. The d+d threshold is at 23.847 MeV, so the state sits within keV of it (the exact E_R above threshold depends on the fit).
  - Proton partial width Γ_p ≈ **40 meV**.
  - Destructive interference with the known broad ⁴He resonance(s) is invoked to explain the energy dependence.
  - With the resonance in the fit, U_e(Zr) comes out at 105 ± 15 eV, close to the 112 eV theory. The resonance removes much of the "excess screening".
- **Why an e⁺e⁻ channel?** [v] A 0⁺ → 0⁺ (ground-state) γ transition is strictly forbidden. If E_R is only eV–keV above threshold, the deuteron width is tiny (strong penetrability suppression) and the nucleon widths are small. E0 internal pair creation (e⁺e⁻) and internal conversion can then dominate the total width.
- **First electron data.** Czerski et al., arXiv:2305.17101 and PRC 109, L021601 (2024), "Indications of electron emission from the deuteron–deuteron threshold resonance". [v] The electron energy spectrum and the electron/proton branching agreed with pair decay of the 0⁺ state. The detector technique (MeV electrons in thin Si) is documented in Measurement 228 (2024) 114392 (arXiv:2312.12446). [v]

### 3.2 Dubey et al., Phys. Rev. X 15, 041004 (published 7 Oct 2025) [v]

Note: this is **University of Szczecin** work, not LBNL or ARPA-E. It was funded within the EU CleanHME programme context. [v/m]

- **Title:** "Experimental Signatures of a New Channel of the Deuteron–Deuteron Reaction at Very Low Energy". Authors: R. Dubey, K. Czerski, Gokul Das H., A. Kowalska, N. Targosz-Ślęczka, M. Kaczmarski, M. Valat. arXiv:2408.07567.
- **Target:** a 0.5 mm thick **ZrD₂ plate** on an ultra-high-vacuum accelerator.
- **Beam:** deuterons down to **E_d = 5 keV** (E_cm = 2.5 keV).
- **Detectors:**
  - Si detectors 0.3–3 mm thick for charged particles (3.02 MeV p, 1.01 MeV t) and for electrons and positrons from internal pair creation (continuum up to **22.8 MeV**, i.e. 23.85 − 1.022 MeV).
  - Large-volume **NaI(Tl) and HPGe** detectors, which recorded **511 keV annihilation radiation and a bremsstrahlung continuum** correlated with the beam.
- **Claims:**
  1. The e⁺e⁻ / proton branching ratio *rises* as beam energy falls.
  2. Interference between the threshold resonance and a broad ⁴He resonance reproduces this trend.
  3. Γ_e⁺e⁻ must be **at least 10× Γ_p**.
  4. **Below about 5 keV, D+D proceeds mainly through the e⁺e⁻ channel.**
- **Rates:** I could not retrieve the absolute count rates or branching values from the paper.

**Assessment (grade C).**

*For the claim:*
- Three independent signatures: the electron continuum, 511 keV lines and bremsstrahlung.
- Consistent across 2022–2025 papers.
- Published in PRX.

*Against, or not yet addressed:*
1. One group proposed the resonance *and* measured the evidence. I know of no independent replication.
2. In 0.3–3 mm Si, multi-MeV electrons deposit only minimum-ionising energy (about 0.39 MeV/mm). The "spectrum" is therefore a weak discriminator against cosmic muons and beam-induced Compton electrons.
3. 511 keV photons have mundane sources (pair production in shielding, β⁺ activity).
4. The A = 4 level scheme (Tilley, Weller & Hale 1992) lists broad states near 20–26 MeV. I know of no ab-initio four-nucleon calculation predicting a narrow 0⁺ at the d+d threshold.
5. If the e⁺e⁻ channel dominated at thermal energies in metals, every "excess heat" claim would come with an intense 511 keV and bremsstrahlung field: about 5×10¹¹ annihilation photons/s per watt. **This mechanism does not provide "heat without radiation".**

**Decisive independent test.** Measure e⁺e⁻ (511–511 keV coincidence plus a high-energy electron telescope) from D+D at E_d = 5–20 keV on a **windowless gas target** with no lattice. The resonance is nuclear, and Czerski claims it appears in gaseous as well as metallic targets.

**Quantitative consequence for design [calc].** A resonance acts as a multiplier K on the S-factor, so it moves the rate linearly, whereas U_e moves it exponentially. Table D (§7.4) shows that K = 10 is worth only +5% in U_eff at 100 eV, and even K = 10⁶ is worth +35%. The resonance changes **which radiation to detect** far more than whether cold rates become measurable.

### 3.3 Thermal D–D fusion in ion tracks (Szczecin) [v]

- **arXiv:2409.02112 (2024).** d(d,p)t on ZrD₂ down to E_cm = 675 eV (E_d = 1.35 keV) on a high-current UHV accelerator. The thick-target yield falls by 7 orders of magnitude and is fitted with U_e = 340 eV. Below E_d ≈ 2.5 keV a **constant yield plateau** appears. The proton energies are higher than beam kinematics predicts, consistent with emission from a near-resting centre of mass. The authors model this as thermal fusion in ion tracks: phonon density locally exceeds the melting temperature.
- **arXiv:2605.27438 (May 2026).** The plateau appears at beam energies as low as 1 keV. Ti and Pd targets, which have different thermal properties and screening energies, "confirmed the thermal spike model" and pointed to the role of enhanced D diffusion in tracks.
- **Assessment: C.** An energy-independent plateau is exactly what beam-contamination artifacts produce: a D₂⁺/D₃⁺ or high-energy tail, scattered or neutral beam, or a small leakage current at the full terminal voltage. Room-background protons can do the same. Ruling these out needs independent beam-energy analysis. If real, it is still *micro-scale hot fusion* (keV deposited in nm³ volumes), not equilibrium cold fusion. It does suggest that **host thermal conductivity and D diffusivity** could be design variables.

---

## 4. ARPA-E LENR Exploratory Topic (2023–2025)

- **Programme** [v]: ARPA-E announced **$10 M for 8 projects** in February 2023. The aim was well-diagnosed experiments to establish or refute LENR nuclear signatures. It included "capability teams" for nuclear radiation diagnostics (University of Michigan) and materials analysis (Robert Duncan's lab, Texas Tech).

| Team | Topic | Tag |
|---|---|---|
| Amphionic LLC | "Nanostructured Pd–ANF (aramid-nanofibre) composites for controlled LENR" ($295,924) | [v] |
| Lawrence Berkeley National Lab (T. Schenkel) | "Quantifying Nuclear Reactions in Metal Hydrides at Low Energies": keV deuterium ions and plasmas on metal hydrides with calibrated neutron, charged-particle and X-ray diagnostics | [m] |
| Texas Tech University | Advanced materials fabrication, characterisation and analysis; detection of nuclear products (capability team, Duncan) | [v] |
| University of Michigan (2 awards) | Deuterium gas cycling experiments with neutron, γ and ion emission measurements; nuclear radiation diagnostics capability team | [v] |
| Stanford University | Project topic not verified | [v] existence |
| MIT | Project topic not verified. MIT researchers (Metzler, Hagelstein et al.) published a review of "known mechanisms that increase nuclear fusion rates in the solid state" (New J. Phys. 2024) | [m] |
| ENG8 | Company-led heat and nuclear test; topic not verified | [m] |

**Results by 2026 (honest status).** I could not verify any peer-reviewed publication from these teams reporting anomalous nuclear products, and none are widely cited as positive. The ARPA-E FY2023 annual report (released August 2025) came up in programme searches, but I could not read it. Pending a direct check of ARPA-E's final programme summary, **treat the programme as a provisional null (N)**.

**The closest peer-reviewed results in this space [m]:**

- **UBC Berlinguette group, Nature 2025**, "Electrochemical loading enhances deuterium fusion rates in a metal target". A benchtop "Thunderbird" reactor fires keV-scale D⁺ plasma ions into a Pd target while a D₂O electrochemical cell loads the back side. Loading raised the conventional D–D neutron rate by about **15%**. The authors attribute this to higher near-surface D density, not to new physics. **C.** It is the first peer-reviewed demonstration that electrochemistry (a "cold" handle) modulates a real nuclear rate.
- **Google programme, Nature 570, 45 (2019)**, "Revisiting the cold case of cold fusion" (Berlinguette, Chiang, Munday, Schenkel, Fork, Koningstein, Trevithick). Electrochemistry, pulsed deuterium plasmas and metal powders were studied. **No anomalous heat or nuclear products. N.**

---

## 5. NASA Glenn "lattice confinement fusion" (LCF)

- **Papers:**
  - Pines et al., "Nuclear fusion reactions in deuterated metals", PRC 101, 044609 (2020): theory and modelling.
  - Steinetz et al., "Novel nuclear reactions observed in bremsstrahlung-irradiated deuterated metals", PRC 101, 044610 (2020): experiment. [m]
- **Mechanism chain:**
  - A 2.9 MeV electron beam produces bremsstrahlung, which photodisintegrates deuterons in ErD₃ or TiD₂ (threshold 2.224 MeV).
  - [calc] The photoneutrons carry at most (2.9 − 2.224)/2 ≈ 0.34 MeV.
  - Elastic n–d scattering transfers up to 8/9 of that, giving deuteron recoils of ≲ 0.3 MeV, mostly tens of keV.
  - These recoils fuse with lattice deuterons, screened by the host electrons. The authors also invoke Oppenheimer–Phillips (d,p) stripping on Er and Ti.
- **Reported observations** [m]:
  - Neutron spectra from a liquid scintillator (EJ-309, with pulse-shape discrimination) and HPGe signatures consistent with 2.45 MeV D–D neutrons plus a higher-energy component.
  - γ activation interpreted as stripping products.
  - Hydrogenated and unirradiated controls showed no such signals.
- **What is "cold" and what is not.** Only the bulk lattice temperature and the fuel density are "cold". The fusing deuterons are **keV–MeV projectiles** produced by MeV photons, so this is beam-target fusion with an internal neutron-driven beam.
- **What screening can do here [calc].** For recoils of 10–300 keV, screening of even 1 keV multiplies σ by only 1.0–2.4. The yield should be close to conventional beam-target estimates; screening matters only in the low-energy tail.
- **Detection difficulty [calc, order of magnitude].** Photoneutrons are ≤0.34 MeV and each makes a few recoils of tens of keV. The thick-target yields in Table C, scaled to ErD₃, give each photoneutron roughly 10⁻¹⁰–10⁻⁸ secondary D–D neutrons. That signal sits on a large photoneutron background, so spectroscopy above about 1 MeV is essential.
- **Grade C.** No independent replication or peer-reviewed follow-up is known to me. NASA technical memoranda and conference papers from 2021 onward exist but were not verified this session.
- **Device relevance:** low. The approach needs an MeV driver, which is out of scope. It is useful only as proof that dense deuterides plus energetic internal deuterons give standard D–D signals.

---

## 6. Schenkel et al. (LBNL) and the plasma/beam line of work

- **J. Appl. Phys. 126, 203302 (2019)**, "Investigation of light ion fusion reactions with plasma discharges" (Schenkel, Persaud, Wang, Seidl et al.; part of the Google programme). [m] Pulsed deuterium plasmas at applied voltages of roughly 0.3–2.5 kV bombarded Pd and Ti. D–D neutrons were measured with ³He proportional counters and scintillators. The reported yields and limits were consistent with conventional beam-target D–D, with **no evidence of anomalous enhancement**. I could not retrieve the exact detection limits.
- **Why kV plasmas sit at the edge of detectability [calc, Table C].** At E_d = 2 keV, 1 mA on PdD gives 10⁻⁴ n/s with no screening and about 0.8 n/s even at U_e = 800 eV. Only sheath voltages ≳ 3–5 kV at mA–A currents give robust neutron counts. For glow-discharge "LENR" claims below 1 kV, any real D–D neutron signal requires U_e ≫ 800 eV or a non-screening mechanism.
- **LBNL's ARPA-E continuation** (keV beams and plasmas on hydrides with in-situ loading) is the right experimental template. Its results were not verified here (§4).

---

## 7. Theoretical rate framework

### 7.1 Equations and constants

- **Cross-section:** σ(E) = [S(E)/E] · exp(−√(E_G/E)), where E_G = 2μc²(παZ₁Z₂)².
  - D+D: μc² = 937.8 MeV, so **E_G = 985.8 keV** [calc].
  - Other pairs: p+D 657 keV, D+T 1,182 keV, p+⁶Li 7.60 MeV, p+⁷Li 7.76 MeV, d+⁶Li 13.3 MeV, p+¹¹B 22.6 MeV.
- **Screening (constant shift):** σ_s(E) = σ_b(E + U_e) ≈ σ_b(E) · exp(πηU_e/E), with 2πη = √(E_G/E).
- **Bound pair:** λ = A · ρ₀ · P.
  - A = S c/(πα μc²) = **1.52×10⁻¹⁶ cm³/s** for S_tot = 109 keV·b [calc].
  - ρ₀ is the pair density where tunnelling starts.
  - P = exp(−√(E_G/U_eff)).
- **Free gas:** rate per deuteron λ_D = ½ n_D ⟨σ_s v⟩, Maxwell-averaged at 300 K. n_D = 6.8×10²² cm⁻³ for PdD.
- **Branches** (S-factors are Bosch–Hale S(0); energies at low E) [m]:

| Channel | Q (MeV) | Products | S(0) (keV·b) | Branching |
|---|---|---|---|---|
| D(d,n)³He | 3.269 | n 2.45 MeV + ³He 0.82 MeV | 53.7 | ~50% |
| D(d,p)T | 4.033 | p 3.02 MeV + T 1.01 MeV | 55.6 | ~50% |
| D(d,γ)⁴He | 23.847 | 23.8 MeV γ | — | ~10⁻⁷ of the p branch at E_cm ~ 10–100 keV |
| D(d,e⁺e⁻)⁴He | 23.847 − 1.022 | e⁺e⁻ continuum | — | Conventionally ≲10⁻⁹; Czerski/Dubey claim dominant below ~5 keV [v] |

  - The n/p ratio is ≈ 0.95–1 at low energy. The mean energy per fusion is 3.66 MeV [calc].
- **Charged-product ranges in Pd** [calc, Bethe estimate]: 3.02 MeV p ≈ **33 µm**; 1.01 MeV t ≈ 7 µm; 0.82 MeV ³He ≈ 2 µm.
- **Other S-factors** [m]: p+D (LUNA, Nature 587, 210, 2020) has S(0) ≈ 0.2 eV·b, about 5×10⁵ below D+D. p+⁷Li has S(0) ≈ 55–60 keV·b but an E_G 8× larger.

### 7.2 Benchmarks: Koonin–Nauenberg, Leggett–Baym, muon catalysis

- **Koonin & Nauenberg (Nature 339, 690, 1989):** λ(D₂ molecule) ≈ **3×10⁻⁶⁴ s⁻¹** per pair [m]. In this framework that corresponds to U_eff = 34 eV with ρ₀ = 1.5×10²⁶ cm⁻³ [calc; matches M0]. They also showed that raising the electron mass 5–10× would be needed to reach the 1989 claims (Jones et al., ~10⁻²³ s⁻¹).
- **Leggett & Baym (PRL 63, 191, 1989):** a rigorous upper bound on barrier penetration for a pair in an equilibrium many-body system, set by bulk properties. Static screening in Pd or Ti cannot approach 10⁻²³ s⁻¹ per pair. [m] The loopholes are **non-equilibrium, dynamic or rare-configuration** effects, which is exactly where §3.3 and §8 point.
- **Muon-catalysed fusion (grade A)** [m]: m_μ = 206.8 mₑ, τ_μ = 2.197 µs.
  - Muonic molecules are about 207× smaller than D₂ (hundreds of fm vs 0.74 Å).
  - Fusion rates: ddμ (J = 1) ≈ 4×10⁸ s⁻¹; dtμ ≈ 10¹² s⁻¹.
  - Up to ~120–150 fusions per muon in D–T.
  - Sticking: ω_s ≈ 0.4–0.9% for dt and ≈ 12% for the dd → ³He branch.
  - It is not net-energetic: a muon costs about 5 GeV, against 17.6 MeV × ~150 ≈ 2.6 GeV returned.
  - [calc] Scaling the D₂ calculation by 207 (ρ₀ × 207³, U_eff × 207 ≈ 7 keV) gives λ ≈ 1.5×10¹² s⁻¹ for s-wave ddμ. This matches measured muonic rates within p-wave suppression. The framework is validated across 76 orders of magnitude. **The only proven "cold fusion" works by a 207× length contraction that no electronic or lattice effect approaches.**

### 7.3 Computed rates vs screening energy at room temperature [calc]

Assumptions:
- kT = 25.9 meV and S_tot = 109 keV·b.
- "Const-shift" means the potential behaves as e²/r − U_e inside the screening radius.
- "Yukawa" means V = (e²/r)e^(−r/λ) with λ = e²/U_e, solved by exact WKB and re-expressed as an equivalent U_eff.
- The free-gas model lets deuterons approach freely; this is an upper-bound picture.

| U_e (eV) | P = e^(−√(E_G/U_e)) | λ per D, free gas, const-shift (s⁻¹) | U_eff if the potential is Yukawa (eV) | λ per D, free gas, Yukawa (s⁻¹) | λ bound pair, ρ₀ = 10²⁵ cm⁻³ (s⁻¹) | 1 cm³ PdD, const-shift: n/s (W) |
|---|---|---|---|---|---|---|
| 25 | 5.8×10⁻⁸⁷ | 2.3×10⁻⁷⁶ | 11 | 1.6×10⁻¹²² | 8.8×10⁻⁷⁸ | 8×10⁻⁵⁴ (9×10⁻⁶⁶ W) |
| 34 (≡ D₂) | 1.1×10⁻⁷⁴ | 4.3×10⁻⁶⁴ | 14 | 3.5×10⁻¹⁰⁴ | 1.7×10⁻⁶⁵ | 1×10⁻⁴¹ |
| 100 | 7.6×10⁻⁴⁴ | 2.8×10⁻³³ | 41 | 2.6×10⁻⁵⁷ | 1.2×10⁻³⁴ | 9×10⁻¹¹ (1×10⁻²² W) |
| 300 | 1.3×10⁻²⁵ | 4.6×10⁻¹⁵ | 123 | 5.0×10⁻²⁹ | 1.9×10⁻¹⁶ | 1.5×10⁸ (0.2 mW) |
| 800 | 5.7×10⁻¹⁶ | 2.1×10⁻⁵ | 331 | 6.9×10⁻¹⁴ | 8.7×10⁻⁷ | 6.8×10¹⁷ (0.8 MW) |
| 2000 | 2.3×10⁻¹⁰ | 8.2 | 841 | 4.9×10⁻⁵ | 0.35 | 2.7×10²³ (3×10¹¹ W) |

- **Inversion against null limits:** 10⁻²³ /D/s ↔ U_eff = **165 eV**; 10⁻²⁵ /D/s ↔ **147 eV** (const-shift, PdD).
- **Local steepness:** d ln λ / d ln U = ½√(E_G/U) = 85, 50, 29 and 18 at 34, 100, 300 and 800 eV.
- **Implication:** the Pd accelerator value (800 eV) cannot describe thermal pairs, even Yukawa-softened (2×10⁹ n/s/cm³, excluded by about 11 orders of magnitude). A 300 eV accelerator value, *if Yukawa-shaped*, is compatible with bulk nulls. A pure-cold signal therefore needs U_eff(thermal) ≈ 150–500 eV on some population of sites (M0 §3). Bulk material is excluded; only rare sites remain.

### 7.4 Beam and plasma regime: where screening is measurable [calc]

**Table B. Enhancement f = σ_b(E+U_e)/σ_b(E) at centre-of-mass energy E.**

| U_e \ E_cm | 0.5 keV | 1 keV | 2.5 keV | 5 keV | 10 keV | 20 keV |
|---|---|---|---|---|---|---|
| 25 eV | 2.9 | 1.5 | 1.10 | 1.04 | 1.01 | 1.00 |
| 100 eV | 48 | 4.3 | 1.5 | 1.15 | 1.05 | 1.02 |
| 300 eV | 1.1×10⁴ | 48 | 3.0 | 1.5 | 1.16 | 1.05 |
| 800 eV | 2.1×10⁷ | 3.0×10³ | 13 | 2.7 | 1.45 | 1.15 |
| 2000 eV | 4.6×10¹⁰ | 5.8×10⁵ | 157 | 8.8 | 2.4 | 1.4 |

**Table C. Thick-target D(d,n) neutron yield from 1 mA D⁺ on PdD₀.₇ (n/s).**

- Stopping power: Lindhard–Scharff electronic stopping ×1.4 (roughly SRIM-like, ±30% absolute) plus ZBL nuclear stopping. The proton yield is about 5% higher.
- The ×(U_e=0) ratios in the last two rows are the U_e = 300 and 800 eV yields relative to the unscreened yield.

| E_d (lab) | 1 keV | 2 keV | 3 keV | 5 keV | 10 keV | 20 keV |
|---|---|---|---|---|---|---|
| U_e = 0 | 2×10⁻¹⁰ | 1.0×10⁻⁴ | 0.034 | 11 | 4.0×10³ | 2.6×10⁵ |
| U_e = 25 eV | 7×10⁻¹⁰ | 1.6×10⁻⁴ | 0.043 | 13 | 4.2×10³ | 2.6×10⁵ |
| U_e = 300 eV | 5×10⁻⁶ | 7.5×10⁻³ | 0.43 | 41 | 6.6×10³ | 3.1×10⁵ |
| U_e = 800 eV | 0.024 | 0.83 | 9.7 | 240 | 1.4×10⁴ | 4.3×10⁵ |
| × (U_e=0) at 300 / 800 eV | 2×10⁴ / 10⁸ | 75 / 8,300 | 13 / 290 | 3.7 / 22 | 1.65 / 3.5 | 1.2 / 1.7 |

- With the M0 neutron threshold of 2×10⁻² reactions/s (³He bank, 14 days, 5σ), every entry at E_d ≥ 3 keV is detectable at 1 mA. For U_e ≥ 300 eV the threshold is reached near 2 keV.
- Charged-particle detection in vacuum is about 10× more sensitive (M0).

**Table D. S-factor multiplier K (e.g. a resonance) expressed as an equivalent U_eff at thermal energy.**

| Base U_eff | K = 10 | K = 10³ | K = 10⁶ | K = 10¹⁰ |
|---|---|---|---|---|
| 34 eV | 35 | 37 | 40 | 45 |
| 100 eV | 105 | 116 | 135 | 170 |
| 300 eV | 326 | 388 | 521 | 838 |

---

## 8. Structural and geometric dependences (the bridge to device design)

| Variable | Reported effect | Evidence | Mechanistic ambiguity | Design lever |
|---|---|---|---|---|
| **Host electronic type** | Metals ≫ insulators and semiconductors (Bochum survey, 58 samples) | B | Stopping power and D mobility differ too | Use conductive hosts: Pd, Pt, Ta, Zr |
| **Specific metal** | Pd, Pt, Co highest (Bochum); PdO > Pd > Fe > Au (Tohoku) | B for ranking, C for values | Surface prep dominates the inter-lab spread | Pd as primary host; Zr/Ti as the Szczecin reference |
| **D concentration** | Small U_e in hydride-forming metals at high x (20 °C); large when x drops at ~200 °C (Bochum 2006). ZrD₂ ≈ 300 eV under UHV (Berlin/Szczecin) | C/B, contradictory | Hydride formation localises electrons, or the analysis assumed the wrong D profile | Measure x in situ (resistance ratio, NRA); scan x rather than maximising it blindly |
| **Temperature** | U_e ∝ T^−½ in Pt, Co (Bochum 2005) | C, disputed | x and profile also change with T | 20–300 °C control, measured jointly with x |
| **Surface contamination and oxides** | C/O layers change fitted U_e strongly (Huke 2008; Kowalska 2023, 2025). PdO ≈ 2× Pd (Kasagi 2002). Au/Pd/PdO heterostructures give large enhancement (Yuki 1998) | B that it matters; C for direction | Oxides can accumulate D at interfaces (a density effect) or change screening | UHV ≤10⁻⁹ mbar; in-situ sputter clean; oxide as an *explicit* variable (split targets) |
| **Lattice defects (vacancies)** | U_e increases with vacancy-type defect concentration (PAS and XRD, Szczecin) | C | Defects trap D, changing the local density | Engineered defect density (ion damage, cold work) characterised by PAS; annealed half as control |
| **Thermal properties and D diffusivity** | Plateau yield differs between Zr, Ti and Pd (thermal-spike model) | C | Possible beam artifact | Include low-conductivity hosts; vary substrate heat sinking |
| **Loading flux (electrochemical)** | +15% conventional D–D rate with back-side electrochemical loading of Pd (Nature 2025) | C, conventional | Density, not screening | Back-side loading membrane geometry (a proven cold handle) |
| **Nanostructure, cracks, grain boundaries** | No credible accelerator measurement known to me. Fracto-emission and "fractofusion" claims (1980s–90s) are D | Gap | — | Worth testing inside the calibrated setup: nanocrystalline vs single-crystal split targets |

**Synthesis.** The only structural variables with more than single-group support are *metallicity* and *surface state*. Every other lever is C or worse. The M0 tail-dominance argument (the rate is set by the rarest, most extreme sites) makes defects, interfaces and oxide/metal boundaries the rational bet. The same rare sites also raise the risk of *artefactual* enhancement (D accumulation, stopping errors). **Every structural comparison must therefore be made in situ against a calibrated conventional-D–D baseline, with D profiles measured.**

---

## 9. Design implications, ranked by expected value

EV = P(credible, measurable signal) × credibility × information value.

1. **Anchor the device on a calibrated keV-deuteron probe of a thin deuterided film. EV highest.**
   - A glow-discharge cathode fall or ion source at 0.5–10 keV guarantees real D–D counts (Table C). This validates every detector in situ.
   - It turns "anomaly" into a measurable departure from conventional beam-target predictions.
   - It works in the 1–3 keV window where screening differences are ×10–10⁴ and robust against normalisation systematics.
   - The film should be at least as thick as the ion range (tens of nm) and **thinner than about 30 µm**, so that 3 MeV protons escape (tritons need ≲7 µm).
   - This probes the lattice; it is not the energy source, so it stays within the cold-lattice scope. The project lead should confirm this framing is acceptable.
2. **Multi-channel detection geometry, co-designed with the target.**
   - Si detectors in vacuum with direct line of sight: p 3.02 MeV, t 1.01 MeV, ³He 0.82 MeV.
   - A back-to-back **511 keV coincidence pair** (LaBr₃ or NaI) plus a plastic/thick-Si **MeV-electron telescope**. This directly tests the Dubey/Czerski e⁺e⁻ claim. If the claim is true, a neutron-only design under-counts the lowest-energy rate by ≥10×.
   - PSD neutron counters (EJ-309 or stilbene) plus a moderated ³He bank.
   - HPGe for activation and 23.8 MeV γ upper limits.
   - Coincidences (p–t, 511–511 keV) buy orders of magnitude in background rejection.
3. **Surface and interface control as first-class geometry.**
   - UHV, in-situ Ar⁺ cleaning, and XPS or NRA access to the same spot.
   - D depth profiling (NRA/ERDA) before and after each run.
   - Oxide/metal (PdO/Pd) and multilayer (Au/Pd/PdO) interfaces as *split-target* variables in the same exposure.
   Without this, any U_e result is uninterpretable at the ±2–3× level (see the table in §2.1).
4. **Defect engineering, with PAS characterisation, and in-run comparisons.** Compare damaged vs annealed and nanocrystalline vs single-crystal halves of one target on a rotating holder or split beam. Include identical H-loaded controls.
5. **Temperature (20–300 °C) and loading control with in-situ x measurement.** This discriminates Debye-like T^−½, thermal-spike, D-density and null explanations. Include back-side electrochemical or gas loading (the proven "cold handle").
6. **Beam-off, purely cold counting runs with the same detectors. Low probability, near-zero marginal cost.** This is the only direct test of equilibrium cold fusion. It needs U_eff(thermal) ≳ 150–500 eV on the active sites (§7.3; M0). Long H-vs-D background comparisons are essential.
7. **Li-in-Pd auxiliary targets (p+⁷Li → 2α, 8.7 MeV α each, neutron-free).** These show the largest reported U_e (keV) and are a sensitive screening probe in beam mode only. Their E_G is 8× that of D+D, so they are irrelevant at thermal energy.
8. **Avoid or deprioritise:**
   - Heat-only observables: 10¹²× less sensitive than particles, per M0.
   - Thick bulk cathodes as the sole sensor: products don't escape and loading is slow.
   - MeV-driven LCF variants: not cold, and the backgrounds swamp the signal.
   - Unshielded 511 keV measurements: cosmic pair production.

**What counts as credible:**
- ≥5σ excesses in **two independent channels** (e.g. 3.02 MeV protons plus PSD neutrons, or 511–511 keV coincidences plus an electron telescope).
- An H-control null.
- Correct scaling with deuteron energy and D profile.
- Reproduction on a second target.

---

## 10. Open questions

1. Why are the measured U_e 3–30× above theory? Answering this needs low-energy stopping powers for H in Pd and Zr measured *on the same targets*, and D profiles tracked during irradiation.
2. Is screening velocity-dependent (dynamic)? That determines whether keV-measured enhancements say anything at all about thermal pairs.
3. Does the 0⁺ threshold resonance exist? A gas-target e⁺e⁻ search at E_d = 5–20 keV by an independent lab (a LUNA-class or university accelerator) would be decisive. If it exists, what are E_R and Γ_tot?
4. Is the sub-2.5 keV "thermal-spike" plateau a beam artifact? This needs a magnetic-analysis/time-of-flight beam check and a gas-target comparison.
5. Do defects raise U_e intrinsically, or only by trapping D (a density effect)? PAS plus NRA on the same spot would separate the two.
6. Are there rare thermal-energy "hot sites" (crack tips, vacancy clusters) with U_eff ≳ 200 eV? There are no data. This is the central unknown that cold-device designs implicitly bet on.
7. ARPA-E final outcomes: fetch the ARPA-E programme summary and the teams' 2024–26 publications. This could not be done this session.
8. Kasagi's Yb value and the JINR/Tomsk ZrD₂ and TiD₂ values need retrieving. The Bochum table values tagged [m] should be checked against EPJA 19, 283 Table 1.

---

## 11. References

*Search links are given where I was not able to confirm a DOI this session.*

**Screening (accelerator)**

- Raiola F. et al., "Enhanced electron screening in d(d,p)t for deuterated metals", Eur. Phys. J. A 19, 283 (2004). https://doi.org/10.1140/epja/i2003-10125-0 [v]
- Raiola F. et al., "Enhanced d(d,p)t fusion reaction in metals", Eur. Phys. J. A 27 s01, 79 (2006). https://doi.org/10.1140/epja/i2006-08-011-0 [v]
- Raiola F. et al., Eur. Phys. J. A 13, 377 (2002) (Ta); Phys. Lett. B 547, 193 (2002); J. Phys. G 31, 1141 (2005) (temperature). https://scholar.google.com/scholar?q=Raiola+electron+screening+d(d,p)t+deuterated+metals
- Czerski K. et al., "Experimental and theoretical screening energies for the ²H(d,p)³H reaction in metallic environments", Eur. Phys. J. A 27 s01, 83 (2006). https://doi.org/10.1140/epja/i2006-08-012-y [v]
- Czerski K. et al., "Enhancement of the electron screening effect for d+d fusion reactions in metallic environments", Europhys. Lett. 54, 449 (2001). https://scholar.google.com/scholar?q=Czerski+2001+Europhysics+Letters+54+449
- Huke A. et al., "Enhancement of deuteron-fusion reactions in metals and experimental implications", Phys. Rev. C 78, 015803 (2008). https://doi.org/10.1103/PhysRevC.78.015803
- Kasagi J. et al., "Strongly enhanced DD fusion reaction in metals observed for keV D⁺ bombardment", J. Phys. Soc. Jpn. 71, 2881 (2002). https://doi.org/10.1143/JPSJ.71.2881
- Yuki H. et al., JETP Lett. 68, 823 (1998). https://scholar.google.com/scholar?q=Yuki+Kasagi+Lipson+1998+JETP+Letters+Pd+PdO
- Greife U. et al., Z. Phys. A 351, 107 (1995). https://scholar.google.com/scholar?q=Greife+1995+d(d,p)t+electron+screening+gas
- Cruz J. et al., Phys. Lett. B 624, 181 (2005). https://scholar.google.com/scholar?q=Cruz+2005+electron+screening+7Li(p,alpha)+environments
- Bystritsky V.M. et al. (JINR), dd reaction with pulsed Hall accelerator. https://scholar.google.com/scholar?q=Bystritsky+dd+reaction+Hall+accelerator+screening+ZrD2
- Assenbaum H.J., Langanke K., Rolfs C., Z. Phys. A 327, 461 (1987). https://scholar.google.com/scholar?q=Assenbaum+Langanke+Rolfs+1987+screening
- Kowalska A. et al., "Crystal lattice defects in deuterated Zr in presence of O and C impurities studied by PAS and XRD for electron screening effect", Materials 16, 6255 (2023). https://doi.org/10.3390/ma16186255 [v]
- Szczecin group, "Electron screening in deuteron–deuteron reactions on a Zr target with oxygen and carbon contamination", Materials 18, 1331 (2025). https://doi.org/10.3390/ma18061331 [v]
- Vesić J., "Electron screening in laboratory nuclear reactions" (review). https://inspirehep.net/files/595aafa9f963e840e0a7ce6969d2a3bd [v]
- "The role of the screening potential in the deuteron–deuteron thermonuclear reaction rates", Nucl. Phys. A (2025). https://www.sciencedirect.com/science/article/abs/pii/S0375947425002738 [v existence]

**Threshold resonance, new channel, thermal spike**

- Dubey R. et al., "Experimental signatures of a new channel of the deuteron–deuteron reaction at very low energy", Phys. Rev. X 15, 041004 (2025). https://doi.org/10.1103/chlp-b215 ; arXiv:2408.07567, https://arxiv.org/abs/2408.07567 [v]
- Czerski K. et al., "Deuteron–deuteron nuclear reactions at extremely low energies", Phys. Rev. C 106, L011601 (2022). https://doi.org/10.1103/PhysRevC.106.L011601 [v]
- Czerski K. et al., "Indications of electron emission from the deuteron–deuteron threshold resonance", Phys. Rev. C 109, L021601 (2024). https://doi.org/10.1103/PhysRevC.109.L021601 ; arXiv:2305.17101 [v]
- Czerski K. et al., Acta Phys. Pol. B 51 (2020), "Deuteron–deuteron reaction cross sections at very low energies". https://s3.cern.ch/inspire-prod-files-6/60a383800667f17c2e74a258230f7c75 [v]
- Czerski K. et al., Europhys. Lett. 113, 22001 (2016). https://doi.org/10.1209/0295-5075/113/22001 [m]
- "High-energy electron measurements with thin Si detectors", Measurement 228, 114392 (2024). https://www.sciencedirect.com/science/article/pii/S026322412400277X ; arXiv:2312.12446 [v]
- Czerski K. et al., "Observation of thermal deuteron–deuteron fusion in ion tracks", arXiv:2409.02112 (2024). https://arxiv.org/abs/2409.02112 [v]
- Czerski K. et al., "Thermal deuteron–deuteron fusion in metallic targets", arXiv:2605.27438 (2026). https://arxiv.org/abs/2605.27438 [v]
- CleanHME publications list. https://cleanhme.eu/?page_id=27 [v]
- Tilley D.R., Weller H.R., Hale G.M., "Energy levels of light nuclei A = 4", Nucl. Phys. A 541, 1 (1992). https://scholar.google.com/scholar?q=Tilley+Weller+Hale+1992+A%3D4

**Programmes, LCF, plasma and beam work**

- ARPA-E LENR selections: https://executivegov.com/2023/02/arpa-e-funds-8-projects-to-break-impasse-in-low-energy-nuclear-reactions-field/ ; https://www.greencarcongress.com/2023/02/20230218-lenr.html ; https://www.ans.org/news/article-4769/arpae-picks-eight-teams-to-proveor-debunklowenergy-nuclear-reactions/ ; https://news.newenergytimes.net/2023/03/01/u-s-department-of-energys-arpa-e-funds-low-energy-nuclear-reaction-research/ [v]
- ARPA-E FY2023 Annual Report (Aug 2025). https://arpa-e.energy.gov/sites/default/files/2025-09/ARPA-E%20FY%202023%20Annual%20Report.pdf [v existence]
- Pines V. et al., Phys. Rev. C 101, 044609 (2020). https://doi.org/10.1103/PhysRevC.101.044609
- Steinetz B.M. et al., Phys. Rev. C 101, 044610 (2020). https://doi.org/10.1103/PhysRevC.101.044610
- Schenkel T. et al., "Investigation of light ion fusion reactions with plasma discharges", J. Appl. Phys. 126, 203302 (2019). https://scholar.google.com/scholar?q=Schenkel+2019+Investigation+of+light+ion+fusion+reactions+with+plasma+discharges
- Berlinguette C.P. et al., "Revisiting the cold case of cold fusion", Nature 570, 45 (2019). https://doi.org/10.1038/s41586-019-1256-6
- UBC (Berlinguette group), "Electrochemical loading enhances deuterium fusion rates in a metal target", Nature (2025). https://scholar.google.com/scholar?q=Electrochemical+loading+enhances+deuterium+fusion+rates+in+a+metal+target [m]
- Metzler F. et al., "Known mechanisms that increase nuclear fusion rates in the solid state", New J. Phys. (2024). https://scholar.google.com/scholar?q=Known+mechanisms+that+increase+nuclear+fusion+rates+in+the+solid+state [m]
- Naranjo B., Gimzewski J.K., Putterman S., "Observation of nuclear fusion driven by a pyroelectric crystal", Nature 434, 1115 (2005). https://doi.org/10.1038/nature03575 (a replicated tabletop D–D neutron source, useful for calibration)

**Theory and benchmarks**

- Koonin S.E., Nauenberg M., "Calculated fusion rates in isotopic hydrogen molecules", Nature 339, 690 (1989). https://doi.org/10.1038/339690a0
- Leggett A.J., Baym G., "Exact upper bound on barrier penetration probabilities in many-body systems: application to 'cold fusion'", Phys. Rev. Lett. 63, 191 (1989). https://doi.org/10.1103/PhysRevLett.63.191
- Jones S.E. et al., Nature 338, 737 (1989). https://doi.org/10.1038/338737a0
- Jones S.E., "Muon-catalysed fusion revisited", Nature 321, 127 (1986). https://doi.org/10.1038/321127a0
- Breunlich W.H., Kammel P., Cohen J.S., Leon M., "Muon-catalyzed fusion", Annu. Rev. Nucl. Part. Sci. 39, 311 (1989). https://doi.org/10.1146/annurev.ns.39.120189.001523
- Bosch H.-S., Hale G.M., "Improved formulas for fusion cross-sections and thermal reactivities", Nucl. Fusion 32, 611 (1992). https://doi.org/10.1088/0029-5515/32/4/I07
- Mossa V. et al. (LUNA), "The baryon density of the Universe from an improved rate of deuterium burning", Nature 587, 210 (2020). https://doi.org/10.1038/s41586-020-2878-4
- Ichimaru S., "Nuclear fusion in dense plasmas", Rev. Mod. Phys. 65, 255 (1993). https://doi.org/10.1103/RevModPhys.65.255
- Project model: `docs/models/M0-rate-budget.md`, `sim/m0_rate_budget.py`.

---

## Appendix: core of the R3 calculation (Python/numpy/scipy)

```python
alpha=1/137.036; mu=1875.612e6/2; c=2.998e10; kT=0.02585; nD=6.8e22
E_G=2*mu*(np.pi*alpha)**2                        # 985.8 keV
S=(53e3+56e3)*1e-24                              # eV cm^2
A=S*c/(np.pi*alpha*mu)                           # 1.52e-16 cm^3/s (bound pair: lam=A*rho0*P)
P=lambda U: np.exp(-np.sqrt(E_G/U))
# free gas, Maxwell-averaged, rate per deuteron
f=lambda E: S/E*np.exp(-np.sqrt(E_G/(E+Ue)))*np.sqrt(2*E/mu)*c*2*np.sqrt(E/np.pi)*kT**-1.5*np.exp(-E/kT)
lam_D=0.5*nD*quad(f,1e-8,60*kT)[0]
# Yukawa mapping: V=(e2/r)exp(-r*Ue/e2); WKB exponent G -> Ue_eff=E_G/G^2
# beam enhancement: f(E)=exp(sqrt(E_G/E)-sqrt(E_G/(E+Ue)))
# thick target: Y(E0)=int_0^E0 x*sigma_s(E/2)/eps_PdDx(E) dE,
#   eps = 1.4*LSS_e(Pd) + ZBL_n(Pd) + x*(LSS_e(D)+ZBL_n(D))
```
