# Red team 0: review of the plan before any design is committed

**Scope.** `00-charter.md`, `01-master-plan.md`, `ADR-001`, `M0` (doc and `sim/m0_rate_budget.py`), R1–R7 and the model briefs, plus the emerging design direction: a cold dual-chamber thin Pd foil, loaded electrochemically from the back, with a vacuum front face facing Si telescopes, neutron and 511 keV detectors, and deuterium-flux modulation ("C3/DFM").
**Date:** 2026-09-29. **Stance:** adversarial. A finding is listed only if it would change a decision.

**Method.** Every M0 number was recomputed independently with scripts kept outside the repo. The engineering checks (membrane stress, exit-face loading, gas load, capacitive coupling, radon plating, ⁴He sensitivity) are first-principles order-of-magnitude calculations; the Appendix gives the formulas so anyone can reproduce them.
Tags:
- **[BK]**: from the reviewer's background knowledge, not re-verified in this session (web access blocked).
- **[calc]**: computed for this review.

Severity:
- **CRITICAL**: the plan as drafted cannot meet the charter's own goal without changing this.
- **MAJOR**: changes a design or ranking decision.
- **minor**: fix the text or code; the decision stands.

---

## 0. Summary: the five findings that matter most

1. **CRITICAL: the detector-facing membrane cannot be highly loaded, carry a large flux and have a clean exit face all at once.**
   - *Why the bulk follows the exit face.* A 10–50 µm Pd foil can carry 1–5 A cm⁻² of D per unit of loading difference Δx across it. The absorbed electrochemical flux is ≤ 50 mA cm⁻², so Δx ≤ 0.02 [calc]. The whole membrane therefore sits at whatever loading the exit face allows.
   - *Clean exit face.* A clean exit face in vacuum (sticking coefficient 10⁻⁶–0.3) sits at x_exit ≈ 10⁻⁵–10⁻² by detailed balance, so the entire foil is α-phase Pd, including the electrolyte face.
   - *What x_exit ≥ 0.9 would take.* The exit face would need to be 10¹⁰–10¹² times less active than clean Pd. That is a seal, and a seal removes the through-flux.
   - *Where the plan contradicts itself.* It assumes all three properties: C3 in the master plan, the segmented "same-flux" skins in the M6 brief, and the bake-out plus Ar⁺ sputter cleaning in R7/R3.
   - *What to do.* Pick two, explicitly, and build the apparatus for that regime (§3.1).
2. **CRITICAL: as drafted, the DFM reproduces none of the better-graded claimed conditions at the face its detectors see.**
   - *Not reproduced:* SRI/ENEA loading ≥ 0.9 held for weeks; SPAWAR co-deposit; Iwamura's entry-face multilayer with 50–110 mA cm⁻² net permeation; Kitamura's 200–300 °C nanopowders; Clean Planet's H₂ at 500–900 °C.
   - *What it does test:* Lipson/NTT-style desorption (C/D) and Czerski's beam-free 511 keV claim (D).
   - *Consequence:* a null would be dismissed as "regime not reached", which is exactly how the 2019 Google null is treated.
   - *Fixes (§3.3):* an 8–12 µm membrane that lets the Si see the entry face; an Au-capped exit; co-deposition on the entry face; a gas-entry Iwamura variant.
3. **MAJOR: ADR-001 optimises the wrong quantity.**
   - *Sensitivity is a weak lever.* By M0's own steepness (d ln λ / d ln U = 35–50), 100× better sensitivity is worth only 10–14 % in U_eff. Under a log-uniform prior spanning 40 decades of anomaly strength, it adds about 5 % detection probability [calc].
   - *Two levers are order-unity:* reproducing the claimed conditions, and running several samples (8 samples at p = 0.2 each give P = 0.83, against 0.20 for one).
   - *Blind spot.* The objective cannot see H2 when the energy goes into the lattice.
   - *Missing hypotheses* (§2.3): the H1/H2 split omits the e⁺e⁻ channel, fracto-driven "hot" D–D, transmutation (H3) and H-involved reactions.
4. **MAJOR: ⁴He is the missing H2 "particle detector", and the DFM is unusually well suited to it.**
   - *Sensitivity.* Furnace or laser extraction of the foil, plus a static, all-metal, getter-cleaned front volume, reaches ~10¹⁰ atoms [BK]. That equals 0.04 J of D+D→⁴He, about 7×10⁵ times the reach of a 10 mW calorimeter over 30 days [calc].
   - *Why the DFM suits it.* The Pd membrane blocks He, which isolates the front volume from the air-equilibrated electrolyte (1.2×10¹⁴ He per 100 mL, R5).
   - *The plan's error.* It demotes this channel to "secondary".
5. **MAJOR: the build will fail mechanically and electrically unless it is redesigned.**
   - *Mechanical.* An annealed 25 µm × 20 mm membrane at 1 atm carries ≈ 250 MPa, against a yield stress of 35–70 MPa [BK], and bulges ≈ 0.5 mm [calc]. α/β cycling under that load ratchets towards rupture, which floods the detector vacuum with electrolyte.
   - *Electrical.* The cathode membrane couples **3–6 MeV-equivalent per volt** into a Si detector 0.5–2 cm away [calc]. Every current-modulation edge, and every 10 mV of bubble noise (≈ 60 keV), lands exactly in phase with the flux lock-in.
   - *Fixes (§3.2, §3.4):* a pressure-balanced 1 atm D₂ front chamber or a support grid; a grounded membrane plus a grounded screen; a flux-blocked twin.

**Other MAJOR items:**
- The charter's 600 K cap makes C5 (Clean Planet, 500–900 °C) impossible to run inside the charter.
- R3, R4 and R7 all recommend a keV probe beam; the charter forbids beams as an energy source. The plan has not ruled on a *diagnostic* beam.
- Nothing in the plan gives an in-situ positive control showing that the Si channel would see D–D protons from this foil.

**M0 verdict.** The physics and code are essentially correct. No numerical error changes a conclusion. Section 1 lists the corrections: the D₂ ρ₀ concept, the 1 fm WKB cutoff, the "300–800 eV does not apply" wording, the matched-control and systematics factors, and the exponent bookkeeping.

---

## 1. M0 audit (`sim/m0_rate_budget.py`, `docs/models/M0-rate-budget.md`)

The code was re-run from a scratch copy: its output is byte-identical to `figs/m0_output.txt`. Note that running the script in place **rewrites tracked files** in `docs/models/figs/`; make OUT configurable.

| # | Item (location) | M0 | Independent value | Verdict | Changes a conclusion? |
|---|---|---|---|---|---|
| 1 | Gamow energy (l. 37) | 985.8 keV | 985.77 keV | Correct | No |
| 2 | S-factors (l. 40–41) | S_n = 55, S_p = 57 keV·b | Bosch–Hale S(0): 53.7 / 55.6 keV·b [BK] | ~2.5 % high. Branching ratio unchanged (0.491). ln A shifts by 0.025 e-folds | No |
| 3 | Bound-pair constant A = S c/(πα μc²) (l. 47) | 1.56×10⁻¹⁶ cm³ s⁻¹ | 1.562×10⁻¹⁶ (1.524×10⁻¹⁶ with Bosch–Hale) | Correct: this is the Jackson/Koonin–Nauenberg form, A = lim σv/C²(η) | No |
| 4 | D₂ calibration ρ₀ (l. 111–114) | ρ₀ = (πx₀²)^−3/2 = 1.5×10²⁶ cm⁻³ → U_eff(D₂) = 34.0 eV | ρ₀ = 1/(√π x₀ · 4πR_e²) = 7.7×10²³ cm⁻³ → U_eff = 36.3 eV, X = 164.9 (not 170.2) | **Conceptual error.** The code treats the D–D relative wavefunction as a 3D oscillator centred at r = 0. It is a radial shell at R_e = 0.74 Å, so ρ₀ is 196× too high. The "calibration" is also a one-parameter fit that reproduces its own input, so it validates nothing | No: the gap shrinks by ~5 e-folds |
| 5 | Koonin–Nauenberg 3×10⁻⁶⁴ s⁻¹ | cited | 3×10⁻⁶⁴ s⁻¹ [BK] | Correct citation | No |
| 6 | Yukawa WKB (l. 58–71): limits, singularity, turning point | Integral from R_n = 1 fm to r_tp; geometric split | Integral from 0, using r = r_tp s² to remove both endpoint singularities: 2G is larger by 1.05 e-folds at every U_e | The turning point and quadrature are correct (agreement to 10⁻² at the same lower limit). The 1 fm cutoff is inconsistent with the point-Coulomb Gamow factor the S-factor is defined against. Corrected table: U_e = 100/300/500/800/1000/2000 → U_eff = **41/120/199/318/397/792** eV (M0: 41/123/206/331/414/841). The mapping is robust to the assumed relative energy: ±3 % for E = 0.025–0.5 eV | No: 2–6 % in U_eff |
| 7 | "Accelerator U_e ~300–800 eV does NOT apply at thermal energies" (l. 169; M0 §3) | Stated as a blanket exclusion | A Yukawa U_e of 300 eV gives a free-gas rate of 1.7×10⁻²⁹ /D/s, **below** the 10⁻²⁵ null. The boundary is U_e ≈ **370 eV** | **Wrong for U_e ≲ 370 eV.** Tohoku's Pd value (310 eV) is allowed; only Bochum-class 500–800 eV is excluded. R3 §7.3 already says this | Not for M0. It matters for M1 scenario (ii), which must run both 310 and 800 eV |
| 8 | Free-gas Maxwell average (l. 75–84) | ⟨E^−½⟩ = 2/√(πkT), with σ = (S/E)·P(E+U) | Full quadrature agrees within 1.3 % | Correct. The prefactor convention matters: S/E (asymptotic 1/k², used by M0) versus Assenbaum's σ_b(E+U) = S/(E+U) gives a null bound of **147 vs 185 eV**. M0's choice is the physical one; state it. "Upper bound" is not rigorous: the rigorous bound is Leggett–Baym | No |
| 9 | Detection threshold (l. 133–135) | Gaussian, known background: 5√(B t) | With the **equal-time matched control the charter demands**: n 2.9×10⁻², p 5.5×10⁻³ s⁻¹. Asimov for protons (b = 24 counts): 3.1×10⁻³. **Neutron systematic floor**: at δB/B = 1–2 % of 0.05 cps, S ≥ 0.05–0.1 fusions s⁻¹ however long the run (R5 §10.1) | Optimistic by ×1.2–5 | No: ≤ 1.6 e-folds, ≤ 4 % in U_eff |
| 10 | Si solid angle (l. 138) | 5 % | R5: 18–31 % for a 12 mm-radius disk at 5–10 mm | Conservative by ×4–6 | No |
| 11 | Neutron dose (l. 146–152) | 400 pSv cm² → 9.6 Sv/h at 1 m | ICRP-74 H*(10) at 2.5 MeV ≈ 416 pSv cm² [BK] → 10.0 Sv/h (bare point source; room return adds tens of %) | Correct | No |
| 12 | Heat equivalence | 1.71×10¹⁰ s⁻¹ per 10 mW (standard branching) | Correct. But the **charter's H2 row** quotes "~10¹⁰ s⁻¹ per 10 mW"; at 23.85 MeV per ⁴He it is **2.6×10⁹** | Charter typo | No. See §2.2 for the real omission: M0 has no ⁴He channel |
| 13 | Exponent budget (l. 171–181; charter "Honest prior" item 2) | Gap 97–118 e-folds; "engineering buys 25–30" | The X-gap mixes ρ₀ = 1.5×10²⁶ (D₂) with 10²⁵ (needed). The true ln-rate gap is **94–115**. The 10⁹ "more sites" lever is already inside the N = 10²¹ figure, so the charter double-counts it. With all 10²¹ pairs active *and* 100× better detection, **89 e-folds remain** | Wording error | No; the conclusion gets stronger |
| 14 | Tail-dominance MC (l. 190–197) | Top 0.1 % of sites carry 16 / 59 / 90 % | **Verified:** the population answer by quadrature is 15.4 / 58.2 / 90.6 %, and the MC is converged (N = 2×10³ to 2×10⁷ changes it by ≤ 6 %) | Correct for a Gaussian, but the number is an *assumption about the tail*. The rate-weighted peak sits at U* = 133 eV (σ = 10) or 181 eV (σ = 20), at 3.3–4.1σ, i.e. a site fraction of 6×10⁻⁴–2×10⁻⁵. A physical cutoff or a heavier tail changes everything | No, but see the reframing below |
| 15 | Steepness ½√(E_G/U) | 70 / 50 / 35 / 29 | Correct | — | — |

**Reframing of the tail result.** Rows 14 and 15 together say that **site quality beats site count**:
- Raising the mean U_eff of a site class by 10 % (100 → 110 eV) multiplies the rate by ×113.
- Matching that with geometry means 113× more visible area.

Processing (defects, loading, surface chemistry) is therefore the dominant lever, and the area a configuration exposes to detectors is secondary.

M0 should also:
- add the ⁴He channel (§2.2);
- stop using "14 d" in M0 but "30 d" in ADR-001.

---

## 2. Objective function (ADR-001)

### 2.1 Sensitivity to a localised anomaly is a weak lever, by M0's own numbers

ADR-001 scores a geometry by "the minimal enhancement factor it would detect at 5σ in 30 days". Three quantitative problems follow.

- **Exponential physics compresses sensitivity.** With d ln λ / d ln U = 35–50 [calc]:
  - 10× better sensitivity is worth only **5–7 %** in the U_eff the device can reach;
  - 100× is worth 10–14 %;
  - 1000× is worth 15–22 %.

  Candidate configurations differ in N_visible × ε by 10–1000×, so the M1 "configuration × site class" matrix will differ between configurations by only ~5–20 % in U_eff. The decision will be almost insensitive to it.
- **Under an honest prior the value is small.** If the anomaly strength K (over the standard rate) is log-uniform over 40 decades (10³⁰–10⁷⁰), 10× / 100× / 1000× better sensitivity adds **2.5 / 5 / 7.5 %** detection probability [calc].
- **The order-unity levers are elsewhere:**
  - *Conditions.* A claimed effect either happens in the configuration or it does not; that factor multiplies everything else.
  - *Sample multiplicity.* SRI, ENEA and Energetics report that only a fraction of vetted cathodes ever work (the hidden variable "M" in R1). With p = 0.2 per sample, 1 / 4 / 8 samples give P(≥ 1 active) = **0.20 / 0.59 / 0.83** [calc].
  - *Saturation.* For every claim with a stated magnitude, the drafted detectors already reach below 1 % of it:

  | Claim | Claimed magnitude | Drafted reach |
  |---|---|---|
  | SRI-scale heat–He | 0.1–1 W → 10¹⁰–10¹¹ He s⁻¹ | about 4×10³ reactions s⁻¹ via ⁴He |
  | Lipson protons | ~10⁻²–1 s⁻¹ [BK] | about 3×10⁻³ s⁻¹ via Si |
  | Iwamura products | 10¹²–10¹⁴ atoms cm⁻² | about 10¹⁰ cm⁻² via ICP-MS |

  Beyond that point more sensitivity buys nothing for testing claims.

### 2.2 H2 is invisible to the chosen primary instruments; ⁴He, not heat, is H2's particle detector

- **What H2 predicts.** The heat–⁴He claims (China Lake: 18/21 runs; SRI M4; ENEA) together with Hagelstein's bound (charged products ≲ 20 keV) describe a channel that deposits 23.85 MeV in the lattice. Si telescopes, ³He counters and 511 keV pairs are blind to it by construction.
- **What ADR-001 keeps for H2.** Only heat, as a "secondary" channel, and in a thin-membrane cell heat is at its weakest.
- **What M0's "heat is 10¹² times less sensitive" comparison actually covers.** It is correct for H1. For H2, the relevant comparison is ⁴He accumulation [calc; limits BK]:

  | ⁴He detection limit | Integrated D+D→⁴He | 30-day mean power | Versus a 10 mW calorimeter |
  |---|---|---|---|
  | 10⁹ atoms (noble-gas-lab blank) | 3.8 mJ | 1.5 nW | ×7×10⁶ |
  | 10¹⁰ atoms | 38 mJ | 15 nW | ×7×10⁵ |
  | 10¹² atoms (conservative engineering) | 3.8 J | 1.5 µW | ×7×10³ |

- **Why the DFM is unusually good for this:**
  1. Pd does not transmit He. The membrane separates the electrolyte, which carries air-equilibrated He and permeates glass and elastomers, from a front volume that can be all-metal and static.
  2. The foil is small (≈ 0.1 g for 25 µm × 3 cm²). Furnace or laser extraction at the end of the run recovers the *retained* He that defeated earlier heat–He tests. Spot extraction (UV laser microprobe [BK]) can even *localise* ⁴He to ~100 µm, which is the H2 analogue of "localised anomaly".
  3. A getter (Zr–V–Fe) removes D₂ by ≥ 10⁶ but not He, so the ⁴He/D₂ mass-4 overlap (m/Δm ≈ 160) is solved without a high-resolution spectrometer at the cell.
- **Conclusion.** ⁴He (and ³He, which needs m/Δm ≥ 510 or getter removal of HD) must be **co-primary**, with a pre-registered foil-extraction protocol and an H₂O-twin foil as the blank.

### 2.3 The H1/H2 split is incomplete

| Hypothesis | Products / signature | In charter? | Sees it in the drafted DFM |
|---|---|---|---|
| H1a: enhanced screening, normal branching | p 3.02, t 1.01, ³He 0.82, n 2.45 MeV | Yes | Si (exit face only), ³He bank |
| H1b: *hot* micro-sources: fracto-emission, desorption-field acceleration (M1 Part B) | The same products, **correlated with cracking or transients**. Conventional physics, so it is a **confound** for any H1 positive | No | Nothing tags cracking. Add an acoustic-emission sensor ($1–3k) |
| H1c: threshold resonance, e⁺e⁻ channel (Czerski) | 511–511 keV coincidences, 3–23 MeV e± and bremsstrahlung, **no p/n**. Not "heat without radiation": ~5×10¹¹ annihilation photons per W (R3) | No (M7 only) | 511 pair, provided the shielding and site are adequate |
| H2: lattice-coupled D+D→⁴He | ⁴He at 2.6×10¹¹ J⁻¹, heat, possibly ≤ 20 keV particles or soft X-rays | Yes | **None as drafted** (§2.2) |
| H3: transmutation (Iwamura Cs→Pr, Sr→Mo), low-Z "fission daughters" (the MIT/ARPA-E focus) | Isotopically anomalous surface species at 10¹²–10¹⁴ cm⁻² | No | None. Needs pre/post ToF-SIMS and HR-ICP-MS, **isotopically tagged targets** (e.g. ⁸⁶Sr → ⁹⁴Mo, per R2) and fiducial-marked spots |
| H4: H-involved reactions (p+d→³He, Ni–H claims, ³He in CNZ) | ³He, heat with H₂ | No | None. Note that the H₂O twin is a "control" only if H is inert, which the Ni–H claimants contest. For Pd–D that is acceptable; state the assumption |

### 2.4 A null is only informative if the claimed conditions were reproduced

The charter calls the design "the most informative null if nothing is there". The informativeness of a null for claim c is

  P(conditions of c reproduced) × P(detect at m_c/100 | conditions).

The first factor is ≈ 0 for every B-grade claim in the drafted DFM (§3.3). Google 2019 teaches the same lesson: a null outside the claimed regime is dismissed.

### 2.5 Proposed replacement (ADR-002)

For design d:

  Score(d) = Σ_c w_c · P(cond_c | d) · [1 − (1 − p_c)^{N_samples}] · P(det(m_c/100) | d, cond_c) · Cred(d)
       + w₀ · S_generic(d)

where:
- c ranges over the claims, each with an evidence-graded weight w_c (B > C > D);
- p_c is the per-sample success rate the claimants report;
- S_generic is the ADR-001 localised-anomaly sensitivity, kept but with a small weight (w₀ ≈ 0.1) for anomalies nobody has claimed;
- P(det) is **capped at 1** once the design reaches 1 % of the claimed magnitude.

Suggested reweighting of the master-plan matrix:

| Criterion | Current weight | Proposed weight |
|---|---|---|
| H1 sensitivity | 0.25 | 0.10, saturating as above |
| Claim-condition reproduction (replaces "evidence base" and absorbs site multiplicity) | 0.10 + 0.20 | 0.30 |
| H2 via ⁴He and calorimetry | 0.15 | 0.20 |
| Independent samples / statistics | — | 0.10 |
| Credibility | 0.20 | 0.20 |
| Safety / cost | 0.10 | 0.10 |

---

## 3. The emerging design direction (C3/DFM)

**Credit where due.** Energy-resolved charged-particle spectroscopy in vacuum is the most artifact-resistant nuclear channel available: intrinsic background ~2×10⁻⁵ cps in-window is plausible, U/Th in Pd contributes ≪ 10⁻⁶ cps, and Si(n,α) from cosmic neutrons ~10⁻⁶ cps [calc]. Electrochemistry as a switchable knob is the right *kind* of control, as Thunderbird showed. The problems below are physical, not stylistic.

### 3.1 Loading at the exit face (CRITICAL)

**Diffusion pins the bulk to the exit face.** Using D_D ≈ 5×10⁻⁷ cm² s⁻¹ (R1; the dilute tracer value, so the β-phase chemical diffusivity only makes this stronger):

| Thickness | Diffusion capacity D·n_Pd/L | Current equivalent |
|---|---|---|
| 10 µm | 3.4×10¹⁹ D cm⁻² s⁻¹ per unit Δx | 5.4 A cm⁻² |
| 25 µm | 1.4×10¹⁹ | 2.2 A cm⁻² |
| 50 µm | 6.8×10¹⁸ | 1.1 A cm⁻² |

Even if all of 50 mA cm⁻² were absorbed, Δx across 25 µm would be 0.023. So **x_entry ≈ x_exit + ≤ 0.02**: the membrane's loading is set by the exit face alone.

**The exit face collapses if it is clean (detailed balance).**
- Desorption equals s · 2Z(P_eq), with Z(1 atm D₂) = 7.7×10²³ cm⁻² s⁻¹.
- α-phase Sieverts law: x_α,max ≈ 0.015 at a plateau pressure of ~0.04 atm [BK].

| Sticking s | J = 10¹⁶ D cm⁻² s⁻¹ | J = 3×10¹⁷ D cm⁻² s⁻¹ |
|---|---|---|
| 0.3 (clean Pd) [BK] | x_exit ≈ 1×10⁻⁵ | 6×10⁻⁵ |
| 10⁻⁴ | 6×10⁻⁴ | 3×10⁻³ |
| 10⁻⁶ | 6×10⁻³ | reaches the plateau (β onset only) |

- Holding x_exit = 0.9 (P_eq ≈ 3×10³ atm, R1) needs s_eff ≲ 2×10⁻¹² at J = 10¹⁶: **10¹¹ times less active than clean Pd.**
- Empirically, contaminated, air-exposed surfaces *do* reach this: PdD₀.₉ foils deload at ~10¹⁵–10¹⁶ D cm⁻² s⁻¹ per face [calc from typical hour-scale deloading, BK]. But such surfaces are uncontrolled and drift as PdO is reduced and carbon builds up.

**Consequences for documents already written:**
- R7 §7.1 and R3 §9 recommend bake-out and in-situ Ar⁺ sputter cleaning. In a beam experiment implantation replenishes the surface D. In the cold DFM, cleaning *maximises* exit desorption and drives the whole foil into α phase. Consistent with this, Thunderbird's electrochemistry raised the rate by only 15 % [BK interpretation].
- **The M6 brief item 6 premise is false.** Sectors with different skins do *not* "share the same deuterium flux". Sectors are mm-scale against a 25 µm thickness, so the flux is locally 1-D, and each skin sets its own x_exit and J. Sector comparisons then conflate chemistry with loading and flux.
- **The M6 brief's "Au overlayer as the non-hydriding null" is inverted.** Au is a strong permeation barrier (it was used for exactly that in NTT 1990), so the Pd under a 20–50 nm Au sector is the *most* highly loaded region. 3 MeV protons from beneath it escape with a few keV of loss.

**Fix: choose the regime explicitly and build for it.**
- **C3-L (loading regime).** Exit face capped with 20–50 nm Au, or an ALD oxide at 5–20 nm, both transparent to MeV protons; use SRI's protocol. The whole foil is then at x ≳ 0.9 with near-zero steady flux. Get "flux" from current steps (McKubre's |dx/dt| term), not through-permeation.
- **C3-F (flux regime).** Clean exit, α-phase or low β-phase Pd, high through-flux. This tests Lipson-type desorption and the Czerski-511 claim. Declare it low-loading.
- **Diagnostics for both.** Measure permeation flux (calibrated pressure rise or RGA) together with 4-wire R/R₀, which gives the average x. Knowing the diffusion capacity, these two yield both faces' loading.
- **Make M3's DFM design chart a gating deliverable** before any C3 hardware is ordered.

### 3.2 Mechanics, pressure differential and embrittlement (MAJOR)

**Stress under 1 atm.** Hencky clamped-membrane solution, E_Pd = 121 GPa [calc]:

| Diameter × thickness | σ_max | Centre bulge |
|---|---|---|
| 10 mm × 25 µm | 155 MPa | 0.18 mm |
| 20 mm × 25 µm | 247 MPa | 0.46 mm |
| 25 mm × 50 µm | 180 MPa | 0.49 mm |
| 20 mm × 10 µm | 455 MPa | 0.62 mm |

- Annealed Pd, which is what the ENEA and SRI recipes produce, yields at ~35–70 MPa; cold-worked Pd at ~200–300 MPa [BK].
- The membrane therefore yields on first pressurisation.
- α↔β cycling (≈ 3.5 % linear misfit) under load causes transformation-induced ratcheting. Industry keeps pure Pd membranes above ~300 °C, or uses Pd–Ag23, for this reason [BK].
- Rupture floods the detector vacuum with LiOD/D₂O, destroying the Si detectors and pumps.

**Fixes, in order of preference:**
1. **Pressure-balance the front with ≈ 1 atm D₂** (§4, item 4).
   - Δp → 0.
   - The exit face is clamped at x ≈ 0.6–0.7 (β phase) instead of ~10⁻⁴.
   - The front becomes a static ⁴He accumulator: bled aliquots go through a getter to the mass spectrometer, with He accounted for by mass balance.
   - 1 atm D₂ is nearly transparent to D–D products: a 3 MeV proton loses ~20 keV cm⁻¹, a 1 MeV triton ~0.1 MeV cm⁻¹, 0.82 MeV ³He ~0.3 MeV cm⁻¹ (rough electron-density scaling [calc]).
   - Si detectors run at 1 atm; Paschen breakdown at 76 Torr·cm is kV-scale [BK].
   - Keep the flight path ≤ 5 mm for ³He.
2. If vacuum is mandatory, add a vacuum-side support grid with 100–200 µm apertures and ≥ 70 % open area: σ falls to 7–21 MPa. Include the grid bars in the escape Monte Carlo (M5), and count their own α emission.
3. Use Pd–Ag23, which has no miscibility gap at room temperature [BK], for the flux regime only (it is not a claimed-active host), and pure annealed Pd only in the pressure-balanced configuration.
4. Add a hardware pressure interlock: a fast gate valve plus Si-bias trip at a ΔP rise of 10⁻³ mbar per second.
5. Pre-screen every membrane for pinholes with H₂/forming-gas leak detection, **not He**, which would contaminate the ⁴He channel.

### 3.3 Does the DFM reproduce any claimed positive condition? (CRITICAL)

| Claim (grade) | Claimed conditions | Drafted DFM | Minimal fix |
|---|---|---|---|
| SRI heat–⁴He (B−/C) | Bulk x ≥ 0.875–0.9, i ≥ 0.1–0.4 A cm⁻², weeks, D flux, ⁴He | Visible face at x ~10⁻⁴ if clean (§3.1); no ⁴He | C3-L: Au-capped exit, **8–12 µm membrane** so Si sees the *entry* face (protons emerge at 2.55 / 2.42 / 2.28 MeV through 8 / 10 / 12 µm [calc]); SRI protocol; foil ⁴He extraction |
| ENEA foils (C) | 50 µm annealed ⟨100⟩ foil, x ≥ 0.95, both faces in electrolyte | Foil type matches; loading and face do not; annealed foil yields at 1 atm | ENEA-processed foil only in C3-L with pressure balance |
| SPAWAR co-deposition (C/D) | Pd from PdCl₂/LiCl onto Au/Ag while D₂ evolves; hours, no incubation | Absent | Co-deposit on the entry face of an 8–12 µm membrane, or the scintillator arm (§4, item 2) |
| Iwamura Pd/CaO transmutation (B, contested) | D₂ gas at 1 atm and 70 °C; multilayer and Cs/Sr on the **entry** face; J = 3–7×10¹⁷ D cm⁻² s⁻¹ (= 48–112 mA cm⁻² *net* permeation [calc]) | M6 puts the multilayer on the **exit** face. The electrolyte dissolves Cs/Sr (their hydroxides are soluble). Electrochemical net permeation is 10¹⁴–10¹⁶ (R2) | C3-G: gas-entry variant with an isotopically tagged ⁸⁶Sr target on the entry face, plus H3 analytics |
| Kitamura / NEDO nanocomposites (B) | 100–200 g Ni-based ZrO₂ composites at 200–300 °C, H₂/D₂ | Absent | Separate sealed ⁴He/³He arm (§4, item 6) |
| Clean Planet Ni/Cu (C) | H₂ loaded at ~250 °C, then heated to 500–900 °C during desorption | Absent; excluded by the ≤ 600 K charter | Charter decision (§3.8) |
| Lipson Pd/PdO:D (C/D) | Electrolytic loading, then desorption in vacuum facing Si/CR-39; ~3 MeV p and 11–16 MeV α [BK] | **Best match** (a C3-F desorption transient with a PdO or Au/Pd/PdO skin) | Keep. Add current-off deloading transients to the protocol. The E detector must be ≥ 1–1.5 mm thick to stop 11–16 MeV α and 14.7 MeV secondary d+³He protons |
| NTT Yamaguchi–Nishioka 1990 (D) | Au on one face; out-diffusion into vacuum | Matches if Au is on the entry side | — |
| Czerski beam-free 511 keV (D) | D-loaded Pd/Zr, no beam, underground | Matches, if the 511 pair is shielded and vetoed and ideally underground | — |

**Bottom line.** As drafted, the DFM tests a combination of α-phase, high-flux, clean-exit Pd that no B-grade claimant has used, plus two C/D-grade claims. The fixes are cheap; they only need to be *decided*.

### 3.4 Electrolysis EMI, microphonics and lock-in confounds (MAJOR)

**Capacitive injection** [calc]:
- Membrane-to-detector capacitance is C ≈ ε₀A/d = 0.13–0.27 pF for 1–3 cm² at 5–20 mm.
- That is 0.8–1.7×10⁶ e⁻ per volt, i.e. **3–6 MeV-equivalent per volt** of membrane potential swing.
- A galvanostatic step of a few volts produces a fake MeV-scale pulse at every modulation edge.
- 10 mV of bubble noise corresponds to 30–60 keV.
- If the detector's facing electrode is biased (e.g. 50 V), 1 µm of membrane vibration at 1 cm injects ~30 keV (microphonics).

**Other confounds, all correlated with current and therefore with the flux lock-in:**
- Joule heating, which drifts detector temperature and gain;
- bubble-driven vibration;
- D₂ partial pressure on the detector side, which causes discharge risk during bursts;
- gauge and RGA filaments, which crack D₂ and emit light and ions.

**Fixes:**
- Make the membrane the **hard ground** of the detector system; drive the anode.
  - Residual I·R across a 25 µm foil (4.4 mΩ/sq): ~4 mV per amp, i.e. ~25 keV, which is below threshold.
- Add a grounded screen (Ni mesh or ≤ 1 µm Al) between membrane and detector. The screen also blocks light and low-energy ions.
- Keep the detector's facing electrode at ground potential.
- Run a **dummy Si detector** behind 200 µm of Al on the same electronics, to see EMI and microphonics but no particles.
- Digitise waveforms; require ΔE–E coincidence.
- Use linear, not switching, supplies.
- Switch gauges and RGA off during counting windows.
- Add a **flux-blocked twin**: an Au-capped D membrane with identical electrochemistry and current waveform. Also run an **orthogonal flux modulation at constant current** (front D₂ pressure steps, or temperature), so any "flux-correlated" signal can be separated from "current-correlated" artifacts.
- Pre-register all of the above in the M7 protocol.

### 3.5 Vacuum side (minor unless noted)

- **D₂ gas load** over 3.1 cm² [calc] is benign for Si:
  - J = 10¹⁶ D cm⁻² s⁻¹ → 6.5×10⁻⁴ mbar L s⁻¹ → 2×10⁻⁶–10⁻⁵ mbar at 70–300 L s⁻¹;
  - J = 3×10¹⁷ → 2×10⁻² mbar L s⁻¹ → ≤ 3×10⁻⁴ mbar.
  - A static NEG-pumped mode (for ⁴He) must absorb 17–1700 mbar·L per 30 days at J = 10¹⁴–10¹⁶, which is feasible only at low flux. This is another reason to prefer the D₂-filled front.
- **D₂ exposure of Si:** no known damage to PIPS at these pressures [BK].
- **Light:** Pd is opaque. Once α/β cycling opens pinholes, cell light and D₂O vapour (~30 mbar) leak through, so use an opaque cell body; thin PTFE transmits light.
- **Laser stimulation** (Letts-type, M6) blinds the Si unless it is gated or behind ≥ 0.5 µm Al.
- **Rupture (MAJOR for equipment):** see the interlock in §3.2.

### 3.6 A new background if the membrane is thinned (minor, MAJOR for C3-L)

- Air-equilibrated electrolyte holds ~1.25 mBq of ²²²Rn per 100 mL. Progeny (²¹⁸Po, ²¹⁴Pb/Bi/Po) plate electrochemically onto the *cathode*: ~4×10⁻⁴ α s⁻¹ heading towards the detector if 30 % plate [calc]. That is **20× the M0 Si window background**.
- A ≥ 25 µm foil stops them. An 8–12 µm foil passes ²¹⁴Po alphas at 4.6–2.8 MeV, straddling the degraded-proton line at 2.3–2.6 MeV.
- **Mitigation:** Rn-free electrolyte (sparged with aged N₂/Ar, sealed), plus a ΔE–E telescope. A 10–15 µm ΔE stops a 2–3 MeV α but passes a 2.3 MeV p.

### 3.7 Controls are mismatched

The H₂O twin at equal *current* does not have equal x or J: the Pd–H plateau is 3–5× lower, and diffusivities and surface kinetics differ. For loading- or flux-dependent hypotheses:
- match the twins on **measured J and R/R₀-inferred x**, by a pre-registered rule;
- add the flux-blocked D twin (§3.4);
- keep the Pt-cathode twin only for the C1 calorimetric cell.

### 3.8 No positive control, and the charter conflicts (MAJOR)

- **No positive control.** Nothing shows *in situ* that the chain (foil → gas or vacuum gap → telescope → cuts) would register D–D protons born in this foil. ²⁴¹Am calibrates energy, not the process. Options:
  - (a) a ≤ 5 keV, ≤ 1 µA D⁺ probe, or a pulsed glow-discharge implant on the front face. R3 §9 item 1, R4 §5 item 1 and R7 §7.1 all rank this highest-EV. It gives real D–D counts from the actual foil, *in-situ* U_e versus loading and flux, and continuity with the only B-grade anomalies (screening, plateau). Beam-off data remain the cold experiment.
  - (b) a pyroelectric D–D source (Naranjo 2005) for the neutron and proton chains.
  - The charter bans beams "as the energy source". Rule explicitly that a diagnostic probe is allowed; below ~5 kV, radiation-machine registration may not be required (R6) [BK: verify locally].
- **The 600 K cap excludes Clean Planet's regime (500–900 °C).** C5 cannot be run inside the charter. The cap has no physics basis: kT = 0.1 eV versus 0.026 eV changes the Gamow factor negligibly next to the keV scale. Either raise the cap for gas-phase samples to ~1200 K while keeping "no accelerated particles above ~100 eV as the driver", or drop C5 explicitly.

---

## 4. Configurations and levers the plan has not considered

| # | Idea | What it buys against the charter's criteria | Quick plausibility [calc/BK] | Verdict |
|---|---|---|---|---|
| 1 | **⁴He/³He as co-primary** (foil extraction + static front volume) | H2 sensitivity ×10⁴–10⁷ over calorimetry; localisation by spot extraction | 10¹⁰ atoms ↔ 38 mJ of D+D→⁴He. The Pd membrane is a He barrier. Noble-gas blanks are ~10⁹ atoms [BK] | **Adopt** in every configuration |
| 2 | **Co-deposition on a scintillator** (ETC "CATHODE") | SPAWAR conditions exactly; ~2π acceptance (≈ 50 %, versus 5–30 %); no incubation; no membrane mechanics | Si is *not* a viable substrate: a cathode film on 1 µm dielectric over Si (3.5 nF cm⁻²) injects **~78 MeV-equivalent per mV**. Optical readout is immune, which is why ETC chose a scintillator. SPAWAR track densities (~10³–10⁴ cm⁻² over days to weeks [BK] → ~10⁻³–10⁻² cm⁻² s⁻¹) would reach 5σ in ≲ 1 day. Risks: Rn-progeny plating (~4×10⁻⁴ α s⁻¹, so use PSD plastic such as EJ-276 or a thin ΔE layer); poor energy resolution. ETC has published no result (R7), which is weak evidence of a null | **Adopt as a cheap parallel arm** ($5–15k) |
| 3 | **Stable deuteride films** (TiD₂, ZrD₂, ErD₂₋₃; 1–10 µm, facing Si in UHV) | No loading problem: x = 2 is retained in vacuum at room temperature. No Δp. D density 1.2–1.5× PdD. Host of the Szczecin plateau and ICCF-26 beam-free 511 claims. Provides the static, high-density reference against the Pd flux cell | Sputter plus deuteriding at ~400 °C costs ~$1–2k per batch. Loses the electrochemical knob, but temperature cycling through TiD₂ gives C6 fracto tests | **Adopt as the reference arm** |
| 4 | **Pressure-balanced gas-front DFM** (1 atm D₂) | Removes the rupture risk, clamps the exit at x ≈ 0.65, makes the front a static ⁴He volume, removes pumping | Products lose ≤ 0.3 MeV over 5 mm (§3.2); Si works at 1 atm | **Adopt as the default C3 build** |
| 5 | **8–12 µm membrane viewing the entry face** | Puts SRI/ENEA/co-deposit conditions in the telescope's view | Protons emerge at 2.3–2.55 MeV. ²¹⁴Po α leak through (§3.6), so it needs a telescope and Rn-free electrolyte. Mechanics need pressure balance | **Adopt** for C3-L |
| 6 | **Sealed PNZ/CNZ nanocomposite ⁴He arm** (200–300 °C, inside the charter) | Tests the best-replicated in-network heat claim (B) at its magnitude | 10 W for 1 week via ⁴He is ~1.6×10¹⁸ He (R2), which is unmissable. The charged-particle view of a powder bed is hopeless: one grain layer (~20 mg cm⁻²) is ~10⁻⁴ of the bed. Si cannot run at 250 °C, but CVD diamond or SiC can [BK], if a particle channel is wanted | **Adopt as the H2 arm**, with calorimetry second |
| 7 | **Parallel multi-membrane array or carousel** | P(≥ 1 active) goes from 0.20 to 0.83 for 8 samples at p = 0.2 | Pd costs ~$100 per membrane (R6); the cost is detector channels | **Adopt** (≥ 4–8) |
| 8 | **Localisation**: double-sided Si strip detectors (DSSD) plus IR microbolometer imaging of the exit face through ZnSe or CaF₂ | Makes "localised anomaly" operational. Separates edges, clamp and grid bars, and contamination spots from true sites. Correlates hot spots (the Iwamura and SPAWAR claims) with particle events | DSSD ~$5–10k plus channels; IR camera ~$5–20k [BK] | **Adopt**, low cost |
| 9 | **Low-threshold channel** (Si drift detector or keV Si) on the exit face | Tests Hagelstein's ≤ 20 keV charged-product bound and Karabut-type soft X-rays: the only direct H2-adjacent particle test | ~$10–20k [BK] | Optional |
| 10 | **Cryogenic stage** (load at RT, cool to 77 K) | Tests Bochum's disputed U_e ∝ T^−½: ×1.97 at 77 K. *If* it applied at sites with 100 eV, the rate would rise ~e²⁹ ≈ 10¹² | The Debye model is physically untenable (R3), and D is immobile at 77 K, so the flux lever is lost. P(real) < 1 % | Cheap add-on to arm 3 only |
| 11 | **Muon-catalysis-adjacent** | Nothing cold | Stopped cosmic μ⁻ ~10⁻⁵ g⁻¹ s⁻¹ [BK estimate]. In Pd they are captured by Pd (Z = 46), not d; in D₂O mostly by O. ddμ fusions ≲ 10⁻⁶ s⁻¹. A muon source is an accelerator | **Reject.** Use only as M1's background floor |
| 12 | **Acoustic-emission sensor** on the membrane | Tags cracking, the H1b confound (real *hot* D–D) | $1–3k | Adopt |

---

## 5. Recommended changes to the plan

**Decisions (lead):**
1. **Supersede ADR-001 with ADR-002** (§2.5):
   - a claim-conditioned objective;
   - sensitivity saturating at 1 % of the claimed magnitude;
   - ⁴He/³He co-primary;
   - sample multiplicity scored;
   - a hypothesis table covering H1a/b/c, H2, H3 and H4.
2. **Split C3 into declared regimes**, each with its own hardware:
   - **C3-L**: Au- or oxide-capped exit, 8–12 µm Pd viewing the entry face, SRI protocol, pressure-balanced D₂ front;
   - **C3-F**: clean exit, flux regime, Lipson and Czerski-511 tests, α/β honestly declared;
   - **C3-G**: gas-entry Iwamura variant with an isotopically tagged target on the entry face.
   - **No C3 hardware purchase before M3 delivers the exit-face chart.**
3. **Charter rulings:**
   - (a) allow a ≤ 5 keV, ≤ 1 µA diagnostic D⁺/H⁺ probe or glow discharge, with beam-off data as the cold measurement;
   - (b) raise the gas-phase temperature cap or drop C5;
   - (c) fix the H2 row to 2.6×10⁹ s⁻¹ per 10 mW;
   - (d) restate the exponent budget: 89 e-folds remain after all engineering.
4. **Add three cheap parallel arms:**
   - co-deposition on a PSD scintillator (§4, item 2);
   - a stable-deuteride reference (item 3);
   - a sealed PNZ/CNZ ⁴He arm (item 6).
5. **Run ≥ 4–8 samples per configuration**, with pre-registered per-sample denominators.

**Model briefs:**
- **M1.**
  - Scenario (ii) must include U_e = 310 eV (allowed by nulls) as well as 800 eV (excluded).
  - Add a "conditions present at the visible face" matrix (x, J, T, site classes) beside the sensitivity matrix; report the latter in U_eff as well as in enhancement factor.
  - Treat fracto-emission as a *confound* class (H1b).
- **M3.**
  - Add the detailed-balance bound of §3.1.
  - Report x_entry, x_exit and J for sticking 10⁻¹²–0.3, and for 1 atm D₂ versus vacuum on the front.
  - Include creep and ratcheting under Δp with α/β cycling, and support-grid shading.
- **M5.**
  - Add the capacitive and microphonic coupling model (§3.4), a Rn-progeny plating background (§3.6), and E-detector thickness ≥ 1–1.5 mm (to stop 16 MeV α and 14.7 MeV p).
  - Add the gas-front (1 atm D₂) geometry; include support-grid transmission.
- **M6.**
  - Drop the "shared flux" premise for segmented skins and correct the Au "null" (§3.1).
  - Specify skins on the **entry** face for Iwamura-type stacks.
- **M7.** Add the flux-blocked twin and the orthogonal constant-current flux modulation to the lock-in protocol. Requiring detector-side ground and a screen is a precondition.
- **M4.** Scope the ⁴He foil-extraction and static-volume system (blank budget, getter capacity, mass balance of bled D₂) as a first-class deliverable, not "guidance".

**M0 fixes (minor):**
- the D₂ ρ₀ shell formula;
- the WKB lower limit at 0 (or state 1 fm and apply it consistently);
- reword "300–800 eV does not apply" to "≳ 370 eV (Yukawa) excluded";
- thresholds with the matched control and systematic floor;
- ln-rate versus X-gap bookkeeping;
- add the ⁴He channel row;
- a configurable output path so running the script does not rewrite tracked figures;
- harmonise 14 d versus 30 d.

---

## 6. Findings register (ranked by decision impact)

| ID | Severity | Finding | Changes a decision? |
|---|---|---|---|
| F1 | CRITICAL | A clean exit face pins the whole membrane in α phase; high loading, high flux and a clean face are mutually exclusive (§3.1) | Yes: C3 split into L/F/G; M3 gating |
| F2 | CRITICAL | The DFM reproduces no B-grade claimed condition at the detector-visible face; its null would be dismissed (§3.3) | Yes: thin membrane, Au cap, co-deposition, gas-entry variants |
| F3 | MAJOR | The ADR-001 objective overweights detector sensitivity (worth 2.5–7.5 % of probability) and ignores conditions and multiplicity (O(1)) (§2.1, §2.4) | Yes: ADR-002 and reweighted matrix |
| F4 | MAJOR | ⁴He missing as H2's primary detector; the DFM is well suited to it (§2.2) | Yes: co-primary, foil extraction |
| F5 | MAJOR | The H1/H2 split omits e⁺e⁻, H1b (fracto confound), H3 and H4 (§2.3) | Yes: added observables (AE sensor, SIMS/ICP-MS with tagged targets, ³He) |
| F6 | MAJOR | Annealed membrane yields at 1 atm (≈ 250 MPa against 35–70 MPa); α/β ratcheting leads to rupture (§3.2) | Yes: pressure-balanced D₂ front or grid; interlock |
| F7 | MAJOR | 3–6 MeV per volt capacitive injection; modulation edges and bubble noise are in phase with the lock-in (§3.4) | Yes: grounding, screen, dummy detector, flux-blocked twin, orthogonal modulation |
| F8 | MAJOR | No in-situ positive control; the charter bars the probe beam recommended by R3, R4 and R7 (§3.8) | Yes: charter ruling |
| F9 | MAJOR | The 600 K cap makes C5 impossible to run (§3.8) | Yes: charter ruling |
| F10 | MAJOR | Single-sample design ignores per-sample success rates (§2.1) | Yes: ≥ 4–8 samples |
| F11 | MAJOR | M6 segmented-skin premise false; Au "null" inverted; Iwamura stack on the wrong face (§3.1, §3.3) | Yes: brief rewrite |
| F12 | minor→MAJOR (thin foil only) | Rn-progeny α penetrate 8–12 µm foils into the proton window (§3.6) | Only if C3-L is adopted |
| F13 | minor | H₂O twin mismatched in x and J (§3.7) | Protocol |
| F14 | minor | M0: D₂ ρ₀ concept; 1 fm cutoff; "300–800 eV" wording; thresholds optimistic ×1.2–5; exponent bookkeeping; H2 row typo; figures overwritten (§1) | No |
| F15 | minor | D₂ gas load, Si exposure to D₂, light: benign with standard practice (§3.5) | No |

---

## Appendix: reproduction notes

- **D₂ shell density:** ρ₀ = [√π x₀ · 4πR_e²]⁻¹, with x₀ = ħc/√(μc²·ħω), ħω = 0.371 eV, R_e = 0.7414 Å.
- **Yukawa WKB:** 2G = 2∫₀^{r_tp} √(2μ(V − E))/ħ dr, with V = (e²/r)e^{−rU_e/e²}. Substitute r = r_tp s² to remove the endpoint singularities.
  - Check: with V = e²/r − U_e, 2G = √(E_G/U_e) exactly (99.286 at 100 eV).
- **Yukawa free-gas null boundary:** ½ n (S c√(2/μc²)) (2/√(πkT)) e^{−2G(U_e, kT)} = 10⁻²⁵ gives U_e ≈ 367 eV.
- **Tail share (population):** ∫_{U > μ+3.09σ} e^{−√(E_G/U)} φ(U) dU / ∫ e^{−√(E_G/U)} φ(U) dU, evaluated in log space on a fine grid.
- **Membrane:** σ_max ≈ 0.423 (E p² a²/t²)^{1/3}; w₀ ≈ 0.662 a (p a/(E t))^{1/3}.
- **Exit face:** J = s · 2Z₁ · P_eq[atm], with Z₁ = P/√(2π m kT) = 7.7×10²³ cm⁻² s⁻¹ at 1 atm D₂; α-phase x = 0.015 √(P_eq/0.04 atm).
- **Diffusion capacity:** D n_Pd/L per unit Δx.
- **Capacitive injection:** Q = ε₀(A/d) ΔV; energy-equivalent = (Q/e) × 3.62 eV.
- **⁴He:** 1 J = 2.62×10¹¹ ⁴He at 23.85 MeV.
- **Residual energies** use power-law range fits to R5's ATIMA ranges in Pd: p, R ∝ E^1.75 with R(3.02) = 31 µm; α, R = 1.06 E^1.33 µm.
