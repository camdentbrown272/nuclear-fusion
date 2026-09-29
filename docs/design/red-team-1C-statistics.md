# Red team 1C: statistics, inference and credibility of iteration-1 (DFM-4)

**Scope.** `docs/design/iteration-1.md` §§1, 4, 9, 10, 13 against ADR-003/004, M5 (`sim/m5_stats.py`), M8, M4, R1/R5/R7 and red-team-0.
**Stance:** a skeptical mainstream nuclear physicist and a trial methodologist. **Date:** 2026-09-29.
**Code:** `sim/rt1c_stats.py` (deterministic; no random numbers) → `docs/design/figs/rt1c_stats.txt`, `rt1c_reach.png`, `rt1c_posterior.png`. Section references §S1–§S8 below point to that output.
**Tags:** [calc] computed here; [BK] reviewer background knowledge, not re-verified (no literature access in this session); [OoM] order-of-magnitude estimate, uncertain ×3–10.

---

## Summary verdict

The design is not pre-registered in any sense a referee would accept. "T1–T4, look-elsewhere-corrected across 4 tests" is a list of intentions:
- T2's detector does not exist yet (M7 pending);
- T3 has no test statistic;
- "Blinded" is one word;
- nothing defines the stopping rule, data-quality cuts, analysis population or failure handling.

The text invites about 1000 analyses (§S1). With a realistic effective trials factor of 100–300, a local 5σ is a **global 3.8–4.0σ**.

The claim rule works against the sample count. The rule requires ≥ 2 membranes, but the design sells P(≥ 1 of 4 active) = 0.59. Under the rule, the probability of a claim even if the effect is real is **0.18** (0.10 per skin). Lot correlation lowers P(≥ 1) to 0.29–0.44.

The primary contrast, Σ A − T-H, is the weakest and least-controlled one available:
- T-H has **one quadrant per skin**, so it needs 44–62 signal counts instead of 6–8;
- D and H cells differ in current, gas evolution, EMI and heating;
- D and H cells also differ in **D-specific cosmic-neutron backgrounds**, which the H twin cannot cancel.

The §10 statement that every channel reaches ≤ 1 % of every claim is false for Tohoku-static (1.7 %) and for heat (up to 20–100 %).

A Bayesian reading is worse still. **An all-null outcome has a Bayes factor of only 1.36 against the effect**, so the claim that "a null is informative for the first time" is not yet earned. A single-membrane 5σ with a realistic 1–5 % per-campaign artifact rate moves a 1 % prior only to 5–21 %.

All of this is fixable on paper, before hardware, at modest cost:
- a frozen, hashed analysis plan with gatekeeping;
- a within-cell flux-on/off primary contrast;
- ≥ 8 membranes over ≥ 2 Pd lots with ≤ 3 skins;
- hardware and software salting, a blind custodian and an independent second analysis.

---

## Ranked findings

| ID | Sev. | Finding | Quantitative evidence | Concrete fix |
|---|---|---|---|---|
| **S-1** | **Critical** | Pre-registration is incomplete. T1 is defined only by an energy window. T2 depends on M7 (pending). T3 and T4 have no statistic, background model or threshold. The claim rule gives no significance for the per-membrane replication. "Look-elsewhere across 4 tests" ignores that T1 alone is ≥ 5 tests (one per skin) | The text invites ≈ 985 analyses (§S1): 5 skins × 5 phases × 2 windows × 3 contrasts, per-quadrant variants, t/³He lines, a lock-in period × lag scan, He windows and melts. At N_eff = 100–300, local 5σ → **global 4.0–3.8σ**; global 5σ needs **local 5.8–6.0σ** | Freeze the skeleton in §A: two primary tests with α-weights and α-propagation, and everything else confirmatory or descriptive. Hash-timestamp it (OSF or a signed git tag) before P1 |
| **S-2** | **Critical** | The claim rule (≥ 2 membranes) is inconsistent with N = 4 at p = 0.2. The design quotes P(≥ 1) = 0.59 | Given that the effect is real: P(≥ 2 of 4) = **0.18**; per skin (3 slots) P(≥ 2 of 3) = 0.10; with lot correlation ρ = 0.3 it is 0.22 but P(≥ 1) falls to 0.44 (§S6). ADR-003's own sample term is overstated by ≈ 3× | ≥ 8 active membranes from ≥ 2 Pd lots, with ≤ 3 skins (Σ_s P(≥ 2) = 0.60 for 5 skins vs 0.87 for 3 and 0.99 for 2; §S6). Drop ED, and drop M if P0b shows x_in < 0.95 (iteration-1 risk 1). Define the membrane as the unit of replication |
| **S-3** | **Critical** | T1's control is under-exposed and mismatched. Σ(A) − T-H per skin sets t_off/t_on = 1/n_q, which contradicts M5 rec. 12 (background exposure ≥ 10× live time) | Signal counts for a median 5σ: **44–62 against T-H**, versus 6.5–11 against a pooled background model (§S2). Per-skin reach worsens from 6×10⁻⁵–2×10⁻⁴ to **4×10⁻⁴–1×10⁻³ fusions s⁻¹**, above 1 % of Lipson's low end (10⁻⁴) | Primary T1 background = a pooled model: blank telescope B, all 16 quadrants of T-Pt and T-H, P0a pre-run, and within-cell D-static (no-flux) blocks, with nuisance scales. T-H becomes an artifact discriminator, not the subtraction |
| **S-4** | **Major** | The D-vs-H contrast is not signed-conservative. Cosmic neutrons on deuterium in the membrane, the 100 µm electrolyte film and the front D₂ produce d(n,np) breakup protons and n–d recoil deuterons (PID leakage). The H₂O twin lacks both. M5 counted only membrane recoils (≤ 0.003/day) | D atoms within product range ≈ 2.4×10²¹, of which 88 % are in the electrolyte film. Breakup protons in the T1 window: **5×10⁻⁴–1.3×10⁻²/day/cell**. Recoil-d leakage 5×10⁻⁵–1.2×10⁻²/day (§S4) [OoM]. Compare the whole modelled background of 0.01–0.03/day | Geant4 or MCNP for all D inventories. Primary contrast = D flux-on vs D static-loaded in the same cell, which cancels any flux-independent D term. Pre-register a D-static block (L-quadrant-like, front pressure raised to null J) of ≥ 25 % of P2 |
| **S-5** | **Major** | Low-count "5σ" is fragile to the background model and to bursts. n_crit is 4–7 counts | If the true background is ×3 the assumed value (M5's own uncertainty), the n_crit event count is only **3.9–4.2σ**; at ×10 it is 1.8–2.9σ (§S3). One unvetoed microdischarge or pile-up cluster supplies n_crit by itself. Tightening local 5σ to 6σ costs almost nothing in counts (5.6 → 5.6 at B ≈ 0.04; §S2), because the credibility risk sits in B, not in α | Profile B over its ±×3 band (use 3B in n_crit). Require **≥ 10 signal events**. Count 10-s clusters as one event and run a separate pre-registered burst test. Require spectral/depth consistency (§A.4) |
| **S-6** | **Major** | §10 "each channel reaches ≤ 1 % of the claimed magnitude" is false | Tohoku-static: 2.2/130 = **1.7 %**. Heat: 4–20 mW vs 0.1–1 W is 0.4–20 %; volume-scaled from the SRI wire (0.024 cm³) to the membrane (0.0047 cm³) it is **2–100 %**. Si per skin in DFM geometry: 6×10⁻⁵–1×10⁻³, up to 10 % of Lipson's low end (§S8). M5's 3–8×10⁻⁵ is for the C3 geometry (whole foil, vacuum, ε = 0.22, 30 d, 100 % live), not the DFM (ε ≈ 0.15 per quadrant, 0.08 for entry-face sources, 21-d P2, ~0.8 live) | Normalise every claim to area or volume and state the scaling. Restate §10 per quadrant and per skin with live time. Mark heat and Tohoku as "partially tested". Do not use "saturated" for them in the ADR-003 score |
| **S-7** | **Major** | The Bayesian value of the campaign is overstated | P(pos \| real) = P_cond × P(≥ k active) × P_det = 0.27 for a single membrane, **0.08 under the claim rule** (P_cond 0.5, P_det 0.9). Bayes factor of an all-null = **1.36** (1.60 at N = 8; 3.1 only if P_cond = 0.9 and N = 8). At a prior of 10⁻², the posterior after a claim-rule positive is 0.45 if the artifact rate is 10⁻³ and 0.89 if it is 10⁻⁴ (§S7, `rt1c_posterior.png`) | State the null's information honestly: an in-regime *upper limit on magnitude*, not a test of existence. Raise P_cond, which is the dominant lever, via P0b/x-verification gates; raise N; push the artifact rate to ≤ 10⁻⁴ with S-8 |
| **S-8** | **Major** | Credibility machinery is missing: no custodian, no salting, no independent analysis, no open data, no hidden He spikes | A skeptic's prior on "artifact" for any single-lab LENR positive is ≥ 10⁻² per campaign [BK: 1989–2019 history]. S-7 shows that a posterior > 0.5 needs ≤ 10⁻⁴ | §B: custodian-held cell map; software salting (0–2 channels, rate drawn from {0, 0.5, 1, 2}× reach); hardware salting via the ADR-004 source boost on a hidden schedule; coded He aliquots with hidden ⁴He spikes; two analysis teams; release of list-mode data and code at unblinding |
| **S-9** | **Major** | The T-H match "on measured x and J" is a garden of forking paths. It takes two knobs with isotope-dependent calibrations: the R/R₀–x curve differs for H and D, and J metering differs. Which x (entry, mean or exit)? Which J? What tolerance? Is the operator in the loop? Matching x and J forces a different current, so current-coupled artifacts (EMI, Joule heating of Si, bubbles) differ between D and H | Current at matched x differs by the isotope factor in the Tafel/isotherm relations (≈ 1.5–3× [BK]). Any current-correlated artifact survives Σ A − T-H | A deterministic controller with targets frozen from P0b (no operator adjustment after P1). Run T-H in two pre-registered blocks, current-matched and x/J-matched. Use T-Pt (same current waveform) as the EMI and heat control in T1 |
| **S-10** | **Major** | ⁴He (T3) has D-specific artifacts that A − T-H does not cancel: a residual D₂⁺ tail at m/z 4 exists only in D cells; each membrane carries its own stock He; tritium in D₂O decays to ³He | D₂O lots carry 10¹–10³ dpm/mL of T (R5 §5). In 12 mL that is 2–200 Bq, i.e. **2–200 ³He atoms/s**, against a ⁴He floor of 43 atoms/s. This biases ³He isotope dilution and fakes the ³He line (H1a/H4). σ_B = 7×10⁵/day is assumed, not measured, and between-cell variance is not modelled | T3 statistic = each cell's excess over **its own** ≥ 5-blank series; D₂-matrix spike for every aliquot class; LSC-assay every D₂O lot and correct ³He for decay; measure between-cell σ in P0a with all cells blank |
| **S-11** | **Major** | Stopping, analysis population and failure handling are undefined. "P2 ≥ 3 wk" is open-ended, which is optional stopping | Unplanned peeking with extension until significance can inflate the type-I error many-fold [BK: Armitage et al. 1969] | Fixed live-time targets (§A.6); interim looks for safety only; if extension is ever allowed, O'Brien–Fleming spending pre-registered. ITT and per-protocol sets; no replacement membranes after P2 starts |
| S-12 | Minor | The "ABBA midpoint swap" is one swap, i.e. AB/BA, not ABBA. Swapping telescopes vents the front volumes, which breaks T3 static windows, admits air He and loads the membrane ΔP | One swap confounds detector with period | Either do a full ABBA at He-window boundaries with post-reseal blanks, or replace swaps with per-telescope characterisation in P0a and P6 (≥ 10× live time each) |
| S-13 | Minor | "Absent in the matched control" has little power at T-H's exposure | A common artifact giving s = 10 A counts over 4 quadrants predicts 2.5 in T-H; P(0) = 0.08, **LR ≈ 12** (§S5) | Quantify: the 95 % UL for T-Pt and pooled T-H per quadrant-exposure must be < 25 % of the A rate |
| S-14 | Minor | Flux lock-in with free period and lag is a hidden scan (≈ 320 trials in §S1) | — | Freeze period (240 s) and a single lag model (τ from P0b permeation transients) |
| S-15 | Minor | T4 reach of 10¹⁰ cm⁻² ⁹⁴Mo is undemonstrated. The DFM grid is **Mo** (if reused in C3-G), and ⁹⁴Zr is a SIMS/ICP isobar | Natural-Mo contamination at 10¹³ cm⁻² already carries 9×10¹¹ ⁹⁴Mo, so seeing 10¹⁰ needs ~1 % ratio precision | No Mo in C3-G; quote the reach as an isotope-ratio precision; add a Zr survey; use coded samples |
| S-16 | Minor | T2 cannot be pre-registered until M7 exists | — | Remove T2 from the primary family unless M7 is frozen before P1 (§A) |

---

## A. Proposed frozen pre-registration skeleton

**A.0 Registry.**
- Freeze this document, the analysis code (run on simulated and P0a data) and the DQ cuts as a signed git tag and an OSF registration **before P1**.
- Amendments are allowed only before unblinding, dated and justified. Unblinding is irreversible.

**A.1 Units and sets.**
- The *membrane* is the unit of replication. Quadrants are sub-units: they share the entry face, electrolyte, current and lot.
- **ITT set:** all membranes that complete P1.
- **Per-protocol set:** ≥ 14 d P2 live time at mean x ≥ 0.90 by R/R₀, calibration frozen in P0b.
- Both sets are reported; the primary analysis uses ITT.

**A.2 Primary family.** Weighted Bonferroni with α-propagation (graphical procedure). Total one-sided α_G = 2.87×10⁻⁷.

| Test | Weight | Data | Statistic | Background model | Local threshold |
|---|---|---|---|---|---|
| **T1 (H1a)** | 0.5 | P2 current-on, **all 16 A quadrants pooled** (one test, not five), PID window 1.41–3.10 MeV, clusters of 10 s counted once | Poisson profile-likelihood ratio, s ≥ 0 | Pooled: B telescope + T-Pt (4 quads) + T-H (4 quads) + P0a pre-run + **A-cell D-static blocks**. Each gets a scale nuisance; global B nuisance log-normal with σ = ln 3; exposure ≥ 10× live time | p ≤ 1.4×10⁻⁷ (**5.13σ**); ≥ 10 net events |
| **T3 (H2)** | 0.5 | Front and headspace ⁴He over all P2–P4 windows, summed per cell and then over A-cells; quadrant melts as a separate pre-specified component | Σ_cells (N_obs − N̂_blank,cell) / σ_tot, where σ_tot has the per-cell blank-series variance plus between-cell variance from P0a; D₂-matrix-spike correction | Per-cell blank series (≥ 3×1 d, ≥ 2×7 d, pre and post) and sibling coupons for the melts | p ≤ 1.4×10⁻⁷ |

- If T1 rejects, its α passes to T3, and vice versa.
- T2 (511–511) enters only if M7 is frozen before P1, with weight taken 0.1 from each of T1 and T3.
- T4 (C3-G) is a separate registered study with its own α. It is not in this family.

**A.3 Confirmatory conditions for a discovery claim.** All must hold; each is pre-registered and one-sided.
1. The primary test rejects at its threshold.
2. **Flux dependence on independent data:** P3 lock-in, one frozen period and lag, in the same channel ≥ 3σ, **and** D-static blocks consistent with background (95 % UL < 25 % of the P2 rate).
3. **Replication:** ≥ 2 membranes, each ≥ 2σ locally, from ≥ 2 lots if available. Heterogeneity χ² across membranes not rejected at p < 0.01. No single cluster or day contributes > 20 % of events.
4. **Signature:**
   - The ΔE–E energy distribution is consistent with protons from 0–15 µm (χ²/KS p > 0.01), using the ⁶LiF-calibrated response.
   - The distribution is inconsistent with the sideband shape at p < 0.01.
   - Rates of t and ³He are consistent with the p rate where kinematically accessible.
5. **Controls:** T-Pt, T-H (4 quads pooled) and B each have a 95 % UL per quadrant-exposure < 25 % of the A rate.
6. **Blind integrity:** salts recovered within ±2σ of their injected rate before the real box is opened. The two analysis teams agree within 1σ.
7. **Neutron consistency statement:** the pre-registered sentence from M5 ("protons below ~10⁻² fusions s⁻¹ appear without detectable neutrons") is evaluated, not used as a veto.

**A.4 Reporting tiers.**
- Global ≥ 3σ on a primary test is "evidence".
- A "claim" requires all of A.3.
- Everything else, including per-skin and per-quadrant results, secondary lines, neutron, heat, AE and 511, is **descriptive**, with FC intervals and no significance language.

**A.5 Null reporting (committed).**
- For every claim c, publish Feldman–Cousins 90 %/95 % upper limits normalised per area and per volume, with the scaling stated, whatever the outcome.
- Submit as a registered report where the journal allows it. Release list-mode data, waveforms and code at unblinding.

**A.6 Stopping.**
- Fixed live-time targets: P2 = 21 d, P3 = 14 d, P4 = 7 d.
- Interim unblinded looks are for safety only (neutron bank at 10× background, pressure).
- Hardware failure:
  - a membrane failing before P2 is replaced under a new ID;
  - after P2 starts, its data are truncated at the last DQ-good period, and no replacement enters the primary analysis.
- No extension conditioned on data. Any extension is decided on DQ grounds only, by the custodian, blind to box counts.

**A.7 Frozen DQ cuts.** Tuned on P0a, blank telescope and sideband data only:
- veto ±1 s around current edges;
- leakage current and ΔE-noise-rms thresholds;
- muon-veto coincidence;
- pulse-shape χ² against pulser templates;
- cluster flag (≥ 2 events in 10 s in one detector);
- He aliquot acceptance (spike recovery within ±10 %).

---

## B. Blinding and salting scheme (fills "Blinded")

1. **Custodian.** One person outside the analysis team owns:
   - the cell → channel map (T-H and T-Pt among A labels);
   - the salt file;
   - the He aliquot codes.
2. **Box.** PID-window events of all telescope channels are hidden. Analysts see sidebands, B and marker windows (¹⁰B α/⁷Li, ⁶LiF t/α), and aggregate live time.
3. **Software salt.** Synthetic ΔE–E events are drawn from the ⁶LiF-validated response for a depth-uniform source. They are injected into 0–2 randomly chosen channels at a rate drawn from {0, 0.5, 1, 2}× the median 5σ reach, with Poissonian timing and a random fraction during current-on. Salts are removed only after the analysis is frozen and salt recovery is reported.
4. **Hardware salt.** The ADR-004 source boost (AmBe/Cf at a partner facility, or cosmic ¹⁰B markers) is toggled on a custodian schedule. It tests the lock-in pipeline end-to-end without a beam.
5. **He.** Aliquots go to the external sector-MS lab under codes. ≥ 20 % are blanks, and ≥ 2 hidden ⁴He spikes of 5–10 σ_B are included.
6. **Two teams.** Independent code bases share only the pre-registration and the DQ list. At least one external member is ideally a known skeptic.

---

## C. Details

### C.1 Trials factor (§S1)
- The M5 "N = 1" and the iteration-1 "4 tests" count intentions, not the analyses the text invites.
- The largest hidden contributors are:
  - per-skin × per-phase × per-contrast T1 variants (150);
  - per-quadrant results needed by the claim rule (160);
  - the lock-in period/lag scan (320).
- The gatekeeping above collapses the discovery family to two tests (local 5.13σ). Everything else is either confirmatory with fixed parameters or descriptive.
- At the actual background levels (≤ 0.3 counts), moving from 5σ to 6σ local costs < 1 signal count. The trials penalty is cheap in counts and expensive in credibility, so there is no reason not to pay it.

### C.2 Sensitivity recomputation (§S2, §S8, `rt1c_reach.png`)
- **Inputs:**
  - ε per quadrant 0.15 (iteration-1 risk 4), or 0.08 for entry-face sources (2.05 MeV p, lower PID acceptance; M5 `U(0-20um)` scaling);
  - background per quadrant = cell B/4 with B = 0.01, 0.03 or 0.09/day;
  - live fraction 0.8 (edge vetoes during the P3 square wave; keep-alive gaps; DAQ);
  - P2 = 21 d.
- **Pooled background model:**
  - per quadrant 1.5–6×10⁻⁵ fusions s⁻¹;
  - per skin (sum over 3–4 quadrants) 6×10⁻⁵–2×10⁻⁴.
- **T-H as control:** per skin 4×10⁻⁴–1×10⁻³.
- **Dead time** is negligible at these rates. The losses are vetoes and downtime (≈ 20 %).
- **Background non-stationarity** (barometric ±0.7 %/hPa on the cosmic terms; solar modulation of a few %) is irrelevant at ≤ 1 count. What matters is the ×3 model uncertainty (§S3) and non-Poisson bursts.
- **He floors.** 0.16 nW corresponds to 43 He s⁻¹, consistent with M8. The floor is systematics-limited (N_min ∝ σ_B·t, so it does not improve with time), and σ_B is an assumed stock-He mobilisation rate. The claim ratio stays ≤ 10⁻⁶ even at 10× σ_B, so the He reach conclusion is robust. The T3 *test*, however, needs the between-cell variance (S-10).
- **Heat.** The 4–20 mW 5σ figure is M4's. Against volume-scaled SRI claims it does not reach 1 %.

### C.3 Matched controls and D-vs-H confounds
| Confound | Direction in A − T-H | Cancelled by |
|---|---|---|
| Current, EMI, Joule heating of Si (current differs at matched x/J) | ± either | T-Pt (same waveform); current-matched T-H block |
| Gas evolution and bubble microphonics (D₂O vs H₂O bubble size and overpotential) | ± | T-Pt; waveform χ² cut |
| d(n,np) protons and n–d recoil d (cosmic) | **+ (fakes D signal)** | D-static blocks in the same cell |
| n–p recoils in H cell (0.006/day) | − (masks) | model it |
| D₂⁺ tail at m/z 4 (He) | **+** | D₂-matrix spike per aliquot class |
| Tritium in D₂O → ³He; electrolytic T enrichment | **+** for the ³He line | LSC assay per lot; decay correction |
| Stock He differing per membrane | ± | per-cell blanks; sibling coupons |
| D₂O vs H₂O heat capacity; thermoneutral 1.527 vs 1.481 V | ± heat | already pre-registered (§6); keep |

The within-cell flux-on vs D-static contrast cancels every row except current/EMI, and T-Pt covers that one. **Make it the primary T1 contrast**; the D-vs-H comparison then tests isotope specificity (A.3.5).

### C.4 Sample count and allocation (§S6)
- **p = 0.2 is not defensible as a per-sample probability in an independent lab.**
  - It is a claimant-conditional rate: ENEA says "a subset", Storms "a minority", and NRL Pd–B gave 7/8 in one group (R1).
  - It is selection-biased: denominators are rarely published (R1 artifact #15).
  - Independent programmes (NHE 1992–97, Google 2019) observed ~0 of many, though Google did not reach the regime.
  - A defensible prior is p ~ 0.05–0.1 conditional on reaching the regime, with a large lot effect (ENEA attributes reproducibility to metallurgy).
- **Correlation matters.** With ρ = 0.3–0.7, P(≥ 1 of 4) at mean 0.2 falls from 0.59 to 0.44–0.29, and adding membranes from the same lot barely helps.
- **Decision.**
  - Maximise the number of membranes and lots before the number of skins.
  - Suggested allocation: **8 active membranes = 2 lots × 4 (or 4 lots × 2), 3 skins (L, F, X)** × 2–3 quadrants each, with the 4th quadrant as a within-membrane repeat of the best-motivated skin.
  - That gives P(≥ 2 of 8) = 0.50 at p = 0.2, and 0.19 at p = 0.1.
  - The marginal cost is ≈ 4 × (telescope $6k + cell hardware ~$5k + skins ~$2k) ≈ $50k, which is cheaper than any other lever on ADR-003's sample term.
- **Staging alternative.** Run the array twice: 4 + 4 membranes, with the second array's allocation frozen before the first is unblinded.

### C.5 Posterior after possible outcomes (§S7)
**All null.**
- Bayes factor 1.36. At a prior of 1 %, the posterior is 0.74 %.
- The informative output is the magnitude limit per claim, conditional on the regime being reached. Report it that way, not as a refutation.

**One membrane positive (5σ local).**
- Posterior 5–21 % from a 1 % prior (artifact rate 5 % / 1 %).
- A skeptic will say: "one cell, one lab, a handful of counts, an unmodelled background; the same as 1989–2019." Correct response: iteration-2 blind replication (§13 already says so), with no claim language.

**Claim-rule positive with S-8 machinery.**
- Posterior 0.45–0.89 from a 1 % prior (artifact rate 10⁻³ / 10⁻⁴).
- It remains single-site. A skeptic will still ask for:
  - an independent site with its own detectors;
  - a dose–response in J (P4 gives only 2 points; pre-register a ≥ 3-level J ladder);
  - absence of the signal with the membrane replaced by an identically processed but unloaded Pd dummy in the same cell;
  - a check that the energy-depth signature is not reproduced by the ⁶LiF marker leaking through a pinhole.

### C.6 What a skeptic will still say after all fixes
1. **"Your 5σ is a Poisson 5σ."** Only the pooled background model with its nuisance terms and the ×3 band answers this, together with the ≥ 10-event floor.
2. **"The electrolyte-film deuterium is a target, not a bystander."** Answered by S-4's simulation and the D-static blocks.
3. **"You tuned the matching."** Answered by the frozen controller and the two-block T-H.
4. **"Where is the heat?"** The pre-registered statement that particle products below the calorimetric floor carry no heat claim must be in the abstract.
5. **"Single site."** Unavoidable in iteration-1. Say so.

---

## Changes requested to iteration-1 (summary)
1. **§9 pre-registration table.** Replace it with §A, and remove "look-elsewhere across 4 tests".
2. **§4.1.**
   - N ≥ 8 membranes, ≥ 2 lots, ≤ 3 skins.
   - State P(≥ 2 active) alongside P(≥ 1).
   - Replace the midpoint swap (S-12).
3. **§4.1 and §9 T1.** Make the primary contrast within-cell flux-on vs D-static, with a pooled background model. T-H becomes an isotope-specificity check; T-Pt becomes a mandatory T1 control.
4. **§10.** Restate reach per quadrant and per skin with live time and ε = 0.08–0.15. Flag Tohoku and heat as > 1 %, and state the claim normalisation.
5. **§5 and §9 T3.** Add the per-cell blank-series statistic, between-cell σ from P0a, a D₂O tritium assay and D₂-matrix spikes.
6. **Protocol.** Add a custodian, salting (software and hardware), hidden He spikes, two analysis teams, data release and fixed durations (§A.6, §B).
7. **M5 follow-up.** Add the D-inventory cosmic background (breakup and recoil-d) for the membrane, electrolyte film and front gas.
