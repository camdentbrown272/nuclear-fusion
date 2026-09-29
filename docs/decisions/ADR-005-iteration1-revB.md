# ADR-005 — Iteration-1 revision B (response to red-team-1A/1C/1D)

**Status:** accepted · **Date:** 2026-09-29 · **Supersedes:** iteration-1 rev A (commit 77887a6) where they conflict
**Inputs:** `docs/design/red-team-1A-physics.md`, `red-team-1C-statistics.md`, `red-team-1D-claim-fidelity.md`.
Red-team-1B (engineering/safety) was cut off by a usage limit before it pushed any output. It is re-run on revision B (§6).

## 1. What the reviews established (accepted as fact)

1. **Quadrant skins with different exit laws cannot coexist on one membrane.**
   - Under any drained exit (bare, oxide or ED at 0.5 bar D₂), the whole 15 µm column sits at x ≈ 0.655, entry face included. Holding x_in = 0.9 against that exit needs ≈ 3 A cm⁻² absorbed, while the cell supplies ≤ 0.15 A cm⁻² (1A-P2, 1D-D1; `m3_common`).
   - Neighbouring L and F quadrants therefore differ by Δx ≈ 0.25. That is a 1.6 % misfit, or GPa-level stress against a 40–150 MPa yield, so the foil will crack along the skin boundaries. This is a built-in H1b confound (1A-P3).
   - A rim van der Pauw gauge on a mixed-regime membrane reports nothing useful (1A-P2, 1D-D2).
2. **The PCTFE cell body and FFKM rim seal swamp both ⁴He channels.** Atmospheric He permeating the polymer adds 10⁸–10⁹ He s⁻¹, against a floor of 43–100 He s⁻¹, and the rate rises with temperature (1A-P1). The rev-A design violated M4 §5.10 and M8.
3. **At 15 µm the ⁶LiF depth calibration fails and entry-face efficiency halves.** Only 3 % of calibration tritons emerge, and entry-face ε_PID = 0.116 over the real acceptance (1A-P4, P6).
4. **The pre-registration was a list of intentions** (1C-S1). The rev-A text invites ~10³ analyses, so a local 5σ is a global 3.8–4.0σ.
5. **The claim rule (≥ 2 membranes) contradicts N = 4.** With N = 4, P(claim | effect real) = 0.18, and one Pd lot means N_eff ≈ 1 (1C-S2, 1D-D3).
6. **The D − H primary contrast is biased.**
   - Electrolyte-borne cosmic recoils are D-specific in D₂O: n–d recoils plus D(n,np) breakup give 0.014 d⁻¹.
   - They are H-specific in H₂O: n–p recoils give 0.036 d⁻¹.
   - Both exceed M5's whole background budget, and neither is seen by the blank telescope (1A-P5, 1C-S4).
7. **Claim-register fidelity of the rev-A array is low.**
   - Under ADR-003 the DFM array scores ≈ 0.42 of the programme's ≈ 2.3, at ~75 % of cost.
   - C3-G with N = 3 (1.10) and C1 ×4 (0.71) each outscore it (1D §4).
   - The DFM "Iwamura" (X) and "Lipson" (F) quadrants are not those experiments (1D-D9, D13).

## 2. Decisions

### 2.1 Membranes: single-regime, grouped by loading state (resolves 1.1)

| Membrane type | Exit face | Loading state (whole foil) | Tests |
|---|---|---|---|
| **H-L** ("closed cathode") | Au ≥ 50 nm, pinholes ≤ 3×10⁴ cm⁻² (H-permeation witness + Cu decoration per sputter batch) | x ≈ x_in ≥ 0.9, J ≈ 0; flux only as \|dx/dt\| from cathodic current steps | SRI/ENEA at a detector-visible location (the closest DFM analogue of a closed SRI rod); NTT-like out-diffusion in P6a |
| **H-M** ("loaded + flux") | Ni ≥ 20 nm (thickness from P0b) | x ≈ 0.91–0.92 with J ~ 10 mA cm⁻² | McKubre (x − x₀)²·flux. **Pre-registered fallback: if P0b shows no loaded-flux window, H-M membranes are built as H-L** |
| **FX** ("flux") | two diagonal quadrant pairs: F (bare) and X (Pd/CaO multilayer) | x ≈ 0.65 everywhere (F and X share the drained state, so there is **no misfit**) | Flux-regime H1 (Lipson/Czerski-type conditions), oxide-interface H1 site test, flux lock-in, P6a desorption into vacuum |

- **ED is dropped.** It has no distinct regime (1A-P12, 1C-S2). X is kept only as a within-membrane comparison with F on FX membranes. **Iwamura claims are removed from the DFM** and tested only by C3-G.
- **Skin count ≤ 3 families (L, M, FX).** The per-skin analyses are descriptive (§2.6).
- **The quadrant ΔE–E telescope and cross septum are retained.** On single-regime membranes the four quadrants give:
  - localisation;
  - within-membrane heterogeneity and burst checks (claim condition A.3.3);
  - an F-vs-X comparison on FX membranes.

### 2.2 Thickness 12 ± 1 µm unloaded, measured per foil to ± 0.3 µm (resolves 1.3)

Computed with `m5_si_telescope` over the real acceptance, 4 mm of 0.5 bar D₂:

| Unloaded thickness (µm) | 8 | 10 | **12** | 15 |
|---|---|---|---|---|
| Entry-face p ε_PID | 0.200 | 0.186 | **0.162** | 0.116 |
| Entry-face p median at ΔE (MeV) | 2.25 | 2.03 | **1.80** | 1.41 |
| ⁶LiF tritons emerging > 0.1 MeV | 72 % | 54 % | **34 %** | 3 % |
| ²¹⁴Po α at normal incidence (MeV) | 4.75 | 3.85 | **2.82** | 0.84 |

- **Why 12 µm.** It is inside M1's 15 ± 3 µm band. It gives a usable ⁶LiF calibration and +40 % entry-face efficiency relative to 15 µm.
- **Why not thinner.** Below 12 µm the gain in efficiency is small and the mechanical margin shrinks.
- **Radon progeny.** ²¹⁴Po α still stop in the 25 µm ΔE (< 5.17 MeV), so PID rejects them (M5). They are not a leakage path.

### 2.3 All-metal, He-tight cell and bonded membrane (resolves 1.2)

- **Cell body.** 316L, PTFE-lined on all electrolyte-wetted surfaces. The liner has no path to air, is baked at 150 °C and D₂-purged before fill, and ≥ 2 blank He windows must pass after the build. PCTFE is removed.
- **Membrane mount.**
  - The membrane is **diffusion-bonded to a corrugated Pt annulus** (0.10 mm, one convolution, ≥ 1 mm radial travel) during the 800–850 °C vacuum anneal. Pd–Pt is a continuous solid solution, so the bond forms under a light dead-weight during the same furnace step.
  - The Pt annulus does not hydride. The rigid seal therefore carries no misfit, and the corrugation absorbs the ≈ 0.45 mm radial growth of the loaded disk (M3's "no rigid clamp on loading Pd").
- **Seals.** The annulus is clamped between the cell flange and the front flange with **Au-wire seals** (CF knife-edge optional). No elastomer anywhere.
- **Reverse-ΔP protection (new).**
  - A 12 µm membrane pushed away from its grid by P_front − P_cell = 50 mbar reaches ≈ 54 MPa (Hencky, a = 10 mm), above the annealed yield.
  - A 20 mbar **reverse-relief burst foil** connects the front to the cell headspace. It is used only in failure. Its He contamination is irrelevant after a failure, and it is logged.
  - A Δp interlock closes the front fill and exhaust valves at ± 40 mbar.

### 2.4 Array and sample count (resolves 1.5)

- **8 active membranes from ≥ 2 screened Pd lots,** run as **two stages of 4** in the same 8-position house.
  - The stage-2 allocation is frozen before stage 1 is unblinded (1C §C.4).
  - Each stage contains, per lot: 1 × H-L; and one of {H-L, H-M, FX}. The totals over both stages are **H-L 4, H-M 2, FX 2**, balanced by lot and position.
- **Controls in every stage:**
  - **T-Pt**: D₂O, Pt membrane of identical geometry, same current waveform;
  - **T-H(L)**: H₂O, H-L type;
  - **T-H(FX)**: H₂O, FX type;
  - **B**: blank telescope, with one quadrant carrying an identical septum and grid (1A-P11).
- **Lot pre-screening** on sibling coupons, before a lot is admitted:
  - loading ≥ 0.95 (D₂O, ≤ 300 mA cm⁻²);
  - swelling;
  - EBSD texture and grain/thickness ratio;
  - ICP-MS impurities.
  - Also request ENEA/SRI archival or processed material (1D-D3).
- **Rim van der Pauw** is valid on single-regime membranes. R(x) is calibrated on sibling coupons for both isotopes and both branches. x is logged at ≥ 1 Hz, and "in-regime hours" (x ≥ 0.90, i ≥ i₀) are reported as the exposure denominator.

### 2.5 Programme re-tiering (resolves 1.7)

| Tier | Arm | Rationale (ADR-003 score, 1D §4) |
|---|---|---|
| **1** | **DFM-8** (§2.1–2.4) | The only instrument that makes H1 nulls informative at a detector-visible location (S_generic), and the most sensitive ⁴He accounting. SRI term ≈ 0.53 with single-regime membranes, 2 lots and 6 weeks |
| **1** | **C3-G ×5**: 3 D₂ + 1 H₂ + 1 no-CaO. **Cs→Pr primary** (Toyota-replicated; Cs added off-site as sample prep), ⁸⁶Sr→⁹⁴Mo secondary. MHI substrate recipe. **No Mo anywhere in C3-G**; Zr survey | Highest single score (≈ 1.10) at the lowest cost |
| **1** | **C1 ×4 + 1 H₂O twin**: 1.0 M LiOD (SRI lineage), 2 screened lots, one cell with an Al additive, one with a SuperWave-type waveform, ≥ 6 weeks, cathodic current steps, Seebeck twin calorimetry, headspace ⁴He, post-run melts | ≈ 0.71; the "not a closed cathode" objection can only be answered here |
| 2 | Nanocomposite / Ni–Cu furnace | High claimed p, but P(cond) depends on claimant-lineage material that is not in hand. It is promoted to Tier 1 only once such material is secured |
| 3 | Co-deposition on PSD scintillator | Cheap, low credibility |

### 2.6 Inference (resolves 1.4, 1.6)

- **The frozen pre-registration is `docs/design/preregistration-v1.md`,** adapted from 1C §A–B. The primary family is T1 (charged particles, pooled over all active membranes and quadrants) and T3 (⁴He), with weighted Bonferroni and α-propagation, global one-sided α = 2.87×10⁻⁷.
- **T1 background model** is pooled: B telescope; T-Pt (carries the D₂O recoils and the current/EMI waveform); P0a pre-run; current-off blocks of FX membranes (≥ 25 % of their P2 live time); and pre-/post-run blocks of H membranes. Any measured n–p term is fitted.
- **T-H cells become isotope-specificity checks** (claim condition A.3.5), not the subtraction. Each runs in two pre-registered blocks, current-matched and x/J-matched, under a deterministic controller frozen after P0b (1C-S9).
- **T2 (511–511)** enters the primary family only if M7 is merged and its analysis is frozen before P1.
- **Blinding:**
  - a custodian holds the cell map, salt file and He codes;
  - software and hardware salting;
  - hidden ⁴He spikes;
  - two independent analysis teams;
  - list-mode data and code are released at unblinding.
- **Claim scaling.** Every claim is normalised per area and per volume. Heat and the Tohoku-static prediction are labelled "partially tested" (1C-S6).

### 2.7 Protocol changes

| Change | Source |
|---|---|
| Cathodic-only modulation on H membranes (e.g. 100 ↔ 500 mA cm⁻² steps). Anodic excursions only on FX (flux lock-in) | 1D-D7 |
| Pre-registered cathodic excursions to 0.5 A cm⁻² (≤ 2 h, flange cooling, thermal limit confirmed on the P0b witness) | 1D-D6 |
| P2 = **42 d fixed** (not data-dependent). In-regime hours reported | 1D-D8 with 1C-S11 |
| P4 warm phase: fill raised to 1.0 bar (cell and front, balanced), **60 °C, 14 d** | 1D-D12 |
| **P6a desorption window**: current off, front pumped through the HR-QMS line, telescopes on, dynamic He logged | 1D-D13 |
| Al additive as a descriptive variant on one C1 cell and one stage-2 H-L DFM; SuperWave-type waveform on one C1 cell | 1D-D5, D10 |
| Certified low-tritium D₂O (≥ 99.9 % D), LSC assay per lot, ³He ingrowth correction, weekly H/D in headspace aliquots | 1A-P9, 1C-S10, 1D-D15 |
| ¹⁰B markers removed from active membranes (kept on B and on the calibration membrane). One ⁶LiF calibration membrane per lot at 12 µm | 1A-P4, P8 |
| ²²⁸Th/²²⁶Ra assay of LiOD and D₂O | 1A-P8 |
| ≥ 2 AE sensors per cell, so events can be localised | 1A-P3 |
| Per-skin He release fractions; X relies on melts | 1A-P12 |
| Loaded thickness (×1.039) used in all transport | 1A-P10 |

## 3. Findings not adopted, or modified

| Finding | Disposition |
|---|---|
| 1C: "3 skins L, F, X" | Modified to L, M (conditional on P0b), FX. The SRI regime carries the most ADR-003 weight, and M is the only DFM flux-at-high-loading state. X survives only inside FX membranes |
| 1D: raise P4 fill to ≥ 1.2 bar | Modified to 1.0 bar at 60 °C. At 1.2 bar the loss-of-front case (≈ 50 MPa at 12 µm) exceeds the annealed yield. 1.0 bar gives ≈ 44 MPa, below the post-transit yield of ~150 MPa |
| 1D: ring anode to avoid Pt redeposit | Not adopted. The M2 planar mesh gives ± 2.2 % uniformity, which a single-regime membrane needs. A post-run XPS check for Pt on every entry face is pre-registered instead |
| 1D: co-deposited Pd-black entry layer | Deferred to iteration 2 ("morphology only") |
| 1A: honeycomb collimator on one quadrant for depth resolution | Deferred. The response-matrix unfold ε(z, E) is pre-registered instead |

## 4. Consequences

- **Cost.** Tier 1 grows to ≈ $0.55M (iteration-1 §12).
- **Cheaper subset.** The minimum claim-weighted subset, in ADR-003 score per dollar, is **C3-G ×5 + C1 ×2 + H₂O twin + one DFM stage with 2 membranes** (≈ $0.28M).
- **Schedule.** Each DFM stage takes ≈ 16 weeks: commissioning, 42 d P2, P3, P4, P6. Two stages take ≈ 8 months.
- **Honest prior.** Unchanged (≲ 3 % for any claim-rule positive). An all-null outcome yields magnitude limits conditional on the regime being reached, with a Bayes factor of order 1.4–1.6 against existence (1C-S7). It is **not** a refutation, and the write-up must say so.

## 5. Traceability

Each change above is carried into `docs/design/iteration-1.md` (rev B) and cites the finding ID.

## 6. Open items

- Red-team-1B re-run on rev B: membrane and bond mechanics, Si in 0.5 bar D₂, thermal limits of the lined cell at 0.5 A cm⁻², He system, safety FMEA, budget.
- M7 merge: 511 layout, Pb/Cu shielding, T2 freeze.
- 2-D misfit model for FX membranes (F vs X at equal x) and at the bonded rim (M3 follow-up).
