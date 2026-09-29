# Pre-registration v1 — DFM-8 (iteration-1 rev B)

**Status:** draft to be frozen. Before P1 of stage 1 it must become a signed git tag `prereg-v1` plus an OSF registration. After that, amendments are allowed only before unblinding, dated and justified, and unblinding is irreversible.
**Basis:** red-team-1C §A–B (adopted with the ADR-005 modifications).

## 1. Units and analysis sets

- **Unit of replication:** the *membrane*. Quadrants are sub-units: they share entry face, electrolyte, current and lot.
- **Active membranes:** 8, in two stages of 4 (ADR-005 §2.4).
- **ITT set:** all membranes that complete P1.
- **Per-protocol set:** membranes with ≥ 14 d of P2 live time at x ≥ 0.90 by rim van der Pauw. The R(x) calibration is frozen from sibling coupons in P0b. For FX membranes the condition is "J within the P0b-predicted band".
- The **primary analysis uses ITT**. Both sets are reported.

## 2. Primary family

Graphical weighted Bonferroni with α-propagation. Global one-sided α_G = 2.87×10⁻⁷ (5σ).

| Test | Weight | Data | Statistic | Background model | Threshold |
|---|---|---|---|---|---|
| **T1 (H1a)** | 0.5 | P2, current-on, **all active membranes and quadrants pooled** (one test). ΔE–E PID window 1.41–3.10 MeV. Clusters of ≥ 2 events within 10 s in one detector count as one event | Poisson profile-likelihood ratio, s ≥ 0 | Pooled: B telescope, T-Pt (4 quadrants), P0a pre-run, FX current-off blocks (≥ 25 % of FX P2 live time), and pre-/post-run blocks of H membranes. One scale nuisance per component. Global log-normal nuisance with σ = ln 3. Background exposure ≥ 10× signal live time | p ≤ 1.43×10⁻⁷ **and ≥ 10 net events** |
| **T3 (H2)** | 0.5 | ⁴He from front volume and headspace, over all P2–P4 windows. Summed per cell, then over the active cells. Quadrant melts form a pre-specified second component | Σ_cells (N_obs − N̂_blank,cell) / σ_tot. σ_tot combines each cell's blank-series variance with the between-cell variance measured in P0a. A D₂-matrix spike correction is applied per aliquot class | Each cell's own blank series (≥ 3 × 1 d and ≥ 2 × 7 d, taken both pre and post); sibling coupons for the melts | p ≤ 1.43×10⁻⁷ |

- **α-propagation:** if T1 rejects, its α passes to T3, and vice versa.
- **T2 (511–511 coincidence, H1c):** enters only if M7 is merged and its analysis is frozen before P1. It would take weight 0.1 from each of T1 and T3.
- **T4 (C3-G transmutation):** a separately registered study with its own α. Its primary endpoint is ¹⁴¹Pr from ¹³³Cs, measured by blinded two-lab ICP-MS and ToF-SIMS on coded samples. The reach is stated as isotope-ratio precision.

## 3. Conditions for a discovery claim

All of the following are required. Each is one-sided and fixed now.

1. **Primary test.** A primary test rejects at its threshold.
2. **Flux or loading dependence, on independent data.**
   - FX membranes: P3 lock-in at the frozen period (240 s) with the frozen lag model (τ from the P0b permeation transient) gives ≥ 3σ in the same channel.
   - H membranes: P3 cathodic step response, with the same frozen parameters.
   - In both cases, the current-off blocks must be consistent with background: 95 % upper limit < 25 % of the P2 rate.
3. **Replication.**
   - ≥ 2 membranes, each ≥ 2σ locally, from ≥ 2 lots.
   - Heterogeneity χ² across membranes not rejected at p < 0.01.
   - No single cluster or day contributes > 20 % of events.
4. **Signature.**
   - The ΔE–E distribution is consistent with protons from 0–12.5 µm depth (response-matrix unfold, ⁶LiF-calibrated; KS p > 0.01).
   - It is inconsistent with the sideband shape (p < 0.01).
   - Triton and ³He rates are consistent with the proton rate where kinematically accessible.
5. **Controls.** T-Pt, T-H(L), T-H(FX) and B each have a 95 % upper limit per quadrant-exposure below 25 % of the active rate.
6. **Blind integrity.** The salts are recovered within ± 2σ before the box is opened, and the two analysis teams agree within 1σ.
7. **Neutron statement.** M5's pre-registered sentence ("protons below ~10⁻² fusions s⁻¹ appear without detectable neutrons") is evaluated. It is not used as a veto.
8. **Heat.** Any heat statement additionally requires the mass-flow cross-check (M4 rec. 12). Particle products below the calorimetric floor carry no heat claim.

## 4. Reporting tiers

- **Evidence:** a primary test at global ≥ 3σ.
- **Claim:** all conditions in §3 are met.
- **Descriptive only:** per-skin, per-quadrant, per-lot and per-variant results (Al additive, SuperWave), plus neutron, heat, AE, triton/³He lines, and 511 if T2 was not frozen. These are reported with Feldman–Cousins intervals and no significance language.

## 5. Null reporting (committed)

- For every claim in the ADR-003 register, publish Feldman–Cousins 90 %/95 % upper limits, normalised per area and per volume with the scaling stated, whatever the outcome.
- State plainly that an all-null outcome is a magnitude limit conditional on reaching the regime. Its Bayes factor against existence is ~1.4–1.6, and it is **not** a refutation.
- Submit as a registered report where the journal allows it. Release list-mode data, waveforms, He raw data and code at unblinding.

## 6. Stopping and failures

- **Fixed durations per stage:** P2 = 42 d, P3 = 14 d, P4 = 14 d, P6a = 2 d.
- **Interim unblinded looks:** safety only (neutron bank at 10× background, pressure and Δp, temperature).
- **Membrane failure:**
  - before P2: replaced under a new ID;
  - after P2 starts: data truncated at the last DQ-good period, and no replacement enters the primary analysis.
- **Extensions:** none conditioned on data. Any extension is decided on DQ grounds only, by the custodian, blind to box counts.
- **Stage 2:** its allocation is frozen and registered as amendment A1 before stage-1 unblinding.

## 7. Frozen data-quality cuts

These are tuned on P0a, blank-telescope and sideband data only.

- Veto ± 1 s around current edges.
- ΔE leakage-current and noise-rms thresholds.
- Muon-veto coincidence.
- Pulse-shape χ² against pulser templates.
- Cluster flag.
- AE-coincident events (± 50 ms) flagged. They are included in T1, and reported separately as the H1b tag.
- He aliquot acceptance: spike recovery within ± 10 %.

## 8. Blinding and salting

1. **Custodian.** One person, outside the analysis team, holds the cell → channel map (T-H and T-Pt among the active labels), the salt file and the He aliquot codes.
2. **Box.** PID-window events of all channels are hidden. Analysts see sidebands, B, the calibration-marker windows and the live time.
3. **Software salt.** Synthetic ΔE–E events are drawn from the ⁶LiF-validated response. They are injected into 0–2 random channels at {0, 0.5, 1, 2} × the median 5σ reach, with Poisson timing and a random current-on fraction.
4. **Hardware salt.** The ADR-004 source boost (at a partner facility) or cosmic ¹⁰B markers, toggled on the custodian's schedule.
5. **He.** Aliquots are coded. ≥ 20 % are blanks, and ≥ 2 hidden ⁴He spikes of 5–10 σ_B are included.
6. **Two analysis teams.** Independent code bases. At least one external member, ideally a known skeptic.

## 9. Pre-computed reach

`sim/design_iter1_reach.py` (output `figs/iter1_reach.txt`) gives the following for T1, P2 of 42 d at 0.8 live fraction, α = 1.43×10⁻⁷, with the ≥ 10-net-event floor:

| Background per cell (d⁻¹) | Single membrane, exit face | Single membrane, entry face | Pooled 8, per membrane (exit / entry) |
|---|---|---|---|
| 0.03 | 4.7×10⁻⁵ | 6.3×10⁻⁵ | 1.0 / 1.4 ×10⁻⁵ |
| 0.09 (3× band) | 5.5×10⁻⁵ | 7.5×10⁻⁵ | 1.6 / 2.2 ×10⁻⁵ |

Units are D+D fusions s⁻¹. The inputs are:
- ε = 0.22 (exit) or 0.162 (entry), × 0.80 grid transmission × 0.9 septum loss;
- a proton branch of 0.5.
