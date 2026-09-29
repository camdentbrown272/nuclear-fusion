# Red team 1D — Claim fidelity (proponent's brief against iteration-1)

**Reviewer stance:** the strongest well-informed LENR proponent (McKubre/SRI, Violante/ENEA, Iwamura, Kitamura/Takahashi, Clean Planet, Szpak/Mosier-Boss, Lipson, Czerski, Storms, Hagelstein). The job is to find every reason a null from `iteration-1.md` would be dismissed as out-of-regime.
**Date:** 2026-09-29 · **Inputs:** charter, ADR-002 rev 2, ADR-003, iteration-1, R1/R2/R7, M2, M3, M6, M8, red-team-0.
**Calculations:** new numbers were computed with `sim/m3_common.Material` (15 µm PdD, 25 °C) and the ADR-003 score formula; the commands are in the appendix. † marks literature values from memory that were not re-verified in this session.

---

## Summary verdict

The DFM-4 does not reproduce the SRI/ENEA regime over most of its entry face, and it cannot show that it does anywhere. The design's own drained-quadrant boundary condition (x_exit ≈ 0.63 on F/X, set by 0.5 bar D₂) combines with 15 µm of Pd to give a diffusion capacity of about 3.8 A cm⁻². That is more than 10× what the electrolysis can supply, so the whole column over F, X and ED, **entry face included**, sits at x ≈ 0.63–0.65 (computed below). Only L and M, 7 of 16 quadrant slots, can be at x ≥ 0.9. M reaches it only with a better surface than P2 will produce. L has no steady flux.

Loading is read by van der Pauw on an un-electrolysed rim that averages over all four regimes, so McKubre's criterion cannot be verified per sample. Every array membrane comes from one commercial, unscreened lot, so the design's "N = 4" is N_eff ≈ 1 for the dominant hidden variable. SRI and ENEA reported that variable directly: behaviour "varies very much with the lots", and lots without the favoured features gave nothing in either lab.

Add 22 °C, i ≤ 0.3 A cm⁻² (below SRI's i₀ = 0.4 for some cathodes), no Al/Si additive, no SuperWave, and a 3-week hold. A proponent would dismiss a DFM null on the SRI claim within one paragraph. The DFM's "Iwamura" and "Lipson" quadrants are not those experiments at all.

Under ADR-003's own scoring, the DFM as specified contributes about 0.4 of the programme's claim-score of about 2.3 while consuming about 75 % of the budget. C3-G (Iwamura) and a multi-cell C1 replica each score higher than the whole Tier-1 array. **Promote C3-G and a 4-cell, 2-lot C1 to Tier 1. Convert ≥ 2 DFM membranes to single-regime (L/M) faces. Screen lots before building.**

---

## 1. Claim-by-claim fidelity (as specified)

Fidelity is P(cond_c | d) in ADR-003, judged at a detector-visible location, 0–1.

| Claim (ADR-003) | Claimed conditions | Iteration-1 provides | Fidelity | Justification |
|---|---|---|---|---|
| **SRI/ENEA heat + ⁴He** (B−/C) | Bulk-average x ≥ 0.90 (best cells 0.95–1.0) by 4-wire R along the cathode. i above i₀ = 0.1–0.4 A cm⁻² (SRI ran to 0.5–1 A cm⁻²). Loading held "several weeks to a month" (R1 §3.1). Non-zero interfacial flux \|dx/dt\| from current steps. Closed Ø1 mm × 30 mm wire (SRI) or closed 50 µm foil (ENEA), vertical. **1.0 M LiOD at SRI**†, 0.1 M at ENEA/F-P. SRI Al/Si additions†. Temperature ≥ 30 °C typical. ENEA-processed or pre-screened lots (<100> texture, grains < 100 µm). The ICCF-14 ENEA/SRI lot replications ran SuperWave (see §5) | 15 µm permeating disk, horizontal, face up, 3.14 cm². 1.0 M LiOD. 200–300 mA cm⁻². 22 °C. One commercial lot, unscreened. No additive. P2 ≥ 3 wk. Current square wave with **anodic** −100 mA cm⁻² half-cycles. Loading read on the rim | **DFM 0.20** (L/M only). **C1 twin 0.40** | Only 7/16 of the entry area can be in regime (§2). x_in ≈ 0.915–0.924 with a "good" surface (M2 §5.1), so the (x−x₀)² factor is 5–25× below SRI's best cells. Below i₀ for x₀ = 0.832/i₀ = 0.4 cathodes. Lot risk is uncontrolled. C1 is closer in geometry but runs at 0.1 M (the F-P/ENEA lineage, not SRI's), N = 1, with no lot control |
| **Iwamura Pd/CaO permeation transmutation** (B, contested) | D₂ gas 1 atm upstream, vacuum downstream, ~70 °C, Pd 25×25×0.1 mm substrate, multilayer + Cs/Sr on the **entry** face, J = 3–7×10¹⁷ D cm⁻² s⁻¹, ~1–2 weeks per sample. MHI substrate prep (anneal, etch†). Cs→Pr is the channel Toyota replicated by ICP-MS (JJAP 2013), Sr→Mo the secondary | **DFM X:** multilayer on the exit face, electrolytic entry, 22 °C, 0.5 bar D₂ downstream, no target element. **C3-G:** 100 µm Pd, 1 atm D₂, 70 °C, multilayer + ⁸⁶Sr on the entry face, J target matches; H₂, no-CaO and static controls; N and Cs not specified | **DFM X 0.03. C3-G 0.60** | The DFM X quadrant reverses the geometry, has no transmuting target, and runs at the wrong temperature, phase and gas side: it is not the Iwamura experiment, yet §1/§10 claim an "Iwamura interface" limit from it. C3-G is close, but drops the replicated Cs→Pr channel for a novel ⁸⁶Sr tag that the claimants never ran |
| **Kitamura/Takahashi CNZ/PNZ** (B in-network) | 100–200 g of claimant-type (Santoku/Technova-lineage†) Pd–Ni or Cu–Ni/ZrO₂, calcined/baked pretreatment, H₂ or D₂ at 0.1–1 MPa, 200–300 °C (up to 450), weeks, oil-flow calorimetry | Tier 3 outline, 100 g, 200–300 °C, flow jacket, inert twin | **0.35** | Right regime, but powder provenance and pretreatment are unspecified, and N = 1 charge. Heat also claimed in H₂, which the ⁴He accounting cannot test |
| **Tohoku / Clean Planet Ni/Cu** (C) | 6×[Cu 2 / Ni 14 nm] on 25×25×0.1 mm Ni, H₂ loaded ~250–270 °C for ~16 h, evacuated, heated to 500–900 °C, repeated cycles, days–weeks | Tier 3 outline, as claimed | **0.50** | Faithful on paper. Radiometric calorimetry is emissivity-limited (R2); sample count unstated |
| **SPAWAR co-deposition** (C/D) | PdCl₂ + LiCl in D₂O on Au/Ag, stepped current, CR-39 in contact, days | Tier 3 PSD-scintillator arm, open cell | **0.60** | Chemistry is faithful. The detector is not CR-39, so a null can be argued away by "the PSD threshold misses the pits' low-energy/heavy component" |
| **Lipson Pd/PdO:D desorption** (C/D) | Au/Pd/PdO heterostructure, electrolytic loading, then **controlled exothermic desorption** facing Si/CR-39; 3 MeV p, 11–20 MeV α claimed | F quadrant: bare Pd exit in 0.5 bar D₂ during electrolysis. PdO and the Au/Pd/PdO skin are rejected (M3, M6). No desorption window | **0.10** | Wrong surface, wrong phase (loading, not desorption), and exit gas instead of vacuum |
| **Czerski beam-free 511 keV** (D) | D-loaded Pd/Zr, no beam, very low (underground) background, 511–511 coincidences (ICCF-26) | LaBr₃ pair at a surface lab, Cu shield, muon veto | **0.30** | Conditions are loose; background is the fidelity issue. M7 is still pending |
| **NTT out-diffusion bursts** (D) | Pd with Au on one face, oxide on the other, gas-loaded, then out-diffusion into **vacuum** | L quadrant Au cap; front at 0.5 bar D₂; no evacuation step | **0.10** | No out-diffusion into vacuum is ever performed |

---

## 2. Is McKubre's criterion met, quadrant by quadrant?

McKubre's criterion has four parts: (i) x ≥ 0.9 bulk-average, (ii) i > i₀, (iii) t ≥ 10 τ_D (in practice weeks), (iv) non-zero flux.

**The 15 µm column follows its exit boundary condition.** The steady flux through 15 µm PdD at 25 °C comes from Φ(x_in) − Φ(x_exit) (`sim/m3_common.Material`):

| x_in → x_exit | J (mA cm⁻² equivalent) |
|---|---|
| 0.92 → 0.63 (design's stated F/X exit) | **3 940** |
| 0.95 → 0.63 | 4 096 |
| 0.92 → 0.90 | 137 |
| 0.95 → 0.90 | 293 |

The cell passes at most 200–300 mA cm⁻² in total, and only part of that is absorbed. Solving for x_in over a quadrant whose exit sits at 0.63 gives:

| Absorbed flux (mA cm⁻²) | 10 | 50 | 150 | 300 |
|---|---|---|---|---|
| x_in | 0.631 | 0.633 | 0.638 | **0.646** |

Lateral coupling is about one membrane thickness (red-team-0 §3.1), so the entry face above F, X and ED sits at **x ≈ 0.64**. That is β-phase minimum, far below x₀. The design says the telescope sees "the electrolyte-facing surface where the classic high-loading claims were made". That holds only over L and M.

| Quadrant | x_entry | x_exit | Steady J | (i) x ≥ 0.9 | (ii) i > i₀ | (iv) flux | Verdict |
|---|---|---|---|---|---|---|---|
| **L** (Au cap) | 0.915–0.924 (good surface, 200–300 mA cm⁻², M2) | same | ≈ 0 | yes, marginal | marginal (≤ 0.3 A cm⁻²) | only during current steps | **Closest SRI analogue** (a closed cathode with transient dx/dt). It is not called that, and P3's anodic half-cycles deload it |
| **M** (Ni 2–20 nm) | needs ≥ 0.95 per M3, but only 0.92 is reachable | ~0.915 with 20 nm Ni at x_in 0.92 (the 0.92→0.90 capacity is 137 mA cm⁻², so J ≈ 10 mA cm⁻² barely tilts the column) | ~10 mA cm⁻² | marginal | marginal | steady through-flux | Can hold x ≈ 0.91–0.92 *with* flux if Ni is ≥ 20 nm. Contrary to §15 risk 1, M does not collapse at x_in = 0.92; it runs at a lower (x − x₀)². Proponents will argue that steady through-flux is not SRI's \|dx/dt\| |
| **F** (bare) | **0.64** | 0.63 (design) | supply-limited | **no** | — | yes | Out of the SRI regime. If M3's saturation-cap branch holds instead, F is loaded but has ~1 mA cm⁻² flux, and the design's F label is wrong. Either way it is not the claimed state |
| **X** (CaO ML) | **0.64** | 0.63 | supply-limited | **no** | — | yes | Out of regime for SRI. For Iwamura it is the wrong face, temperature and gas |
| **ED** (1 µm e-dep Pd) | **0.64** (behaves like F) | ≈ 0.63 | supply-limited | **no** | — | yes | Out of regime |

- **Criterion (iii)** is trivially met in charge: 10⁷ C per mol Pd (Cravens–Letts) is 5.3 kC for a 4.7 mm³ membrane, about 2 h at 0.8 A. It is not met in calendar time. P2 is 3 weeks, while McKubre gives "several weeks to a month" and F-P runs went to months. R1 §3.2 notes that the incubation is not diffusion-limited, so thinness buys nothing.
- **The loading gauge is invalid.** The van der Pauw contacts sit on the Ø26 rim, outside the Ø20 wetted area and under the seal. The rim loads only laterally, so it averages over the L/M/F/X column states. R/R₀ is double-valued, and the rim mixes a 0.92 region with 0.64 regions. The measured "mean x" therefore does not test any quadrant against x ≥ 0.9. **A proponent will say, correctly, that the design never demonstrated the loading at the location where it reports a null.**
- **M2/M3 numbers, critically.** The "good" surface class (x = 0.90 at 86 mA cm⁻²) is fitted to ENEA-type *vacuum-annealed + etched* samples (M2 §4). It is calibration, not prediction. A single commercial lot, anodised by P3 half-cycles and exposed to Pt redeposition from a mesh 4 mm above a horizontal face-up disk, will more likely behave as "typical" (x ≈ 0.87 at 300 mA cm⁻², M2 §5.1). In that case **no quadrant reaches 0.9.** The only quantitative lever M2 offers is a recombination poison (P = 10 → +0.065). The design omits it.

---

## 3. Ranked findings

| ID | Severity | Finding | Evidence | Concrete fix |
|---|---|---|---|---|
| D1 | **Critical** | Entry face over F/X/ED is at x ≈ 0.64: 9 of 16 slots are outside the SRI regime, contradicting §0/§1/§3.4 | §2 table; M3 Kirchhoff capacity 3.9 A cm⁻² vs ≤ 0.3 A cm⁻² supply | Make ≥ 2 of 4 active membranes **single-regime L or M over the whole face**. Keep quadrant membranes for H1 site comparison only, and drop "SRI at the entry face" claims for F/X/ED |
| D2 | **Critical** | Loading is not verifiable where the null is reported: rim vdP averages mixed regimes | §3.2 "4-point van der Pauw on the membrane rim" | Single-regime membranes (D1) make rim R meaningful. Also calibrate R(x) on sibling coupons, including the α→β→high-x branch in D and H. Log x per membrane at ≥ 1 Hz; pre-register "in-regime hours" (x ≥ 0.9, i ≥ i₀) as the exposure denominator |
| D3 | **Critical** | One commercial lot per array, so N_eff = 1 for the dominant hidden variable. ADR-003's [1 − (1 − p)^N] is overstated: 0.59 becomes 0.20 | McKubre, ICCF-14 ENEA/SRI replications: L14 gave ~80 % excess at both labs; L17 13 %/500 %; L19 43 %/100 %; "behaviour varies very much with the lots", and other lots gave none at either lab (McKubre, *One perspective…*). R1 G5–G8 | **≥ 2 lots, ideally 3**, balanced across positions. **Pre-screen** each lot on sibling coupons: x ≥ 0.95 in D₂O (or H/Pd ≥ 0.95) at ≤ 300 mA cm⁻², low swelling, EBSD <100> texture, ICP-MS impurities. Ask ENEA/SRI for archival or processed lot material. Report per-lot denominators |
| D4 | **Major** | Allocation: under ADR-003, C3-G (0.54 → 1.10 with N = 3) and a 4-cell C1 (0.71) each outscore the whole DFM claim-term (~0.42) at ~15–25 % of its cost | §4 table | Promote **C3-G (N ≥ 3, Cs→Pr primary)** and **C1 ×4 (2 lots)** to Tier 1. Fund by cutting DFM Tier 1b to Tier 1a (4 → 2 membranes plus twins) if the budget is fixed |
| D5 | **Major** | x_in only 0.915–0.924 even for a "good" surface; (x − x₀)² is 5–25× below SRI's best cells. No additive, although SRI used Al (and Si from glass)† | M2 §5.1, §6.3 rec. 13; R1 table (SRI "Al/Si additives") | Add an SRI-type Al additive (as LiAlO₂, order 10² ppm†) to half the D₂O cells, and the same to the H₂O twin. Treat a thiourea-class poison as a pre-registered variant. Target x ≥ 0.95 |
| D6 | **Major** | Current density ≤ 0.3 A cm⁻² sits at or below SRI i₀ (0.4 A cm⁻² in the x₀ = 0.832 fit); SRI ran to 0.5–1 A cm⁻² | R1 §3.1 | Pre-register **cathodic-only excursions to 0.5–0.8 A cm⁻²** (hours) on L/M membranes, with active flange cooling. The M2 ΔT bound scales ~i², ≈ 45 K at 0.6 A cm⁻² conduction-only, so M4 must confirm |
| D7 | **Major** | P3 square wave +300/−100 mA cm⁻² anodically deloads a 15 µm membrane every half-cycle (τ_D ≈ 2 s). Repeated deloading degrades maximum x (Storms; R1 G7). SRI's flux came from cathodic steps | iteration-1 §9 P3; M3 §5.2 transients | Use cathodic-only steps (e.g. 100 ↔ 500 mA cm⁻²). Put anodic excursions only on the site-factory membrane |
| D8 | **Major** | Duration: P2 is 3 weeks and the whole electrolysis is ~6 weeks, against "several weeks to a month" of *in-regime* hold before initiation (SRI) and months (F-P) | R1 §3.1–3.2; Cravens–Letts criterion (ii) | P2 ≥ 6–8 weeks on L/M membranes, with a pre-registered stopping rule on in-regime hours (≥ 1000 h at x ≥ 0.9), not calendar time |
| D9 | **Major** | The Iwamura test is compromised: the DFM X quadrant is not an Iwamura experiment, yet §1/§10 claim limits from it. C3-G replaces the replicated Cs→Pr channel with a novel ⁸⁶Sr tag | R2 #12–13 (Toyota: Pr ≤ 2×10¹¹ → 1.6×10¹² cm⁻², 3/3 runs, H₂ null) | C3-G primary: **Cs→Pr** (¹³³Cs and ¹⁴¹Pr are monoisotopic, and the Pr background is measurable pre-run), with Toyota-style Cs ion-implantation done off-site as sample prep. ⁸⁶Sr→⁹⁴Mo becomes the secondary. N ≥ 3 D₂ + 1 H₂ + 1 no-CaO. Follow the MHI substrate recipe. Delete Iwamura claims from the DFM X row |
| D10 | **Major** | The ENEA/SRI lot replications that proponents cite ran Dardik SuperWave current; the design has DC + square wave | McKubre ICCF-14 (ENEA/SRI replication of Energetics); R1 table (Energetics) | Run a SuperWave-type waveform on one C1 cell and one L-regime DFM, with ⟨V·I⟩ at ≥ 100 kS/s (already specified) |
| D11 | **Major** | C1 "SRI replica" uses 0.1 M LiOD for "lineage", but SRI's lineage is 1.0 M†. N = 1, no lot control, no incubation spec | iteration-1 §7.2; R1 table row SRI | C1 at **1.0 M LiOD** (+ Al variant); 4 cells over 2 screened lots; ≥ 6 weeks; cathodic current steps; 1 H₂O twin |
| D12 | **Major** | P4 warm phase at 70–80 °C in a 0.5 bar-abs cell. D₂O vapour pressure at 80 °C is ≈ 0.43–0.47 bar†, so the cell is near boiling and the ΔP ≥ 10 % trip, O₂ monitor and recombiner are confounded by vapour | iteration-1 §3.2 cell gas; §9 P4 | Raise the cell and front fill to ≥ 1.2 bar abs before P4 (the front already steps to 1.0 bar in P3), or cap P4 at 55 °C. Extend P4 to ≥ 2 weeks: Storms/F-P report a positive temperature coefficient of heat, so one week is a token |
| D13 | **Major** | No Lipson protocol: Lipson's signal is during **controlled exothermic desorption** from Pd/PdO:D, while F is a loading-phase face in 0.5 bar D₂ | Lipson et al., ICCF-12 (Pd/PdO:Dₓ, SSB + CR-39, 3 MeV p and 11–20 MeV α) | Add a P6a **desorption window**: current off, front pumped to vacuum (the worst-case ΔP of 0.5 bar gives 24 MPa, which the design already accepts), telescopes on, over the full deload for all A-cells and T-H. It costs no hardware. It also gives NTT out-diffusion into vacuum at the L (Au) quadrant |
| D14 | Minor | Horizontal, face-up entry disk under a Pt mesh 4 mm away collects Pt redeposit and particulates, which lowers x (R1 G11, §5.14); SRI/ENEA cathodes were vertical | M2 rec. 7 | Use the M2 alternative ring anode (Ø10 disk, g = 7 mm) on the L/M membranes, or keep anode i ≤ 20 mA cm⁻². Post-run XPS for Pt on every entry face |
| D15 | Minor | D₂O isotopic purity and in-run H ingress are not specified; H > 1 % is reported to suppress heat (R1 G19). The H₂O twin shares the house and manifold | R1 G19 | Spec ≥ 99.9 % D. Measure H/D in headspace aliquots weekly (the HR-QMS already exists) |
| D16 | Minor | 15 µm foil with 20–50 µm grains gives a bamboo microstructure (through-thickness boundaries), unlike ENEA's 50 µm/<100 µm-grain foils; texture after rolling to 15 µm is unknown | iteration-1 §3.2; R1 G5–G6 | EBSD each lot at thickness. Report texture and grain/thickness ratio as covariates |
| D17 | Minor | Magnitude comparison is done against whole-cell W claims. A proponent will rescale by active volume (DFM 4.7 mm³ × 7/16 in regime vs SRI 24 mm³): ~0.1× | geometry | Pre-register the claim scaling (per cm² and per cm³). Sensitivity is still saturated (⁴He at 10⁻⁹), so this is presentation only |
| D18 | Minor | Co-deposition was moved out of the closed cell (Cl₂ poison). The telescope-visible entry face therefore never carries SPAWAR/Letts-type dendritic Pd | iteration-1 §8.1 | Optional: one membrane with an ex-situ Pd-black entry layer (co-deposited from PdCl₂/LiCl in D₂O, rinsed Cl-free) is visible at 2.05 MeV through 15 µm. Pre-register it as "morphology-only", not a SPAWAR test |

---

## 4. Allocation under ADR-003 (Q3)

Score term per claim: w·P(cond)·[1 − (1 − p)^N]·min(1, P_det)·Cred. P_det = 1 everywhere (§10 of the design shows each channel ≤ 1 % of the claimed magnitude). Weights: B = 3, B−/C = 2.5, C = 2, C/D = 1.5, D = 1. The claimed per-sample p is 0.5 where claimants report most runs positive (Toyota 3/3; the Kobe and Tohoku in-network series) and 0.2 otherwise. P(cond) and Cred are this reviewer's estimates from §1, and they are the inputs most open to dispute.

| Claim | Arm | w | P(cond) | p | N | 1−(1−p)^N | Cred | Score |
|---|---|---|---|---|---|---|---|---|
| SRI/ENEA | DFM as specified (1 lot → N_eff = 1) | 2.5 | 0.20 | 0.2 | 1 | 0.20 | 0.9 | **0.09** |
| SRI/ENEA | DFM as specified, if lots were independent | 2.5 | 0.20 | 0.2 | 4 | 0.59 | 0.9 | 0.27 |
| SRI/ENEA | DFM single-regime L/M, 2 screened lots, 6 wk, Al | 2.5 | 0.40 | 0.2 | 4 | 0.59 | 0.9 | **0.53** |
| SRI/ENEA | C1 as specified (1 cell, 0.1 M) | 2.5 | 0.40 | 0.2 | 1 | 0.20 | 0.8 | 0.16 |
| SRI/ENEA | C1 ×4, 1.0 M ± Al, 2 lots, SuperWave, ≥ 6 wk | 2.5 | 0.60 | 0.2 | 4 | 0.59 | 0.8 | **0.71** |
| Iwamura | DFM X quadrant | 3 | 0.03 | 0.5 | 3 | 0.88 | 0.9 | 0.07 |
| Iwamura | C3-G as specified (1 D₂ cell) | 3 | 0.60 | 0.5 | 1 | 0.50 | 0.6 | 0.54 |
| Iwamura | C3-G, Cs→Pr + ⁸⁶Sr, N = 3 | 3 | 0.70 | 0.5 | 3 | 0.88 | 0.6 | **1.10** |
| Kitamura | nanocomposite furnace, 1 charge | 3 | 0.35 | 0.5 | 1 | 0.50 | 0.7 | 0.37 |
| Tohoku Ni/Cu | furnace variant, 2 samples | 2 | 0.50 | 0.5 | 2 | 0.75 | 0.6 | 0.45 |
| SPAWAR | PSD arm, N = 4 | 1.5 | 0.60 | 0.2 | 4 | 0.59 | 0.6 | 0.32 |
| Lipson | DFM F quadrant (no desorption) | 1.5 | 0.10 | 0.2 | 3 | 0.49 | 0.9 | 0.07 |
| Lipson | DFM + P6a desorption window | 1.5 | 0.30–0.40 | 0.2 | 4 | 0.59 | 0.9 | 0.24–0.32 |
| Czerski 511 | DFM LaBr₃ pair | 1 | 0.30 | 0.2 | 4 | 0.59 | 0.8 | 0.14 |
| NTT | DFM L quadrant, 0.5 bar front | 1 | 0.10 | 0.2 | 4 | 0.59 | 0.8 | 0.05 |

**Reading:**
- **DFM Tier 1 as specified** scores Σ ≈ 0.09 + 0.07 + 0.07 + 0.14 + 0.05 = **0.42**, at ≈ $226–364k.
- **Tiers 2–3 as specified** score Σ ≈ 0.54 + 0.16 + 0.37 + 0.45 + 0.32 = **1.84**, at ≈ $135–165k.
- The claim-weighted objective the project adopted therefore says the Tier-1 geometry is the *least* productive arm per dollar. Its unique value lies in H1 spectroscopy (S_generic, weight 0.1) and in making H1 nulls informative. That is legitimate, but it is not what ADR-003 optimises.
- **Recommended re-tiering:**
  - Tier 1 = C3-G (N ≥ 3, Cs→Pr) + C1 ×4 (2 lots, 1.0 M) + DFM Tier 1a with ≥ 2 single-regime membranes and the D13 desorption window.
  - Tier 2 = the nanocomposite/Ni-Cu furnace. It is promoted from Tier 3 because W-scale claims with p ≈ 0.5 are the cheapest high-score term, but it stays at Tier 2 until claimant-lineage powder or multilayer samples are secured. Without them P(cond) drops to ≲ 0.15.
  - The ⁴He manifold is shared, so the marginal cost of C1 ×4 and C3-G ×3 is dominated by cells and external analyses, not new detectors.
- **Should C1 or the nanocomposite arm be Tier 1?** C1: yes (0.71 vs DFM-SRI 0.09–0.53). Nanocomposite: not yet. Its score depends almost entirely on material provenance, which is not in hand, and one furnace charge gives N = 1.

---

## 5. Minimal fidelity-raising changes that keep detector visibility (Q4)

| Change | Raises | Visibility cost |
|---|---|---|
| **Single-regime membranes** (≥ 2 of 4 all-L or all-M; Ni ≥ 20 nm for M) | SRI P(cond) 0.2 → ~0.4; makes rim R valid | Loses the in-membrane sector control on those cells only (M6 already showed sectors do not share conditions) |
| **Lot pre-screen, ≥ 2 lots, ENEA/SRI archival material if obtainable** | N_eff 1 → 2–4 | None |
| **Al additive (± thiourea variant), cathodic-only steps, 0.5–0.8 A cm⁻² excursions** | x_in toward 0.95; i above i₀ | Thermal (M4 check); none to detectors |
| **P2 ≥ 6–8 wk on in-regime hours; P4 ≥ 2 wk at ≥ 1.2 bar fill** | Incubation, temperature | Schedule only |
| **P6a desorption window into vacuum** | Lipson, NTT | None (uses existing telescopes) |
| **C3-G: Cs→Pr primary, N ≥ 3** | Iwamura P(cond) 0.6 → 0.7, N 1 → 3 | None |
| **C1 ×4 at 1.0 M in the neutron bank** | SRI P(cond) 0.6 at N = 4 | C1 is blind to charged particles (accepted; the DFM covers H1) |
| *Rejected:* a thicker DFM cathode | — | Kills entry-face visibility (3 MeV p range in Pd is 31 µm, ADR-002) with no gain in x (thinness is not what limits x; the surface is) |
| *Rejected:* codeposition on a closed-cell entry face | — | Cl₂ poisons the recombiner (design §8.1 is right) |

---

## 6. What proponents will say about a null, and what can be pre-empted (Q5)

| Objection | Valid against iteration-1? | Pre-emptable? |
|---|---|---|
| "You never showed x ≥ 0.9 where you looked" | **Yes** (D1, D2) | Yes: single-regime membranes, coupon-calibrated R(x), in-regime-hour denominators |
| "Wrong palladium: one commercial lot, unscreened" | **Yes** (D3) | Yes: ≥ 2 screened lots, archival lot request, per-lot denominators |
| "Current below threshold; anodic cycling destroyed loading" | Yes (D6, D7) | Yes |
| "Too short; incubation takes weeks to months" | Yes (D8) | Partly: ≥ 6–8 weeks in regime. "Months" can always be demanded, so pre-register the stopping rule and cite SRI's own "weeks to a month" |
| "Clean PTFE cell lacks the Al/Si/B surface chemistry of glass cells" | Yes (D5) | Yes: Al variant (glass itself is barred by He permeation, M8) |
| "No trigger (SuperWave/laser)" | Partly (D10) | SuperWave: yes, cheaply. Laser/THz: M6 excluded it on physics grounds; say so explicitly in the pre-registration |
| "A permeating 15 µm membrane is not a cathode; D leaks through the sieve" | Yes, for F/X/ED | Only by C1 ×4 (closed wires). The L membrane is closed, so the DFM can answer this for L only |
| "Wrong temperature" | Partly (D12) | Yes: longer, pressurised warm phase |
| "Iwamura tested on the wrong face / wrong element" | Yes (D9) | Yes: C3-G with Cs→Pr |
| "Transmutation null from a single sample" | Yes | Yes: N ≥ 3 |
| "Kitamura/Clean Planet materials not ours" | Will be, unless sourced | Only by obtaining claimant-lineage material or published recipes with verification (particle size, calcination, XRD) |
| "He retained in the Pd" | No: melts and release calibration (M8) are specified | Already pre-empted |
| "Detector insensitive" | No (§10: ≤ 10⁻² of every claim) | Already pre-empted |
| "H contamination suppressed the effect" | Possibly (D15) | Yes: D purity spec plus weekly H/D |

**Residual that cannot be pre-empted:** the "hidden variable" defence (SRI's cell constant M, ENEA's "material with the right features remains an open problem"). The only answer is denominators across lots and cells, published whatever the outcome. The design should commit to that in the pre-registration.

---

## Sources

- McKubre, "Cold Fusion (LENR): One Perspective on the State of the Science" (ICCF-15). ENEA/SRI lot replications L14/L17/L19 and lot dependence: https://lenr-canr.org/acrobat/McKubreMCHcoldfusionb.pdf
- Cravens & Letts, "The Enabling Criteria of Electrochemical Heat" (ICCF-10). Criteria (i)–(iv) and 10⁷ C per mol Pd: https://www.lenr-canr.org/acrobat/CravensDtheenablin.pdf
- Violante et al., ICCF-15 proceedings (material features, reproducibility "open problem"): https://lenr-canr.org/acrobat/ViolanteVproceeding.pdf
- Current Science 108 (2015) LENR special section (loading x > 0.85–0.9 as whole-cathode average): https://www.currentscience.ac.in/Volumes/108/04/0495.pdf
- Iwamura et al., JCMNS 10 (2013) 63 (permeation transmutation): https://jcmns.org/article/72216-recent-advances-in-deuterium-permeation-transmutation-experiments/attachment/150206.pdf ; Toyota replication (Hioki et al., JJAP 2013) summary: https://news.newenergytimes.net/2013/10/22/journal-publishes-toyotas-independent-replication-of-mitsubishi-lenr-transmutation/
- Lipson et al., "Reproducible nuclear emissions from Pd/PdO:Dₓ heterostructure during controlled exothermic deuterium desorption" (ICCF-12): https://ui.adsabs.harvard.edu/abs/2006cmns...12..293L/abstract ; Au/Pd/PdO:D: https://www.osti.gov/etdeweb/biblio/20845769
- In-repo: R1 §3.1–3.2, table rows SRI/ENEA/Energetics; R2 #6, #12–16; M2 §5.1, §5.7, §6.3; M3 §5.2 and Table 5.4; M6 summary; M8 summary.

## Appendix: reproduction

```
cd sim && python3 -c "
import m3_common as m; M=m.Material(); L=15e-6
mA=lambda J: J*1.602e-19/1e4*1e3
print(mA((M.Phi(0.92)-M.Phi(0.63))/L))          # 3940 mA/cm2 capacity
for i in (10,50,150,300):
    J=i*1e-3/1.602e-19*1e4; print(i, M.x_of_Phi(M.Phi(0.63)+J*L))   # x_in over a drained quadrant
"
```
