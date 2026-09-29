# Red-team 1A — Nuclear / condensed-matter physics review of iteration-1 (DFM-4)

**Scope:** `docs/design/iteration-1.md` against the charter, ADR-002 rev 2, ADR-003, ADR-004 and M0–M8.
**Reproduce:** `python3 sim/rt1a_physics.py` writes `docs/design/figs/rt1a_physics.txt` and `figs/rt1a_tomography.png`. The script reuses `m5_stopping`, `m5_si_telescope` (geometry, deposit, PID) and `m3_common` (isotherm, Kirchhoff potential).

## Summary verdict

At normal incidence the charged-particle arithmetic is correct: 3.015 / 2.54 / 2.00 MeV for exit / mid / entry protons through loaded PdD₀.₉ plus 4 mm of 0.5 bar D₂. The physics built around those numbers has three serious problems.

1. **The ⁴He channel is compromised by the cell's own materials.** A PCTFE cell body and an FFKM rim seal put polymer He permeation and outgassing 10⁶–10⁷× above the quoted 0.16–0.39 nW floors in the headspace channel, and 10²–10⁵× above them in the front channel through cell→front crossover. The effect is temperature-dependent, so the P4 warm phase manufactures a correlated heat/He artifact. M4 §5.10 already forbids exactly this.
2. **The "SRI conditions at the entry face" premise holds only under the L and M quadrants.** A bare or oxide exit at 0.5 bar D₂ needs ≈3 A cm⁻² of *absorbed* current to keep x_in = 0.9 across 15 µm. The cell supplies ≤0.15 A cm⁻², so the whole foil under F, X and ED (10 of the 16 rotation slots) sits at x ≈ 0.655, entry face included. The same quadrant structure imposes a 1.6 % in-plane misfit between neighbouring quadrants. That is GPa-level stress against a 40–150 MPa yield, so the foil will plastically deform and crack at the skin boundaries, which is a built-in H1b confound.
3. **The chosen 15 µm thickness breaks the ADR-004 ⁶LiF depth calibration and halves entry-face proton efficiency** relative to what §3.4 implies.

Electrolyte-borne cosmic recoils (n–p in the H₂O twin; n–d and D(n,np) in D₂O) are missing from the background budget and unbalance the D − H control. None of these makes the experiment unbuildable. Findings 1–2 must be fixed before the design is frozen, otherwise H2 and the SRI arm of the claim register score ≈ 0 under ADR-003.

## Ranked findings

| ID | Sev. | Finding | Quantitative evidence | Concrete fix |
|---|---|---|---|---|
| **P1** | **Critical** | Polymer He sources: the PCTFE cell body and the FFKM rim seal swamp both ⁴He channels | Air→cell through PCTFE (71 cm² × 5 mm, 1–30 Barrer): **1.5×10⁸–4.5×10⁹ He/s**. FFKM seal if air-exposed: 1.7×10⁸ /s. Stored air-He in ~36 cm³ of polymer: 1×10¹⁴ atoms (~2×10⁸ /s over ~5 d). Cell→front crossover via FFKM: 10⁴–10⁷ /s. **Floor: 43–102 He/s.** The Kr tracer does not bound He permeation through polymers | All-metal cell: PTFE-*lined* 316L, as in M4 §5.10. Replace the FFKM with a metal seal to a non-hydriding rim (P1-fix in §P1). Bake and pump-purge polymer liners, and run ≥2 blank windows after the fix |
| **P2** | **Critical** | Under F/X/ED the entry face is at x ≈ 0.655, not ≥ 0.9. SRI/ENEA conditions are visible only under L (4 slots) and M (3 slots, conditional). The rim van der Pauw reads a meaningless mean | M3 potential at 295 K: holding x_in = 0.90 against x_exit = 0.652 over 15.6 µm needs **3.0 A cm⁻² absorbed**. At i = 300 mA cm⁻² and η = 0.5, x_in over F = **0.662** | Re-state P(cond_SRI) per slot. Put the L regime (with current modulation, which *is* SRI's no-through-flux geometry) on ≥2 quadrants per membrane. Add per-quadrant loading sensing (4-probe strips along each quadrant, or post-run XRD a-spacing). Match T-H per quadrant regime, not on the rim reading |
| **P3** | **Major** | Quadrant-to-quadrant misfit | Δx = 0.90 − 0.655 → Vegard strain 0.0652·Δx = **1.6 %**. Biaxial stress ≈ E/(1−ν)·ε/2 ≈ **1.6 GPa**, versus 40–150 MPa yield. Excess length gives ~0.4 mm wrinkles at λ ≈ 5 mm. Every P3 modulation cycles it | Model 2-D misfit in M3 (quadrant pattern plus Au-capped cross web). Pre-register AE event *positions* (≥2 AE sensors to localise) and exclude boundary-correlated events from H1a. Or group L-type and flux-type skins on different membranes |
| **P4** | **Major** | 15 µm breaks the ADR-004 ⁶LiF entry-face triton calibration, which was specified for 8–12 µm | 2.73 MeV t through 15.58 µm (loaded): **3 % of accepted tritons emerge**, median 0 MeV (≤0.16 MeV at normal incidence). At 12 µm unloaded they would emerge at 0.87 MeV | One calibration membrane per lot at 10 µm, or a locally thinned (≤10 µm) window with the ⁶LiF, measured by profilometry. Otherwise drop the claim that the tomography is calibrated through the real foil |
| **P5** | **Major** | Electrolyte-borne cosmic recoils are not in M5's budget and unbalance T1 (A − T-H) | Unshielded, PID window per cell: **n–p in the H₂O twin 0.036 /d** (peak window 0.013); D₂O n–d misID 0.006 /d; **D(n,np) breakup protons 0.008 /d**. Compare M5's total of 0.028 /d. The D-specific part is 0.014 /d (≈1.7 counts over 4 cells × 30 d); the H-specific part is 2.6× larger | Primary T1 background = the **same cell current-off and the D₂O/Pt twin** (both carry the D₂O recoils). T-H becomes secondary, with a modelled n–p term. Add these terms to M5 and fit them with the measured fast-n flux |
| **P6** | **Major** | Entry-face efficiency and tomography over the real acceptance | Accepted direction cosines: median 0.68, 10th percentile 0.37. Entry-face p: median **1.41 MeV** at the ΔE, which is the PID lower edge. ε_PID (δ source): 0.220 at the exit face, 0.203 mid-foil, **0.116 at the entry face**, 0.084 at the +2 µm tolerance. Energy does not map one-to-one to depth: z = 4 µm spans 2.34–2.76 MeV | Quote the reach per depth: entry face ≈ 2× worse (6×10⁻⁵–1.5×10⁻⁴ fusions/s). Measure each foil's thickness to ±0.3 µm. Pre-register a response-matrix unfold ε(z, E) instead of reading "energy = depth". An optional ≤30° honeycomb collimator on one quadrant gives real depth resolution |
| **P7** | **Major** | The regime claims for F and L depend on unmeasured surface kinetics | The F "clamp" holds only if gas exchange ≫ J. M3's bare-Pd lower bound (J_sat ≥ 1 mA cm⁻² at 25 °C) allows x_exit anywhere in [0.652, x_in]. L needs J ≈ 0: with 50 nm pinholes each drains ~3×10¹¹ D/s, so the Au cap must have ≤ 3×10⁴ pinholes/cm² to keep J < 1.6 mA cm⁻² | Run the P0b witness **at 0.5 bar D₂** (not vacuum) and include an L-type witness. Specify Au ≥ 50 nm and measure pinhole density (Cu-decoration or H-permeation witness) |
| P8 | Minor | Rn/Tn progeny and ¹⁰B markers land in the ΔE-only triton and ³He windows (secondary observables) | ²¹²Po α plated on the entry face: 0.7 % per decay in the t window (0.85–1.10 MeV) and 0.8 % in the ³He window (0.55–0.85). ²¹⁴Po: 1 % at the −2 µm tolerance. Sparging removes ²²²Rn, not ²²⁰Rn/²¹²Pb from Th in the salt or polymer. ¹⁰B marker ⁷Li (0.84 MeV) → ≈0.4–0.7 MeV at the ΔE (Li stopping estimated): 0.04–0.5 /d, and ×10³–10⁴ with the source boost | Remove the ¹⁰B dot from the telescope view during physics windows, or run the boost only in calibration blocks. Assay LiOD/D₂O for ²²⁸Th/²²⁶Ra. Treat the t/³He windows as non-primary unless they are backed by E-coincidence |
| P9 | Minor | D₂O tritium is unverified; M8 assumes T/D = 1.4×10⁻¹⁵ (≈150 Bq/kg) | Commercial D₂O at 10⁵–10⁷ Bq/kg gives ³He ingrowth of **10⁸–10¹⁰ /day** in 13 g. That compares with the ³He isotope-dilution spike (10⁹–10¹⁰ per aliquot) and dominates any H4 ³He test | Buy certified low-tritium D₂O. Assay every lot by LSC. Size the ³He spike ≥100× the ingrowth, or use ⁴He-pipette calibration |
| P10 | Minor | Internal numeric inconsistencies | x_exit(0.5 bar, 295 K): design 0.63, M3 **0.652**. P3 swing: 0.61↔0.67 in the design, M3 0.629↔0.670. ²¹⁴Po exit energy: design 1.26 MeV; recomputed 1.15 MeV (15 µm), 0.72 MeV (loaded 15.58 µm), 2.46 MeV (at −2 µm). Thickness 15 µm vs ADR-002 rev 2's 8–12 µm. §5's He-window rule (≤10 mA cm⁻²) vs P2, where F/X/ED permeate ~50–150 mA cm⁻²-equivalent | Fix the numbers. Use loaded thickness (×1.039) in all transport. State that the P2 He windows rely on the Pd–Ag path and not the M8 NEG sizing rule |
| P11 | Minor | Local (n,p)/(n,α) sources and blank-telescope blind spots | The Cu septum (≈6 cm² facing Si) and Mo grid give an (n,xp) estimate of ~0.01 /d in the PID window. Pd(n,xp) ~6×10⁻⁴ /d. The blank telescope (under 100 µm Al) sees none of the electrolyte, septum or grid terms. Muon-spallation neutrons from 5 cm Cu are not modelled | Add these to M5. Make one quadrant of B carry an identical septum. Measure fast-n inside the house |
| P12 | Minor | He release for X/ED is far below the quoted 0.85 | M8: a uniform source over δ gives f ≈ λ/δ. For ED (1 µm) with λ = 3–300 nm, **f = 0.003–0.3**. X (140 nm) is 0.02–1 | Quote floors per skin. Rely on quadrant melts for X/ED |

**Checked and negligible** (asked for explicitly):

| Candidate | Result |
|---|---|
| ⁶Li(n,t) tritons from 1 M LiOD reaching the telescope | 2.73 MeV t cannot cross 15.6 µm except from the first ~1 µm at normal incidence. Rate ≲10⁻⁹ s⁻¹, and PID rejects tritons anyway |
| ⁶Li(n,α)t ⁴He in the electrolyte | 5×10⁻⁴–5×10⁻³ /s (45–450 /d), against σ_B = 7×10⁵ /d |
| ¹⁰B marker ⁴He | ~1 /d (10³–10⁴ /d with the boost), ≪ σ_B |
| Radiogenic He in Pd, Au, Mo, Cu | < 10³ /d produced; not released at room temperature |
| n–d recoils in the 0.5 bar front gas | 2×10⁻⁴ /d |
| Energy loss in the D₂ gap | ≤ 6 keV for p. ε_PID changes < 2 % versus vacuum. The design's "gas is harmless for protons" holds |
| α misidentification | ²¹²Po α exits the foil at 3.1–3.4 MeV and stops in the 25 µm ΔE (limit 5.17 MeV), so PID rejects it (M5) |

## Details

### P1 — Polymer He paths (Critical)
- **Where the polymers are.** §3.2 specifies a PCTFE cell body with a 316L lid, and an FFKM "floating" seal between the electrolyte and the membrane rim. The front volume is all-metal, but its only barrier to the cell at the rim is that FFKM seal.
- **What M4 already says.** M4 §5.10 quantifies a single Viton O-ring at 0.8–1.8 mW-equivalent and allows PTFE only as a liner with no path to air.
- **Permeation estimate.** Take PCTFE He permeability as 1–30 Barrer. This is an order-of-magnitude range, but even 1 Barrer is ~100× Pyrex (Pyrex back-derived from M4 is ≈0.01 Barrer). Permeation then gives 10⁸–10⁹ He/s into the headspace, i.e. 0.4–4 mW-equivalent. §10's "0.4–1.6 nW" headspace/melt reach is wrong by ~10⁶ for the headspace part. The melt part is unaffected.
- **The front channel is not protected either.**
  - Any cell-side He inventory (10¹¹–10¹⁴ atoms from polymer outgassing, dissolved air, or incomplete sparging) sets a He partial pressure of 10⁻¹⁰–10⁻⁶ bar at the FFKM.
  - Crossover to the front is then 10⁴–10⁷ /s, versus the 43 /s floor.
  - The 1 % Kr tracer cannot detect this. Fluoroelastomer and fluoropolymer permeabilities for Kr are ≥10× lower than for He, so a clean Kr reading does not bound He.
- **The artifact is temperature-correlated.** Permeability rises by ~×2–3 per 20–30 K, so the P4 70–80 °C phase produces a **He rise that is correlated with temperature and therefore with current and heat**. That is the worst possible artifact for H2.
- **Fix (P1-fix).**
  - PTFE-lined 316L cell.
  - Membrane diffusion-bonded or brazed to a **non-hydriding annulus** (Pt, or Au-plated Cu, ~0.1 mm) that is then sealed with a Au-wire or CF knife-edge. The annulus does not load, so the rigid seal carries no misfit, which meets M3's "no rigid clamp on loading Pd" rule. Radial compliance comes from a thin corrugated annulus.

### P2 — Entry-face loading under the flux quadrants (Critical)
- **Result.** From `m3_common.Material(T = 295 K)`, x_eq(0.5 bar) = 0.652.
  - Φ(0.90) − Φ(0.652) over 15.58 µm gives J = 1.9×10²³ D m⁻² s⁻¹ (≈3 A cm⁻²). M6 §5 found ~1 A cm⁻² for 50 µm, so this is consistent.
  - The cell can absorb at most η·i ≈ 60–150 mA cm⁻².
  - Under a clamped exit the entry face therefore follows x_in = Φ⁻¹(Φ(x_exit) + ηiL) ≈ **0.655–0.662**.
- **The fallback in §15.1 is not what happens.** The design lists "M collapses into L or F" as an open risk. The stronger statement is that F, X and ED are *never* at SRI loading anywhere in the foil.
- **L-type slots reproduce SRI's geometry better than M.** SRI's cathodes are solid rods with in/out flux through one face and no through-flux. L plus the P3 current square wave is exactly that geometry at x ≥ 0.9 (given M2's x_in). L therefore deserves more slots.
- **Lateral coupling.** The first-eigenmode decay length is L/π to 2L/π = 5–10 µm, which confirms that quadrants are independent. It also means a single rim van der Pauw cannot report any quadrant's x.

### P3 — Misfit between quadrants (Major)
- **Magnitude.** L quadrants at 0.90 and F/X/ED at ~0.655 differ by Δx ≈ 0.245 → 1.6 % linear misfit (M3 η = 0.0652). Elastic accommodation would need ~1.6 GPa; annealed Pd yields at 40 MPa and 150 MPa after one transit.
- **Consequences:**
  - plastic flow and dislocation injection concentrated along the 1 mm Au cross web and quadrant edges;
  - out-of-plane wrinkling of the L quadrants (≈0.4 mm amplitude);
  - fatigue cracking under P3/P5 cycling.
- **Why it matters for H1.** Cracks at a skin boundary are the H1b confound, located exactly where the skin comparison is made. They are also leak paths for the cell→front seal.
- **Status of the model.** M3's mechanics are 1-D and uniform-x, so this is unmodelled.

### P4 — ⁶LiF positive control (Major)
- ADR-004 item 2 assumes an 8–12 µm membrane.
- At 15 µm (15.58 µm loaded) the 2.73 MeV t leaves the foil at 0.16–0.30 MeV only near normal incidence. Over the acceptance, 97 % of tritons stop or fall below the 0.1 MeV threshold.
- The foil-integrated depth calibration therefore does not exist for the as-designed membrane.

### P5 — Electrolyte recoils (Major)
- **Method.** MC of cosmic fast neutrons (M5 spectrum, 6.1×10⁻³ cm⁻² s⁻¹ above 1 MeV, unshielded) producing recoils uniformly in 2 mm of liquid above the Ø20 entry face. Recoils are transported through liquid → 15.58 µm PdD₀.₉ → gas, then through M5's `deposit` and `proton_cut`.
- **Cross-sections.** n–d elastic uses Gammel n–p, as in M5. Breakup σ(n,np) uses a threshold at 3.34 MeV, ~0.2 b at 10 MeV and 0.13–0.18 b above that (±50 %).
- **Size.** Multiply the table values by M5's shield factor 0.6 for the house. The H₂O twin's recoil protons (0.022 /d shielded) exceed M5's entire D-cell budget.
- **Consequences:**
  - The equal-time D − H test is biased negative, which weakens limits. The D₂O terms are D-specific and biased positive (≈1–2 counts over the 30-day, 4-cell sum, against a 5σ threshold of ~8).
  - Neither term is visible to the blank telescope.
  - The fix is to change the reference: current-off of the same cell and the D₂O/Pt twin both carry the same recoils.

### P6 — Acceptance smearing (Major)
- **Energy distributions.** See `figs/rt1a_tomography.png`.

  | Proton origin | E at ΔE, 10–90 % (MeV) |
  |---|---|
  | Exit face | 3.007–3.015 |
  | z = 4 µm | 2.34–2.76 |
  | Mid-foil | 1.56–2.51 |
  | Entry face | 0–1.93 (half below the 1.41 MeV PID edge) |

- **Consequences for §3.4:**
  - The table is the normal-incidence limit.
  - The entry-face ε is 0.116, not ~0.22. The M5 reach of 3–8×10⁻⁵ fusions/s applies only to the exit face.
  - ±2 µm thickness tolerance moves the entry-face ε between ≈0.15 and 0.084.

### P7 — Unverified regimes (Major)
- **F.** The clamp at x_eq assumes exchange ≫ J. M3 could not bound the room-temperature exit law to better than 10³–10⁴. If the bare-Pd exit is kinetics-limited near its 1 mA cm⁻² lower bound, F carries little flux and drifts toward L. P0b must therefore be run against 0.5 bar D₂ with a front-pressure-rise flux measurement per quadrant. Masking the other three quadrants is one way.
- **L.** Disk-sink flux per pinhole is 4aDnΔx ≈ 3×10¹¹ D/s (a = 50 nm, D = 10⁻⁶ cm² s⁻¹, Δx = 0.25). Typical 20 nm sputtered Au with 10⁴–10⁶ pinholes/cm² turns L into a 0.5–50 mA cm⁻² flux sector.

### P8–P12
See the table. The numbers come from §1b and §4 of `figs/rt1a_physics.txt`.

## Q6 — The single physics change that most raises P(detect) under ADR-003

**Make the cell and the membrane rim He-tight and polymer-free (P1-fix).**

Under ADR-003, H1 sensitivity is already saturated at ≤1 % of every claimed magnitude, so better silicon buys nothing. H2 carries 0.20 of the weight and is the only quantitative test of the SRI heat–⁴He claim. As designed, that channel's background (10⁸–10⁹ He/s, temperature-correlated) sits at or above 1 % of the SRI claim scaled to one L quadrant (≈2–5×10⁸ He/s after 0.6 release). Its credibility factor is also ≈0 because of the P4 warm-phase correlation.

The fix is cheap, costing only an all-metal cell and a bonded Pt annulus. It moves the headspace and front floors back to M8's 43–100 He/s, which is 10⁶× better.

Runner-up: re-allocate slots so that ≥2 of 4 quadrants per membrane are L-regime under current modulation (P2). That raises P(cond_SRI at a visible location) from 4/16 to ≥8/16 slots.
