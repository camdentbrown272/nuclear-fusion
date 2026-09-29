# Red-team 1B — Engineering, buildability and safety review of iteration-1 rev B

**Scope:** `docs/design/iteration-1.md` (rev B), `docs/decisions/ADR-005-iteration1-revB.md`, with the M2/M3/M4/M5/M8 models and R6.
**Numbers:** `sim/rt1b_checks.py` → `docs/design/figs/rt1b_checks.txt`, `rt1b_checks.png`. Values tagged [est.] are engineering estimates that were not taken from a source.
**Date:** 2026-09-29

## Summary verdict

**Rev B cannot be built and run as drawn.**

**The fatal problem is the gas balance of the closed cell, which no document has worked out.**
- Every deuterium atom that enters the membrane leaves an excess O₂ at the anode. The recombiner then burns headspace D₂ to remove that O₂.
- **Loading.** Bringing one 12 µm foil to x = 0.9 burns 69 % of the 306 µmol headspace D₂. Cell pressure falls from 500 to ≈155 mbar. The reverse-relief foil (20 mbar) bursts during P1 in every cell.
- **FX and H-M cells.** Permeation empties the headspace D₂ in 2–31 min. After that, O₂ builds up **while the pressure falls**, which the R6 "ΔP ≥ +10 %" trip cannot see. The same cells lose 12–177 g of D₂O over P2, from a 13 g inventory.
- **No make-up path exists.** The design has no He-free D₂ make-up path to the cell. Even an ideal Pd-purified make-up supply (He fraction 10⁻¹⁵) puts 10²–10³ He s⁻¹ into the headspace channel, above its 43–100 He s⁻¹ floor.

**The pressure-protection hierarchy is inverted.**
- The cell burst disk (3.5 bar) sits far above what the membrane and grid can take (≈0.5 bar ΔP).
- The 0.10 mm Mo grid reaches 258–516 MPa in the specified loss-of-front cases, and 1–2 GPa if its "ribs" are only in-plane.
- The P3 FX protocol breaks the membrane in two ways:
  - it pulls 1.4× the foil's D inventory per anodic half-cycle, over 5,040 cycles;
  - its 0.2↔1.0 bar front steps conflict with the ±40 mbar interlock.

**Secondary but real:**
- The bonded Pt annulus shunts the rim van der Pauw gauge. Loading sensitivity drops from 80 % to 14 %, so 1 K of drift is worth 3 % of the whole loading signal.
- "22 ± 1 °C" can only be a jacket set-point. The electrolyte runs 3–8 K above it at 200–300 mA cm⁻² and 11–19 K above at 500 mA cm⁻².
- The house needs internal water cooling: +64 K otherwise.
- The He floors are targets, not credible design values, until P0a blanks show them.
- The budget omits or underprices valves, gauges, sector-MS samples, custom Si and ³He tubes. A realistic DFM-8 costs ≈ $0.8–1.5M, not $0.415M. The 8-month schedule leaves out 6–10 months of procurement and development.

**Every Critical item below has a concrete fix.** None requires leaving the COLD envelope.

## Ranked findings

| ID | Sev | Finding | Quantitative evidence | Fix |
|---|---|---|---|---|
| **B1** | **Critical** | Closed-cell D balance. Loading and permeation consume headspace D₂ through the recombiner. There is no make-up path, O₂ enriches at falling pressure, and D₂O runs out | Loading to x = 0.9: 422 µmol D ⇒ 211 µmol D₂ burned = 69 % of the headspace; P_cell 500 → 155 mbar (ΔP −345 mbar vs foil −20). H-M 10 / FX 50 / FX 150 mA cm⁻²: headspace D₂ gone in 31 / 6.3 / 2.1 min, dP/dt −0.4 / −2 / −6 mbar s⁻¹, then O₂ at 0.08–1.2 µmol s⁻¹. D₂O loss 0.28 / 1.4 / 4.2 g d⁻¹ vs 13 g inventory | (1) Active cell-side pressure control with **He-free D₂ make-up**: recycle the front's Pd–Ag exhaust into the headspace through a metal-bellows pump, or use a dedicated electrolytic D₂ generator permeating through Pd. (2) A **two-sided** trip: low pressure at −10 % of fill plus a direct O₂ measurement (QMS m/z 32/36 on a headspace sniff, or a sealed electrochemical O₂ cell). (3) Liquid-level sensor. (4) Declare the headspace ⁴He channel **descriptive only on FX/H-M**: make-up adds ≥ 10²–10³ He s⁻¹ even at a 10⁻¹⁵ He fraction. The front channel is unaffected |
| **B2** | **Critical** | Overpressure hierarchy inverted: the membrane and grid are the weakest boundary, but the only cell relief is 3.5 bar | Loss-of-front at 0.5 bar is already 28 MPa (membrane) and 258 MPa (grid). An isochoric deflagration of the headspace reaches 8–10 × P₀ = 4–5 bar (R6 §3), far beyond both. The ruptured membrane then passes gas and electrolyte to biased Si and a 350 °C Pd–Ag element, which acts as an H₂/O₂ catalyst | Add **forward relief** from the cell to a vent line at ΔP ≤ +150 mbar (sized from the grid rating after B3). Add a front burst disk to vent. A Δp or m/z 20 trip must **isolate and de-energise the Pd–Ag element and cut cell current**. Rate each boundary in a written pressure-hierarchy table |
| **B3** | **Critical** | The Mo grid (0.10 mm, open fraction 0.80) is ~2× undersized for the loss-of-front cases the design itself invokes (P6a pumps the front for 2 d) | Ligament efficiency 0.09. σ = 15 / 258 / 516 MPa at +30 mbar / 0.5 bar / 1.0 bar with 5 mm cells, or 62 / 1031 / 2062 MPa if the ribs are in-plane only (10 mm cells on the septum). Mo yield: 550–700 MPa stress-relieved, 350–450 MPa recrystallised; DBTT ≈ RT [est.] | Use the M3 spec: a 0.3 mm sheet (40 MPa at 1 atm), or ribs **≥ 1 mm deep** under the existing septum shadow and at the 5 mm pitch. Proof-test each grid at 1.5 bar. Forbid recrystallising heat treatment of the grid |
| **B4** | **Critical** | The P3 FX protocol is mechanically incompatible with the mount | (a) Front steps 0.2↔1.0 bar with the cell at 0.5 bar give 500 mbar reverse ΔP on a lifted membrane: 251 MPa > UTS 170. If balanced, both volumes must track to ±20 mbar over an 800 mbar step (2.5 %) while B1 is also acting. (b) The −100 mA cm⁻² anodic half-cycle (120 s) pulls 12 C cm⁻², 1.4× the foil inventory at x = 0.65. The entry face therefore crosses α/β on each of 5,040 cycles, against an M3 perforation life of 12–180 transits | (a) Specify **balanced** ramps ≤ 1 mbar s⁻¹ run by one controller that commands both volumes, with the trip referenced to ΔP. (b) Cap the anodic charge at ≤ 10 % of inventory per half-cycle (e.g. −100 mA cm⁻² for ≤ 8 s, or −10 mA cm⁻²). Pre-register an `m3_cycling` run showing x_entry ≥ 0.60 throughout |
| **B5** | **Critical** | The Δp interlock logic cannot keep FX/H-M cells intact through a routine fault. Its bands overlap the operating window | Operating band +10…+50 mbar vs trip ±40. Thermal: 3.1 mbar K⁻¹ (gas 1.7 + D₂O vapour 1.4), so a 500 mA cm⁻² step moves ΔP ≈ 31 mbar. With the Pd–Ag exhaust stopped (heater or mains loss), the front rises 2.0 / 0.67 / 0.27 / 0.13 mbar s⁻¹ (FX 150 / FX 50 / keep-alive 20 / H-M), and the relief foil bursts in **25 / 75 / 188 / 376 s**. The trip action, "close valves", makes an FX front rise faster | Trip action for all cells: **galvanostat to 0 on FX/H-M**; open circuit leaves FX at x ≈ 0.65, its operating state, and H-L drops to keep-alive. Put the Pd–Ag heaters, the pressure controllers and the PLC on the UPS. Re-derive the bands: operate +30 ± 10, alarm ±15, trip at +5 / +60 mbar. Add cell-temperature feed-forward to the ΔP set-point |
| B6 | Major | The rim van der Pauw gauge sits over the bonded Pt annulus, which carries most of the current | 2-D finite differences (contacts at Ø21): ρ_Pd × 1.8 changes R_vdP by **14 %** with the annulus vs 80 % for Pd alone. The Pt TCR (0.39 % K⁻¹) makes 1 K equal to 3 % of the full loading signal. Near x = 0.9, Δx = 0.01 ≈ 0.5 K. The "x ≥ 0.90 in-regime hours" gate is therefore unmeasurable | Keep the Pd–Pt bond **outside** a Pd-only rim of ≥ 1.5 mm and put the contacts there, or add 4-wire Pd tabs that extend beyond the annulus. Calibrate R(x) on **bonded** sibling assemblies. Thermostat the rim to ±0.05 K and log it |
| B7 | Major | A diffusion bond at 850 °C/1 h under a "light dead-weight" is not credible as a He-tight, fatigue-resistant joint on 12 µm foil | Pd self-diffusion at 1123 K: D ≈ 9×10⁻¹⁸ m² s⁻¹, so √(Dt) ≈ 0.2 µm in 1 h. That is below the asperity heights of rolled foil (0.1–1 µm). Diffusion bonding needs ~5–20 MPa contact pressure [est.]; a dead-weight gives kPa. The bond then sees ≈ 3.3 % hydride misfit (≫ yield) once per full load/deload, and the bond edge is a stress singularity. The front gas (0.5 bar) loads the rim to x ≈ 0.65 whatever the PTFE does | Bond in a press with ceramic platens at 5–10 MPa, or add a sputtered Pt/Au interlayer and thermocompression. **Qualification before any cell:** ≥ 10 coupons; He leak ≤ 10⁻¹¹ mbar L s⁻¹ before and after 3 full D load/deload cycles; peel and cross-section. Budget 6–8 weeks and $15–25k. Keep the fallback (iteration-1 §13) ready |
| B8 | Major | Thermal: 22 ± 1 °C is achievable only as a jacket set-point, and the heat path is not in the design. The house has no heat sink | G_liner = 0.60 W K⁻¹ (24 cm² of 1 mm PTFE). ΔT = 3–8 K at 200–300 mA cm⁻² and **11–19 K at 500 mA cm⁻²** (V_cell 4.1–7.3 V with 27–60 % void; Q 6.5–11 W). τ ≈ 100 s. An unjacketed body (free convection, 0.10 W K⁻¹) would run +36–80 K. The house carries ≈ 230 W inside 20 cm of HDPE (3.6 W K⁻¹), so +64 K. The P3 steps put a current-locked temperature (and Si leakage) swing onto the lock-in channel | Water-jacket each body. Put an electrolyte RTD in a PTFE-sheathed thermowell and servo with I·V feed-forward. Run a chilled loop inside the house. Thermostat the Si mounts separately and log leakage current as a pre-registered nuisance covariate. Pre-register the electrolyte-temperature definition. For Seebeck cells, accept a floating electrolyte temperature and report it |
| B9 | Major | The He floor (0.16–0.39 nW = 43–100 He s⁻¹) depends on leak and valve performance that the build cannot verify before P0a | Air in-leak: 10⁻¹² std-He total gives 129 He s⁻¹. Keeping the leak ≤ 10 % of the floor needs ≤ 3×10⁻¹⁴ mbar L s⁻¹ over all seals, below He-leak-detector reach (Ar route only). One closed valve seat facing a manifold at 10⁻⁶ mbar He after a spike gives ≈ 2,500 He s⁻¹. The headspace side alone carries lid, anode feedthrough, P, T, burst disk, sample valve and relief foil. The ≤ 10-seal budget is already spent | **Pumped (or Ar-guarded) interspace on every demountable seal.** Put a 3-valve pumped-interspace isolation between the manifold and each cell. Do a seal-count audit per volume. Set a P0a acceptance rule: floor ≤ 3× target, otherwise the cell runs as descriptive. Pre-register the fallback reach at a 10× floor, which is still ~10⁻⁷ of the claims |
| B10 | Major | Si in 0.5 bar D₂: margins are thin and there are no long-term data | Paschen minimum ≈ 300–330 V at pd ≈ 2.5 Torr cm, so 150 V gives a ~2× margin, eroded by bond-wire edge fields and by Kr (Penning). P6a pump-down crosses the minimum. The 25 µm ΔE as a diaphragm (a = 12 mm) reaches 13 / 61 / 178 MPa at 10 / 100 / 500 mbar across it. ΔE capacitance is 624 pF per quadrant, so 7–11 keV rms noise (acceptable). At 60 °C (P4), leakage rises ~30× | Interlock bias off for 0.05–50 mbar and during any ramp > 5 mbar s⁻¹. Vent the ΔE–E gap (≥ 2 mm² per quadrant). Soak a sacrificial telescope in D₂ from P0a onward. **Fallback:** Si in a separate vacuum pocket behind a **2 µm Ti window on an aligned twin grid** (≈ 92 MPa at 0.5 bar vs Ti yield ≥ 300; ΔE_p ≈ 60–90 keV). This also protects the Si from B11. The vacuum front for H-L stays as the second option |
| B11 | Major | Membrane rupture sends LiOD onto biased Si faster than any interlock can act | A 1 mm hole at 34 mbar (30 + 4 mm liquid head) passes ≈ 2 mL s⁻¹ and reaches the ΔE in < 1 s. LiOD etches Si, shorts the HV and destroys the telescope. Headspace gas (O₂ ≤ 3 %, more in the B1 case) reaches the hot Pd–Ag | Current-limit the HV (≤ 1 µA, µs trip). Budget ≥ 2 spare telescopes. Adopt the window fallback in B10. Isolate and cool the Pd–Ag on an m/z 20 or Δp trip |
| B12 | Major | Budget is missing or mispriced | See §6: all-metal valves ≈ 60 × $1.5–3k; 24 capacitance gauges; external sector-MS ≈ 560 samples × $0.5–2k; custom epoxy-free quadrant telescopes $10–15k each plus NRE; new ³He tubes $3–8k each; ≈ 60 digitiser and preamp channels; chillers and jackets | Re-baseline DFM-8 at ≈ $0.8–1.5M and Tier 1 at ≈ $1.0–1.7M, excluding labour (3–4 FTE × 1.5 y). Keep HR-QMS as primary and use the sector MS only on a pre-registered subset (≤ 60 samples) |
| B13 | Major | Schedule: "two stages ≈ 8 months" is run time only | Custom Si lead time 16–26 wk; bond development 6–8 wk (B7); He-system commissioning; stage-2 rebuild needs bonding + bake + ≥ 3 × 1 d + 2 × 7 d blanks = ≥ 17 d per cell. P0b fixes the Ni thickness **after** H-M (stage-1 A3) must already be installed for P0a blanks | Plan 15–20 months to the end of stage 2. Run P0b on a separate bench **before** stage-1 assembly, or build stage-1 A3 as H-L |
| B14 | Minor | Transformation plasticity when α↔β transits happen under load | ε_tp ≈ 1.2 % (P1 at +30 mbar); 1.6–5.8 % (P6a at 0.5 bar, hardened/annealed); 2.4 % (after P4 at 1.0 bar). The foil dimples into the 1 mm holes and thins at hole edges | Hold ΔP ≤ +10 mbar during P1 and P6a transits. In P6a, lower the cell to ≤ 50 mbar before pumping the front |
| B15 | Minor | Pt corrugation and the PTFE-liner edge | Guided-leg strain 0.75 % at 0.45 mm travel, far above annealed-Pt yield, so it ratchets on each full cycle (few cycles, acceptable). As a diaphragm it reaches 22 / 45 MPa at 0.5 / 1.0 bar, against a yield of 35–50. The PTFE-vs-316L CTE gives ≈ 0.16 mm liner movement over the 38 K of P4: the liner lifts off the foil, electrolyte creeps to the rim and bond, and the Pt becomes cathode area | Use 0.15 mm Pt, or two convolutions. Spring-load the liner onto the foil (Au-plated wave spring). Accept the Pt rim as cathode area and correct j |
| B16 | Minor | P4 at 60 °C | D₂O vapour ≈ 0.17 bar plus 13 % gas expansion trips the "ΔP ≥ 10 % of fill" rule. The plateau at 60 °C is 0.20 bar, so FX drops to x ≈ 0.59 | Re-zero the trip at temperature with a pre-registered ramp (≤ 5 K h⁻¹) and an O₂ sensor (B1) covering the re-zero window. Keep the Si at ≤ 30 °C |
| B17 | Minor | Burst foil at 20 mbar reverse | Commercial low-pressure disks have ±25–50 % tolerance. A foil that bursts at 30 mbar gives 38 MPa on the lifted membrane, ≈ the annealed yield. The foil sees +30…+50 mbar forward continuously | Reverse-acting disk with vacuum support, tolerance ≤ ±5 mbar, qualified by lot. Log its He contribution as a seal (B9) |
| B18 | Minor | He from the electrolyte and PTFE | Air-saturated D₂O holds ≈ 4.6×10⁻⁸ cm³STP g⁻¹, so ≈ 1.6×10¹³ He in 13 g. This must be stripped by ~10⁻⁷ before the first window. PTFE He diffusion time over 1 mm is ~10³ s, so the 150 °C bake empties it. Neither drives the floor if blanks are taken after the fill | Sparge ≥ 30 min, then take the P0a blanks **after** the fill (not on a dry cell) |

## FMEA (top 10)

S = severity, O = occurrence, D = detection difficulty (1–10); RPN = S·O·D.

| # | Failure mode | Cause | Effect | S | O | D | RPN | Mitigation (see finding) |
|---|---|---|---|---|---|---|---|---|
| 1 | O₂-rich headspace at falling pressure | Permeation or loading consumes D₂ via the recombiner; no make-up (B1) | Flammable mixture in a 15 mL headspace above a 12 µm wall; deflagration ruptures the membrane | 9 | 9 | 7 | **567** | D₂ make-up, two-sided trip, O₂ measurement |
| 2 | Relief foil burst or membrane lift in P1 | Cell pressure falls 345 mbar during loading (B1) | Loss of both He channels for the cell for the stage | 7 | 10 | 3 | 210 | Make-up and ΔP servo during P1 |
| 3 | FX cycling fatigue and perforation | Anodic half-cycle > inventory; 5,040 cycles (B4) | Pinholes, electrolyte to Si, H1b cracking confound | 7 | 7 | 5 | 245 | Cap the anodic charge; simulate before freeze |
| 4 | Spurious ⁴He via valve seats, spikes, leaks | Manifold spikes, 10⁻¹² class seals (B9) | False T3 positive (scientific severity) | 8 | 5 | 5 | 200 | Pumped interspaces, blind spikes, Ar/Kr tracers |
| 5 | Headspace over-pressure or deflagration | Recombiner flooding, hot spot, or B1 | Membrane and grid fail long before the 3.5 bar disk; hot gas and LiOD to Si and Pd–Ag (B2) | 10 | 3 | 6 | 180 | Forward relief ≤ +150 mbar; Pd–Ag isolation |
| 6 | Grid yield on loss-of-front or P6a | 0.10 mm grid, 258–2062 MPa (B3) | Membrane lift, tear, rupture | 7 | 5 | 5 | 175 | 0.3 mm sheet or deep ribs; proof test |
| 7 | FX front flood on power or heater loss | Pd–Ag exhaust stops; foil bursts in 25–75 s (B5) | Stage loss for the cell; D₂ into the headspace | 7 | 6 | 4 | 168 | UPS on the Pd–Ag and controllers; current-off on trip |
| 8 | Membrane rupture onto biased Si | Bond voids, pinholes, fatigue (B7, B11) | Telescope destroyed, HV short, LiOD in front | 7 | 5 | 4 | 140 | HV current limit, window fallback, spares |
| 9 | Unattended power loss over weeks | UPS sized for galvanostats only | Loss of chillers and heaters, ΔP excursions, dry-out | 6 | 6 | 4 | 144 | UPS scope: Pd–Ag, PLC, gauges, chiller controls; hardwired trips; remote paging |
| 10 | Thermal excursion at 500 mA cm⁻² | Jacket or chiller failure; void-fraction V_cell rise (B8) | +11–19 K with a jacket, +36–80 K without; D₂O dry-out; recombiner > 180 °C | 8 | 3 | 3 | 72 | Hardware over-temperature trip at +5 K; limit excursions to jacketed cells |

Below the top 10:
- ³He tubes at 1–1.5 kV in a closed HDPE house carrying D₂ lines. Without make-up the house holds ≈0.2 L STP of D₂ in a 40 L cavity (≤ 0.5 %). With make-up there is a supply line: fit an excess-flow valve and restrictor outside the house plus an H₂ sensor inside.
- LiOD handling (caustic). Fill under Ar in a glove bag, with a face shield.
- Loaded-Pd stored energy (0.06 kJ) is negligible.

## Details

### 1. Membrane and mount (B3, B4, B6, B7, B14, B15, B17)

- **Membrane stresses are internally consistent.** The ADR-005 figures (28 / 44 / 54 MPa) scale correctly from M3.
- **The problems sit in what carries those stresses.**
  - **Grid.** M3 recommended a 0.3 mm sheet on ribs. Rev B shrank it to 0.10 mm and raised the open fraction to 0.80. Plate stress scales as 1/(t²η), so the grid moved from 40 MPa to 516 MPa at 1 atm.
  - **Ribs.** "0.5 ribs at 5 pitch" is only meaningful if the ribs are deep. As drawn they are wider webs in the same sheet, which leaves the cross septum as the only support.
- **Bond.** The Pd–Pt bond is elegant on paper: no hydride forms on the Pt side. But:
  - a dead-weight cannot close asperities;
  - the Pd *on* the bond still hydrides from the front gas;
  - the bonded Pd/Pt bilayer is plastically cycled at every full load/deload.
- **Low-cycle fatigue.** The few full cycles in the protocol (P1, P6a, possibly P4) are acceptable. The 5,040 FX anodic cycles are not, if they cross α/β at the entry face (§E of the output).
- **Pinholes and leaks.** Pinholes in 12 µm rolled foil, and bond voids, are cell→front leak paths. Kr tracing and the m/z 20 trip detect them; nothing prevents them. Hence the per-foil He leak test before bonding and the bond qualification in B7.

### 2. Si telescope in D₂ (B10, B11)

- **Discharge margin.** 150 V bias sits ~2× below the Paschen minimum. That margin is adequate at 0.5 bar steady state. It is not adequate during P6a pump-down or with Penning admixtures, so bias must be interlocked on pressure.
- **Hydrogen effects.** Molecular D₂ at RT does not dissociate on SiO₂/Al, so B–H passivation and interface-trap changes should be small. There are, however, **no months-long data** for PIPS detectors in D₂. A sacrificial soak telescope costs one device and removes the unknown.
- **Recommended fallback.** The thin-window vacuum pocket (B10) beats "vacuum front for H-L":
  - it keeps the front D₂ ambient for all membrane types;
  - it isolates the Si from rupture;
  - it avoids 28 MPa sustained on the membrane.
  The cost is ≈60–90 keV of proton energy and a co-aligned second grid.

### 3. Cell (B1, B8, B15, B16, B18)

- **Geometry.** The 12 mL electrolyte forms a 38 mm column in a 20 mm bore. Liner conduction (0.60 W K⁻¹) is the heat path, and the electrolyte time constant is ≈100 s.
- **Void fraction in the 4 mm gap.** At 0.5 bar the gas volume doubles relative to 1 bar. Void fraction is 0.16–0.60 at 300 mA cm⁻², depending on bubble size (50–100 µm, Stokes rise). Cell voltage is 3.1–5.2 V at 300 mA cm⁻² and 4.1–7.3 V at 500 mA cm⁻².
  - SELV holds.
  - The "≤ 3.2 V at 200 mA cm⁻²" line holds only for ≥ 100 µm bubbles.
  - Measure V_cell on the P0b witness.
- **Recombiner.** Load is 0.96 / 1.44 / 2.40 W at 200 / 300 / 500 mA cm⁻². The headspace turns over every ≈40 s at 300 mA cm⁻². A 1 % recombiner lag gives +0.12 mbar s⁻¹, so the ΔP servo in B1/B5 is needed even for H-L cells.
- **Carbonate.** Negligible once sealed. Fill under Ar.

### 4. He system (B1, B9, B18)

- **The Pd–Ag element.** It is He-tight as a solid. Its risks are pinholes, braze joints and the upstream He content. With the exhaust facing vacuum, a pinhole leaks outward, which is harmless. The fill direction needs a Pd–Ag-purified source.
- **The make-up flow B1 requires** (0.3–4.7 L STP d⁻¹) makes the headspace channel of FX/H-M cells unusable at the stated floor, even at a 10⁻¹⁵ He fraction.
- **Shared manifold.** It needs pumped-interspace isolation per cell, as M8 already specifies for pipettes.
- **Floor credibility.** The M8 floor (σ_B ≈ 7×10⁵ atoms d⁻¹) is about 10× below typical cold line blanks reported by noble-gas labs (~10⁻¹² cm³STP ≈ 3×10⁷ atoms per extraction) [est.]. Treat 0.16–0.39 nW as a target to be earned in P0a. Pre-register the 10× floor as the planning case.
- **Bake.** With Si inside, bakes are limited to ≤ 100 °C [est.], as M8 notes. This does not limit He, but it leaves more H₂O/D₂ load for the getters.

### 5. Safety

The FMEA is above.

- **The single most important safety change is B1.** The rev-B trip philosophy, inherited from R6, assumes stoichiometric gas production. A permeating cathode breaks that assumption and turns pressure *loss* into the hazard signal.
- **B2 inverts the relief hierarchy.** The membrane is the de facto burst disk, and it bursts toward biased electronics and a hot catalyst.
- **P6a (front pumped).** It is the designed loss-of-front case. Do it only after B3 is fixed, and with the cell brought down to ≤ 50 mbar (B14).
- **P4 (1.0 bar, 60 °C).** It needs the trip re-zeroed and O₂ monitoring (B16).
- **UPS scope must include:**
  - Pd–Ag heaters;
  - pressure controllers;
  - the PLC;
  - chiller control.
- **Hardwired trips.** Every trip must act in hardware, per R6 #16.

### 6. Budget and schedule (B12, B13)

| Item | Rev-B line | Realistic [est.] |
|---|---|---|
| All-metal valves (≈ 6 per cell × 8 + manifold ≈ 60) | inside $30k "gas and He" | $90–180k |
| Capacitance and Δp gauges (3 per cell × 8) | same | $50–70k |
| Pd–Ag elements (9), NEGs, ion/turbo pumps for the manifold | same | $30–50k |
| Make-up/recycle loop (B1), per cell | none | $5–10k × 8 |
| Si telescopes, custom epoxy-free quadrant ΔE + E, NRE, 2 spares | $48k | $110–160k |
| Digitisers and low-noise preamps (≈ 60 channels) | $35k | $80–120k |
| ³He bank, 24 new 4-atm tubes plus electronics | $65k | $100–200k |
| Jackets, chillers, house cooling | none | $25–45k |
| External sector-MS (7 He cells × 2 channels × ~20 windows × 2 stages ≈ 560 samples) | inside $45k | $0.28–1.1M if all go to sector MS; ≈ $60k if limited to a pre-registered subset |
| Bond development and qualification (B7) | none | $15–25k |
| **DFM-8** | **$415k** | **≈ $0.8–1.5M** |

- **Labour** is excluded throughout, but it is the dominant cost: 3–4 FTE for ≈1.5 years.
- **The $0.28M minimum subset** carries the same valve, gauge and He-analysis underpricing. It is realistic at ≈ $0.45–0.6M.
- **Schedule.** The 8 months counts run time only. Procurement, bond development and He commissioning add 6–10 months, so plan **15–20 months** to the end of stage 2.
- **P0b sequencing** must move before stage-1 assembly (B13).

### 7. What was checked and found adequate

- **Stored energy.** 0.06 kJ in the membrane and ≤ 50 J in a stoichiometric 15 mL at 0.5 bar.
- **SELV** at ≤ 7.3 V.
- **Pd–Ag exhaust area** for FX 150: ≈ 1.2 cm², which is feasible.
- **Recombiner heat** through the lid.
- **ΔE noise** of 7–11 keV rms against a 0.3–0.6 MeV proton ΔE deposit.
- **Corrugation fatigue** under H-cell P3 steps (Δx ≈ 0.03 ⇒ ≈ 15 µm radial, elastic).
- **Loss-of-front membrane stresses** at 0.5 / 1.0 bar (28 / 44 MPa, once hardened by the first transit).
