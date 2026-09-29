# ADR-007 — Iteration-1 revision C (response to red-team-1B, engineering and safety)

**Status:** accepted · **Date:** 2026-09-29 · **Input:** [`red-team-1B-engineering.md`](../design/red-team-1B-engineering.md) (numbers in `sim/rt1b_checks.py`)
**Amends:** iteration-1 rev B (ADR-005) and the γ station and modulation decisions (ADR-006).

## 1. Findings accepted as fact

1. **The closed cell has no deuterium balance (B1).** Every D atom that enters or crosses the membrane leaves an unpaired O₂. The recombiner then burns headspace D₂:
   - loading one foil to x = 0.9 consumes 69 % of the headspace D₂;
   - FX and H-M cells empty their headspace in 2–31 min;
   - after that, O₂ accumulates while pressure *falls*, which a "ΔP ≥ +10 %" trip cannot see.

   This is the most important engineering error in rev B, and it traces back to R6's stoichiometric-gas assumption.
2. **The overpressure hierarchy is inverted (B2).** The membrane and grid (≈ 0.5 bar) are the weakest boundary, while the only relief sits at 3.5 bar.
3. **The 0.10 mm grid is under-rated (B3).** At open fraction 0.80 it reaches 258–516 MPa in the loss-of-front cases the design invokes.
4. **The P3 FX protocol breaks the mount (B4).** The −100 mA cm⁻² anodic half-cycles pull 1.4× the foil's inventory, over 5,040 α/β crossings.
5. **The interlock bands overlap the operating window (B5).** In addition, "close valves" makes an FX front flood faster.
6. **The Pt annulus shunts the rim gauge (B6).** Loading sensitivity falls from 80 % to 14 %, and 1 K of drift equals 3 % of the whole signal.
7. **A dead-weight diffusion bond on rolled 12 µm foil is not credible (B7).**
8. **The thermal path is missing (B8).** The electrolyte runs 3–8 K above the jacket (11–19 K at 500 mA cm⁻²), and the house needs internal cooling.
9. **He floors, budget and schedule were optimistic (B9, B12, B13).**
   - Treat the floors as targets to be earned in P0a.
   - DFM-8 costs ≈ $0.8–1.5M.
   - Plan 15–20 months to the end of stage 2.

## 2. Decisions

### 2.1 Deuterium balance and pressure control (B1, B5)
- **Closed D₂ recycle loop per FX/H-M cell.** The permeate leaves the front through the Pd–Ag element, passes a thermal mass-flow meter (the ADR-006 flux regressor) and a metal-bellows pump, then a buffer, then a **second Pd–Ag element** into the cell headspace under a pressure controller. Permeated D is returned to the cell and recombines with the anode O₂. The D₂O inventory stays constant.
- **Make-up for H-L cells and P1.** An external D₂ supply, purified through Pd–Ag, tops up the headspace during P1 loading and after any deload. He blanks start only after P1.
- **Two-sided pressure trip:**
  - low at −10 % of fill, high at +10 %;
  - **O₂ inferred continuously** from the D balance (make-up and recycle flow vs coulometric absorption);
  - O₂ measured by HR-QMS sniff every 6 h (aliquot ≤ 0.1 % of headspace, volume logged and included in He accounting);
  - liquid-level sensor (conductance pins).
- **Headspace ⁴He channel on FX and H-M cells is descriptive by default.** It enters T3 only if P0a blanks, taken *with the loop running*, meet ≤ 3× the target floor. This data-quality gate is fixed before P1. The front channel is unaffected.
- **ΔP bands** (cell − front): operate +30 ± 10 mbar, alarm ± 15, trip at +5 and +60 mbar, with cell-temperature feed-forward on the set-point.
- **Trip action.** The galvanostat drops to 0 on FX and H-M (open circuit leaves FX at its operating state, x ≈ 0.65); H-L drops to keep-alive. The Pd–Ag heaters stay on, so the front keeps exhausting. On an m/z 20 or rupture signature the Pd–Ag element is isolated and cooled instead.
- **UPS scope:** galvanostats, Pd–Ag heaters, pressure controllers, PLC, gauges and chiller control. Every trip acts in hardware (R6 #16), with remote paging.

### 2.2 Pressure hierarchy (B2, B3, B14, B17)

| Boundary | Rating | Protection |
|---|---|---|
| Membrane on grid, forward (cell > front) | 1.0 bar at 44 MPa (post-transit yield ~150 MPa) | grid (below) |
| **Grid**: stress-relieved Mo, **0.30 mm**, 1.0 mm A/F hex holes, 0.15 mm webs, carried by the 4 mm cross septum as a deep rib (10 mm spans). Recrystallising heat treatment forbidden | ≈ 160 MPa at 1.0 bar, against a yield of 550–700 (stress-relieved); **proof-tested at 1.5 bar** per grid | forward differential relief |
| **Forward differential relief**, cell → front | opens at +150 mbar | protects the membrane above all |
| Membrane lifted, reverse (front > cell) | 50 mbar ≈ 54 MPa | **reverse-acting burst disk, vacuum-supported, 20 ± 5 mbar, qualified per lot**, front → headspace |
| Front volume | 3.5 bar | front burst disk to a vent volume |
| Cell headspace | 3.5 bar (deflagration of 15 mL at 0.5 bar: 4–5 bar peak) | burst disk; forward relief acts first |

- **P1 and P6a transits** are held at ΔP ≤ +10 mbar (B14).
- **P6a:** the cell is lowered to ≤ 50 mbar *before* the front is pumped.
- **Efficiency cost of the thicker grid.** Grid transmission over the real acceptance falls from 0.70 (0.10 mm) to **0.43** (0.30 mm). Rev B had used a flat 0.80. Per-cell ε drops by ×0.62. This is accepted, because sensitivity is already saturated for all but two claims (ADR-003; 1C-S6). The reach is recomputed in §10 of the design.

### 2.3 FX protocol (B4)
- Anodic excursions are capped at ≤ 10 % of the foil inventory per half-cycle: −10 mA cm⁻², or −100 mA cm⁻² for ≤ 8 s. They remain descriptive.
- An `m3_cycling` simulation showing x_entry ≥ 0.60 throughout is pre-registered before the freeze.
- Front-pressure steps are **balanced ramps ≤ 1 mbar s⁻¹** commanded on both volumes by one controller, with the trip referenced to ΔP. The range is narrowed to **0.2 ↔ 0.6 bar**.

### 2.4 Mount and gauge (B6, B7, B15)
- **Membrane blank Ø25.** The Pd-only rim extends from Ø20 to Ø23 and carries the four van der Pauw contacts at Ø22. The Pd–Pt bond runs from Ø23 to Ø25.
- **R(x) calibration** is done on **bonded** sibling assemblies. The rim is thermostatted to ± 0.05 K, and its temperature is logged.
- **Bond.**
  - Made in a press with ceramic platens at 5–10 MPa, 850 °C, 1 h, **with a sputtered Pt interlayer** (thermocompression).
  - **Qualification before any cell:** ≥ 10 coupons; He leak ≤ 10⁻¹¹ mbar L s⁻¹ before and after 3 full D load/deload cycles; peel and cross-section.
  - Budget 6–8 weeks and $15–25k.
  - Fallback: iteration-1 §13.
- **Pt annulus:** 0.15 mm, ID 23, OD 36, **two convolutions**.
- **Liner.** The PTFE liner is spring-loaded onto the foil with an Au-plated wave spring. Any wetted Pt rim is counted as cathode area in j.

### 2.5 Thermal (B8, B16)
- **Cell and house cooling.** Every cell body gets a water jacket, and a chilled loop runs inside the house (≈ 230 W).
- **Temperature control.** An electrolyte RTD sits in a PTFE-sheathed thermowell, servoed with I·V feed-forward. **"22 ± 1 °C" now means the electrolyte**, as measured by that RTD.
- **Si.** Si mounts are thermostatted separately. Leakage current is a pre-registered nuisance covariate.
- **P4:** fill stays at **0.5 bar**, temperature is capped at **55 °C**, ramp ≤ 5 K h⁻¹. The trip is re-zeroed at temperature, and O₂ is inferred during the re-zero window. Si stays ≤ 30 °C. This replaces ADR-005's 1.0 bar / 60 °C, which the grid and relief hierarchy do not need to carry.

### 2.6 Si protection (B10, B11)
- **Baseline stays Si in 0.5 bar D₂**, with these protections:
  - bias interlocked off between 0.05 and 50 mbar and during ramps > 5 mbar s⁻¹;
  - ΔE–E gap vented (≥ 2 mm² per quadrant);
  - HV current-limited (≤ 1 µA, µs trip);
  - ≥ 2 spare telescopes;
  - one **sacrificial telescope soaked in D₂ from P0a**.
- **Qualified fallback:** a vacuum pocket behind a 1–2 µm window on a grid aligned with the membrane grid. Use SiNₓ or Ni, not Ti, which hydrides in D₂. It is adopted if the soak telescope degrades.

### 2.7 He system (B9, B18)
- **Seals.** Every demountable seal gets a pumped or Ar-guarded interspace. Each cell connects to the manifold through a three-valve, pumped-interspace isolation. The seal count is audited per volume.
- **P0a acceptance.** Floor ≤ 3× target, otherwise the cell's He is descriptive. The **10× floor is the pre-registered planning case**: its reach is still ~10⁻⁷ of the claims.
- **Mass spectrometry.** HR-QMS is primary. The external sector MS is used on a pre-registered subset of ≤ 60 samples.
- **Blanks.** Sparge ≥ 30 min, and take the P0a blanks *after* the fill.

### 2.8 Budget and schedule (B12, B13)
- **Re-baselined, excluding labour:**
  - DFM-8 ≈ **$0.8–1.5M**;
  - γ station ≈ $60–90k;
  - C3-G ≈ $80k;
  - C1 ≈ $70k.
- **Tier 1 ≈ $1.0–1.7M.** Labour is 3–4 FTE for ≈ 1.5 years.
- **Minimum claim-weighted subset ≈ $0.45–0.6M.**
- **Schedule:** 15–20 months to the end of stage 2.
- **Sequencing.** P0b (witness membranes, exit law, τ_eff, thermal and V_cell at 500 mA cm⁻²) runs on a **separate bench before stage-1 assembly**. If it is late, stage-1 A3 is built as H-L.

## 3. Not adopted or modified
| Finding | Disposition |
|---|---|
| B3 option "ribs ≥ 1 mm deep at 5 mm pitch under a 0.10 mm sheet" | Rejected. 5 mm cells still give 258–516 MPa at 0.5–1.0 bar for 0.10 mm, and deep ribs add shadow. The 0.30 mm sheet on the septum rib is simpler and has margin |
| B10 fallback "2 µm Ti window" | Modified to SiNₓ or Ni, because Ti hydrides in D₂ |
| B1 "headspace channel descriptive on FX/H-M" | Modified to "descriptive unless P0a blanks with the loop running pass ≤ 3× target" |

## 4. Consequences
- **Hardware** grows by one recycle loop per FX/H-M cell (4 of 8 active cells, plus T-H(FX)) and one make-up line per cell.
- **Detection efficiency** falls by ×0.62. T1 pooled reach becomes ≈ 2.0–4.2×10⁻⁵ fusions s⁻¹ per membrane (see design §10).
- **The ADR-003 ranking does not change.** C3-G and C1 now look even better per dollar next to a $0.8–1.5M DFM-8. The minimum subset (C3-G ×5 + C1 ×2 + twin + one 2-membrane DFM stage) is the rational first purchase.
