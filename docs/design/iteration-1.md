# Iteration 1 — Detector-Facing Membrane Array ("DFM-4")

**Status:** iteration-1 specification (pre-red-team-1) · **Date:** 2026-09-29
**Decision basis:** ADR-002 rev 2, ADR-003, ADR-004; models M0–M8; research R1–R7; red-team-0.
**Lead designer:** Claude. Every parameter below cites its source in §14.

---

## 0. The configuration in one paragraph

The configuration has seven parts:

1. **Membrane.** A **15 µm annealed palladium membrane, 20 mm active diameter**, is the wall between two chambers. Both chambers hold D₂ at 0.5 bar, so the membrane carries almost no pressure load.
2. **Entry face (top).** This face is the cathode of a small **sealed electrolytic cell** (1.0 M LiOD in D₂O, Pt mesh anode 4 mm above). The cell drives the face to D/Pd ≈ 0.9 and pushes deuterium through the foil.
3. **Exit face (bottom).** This face is divided into **four quadrants with different nanometre skins**:
   - a 20–50 nm Au cap (high-loading, no-flux regime);
   - a thin Ni barrier (high loading *and* flux);
   - bare Pd (maximum flux);
   - a Pd/CaO multilayer (oxide-interface regime).
4. **Charged-particle detector.** A **quadrant-matched silicon ΔE–E telescope** sits 4 mm below the exit face, with a thin cross-shaped septum so that each quadrant has its own detector. Because the foil is thinner than the range of a 3 MeV proton, the telescope sees reactions **anywhere in the foil**, including the electrolyte-facing surface where the classic high-loading claims were made. The proton energy tells which depth a reaction came from.
5. **Helium-4.** Each cell has its own sealed front volume (≈30 cm³) for **static ⁴He accumulation**. The only D₂ path into or out of the volume is a He-tight Pd–Ag element, and gas is analysed as getter-cleaned aliquots. After the run each quadrant is melted to recover retained He.
6. **Array.** Six such cells form the array: four active D₂O cells, one light-water twin, and one flux-blocked D₂O/Pt twin. They share a low-background house with a ³He neutron bank, a 511 keV coincidence pair, a muon veto and a blank telescope.
7. **Drive.** The primary variable is **deuterium flux**. It is modulated two independent ways: by cell current, and at constant current by front-pressure steps.

![DFM cell cross-section](figs/iter1_dfm_cell.png)

## 1. What this configuration is optimised for, and the honest odds

- **Objective (ADR-003).** Maximise Σ over claims c of [weight × P(claimed conditions reproduced where a detector can see) × P(≥1 of N samples active) × min(1, P(detect ≤ 1 % of the claimed magnitude)) × credibility].
- **Why a membrane array.** It is the only configuration found that puts all of the following in one instrument:
  - *every* graded claim regime at a detector-visible location: SRI/ENEA loading at the entry face; Iwamura/Lipson/Czerski flux regimes at the exit face;
  - ≥ 4 independent samples;
  - three orthogonal nuclear channels (charged-particle spectroscopy, ⁴He, 511 keV) with internal controls.
- **Physics prior (M0 rev 2).** Standard physics falls short of a detectable rate by **89–115 e-folds**. Static screening theory gives U_e,eff ≈ 9–13 eV at every lattice site (M1). Prior bulk nulls already exclude accelerator-scale screening acting statically in bulk, vacancy, grain-boundary and surface sites (M1).
- **So a positive result needs new physics.** My estimate of P(any pre-registered 5σ positive that survives all controls) is **≲ 3 %**. The value of the design is that **a null is informative for the first time**. It would set in-regime upper limits on the SRI heat–⁴He, Iwamura interface, Lipson, and Czerski beam-free-511 claims, at ≤ 10⁻⁵–10⁻⁷ of their claimed magnitudes (§10).

## 2. System architecture (tiers)

| Tier | Element | Purpose | Status |
|---|---|---|---|
| **1 (core)** | **DFM array**: 4 active D₂O cells + H₂O twin + D₂O/Pt flux-blocked twin + blank telescope, in one shielded house | H1a/H1c/H2/H3 with internal controls; the geometric answer | specified here |
| 2 | **C3-G gas-entry cells** (Iwamura variant, D₂ + H₂ control), in the same house | Iwamura transmutation (grade B, contested) at its claimed conditions, with isotope-tagged targets | §7.1 |
| 2 | **C1 SRI/ENEA replica twin** (Ø1 mm × 30 mm wire, Seebeck twin calorimeter, headspace ⁴He) | Highest-fidelity reproduction of the heat–⁴He claim | §7.2 |
| 3 | Co-deposition on PSD scintillator (open cell) | SPAWAR conditions at ~2π acceptance, cheap | §8.1 |
| 3 | Sealed CNZ/PNZ nanocomposite ⁴He furnace (200–300 °C; Ni/Cu film variant to 900 °C) | Kitamura / Clean Planet heat claims with ash accounting | §8.2 |

Tier 1 is the deliverable geometry. Tiers 2–3 raise Σ_c at modest cost and are specified at lower detail.

## 3. The DFM cell (Tier 1) — geometry

### 3.1 Cross-section (not to scale; dimensions in mm unless noted)

```
                 Pt/Ir feedthrough (anode, +)       pressure transducer, burst disk (3.5 bar), fill/He-sample valve
                          |                                   |
        +=================|===================================|=================+   316L lid, Cu-gasket CF seal
        |   headspace <=15 mL, D2 0.50 bar abs, recombiner (PTFE-bonded Pt/Al2O3,  |
        |   >=20 mm above liquid, behind baffle, flame-arrestor mesh, thermocouple) |
        |                                                                            |
        |   +------------------------ PCTFE tube, bore 20.00 +/- 0.05 -----------+   |
        |   |   1.0 M LiOD / D2O  (~12 mL, Rn-free: sparged with Pd-Ag-purified D2)|  |
        |   |                                                                    |  |
        |   |   ======== Pt mesh anode (>=50% open, pitch <=1) ========  g = 4.0 +/- 0.2
        |   |                                                                    |  |
        |   |   ~~~~~~~~ ENTRY FACE (cathode, ground) — annealed, etched Pd ~~~~~|  |
        +---+====================================================================+--+  <- floating seal (FFKM), radial travel >= 1
            |  Pd MEMBRANE 15 +/- 2 um, active dia 20.0, blank dia 26            |      4-point van der Pauw contacts on rim (D/Pd)
            |  EXIT FACE: 4 quadrant skins (Au | Ni | bare | Pd/CaO), Au-capped   |
            |  1.0-mm cross web, 10B marker dot at centre                         |
        +---o--------------------------------------------------------------------o--+
        |   Mo support/catch grid: 0.10 thick, 1.0 hex holes, 0.10 webs (open 0.80),|  <- grounded; carries the
        |   on 0.5-mm ribs along the quadrant cross + 5-mm cells                    |     membrane only if dp flips
        |   Cross septum, 4.0 tall, 0.3 thick, low-alpha Si/electroformed Cu (ground)|
        |   [ dE 25 um, 600 mm2, 4 quadrants ] <- 4.0 +/- 1.0 below exit face       |
        |   [ E 500 um, 700 mm2 ]              <- 1.5 below dE                      |
        |   FRONT VOLUME ~30 cm3, D2 0.50 bar abs (only via hot Pd-Ag element),    |
        |   vacuum-fired 316L, <=10 metal seals, no glass/elastomer/epoxy          |
        +-----+----------------------+-----------------+----------------------------+
              |                      |                 |
      Pd-Ag element (350 C)   all-metal valve to    capacitance gauge (no ion gauge)
      = only D2 in/out path   He analysis manifold
```

### 3.2 Dimension and material table

| Item | Specification | Source |
|---|---|---|
| Membrane material | Pd ≥ 99.95 %, one lot per array (plus a sibling coupon per lot for the He blank) | R1 §4 G8, M8 §6.13 |
| Membrane thickness | **15 ± 2 µm** | M1 rec. 3; M5 (PID depth 16–23 µm); red-team-0 §3.3 |
| Membrane blank / active diameter | Ø26 blank, **Ø20.0 wetted/active** | M2 rec. 6 |
| Heat treatment | Vacuum anneal 800–850 °C, 1 h, ≤ 10⁻⁵ mbar; grains 20–50 µm; light etch (dilute aqua regia, 30 s) on the entry face only | M3 rec. 5; R6 recipes; R1 G5 |
| Mounting | **Floating** seal (FFKM, LiOD-compatible) on the electrolyte side, radial travel ≥ 1 mm; no rigid clamp; masked annulus ≤ 1 mm | M3 rec. 5 |
| Entry face | Bare annealed Pd, horizontal, wetted face up | M2 rec. 5 |
| Anode | Planar Pt mesh (≥ 50 % open, pitch ≤ 1 mm), **g = 4.0 ± 0.2 mm**, parallel within 0.2 mm | M2 rec. 7 |
| Electrolyte | **1.0 M LiOD in D₂O**, ~12 mL (≤ 20 mL for the photoneutron background); Rn-free by sparging with Pd–Ag-purified D₂; sealed | M2 rec. 7; M1 rec. 7; red-team-0 §3.6 |
| Cell body | PCTFE tube, bore 20.00 ± 0.05 mm; 316L lid with Cu-gasket seal; headspace ≤ 15 mL | M2 rec. 6; R6 constraint 1 |
| Cell gas | D₂ **0.50 ± 0.05 bar abs** prefill; trip at ΔP ≥ 10 % of fill (O₂ ≤ 3 %) | R6 constraint 3 |
| Recombiner | Hydrophobic Pt/Al₂O₃, ≥ 20 mm above the liquid, behind a baffle, flame-arrestor mesh, thermocouple | R6 constraint 4; M2 rec. 4 |
| Exit-face skins | Quadrants per §3.3; 1.0 mm Au-capped cross web; skin fabrication order: masked Ni → masked CaO/Pd sputter → Au (room temperature); no anneal after Au | M3 rec. 3; M6 rec. 6 |
| Support/catch grid | Photo-etched Mo, 0.10 mm, 1.0 mm hex holes, 0.10 mm webs (open 0.80), on 0.5 mm ribs under the cross web and at a 5 mm pitch; hole edges radiused ≥ 50 µm | M3 rec. 4, adapted to the balanced ΔP |
| Pressure balance | P_cell − P_front = **+30 ± 20 mbar** (membrane pressed onto the grid). Membrane stress at 50 mbar ≈ 5 MPa; worst case (front lost, 0.5 bar) ≈ 24 MPa, below the 40 MPa annealed yield and the ~150 MPa after the first transit | M3 §5 scaling σ ∝ ΔP^(2/3) |
| Telescope | ΔE 25 ± 2 µm, 600 mm², 4 quadrants; E 500 µm, ≥ 700 mm²; ΔE–E gap 1.5 mm; ΔE at **4 ± 1 mm** from the exit face; ceramic, epoxy-free mounts; ULTRA-AS low-α class; bias ≤ 150 V (below the D₂ Paschen minimum of ~270 V) | M5 rec. 1; M8 §6.2 |
| Septum | Cross, 4.0 mm tall, 0.3 mm thick, low-α Si or electroformed Cu, grounded, aligned to the quadrant web | this document (quadrant isolation; M6 cross-talk ≤ 1 % requirement) |
| Front volume | ~30 cm³, vacuum-fired 316L(N), ≤ 10 metal seals, alumina-brazed feedthroughs, no glass, elastomer, epoxy or ion gauge | M8 rec. 1–2, 5 |
| Front gas | D₂ 0.50 bar via a **hot Pd–23Ag element** (350 °C) as the only fill and exhaust path; 1 % Kr in the *cell* fill gas as a leak tracer | M8 rec. 4, 10 |
| Electrical | Membrane = hard ground of the detector system; floating linear galvanostat drives the anode; dummy Si behind 200 µm Al on the same electronics; waveform digitisation ≥ 100 MS/s | red-team-0 §3.4 |
| Loading gauge | 4-point van der Pauw on the membrane rim → R/R₀ → mean x | M2, M3 |
| Other sensors | Acoustic-emission piezo on the flange (H1b crack tag); membrane-flange RTD; cell pressure; front capacitance gauge | red-team-0 §4.12 |

### 3.3 Exit-face skins: regimes by quadrant

With a thin membrane, each quadrant sets its own local (x, J). Lateral coupling reaches only about one membrane thickness (red-team-0 §3.1; M6 §5.7).

| Code | Skin (exit face) | Regime | Expected local state | Tests |
|---|---|---|---|---|
| **L** | Au 20–50 nm, continuous | loading | x ≈ x_in ≥ 0.9 through the foil; J ≈ 0 | SRI/ENEA static loading; surface-vs-bulk discriminator; He retained (measured by melt) |
| **M** | Ni 2–20 nm e-beam (thickness chosen after the witness test, §9 P0b), or a photolithographic Au mask with open fraction φ (pitch ≤ 10 µm) | loading + flux | x_exit 0.92–0.945, J ~ 13–86 mA cm⁻² equivalent, **if x_in ≥ 0.95** | McKubre's (x − x₀)²·flux product |
| **F** | bare annealed Pd | flux | x_exit ≈ 0.63 (clamped by 0.5 bar D₂); J maximal | Lipson / NTT / Czerski-511 flux regime; the flux lock-in segment |
| **X** | Pd 40 nm / [CaO 2 nm / Pd 18 nm]×5 | oxide interface + flux | behaves like F (no barrier function), with irreducible oxide/metal interfaces within 140 nm | Iwamura-type interfaces (reversed geometry); stable-oxide analogue of Kasagi's PdO |

**Rejected skins:**
- **PdO:** permeating D reduces it in 0.2–100 s (M3), so it cannot hold the regime M1 wanted.
- **Ni/Cu stack:** its claims apply at 250–900 °C (M6).
- **Lipson back-face Au heterostructure:** it starves the sector (M6).

**Additional variant:** a 1 µm electrodeposited, defect-rich Pd layer (**ED**) on the exit face. It raises vacancy density within the escape depth, as M1 rec. 4 asks. Room-temperature electrolysis cannot create superabundant vacancies (M3), so they must be deposited.

**Rotation plan** (4 active membranes × 4 quadrants = 16 slots; position is balanced against skin):

| Membrane | Q1 (NE) | Q2 (NW) | Q3 (SW) | Q4 (SE) |
|---|---|---|---|---|
| A1 | L | M | F | X |
| A2 | F | L | X | ED |
| A3 | X | ED | L | M |
| A4 | M | F | ED | L |

The skins appear L ×4, F ×3, M ×3, X ×3, ED ×3. The H₂O twin and the D₂O/Pt twin carry the A1 layout.

![Exit-face skins](figs/iter1_skins.png)

Drawings are generated by `sim/design_iter1_drawing.py`.

### 3.4 What the telescope sees (depth tomography)

| Origin of a 3.02 MeV proton | Energy at the ΔE (15 µm PdD₀.₉ + 0.5 bar D₂, 4 mm) | Counted in |
|---|---|---|
| Exit face (skins, top ≤ 1 µm) | 3.00 MeV | peak window 2.60–3.10 MeV |
| Mid-foil (7.5 µm) | ~2.55 MeV | PID window 1.41–3.10 MeV |
| **Entry face (electrolyte side)** | **2.05 MeV** | PID window |

- Tritons (1.01 MeV) escape only from the top ~4 µm; ³He (0.82 MeV) from ≤ 0.9 µm, arriving at ≈ 0.70 MeV (ΔE-only).
- ²²²Rn/²¹⁸Po α from the electrolyte are stopped by 15 µm. ²¹⁴Po α emerge at 1.26 MeV and are rejected by PID.
- Computed with `sim/m5_stopping.py` (validated against NIST PSTAR/ASTAR).

## 4. Array and detector house (Tier 1)

### 4.1 Positions

| Position | Cell | Role |
|---|---|---|
| A1–A4 | DFM, D₂O, rotation above | **active samples** (N = 4 → P(≥ 1 active) = 0.59 at p = 0.2) |
| T-H | DFM, H₂O/LiOH, A1 skins, **matched on measured x and J**, not current | charged-particle and ⁴He isotope control |
| T-Pt | DFM geometry, Pt membrane (non-absorbing; flux-blocked), identical D₂O volume and current waveform | neutron control (cancels D(n,2n)/D(γ,n) in D₂O), EMI/thermal/electrolysis-artifact control, cell→front leak control |
| B | Blank telescope under 100 µm Al, same electronics | internal, cosmic and EMI background |

A midpoint ABBA swap of telescopes between A-cells and T-H is pre-registered (M5 rec. 6).

### 4.2 House layout (plan view, schematic)

```
   +------------------------------------------------------------------+  20 cm borated HDPE (30 cm if floor allows)
   |  1 mm Cd                                                          |
   |  +------------------------------------------------------------+   |
   |  |  3He tubes (24 x 1" x 40 cm, 4 atm) in 4.5 cm HDPE front / |   |
   |  |  6 cm behind, lining 4 walls                                |   |
   |  |   +-----------------------------------------------------+   |   |
   |  |   |  5 cm OFHC Cu inner shield (no Pb near the bank)     |   |   |
   |  |   |    [A1]   [A2]   [A3]      <- 90 mm pitch            |   |   |
   |  |   |    [A4]   [T-H]  [T-Pt]   [B]                        |   |   |
   |  |   |  LaBr3 pair (2"x2") back-to-back across the array,   |   |   |
   |  |   |  on a rail so it can be positioned at any cell        |   |   |
   |  |   +-----------------------------------------------------+   |   |
   |  +------------------------------------------------------------+   |
   |  plastic-scintillator muon veto on top and on two sides            |
   +------------------------------------------------------------------+
   Inner cavity ~ 40 x 30 x 30 cm; expected 3He-bank efficiency 0.12-0.18 (M5 scaling from 0.23 at r = 8 cm)
```

**Provisional γ/e± layout (to be finalised by M7):**
- A back-to-back LaBr₃(Ce) or CeBr₃ pair gives 511–511 keV coincidences, target ε ≥ 2 % and background ≤ 2×10⁻³ cps (M1 rec. 6).
- A 3–25 MeV window covers the e± pair/bremsstrahlung sum.
- An inner Cu shield (not Pb) handles the 2.614 MeV ²⁰⁸Tl line. Pb stays out, per M5's warning about muon-induced neutrons.
- ²²Na calibration.

## 5. Helium system (Tier 1, co-primary H2 channel)

- **Two He channels per cell** (M8):
  - **Front volume** for exit-face reactions. Release fraction 0.85 for ≤ 10 nm skins; the Au-capped L quadrant retains its He.
  - **Cell headspace** for entry-face reactions. The exit chamber cannot see these (release fraction 0 into the front).
- **Shared analysis manifold** (0.1–0.2 L, 250 °C bakeable):
  - room-temperature getter plus an optional hot getter;
  - an **HR-QMS (R ≥ 500)**, never a unit-resolution RGA for ⁴He;
  - ⁴He and ³He 3-valve pipettes;
  - capacitance gauge;
  - two all-metal sample ports for an external magnetic-sector MS with ³He isotope dilution (the primary result).
- **Cycle.** Static windows of 1 d (at ≤ 10 mA cm⁻² equivalent permeation) or 7 d (≤ 1 mA cm⁻²). Aliquots at 1, 2, 4 and 7 d, with ≥ 2 flux on/off cycles per window.
- **Floors (5σ, full release):** 0.39 nW (1 d) and 0.16 nW (7–30 d) at 23.85 MeV/He. The floor is set by background stability (σ_B ≈ 7×10⁵ atoms/day).
- **Post-run.** Laser-cut the quadrants and melt each separately: sector-resolved He, detection limit ≈ 0.4–1.6 nW-equivalent. Also melt a sibling unexposed coupon per lot.
- **Release calibration.** ³He implanted off-site into designated sectors of a sibling membrane (sample preparation, not a device beam; ADR-004) measures f_exit(depth) in real loaded, cycled PdDₓ.
- **Hygiene:**
  - vent with LN₂ boil-off N₂ only;
  - no He leak testing after the final bake (use Ar/Kr);
  - no He cylinders in the room;
  - log room He.

## 6. Calorimetry (Tier 1 subset + Tier 2)

- **Instrumented cells:** A1 and T-H each sit in a **Seebeck envelope with a ≥ 10 mm OFHC Cu isothermal shell** (SEEB1\*). A ceramic thermal break (≤ 0.01 W/K) separates each from its front volume, and a rim heater calibrates membrane heat.
- **Performance:** 5σ ≈ 4–20 mW (M4 rec. 9).
- **Power metrology:** simultaneous V and I at ≥ 100 kS/s; never ⟨V⟩⟨I⟩, which gives +13 % false excess at m = 0.5.
- **Energy bookkeeping:** permeation enthalpy (1.527 V × I_perm) and loading enthalpy (1.347 V × χI) are pre-registered.
- **Role.** Calorimetry is **secondary** to ⁴He for H2, which is ~10⁷× more sensitive. It is kept for comparability with the legacy heat claims.

## 7. Tier-2 arms

### 7.1 C3-G gas-entry cells (Iwamura variant)

- **Membrane.** Same cell body and telescope. The entry module is D₂ gas (1 atm, 70 °C) instead of electrolyte.
- **Stack.** Pd 100 µm substrate carrying Pd 40 nm / [CaO 2 nm / Pd 18 nm]×5 on the **entry** face, plus an isotope-tagged target: ⁸⁶Sr-enriched Sr predicts ⁹⁴Mo under the claimed ΔA = +8 rule.
- **Geometry.** The exit side is at low pressure with an RGA to log J. Target J = 3–7×10¹⁷ D cm⁻² s⁻¹. The telescope sits on the **entry** side, looking through 4 mm of 1 atm D₂.
- **Controls.** H₂ twin; a no-CaO stack; static D₂ with no flux.
- **Analysis.** Blinded pre/post ToF-SIMS and HR-ICP-MS at two external labs, with a pre-run contamination survey.

### 7.2 C1 SRI/ENEA replica twin

- **Cathode and anode.** Pd wire Ø1.00 × 30 mm between PTFE end plates; Pt helix at R = 10 mm (M2 §6.1).
- **Electrolyte.** 0.1 M LiOD (lineage), 60 mL in an all-metal 316L cell.
- **Calorimeter.** SEEB1\* twin (D₂O/H₂O) in one bath (M4).
- **He.** Headspace ⁴He, with the cathode melted after the run.
- **Other channels.** Placed in the neutron bank; CR-39 channel C (6 µm Mylar + 55 µm Al) clamped against the wire with ≤ 10 µm of electrolyte film (M5 rec. 9).
- **Protocol.** SRI: ≥ 3 weeks at x ≥ 0.9, current steps for |dx/dt|.

## 8. Tier-3 arms (outline)

### 8.1 Co-deposition on a PSD scintillator (SPAWAR conditions)
- Au film (50 nm) on EJ-276 PSD plastic, read by a PMT through a light guide.
- Pd/D co-deposited from PdCl₂ + LiCl in D₂O, in an open, vented cell. An H₂O twin runs in parallel.
- Cost ~$5–15k (red-team-0 §4.2).
- The optical readout is immune to the capacitive pickup that rules out Si under a cathode.

### 8.2 Sealed nanocomposite ⁴He furnace (Kitamura/Takahashi; Clean Planet variant)
- 100 g CNZ/PNZ (Cu–Ni or Pd–Ni in ZrO₂) at 200–300 °C in an all-metal reactor. A black, water-cooled, flow-read jacket gives ±0.1–0.2 % (M4 rec. 10).
- An inert-bed twin runs alongside. ⁴He/³He is measured by the same manifold (³He at m/Δm ≥ 510).
- The Ni/Cu multilayer variant (6×[Cu 2 / Ni 14 nm] on 25×25×0.1 mm Ni) is H₂-loaded at 250 °C, then desorbed while ramping to 900 °C.

## 9. Protocol (pre-registered)

| Phase | Duration | Action | Key observables |
|---|---|---|---|
| **P0a Commissioning** | 3–4 wk | ²⁴¹Am + pulser (Si), ²²Na (511), ⁴He/³He spikes, blank He accumulations (≥ 3 × 1 d, ≥ 2 × 7 d), calorimeter calibrations (M4 §6.3), background runs ≥ 10× planned live time | backgrounds, efficiencies |
| **P0b Witness membrane** | 1 wk | Bare-Pd witness DFM: measure J(i) and x at 22 °C (M3 rec. 3). Choose the M-skin (Ni thickness or Au-mask φ) from the result | exit-face law |
| P1 Loading | 1 d | One α→β transit at 5 mA cm⁻², then 20 → 50 → 100 → 300 mA cm⁻² in 30 min steps, tracked by R/R₀ (M3 Protocol A) | x(t) |
| **P2 High-load hold** | ≥ 3 wk | 200–300 mA cm⁻², 22 ± 1 °C, keep-alive ≥ 20 mA cm⁻² on UPS; ⁴He static windows | primary H1/H2 data |
| P3 Flux modulation | 2 wk | (a) current square wave +300/−100 mA cm⁻², period 20–240 s (Protocol B1); (b) **orthogonal**: constant current, front-and-cell pressure steps 0.2 ↔ 1.0 bar (x_exit 0.61 ↔ 0.67 on F/X) | flux-correlated counts |
| P4 Warm phase | 1 wk | Cell at 70–80 °C (flux over a wider J range; M1 rec. 9, limited by the electrolyte) | J-dependence |
| P5 Site factory (optional) | ≤ 10 cycles | α/β cycling (Protocol B2) **on A4 only**, with leak interlock; then repeat P2 | crack- and vacancy-correlated counts (AE-tagged) |
| P6 End | 1 wk | Deload; final He windows; section and melt quadrants; SIMS/ICP-MS; AFM PSD/EBSD covariates | ash |

### Primary pre-registered tests (one per hypothesis; look-elsewhere-corrected across 4 tests)

| Test | Hypothesis | Statistic |
|---|---|---|
| T1 | H1a | Si PID-window counts, Σ(A-cells) − T-H, per skin type, current-on. Blinded; background exposure ≥ 10× live time |
| T2 | H1c | 511–511 coincidence rate, A-cells − T-Pt |
| T3 | H2 | Front-volume and headspace ⁴He, A-cells − T-H, plus quadrant melts |
| T4 | H3 | Tagged-target products in C3-G versus H₂ and no-CaO controls |

**Secondary observables:** neutron bank rate (A − T-Pt), heat (A1 − T-H), triton and ³He lines, AE correlation, flux lock-in phase.

**Stopping and claim rules:**
- A positive requires 5σ in a primary test **and** replication in ≥ 2 membranes **and** the expected depth/energy signature **and** absence in the matched control.
- The pre-registered statement "Protons below ~10⁻² fusions s⁻¹ appear without detectable neutrons" holds (M5).
- An excess-heat claim additionally requires the independent mass-flow cross-check (M4 rec. 12).

## 10. Predicted sensitivities versus claims (30 days, 5σ)

| Channel / claim | This design reaches | Claimed magnitude | Ratio (reach / claim) |
|---|---|---|---|
| Si PID window, per cell (H1a) | 3–8×10⁻⁵ fusions s⁻¹ (M5) | Lipson ~10⁻²–1 s⁻¹ | ≤ 10⁻² — saturated |
| Si, PdO/oxide-interface prediction (M1) | threshold 2.2 events/day | Tohoku-static: ~130/day | 0.02 |
| ⁴He front (H2, exit skins) | 0.16–0.39 nW | SRI-scale 0.1–1 W | 10⁻⁹ |
| ⁴He headspace + melt (H2, entry/bulk) | ~0.4–1.6 nW-equivalent (melt) | 0.1–1 W | 10⁻⁹ |
| Heat (A1, C1) | 4–20 mW (DFM), 7–27 mW (C1) | 0.1–1 W | 10⁻¹–10⁻² |
| 511–511 (H1c) | M7 (pending) | not stated (≈10⁻³ s⁻¹ assumed) | — |
| Neutron bank | ~2×10⁻² fusions s⁻¹ (systematics 0.05–0.1) | — | secondary |
| C3-G products (H3) | ~10¹⁰ atoms cm⁻² | 10¹²–10¹⁴ cm⁻² | 10⁻²–10⁻⁴ |

Each channel reaches ≤ 1 % of the claimed magnitude where a magnitude exists. Under ADR-003, detector sensitivity is therefore saturated. The remaining levers are condition fidelity and sample count.

## 11. Safety envelope (R6 constraints → design)

| R6 constraint | Implementation |
|---|---|
| Headspace ≤ 50 mL (closed) | ≤ 15 mL per DFM cell |
| Closed-cell O₂ control | D₂ prefill 0.5 bar; trip at ΔP ≥ 10 % of fill; recombiner thermocouple |
| Recombiner placement | ≥ 20 mm above liquid, behind a baffle, flame arrestor |
| O₂-clean parts | PCTFE/PTFE, Pt, Pd, 316L, Mo, Cu only; no oils |
| Cathode ≤ 0.3 cm³ | Membrane 4.7×10⁻³ cm³ (stored D ≈ 0.07 kJ) |
| SELV ≤ 60 V | 3.2 V cell at 200 mA cm⁻² (M2) |
| Monitoring and trips | pressure, temperature, level, H₂ 1 % room sensor, neutron bank at 10× background, Δp interlock (Si bias off, valves close), front pressure-rise or m/z 20 leak interlock |
| Radiation | exempt check sources only (²⁴¹Am, ²²Na, ¹³⁷Cs); licensed neutron sources used only at a partner facility |

## 12. Bill of materials and budget (estimates; quotes required)

| Block | Tier 1a (minimum credible: A1, A2, T-H, T-Pt, B) | Tier 1b (full: A1–A4, T-H, T-Pt, B) |
|---|---|---|
| Pd foil, skins (Au/Ni/CaO sputter, ED Pd), masks, Mo grids, septa | $12k | $20k |
| Si telescopes (quadrant ΔE + E, low-α) | 5 × $6k = $30k | 7 × $6k = $42k |
| Digitisers and HV (≈ 5 channels per telescope + veto/LaBr₃/³He) | $20k | $32k |
| Cell and front hardware (316L CF, PCTFE, Pt mesh, feedthroughs, gauges) | $12k | $20k |
| Gas and He system (Pd–Ag elements, D₂, manifold, getters, turbo, pipettes) | $25k | $30k |
| HR-QMS (or external sector-MS service contract) | $15k service | $55k instrument |
| Neutron bank (³He tubes, HDPE, borated PE, Cd) | 12 tubes: $35k | 24 tubes: $65k |
| LaBr₃/CeBr₃ pair, Cu shield, plastic veto | $30k | $30k |
| Galvanostats (linear, floating), UPS | $8k | $14k |
| Seebeck calorimetry (A1, T-H) | $15k | $15k |
| Consumables (D₂O, LiOD, ¹⁰B, ⁶LiF, CR-39) | $4k | $6k |
| External analyses (sector-MS He melts, SIMS/ICP-MS, AFM/EBSD) | $20k | $35k |
| **Total (hardware + analyses, excluding labour)** | **≈ $226k** | **≈ $364k** |
| Tier 2 (C3-G pair + C1 twin) | | +$60–90k |
| Tier 3 (scintillator arm; nanocomposite furnace) | | +$15k / +$60k |

R6's recommended $50k tier cannot buy Tier 1. The cheapest credible subset is **one DFM D₂O cell, one H₂O twin and a blank telescope, with external He analysis and no neutron bank or LaBr₃ (~$95k)**. It keeps the primary H1 channel (T1) and the He channel (T3) and gives up T2, the neutron cross-check and multiplicity.

## 13. Go / no-go logic for iteration 2

| Outcome after P6 | Iteration 2 |
|---|---|
| Any primary test positive under the claim rules | Blind independent replication at a second site with the same drawings. Add DSSD localisation and IR imaging of the exit face. Scale the number of membranes of the positive skin type |
| T1/T3 null, C3-G positive | Pivot to the gas-permeation geometry. Vary the tag isotopes and the multilayer period |
| All null | Publish the in-regime limits (§10). Then redirect to the highest remaining ADR-003 weight: the Tier-3 nanocomposite ⁴He arm at scale, and a C1 SRI replica using *SRI's own* cathode lots if obtainable. Stop the membrane line unless a new claim appears |
| Hardware failure (membrane rupture, leak) | Fix per the §11 interlocks. Consider 25 µm Pd with a finer grid (M3 baseline) and accept that the entry face becomes invisible |

## 14. Traceability

| Parameter | Value | Governing source |
|---|---|---|
| Objective | claim-conditioned, multi-sample | ADR-003; red-team-0 §2 |
| Platform | DFM | ADR-002 rev 2; R7 §7; R5 §4.4 |
| Thickness | 15 µm | M1 rec. 3; M5 useful depth; red-team-0 §3.3 (M3 prefers 25 µm under vacuum loading; pressure balance removes that driver) |
| Active diameter | 20 mm | M2 rec. 6; M5 rec. 2 |
| Anode gap, electrolyte | 4.0 mm, 1.0 M LiOD | M2 rec. 7 |
| Front gas | 0.5 bar D₂ via Pd–Ag | red-team-0 §3.2; M6 rec. 1; M8 rec. 4; triton/³He loss from `sim/m5_stopping.py` |
| Skins | L/M/F/X/ED; PdO rejected | M3 rec. 3; M6 §6; M1 rec. 2 and 4 |
| Telescope | 25/500 µm, 4 mm, quadrants | M5 rec. 1 |
| Controls | T-H matched on x and J; T-Pt for neutrons | M3 rec. 8; M1 rec. 7; M5 rec. 6 |
| He | two channels, HR-QMS, melts | M8 |
| Calorimetry | SEEB1\* on A1 and T-H | M4 |
| Protocols | A, B1, B2 (≤ 10 cycles) | M3 rec. 6–7 |
| Positive controls | ¹⁰B/⁶LiF markers, exempt sources, off-site calibration | ADR-004 |
| Stimuli excluded | THz beats, ultrasound, nanogap/hot-spot arrays | M6 §6; M1 rec. 8; R4 |

## 15. Open risks (for red-team-1)

1. **Entry-face loading ceiling.** M-regime needs x_in ≥ 0.95, which M2 says requires an "exceptional" (recombination-poisoned) surface. The fallback is that M collapses into L or F. A thiourea-class poison would have to be added to the H₂O twin too.
2. **Exit-face law uncertainty** (10³–10⁴× across the literature; M3). This is mitigated by P0b.
3. **Si in 0.5 bar D₂.** Long-term leakage current and noise are not demonstrated. Fallback: detector chamber at ≤ 10⁻⁴ mbar and front ΔP carried by the grid (M3 baseline: 25 µm on a 1.0 mm grid).
4. **Grid shadowing** (open fraction 0.80 at 0.10 mm thickness) is not yet in M5's efficiency Monte Carlo; expect ε per quadrant ≈ 0.15.
5. **Septum edge effects** (partial collection near quadrant boundaries). The M5-style aperture mask sits 1 mm inside the edges.
6. **γ/e± layout pending M7.**
7. **Cost.** Tier 1b ≈ $364k exceeds hobby scale. Tier 1a or the $95k minimum subset trade multiplicity for cost.
