# ADR-006 — γ/e⁺e⁻ channel as a separate station; modulation protocol from M7

**Status:** accepted · **Date:** 2026-09-29 · **Basis:** [M7](../models/M7-gamma-epem-detection.md), M5, ADR-005. It amends iteration-1 rev B §2, §4, §9, §10 and §12, and pre-registration v1 §2–§3.

## Context
M7 overturned the provisional γ layout in rev B, a LaBr₃ pair on a rail inside the neutron house, on four counts.

1. **The 511 keV source is extended.** Only 0.2–3.8 % of the IPC positrons annihilate in the Pd foil. The rest annihilate within ~15 cm (steel, electrolyte, detectors, moderator), and 14 % annihilate in flight. A back-to-back pair aimed at one foil gains nothing from pointing at it.
2. **LaBr₃ is the wrong crystal.** BGO gives 4× NaI and 6× HPGe in 511–511 efficiency per volume. The best option overall is a **well calorimeter**: an 8″×8″ NaI well around a mini permeation cell reaches 68 % efficiency in 12–30 MeV. Its MDA is **5×10⁻³ pairs s⁻¹ at sea level** and **3×10⁻⁴ at ≥ 30 m w.e.**, against 9×10⁻³ for the BGO pair. It also costs less ($20–35k vs $55–90k).
3. **γ shielding needs 10 cm of Pb** (15 cm underground). M5 forbids Pb near the ³He bank because muon-induced neutrons from Pb fake the neutron channel. The two cannot share a house.
4. **The modulation period is set by electrochemistry, not diffusion.** In a 12 µm foil diffusion takes seconds. The electrochemical and surface relaxation τ_s sets the lock-in efficiency: at P = 2 h, η = 0.96 for τ_s = 1 min but 0.24 for τ_s = 30 min. A heat or field change of 5 % between on and off fakes 4.5σ unless current is **steered to an auxiliary cathode** at constant total cell current.

## Decision

### 1. γ station: a separate Tier-1 stand
- **Detector.** An 8″×8″ NaI(Tl) well calorimeter.
- **Cells.** A **mini permeation cell** sits inside the well. It uses the same 12 µm Pd, electrolyte and anode as the DFM, and runs as an FX-type membrane with its exit face on a small UHV volume. An RGA (m/z 3, 4) measures the exit flux. There is no Si and no He accounting, so the cell can be a simpler vented design with a recombiner.
- **Isotope factorial.** Two identical mini cells, one D₂O and one H₂O, alternate in the well in **5-day blocks**.
- **Shield, inside out.** 5 cm 5 %-borated HDPE, then 10 cm Pb (15 cm if a ≥ 30 m w.e. site is available), then a 5 cm plastic veto on 5 faces with a 20 µs extended window. N₂ purge. The stand is ≥ 2 m from the ³He house, and there is no Pb in the house.
- **Channels:**
  - summed deposit in 12–30 MeV, which sits above every thermal-capture γ line in the setup;
  - 511 keV coincidences within the well segments;
  - MIP tag.
- **Calibration.** Exempt ²²Na (0.1 µCi), plus cosmic muons and Michel electrons above 12 MeV.
- **Blank.** 14 days with an Au foil at the same current, to fix the hadronic background normalisation (×3 uncertainty).
- **Site.** Sea level by default. A ≥ 30 m w.e. site is preferred if one is accessible, since it gives ~17× lower MDA.

### 2. The neutron house loses its γ detectors
The LaBr₃ pair and its rail are removed. The house keeps its Cu inner shield, ³He bank, borated HDPE and a muon veto for the Si channel. The Si telescopes still record MIPs (13 % of IPC events), and the ΔE–E classification must keep them out of the proton window (≤ 10⁻⁴ leakage, M7).

### 3. Modulation protocol, for all DFM cells and the γ station
- **P0b** measures the step response of each membrane type, giving τ_eff.
- **Lock-in period** P = 12 τ_eff, clipped to 1–6 h (default 2 h), with randomised block order. This **replaces the 240 s period** in rev B P3 and pre-registration §3.2.
- **Current steering.** An **auxiliary Pt cathode** (ring, in the same electrolyte) takes the current during "off" half-periods, so total cell current, heat and field stay constant to < 1 %. Anodic half-cycles (+300/−100 on FX) remain a separate descriptive test.
- **Single pre-registered regressor per membrane type:**
  - FX: exit flux, from a thermal mass-flow meter on the Pd–Ag exhaust (DFM) or the RGA (γ station);
  - H-L and H-M: |dx/dt| from the rim van der Pauw.

### 4. Pre-registration
M7 is now merged, so **T2 is frozen now**. It enters the primary family with weights T1 0.4, T3 0.4, T2 0.2.
- **T2 statistic:** a joint on/off likelihood of the γ-station 12–30 MeV and 511-coincidence rates, D₂O vs H₂O blocks, with the MC-fixed ratio between the two channels.
- **Local thresholds:** T1 and T3 at p ≤ 1.15×10⁻⁷ (5.18σ); T2 at p ≤ 5.7×10⁻⁸ (5.30σ).

### 5. Budget
The LaBr₃ pair, Cu block and rail are removed (−$30k); the house veto stays (+$15k). The γ station adds:
- NaI well with PMT: $20–35k;
- ~1 t Pb: $8–12k;
- veto: $12k;
- 2 mini cells and RGA: $10k.

That is ≈ $60k, bringing Tier 1 to **≈ $0.6M**.

## Consequences
- **Coverage.** The Czerski claim is now tested where a proponent would accept it: a D-loaded, D-permeated Pd foil inside a near-4π calorimeter with an isotope factorial. It no longer depends on a pair of crystals viewing an 8-cell array.
- **Priority between channels.** For surface sites, the proton channel beats the γ channel unless e⁺e⁻/p ≳ 3 (M7 §7). For sites deeper than the proton range, only the γ channel sees anything. The γ station therefore covers the bulk.
- **Hardware changes.** Every DFM cell gains an auxiliary Pt cathode and a mass-flow meter on the Pd–Ag exhaust. The CAD model will be updated with the red-team-1B changes.
