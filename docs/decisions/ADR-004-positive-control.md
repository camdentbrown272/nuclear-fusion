# ADR-004 — Positive controls without beams

**Status:** accepted · **Date:** 2026-09-29 · **Basis:** red-team-0 §3.8 (F8), charter scope

## Context
Red-team-0 notes that nothing in the design shows *in situ* that the chain would register D–D-like products born in the actual foil: foil → gap → telescope → cuts → lock-in analysis. R3, R4 and R7 recommend a keV D⁺ probe beam for this. The charter forbids any beam, including diagnostic ones, because the project is strictly cold. We need positive controls that respect that constraint.

## Decision
Use **neutron-capture depth markers**: cold nuclear reactions at a known location with known product energies.
1. **¹⁰B marker dots.** 1 µm ¹⁰B-enriched boron, or ¹⁰B₄C, as 1–2 mm dots on a masked region of the exit face (outside the active sectors) and on the grounded screen.
   - Thermal neutrons give ¹⁰B(n,α)⁷Li: α at 1.47/1.78 MeV and ⁷Li at 0.84/1.01 MeV. σ_th is 3840 b.
   - A 1 µm ¹⁰B layer captures ~5 % of incident thermal neutrons.
   - At sea level the cosmic thermal flux is ~4×10⁻³ n cm⁻² s⁻¹ outdoors [BK: Gordon et al. 2004], rising to ~10⁻² inside the HDPE moderator. That gives **~20–40 captures per cm² per day, or ~5–15 detected particles per cm² per day at Ω/4π ≈ 0.2–0.3, with no radioactive source at all**.
   - The alphas (1.47/1.78 MeV) and ⁷Li ions (0.84/1.01 MeV) fall below the degraded-proton window (≥ 2.2 MeV), so the markers do not contaminate the H1 window.
   - The rate is cross-checked by the ³He tubes, which see the same thermal flux.
2. **⁶LiF marker on the entry side** of one membrane per array (electrolyte-insoluble variant: ⁶LiF under a 20 nm Pd cap).
   - ⁶Li(n,t)α gives t at 2.73 MeV and α at 2.05 MeV. The tritons cross an 8–12 µm Pd membrane and emerge at a predicted energy.
   - This **calibrates the depth-to-energy tomography through the real foil**.
3. **Source-driven boost (optional).** Placing a moderated AmBe or ²⁵²Cf source next to the house (if an institutional licence exists) raises the marker rates 10³–10⁴×. Toggling the source on and off is a **positive control for the lock-in pipeline**: a known modulated nuclear signal through the same electronics and cuts.
4. **Exempt-quantity check sources:**
   - ²⁴¹Am (smoke-detector class) for Si energy scale;
   - ²²Na (≤ 10 µCi exempt) for 511–511 keV coincidence efficiency;
   - ¹³⁷Cs/⁶⁰Co for γ energy scale.
5. **⁴He spike standard** (calibrated leak) and a ³He isotope-dilution spike for the helium channel.

## Rejected
- keV D⁺ probe or glow-discharge implant: excluded by the charter's cold-only scope, despite the higher information value. Recorded here so reviewers know it was a deliberate choice.
- Pyroelectric D–D source: this is beam-target fusion by another name, so it is also excluded.
