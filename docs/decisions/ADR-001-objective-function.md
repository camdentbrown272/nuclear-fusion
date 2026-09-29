# ADR-001 — Objective function: sensitivity to a localised anomaly, not excess heat

**Status:** accepted · **Date:** 2026-09-29 · **Basis:** M0, charter

## Context
Most LENR experiments since 1989 were built around excess heat in bulk Pd cathodes. M0 shows that:

1. For standard D+D branching, heat is ~8×10¹¹ times less sensitive than neutron detection. One watt of standard D+D would be lethal (~10 Sv/h at 1 m).
2. The reaction rate scales locally as U_e,eff^35–50. The total is dominated by the rarest high-enhancement sites; the top 0.1 % of sites can carry 60–90 % of the rate.
3. The gap between known physics and a detectable rate is ~100 e-folds. Engineering buys ~25–30. A positive result therefore requires a real anomaly localised somewhere in the material.

## Decision
The design objective is to **maximise the probability of detecting a localised anomalous enhancement, wherever the best-supported claims and theories place it, at the lowest possible anomaly strength, with controls that make a detection credible.**

Concretely, geometry is scored by the minimal enhancement factor it would detect at 5σ in 30 days, for each candidate site class (the M1 matrix). Heat is kept as a secondary channel only to test H2 (heat + ⁴He correlation).

## Consequences
- Active material must be **thin and interface/defect-rich**, and located **within the charged-particle escape depth** facing spectrometers.
- **Charged-particle spectroscopy in vacuum** (energy-resolved, near-zero background) is the primary H1 channel. Neutrons are secondary, covering the whole volume.
- An identical light-hydrogen **control** is part of the geometry.
- Bulk-heat-only designs (a classic Fleischmann–Pons Dewar with a thick rod) are rejected as primary instruments.

## Alternatives rejected
- *Maximise excess heat*: insensitive under H1, and artifact-dominated historically.
- *Maximise total deuterium inventory* (large cathodes): the rate sits in rare sites, and bulk volume adds background and slow loading without adding visible sites.
