# 00 — Project Charter

**Owner of all design decisions:** Claude (lead designer). Every decision is recorded in `docs/decisions/` with its reasoning so any reviewer can audit it.

## Mission

Design the geometric configuration of a **cold** condensed-matter nuclear experiment that maximises the probability of a **measurable, credible** nuclear signal.

- **Cold** means:
  - solid-state samples at ≤ 1200 K (thermal energies ≤ 0.1 eV, versus ~10 keV for hot fusion);
  - no plasma, no fusor, no glow discharge;
  - **no ion beam or accelerator of any kind, including "diagnostic" probe beams.**

  Allowed drivers: electrochemistry, gas pressure/permeation, temperature cycling, current/field modulation, low-power optical/acoustic stimulation, magnetic fields. Sealed radioactive calibration sources and cosmic-ray neutrons may be used as **positive controls** (ADR-004). *Amended 2026-09-29 after red-team-0 §3.8: the cap was raised from 600 K so that the gas-phase Ni/Cu and nanocomposite claims can be tested; the beam ban was made explicit.*
- **Measurable** means: a pre-registered observable exceeds background by ≥5σ, with a matched control that does not show it.
- **Credible** means: the result survives the artifact catalogue in `docs/research/R5-detection-and-artifacts.md` and could be handed to a skeptical nuclear physicist without embarrassment.

## Honest prior (stated once, then we work)

1. No cold-fusion/LENR claim since 1989 has met the reproducibility standard of mainstream physics. Several careful programmes (NHE Japan 1992–98, Google/Nature 2019, most of ARPA-E 2022–25) reported nulls or ambiguities.
2. First-principles numbers (`docs/models/M0-rate-budget.md`) show that ordinary physics falls short of a neutron-detectable D–D rate by **94–115 e-folds** in rate (41–50 decades). Even with every deuteron in a 0.5 g cathode active and 100× better detectors, **89 e-folds remain**. **Geometry alone cannot close the gap.** A positive result requires an anomalous enhancement that is not in standard physics.
3. Therefore the rational definition of "best geometry" is (revised in ADR-003):

   **the configuration that reproduces, at a detector-visible location, the conditions under which the best-graded claims report anomalies, in enough independent samples, with detectors that reach ≤ 1 % of each claimed magnitude and controls that make a detection credible.**

   A null from such a design is informative. A null from outside the claimed regime would be dismissed, as the 2019 Google study was.

This framing is not a hedge; it is the optimisation target. Designs that maximise "excess heat" alone are rejected because heat is ~10¹² times less sensitive than particle detection for standard D–D and is the most artifact-prone observable in the field's history.

## Two hypotheses the design must serve simultaneously

| | H1: enhanced conventional D–D | H2: anomalous "LENR" channel |
|---|---|---|
| Products | n (2.45 MeV) + ³He; p (3.02 MeV) + T | ⁴He + heat, no strong γ/n |
| Most sensitive signature | charged-particle spectroscopy, neutrons | **⁴He** (extraction and static accumulation), then heat |
| Reactions/s for 5σ | ~10⁻³–10⁻² s⁻¹ | 2.6×10⁹ s⁻¹ per 10 mW (calorimetry); ~4×10³ s⁻¹ via ⁴He collection |
| Mainstream compatibility | Compatible in kind (screening/threshold resonance) | Requires new physics |

A single geometry that is blind to either hypothesis is inferior to one that tests both.

The full hypothesis set (ADR-003) also includes:
- **H1b**: crack- or transient-driven conventional fusion, a confound tagged by acoustic emission;
- **H1c**: the e⁺e⁻ threshold-resonance channel, seen as 511 keV coincidences;
- **H3**: transmutation, tested with isotopically tagged targets;
- **H4**: reactions involving light hydrogen.

## Principles

1. **Quantify before choosing.** Every geometric choice cites a model in `docs/models/` or a source in `docs/research/`.
2. **Controls are geometry.** The H₂O/light-hydrogen twin, the Pt-cathode twin, and the detector blank are part of the configuration, not an afterthought.
3. **Energy spectra over totals.** A 3.02 MeV proton peak is worth more than any calorimetric excess.
4. **Tail, not mean; quality over count.** Rates scale as ~U_e^35–50, so the rarest high-enhancement sites dominate. Processing that raises site quality beats area. Geometry must place high-quality sites where products can escape to detectors.
5. **Reproduce claimed conditions, in multiple samples.** Aim for ≥ 4–8 independent samples per configuration (ADR-003).
6. **Safety is a hard constraint** (`docs/research/R6-safety-materials-build.md`).
7. **Pre-register** observables, analysis windows and stopping rules before data are taken.

## Deliverables

- `docs/research/R1…R7` — literature base with evidence grades.
- `docs/models/M0…M6` — quantitative models (code in `sim/`).
- `docs/decisions/ADR-*.md` — decision records.
- `docs/design/iteration-1.md` — the first full geometric configuration: drawings (ASCII/figures), dimensions, materials, detectors, controls, protocol, predicted sensitivity, and go/no-go criteria for iteration 2.
- `docs/design/red-team-*.md` — adversarial reviews and responses.
