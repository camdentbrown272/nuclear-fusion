# 00 — Project Charter

**Owner of all design decisions:** Claude (lead designer). Every decision is recorded in `docs/decisions/` with its reasoning so any reviewer can audit it.

## Mission

Design the geometric configuration of a **cold** condensed-matter nuclear experiment that maximises the probability of a **measurable, credible** nuclear signal.

- **Cold** means: bulk apparatus temperature ≤ ~600 K; no plasma confinement, no fusor, no ion beam or accelerator as the energy source, no thermonuclear conditions. Allowed drivers: electrochemistry, gas pressure/permeation, temperature cycling, current/field modulation, low-power optical/acoustic stimulation, magnetic fields.
- **Measurable** means: a pre-registered observable exceeds background by ≥5σ, with a matched control that does not show it.
- **Credible** means: the result survives the artifact catalogue in `docs/research/R5-detection-and-artifacts.md` and could be handed to a skeptical nuclear physicist without embarrassment.

## Honest prior (stated once, then we work)

1. No cold-fusion/LENR claim since 1989 has met the reproducibility standard of mainstream physics. Several careful programmes (NHE Japan 1992–98, Google/Nature 2019, most of ARPA-E 2022–25) reported nulls or ambiguities.
2. First-principles numbers (`docs/models/M0-rate-budget.md`) show that ordinary physics falls short of a neutron-detectable D–D rate by roughly **100 e-folds of tunnelling exponent** (~43 decades). Geometry and engineering can buy perhaps 25–30 e-folds (more sites, better detectors, longer runs). **Geometry alone cannot close the gap.** A positive result requires an anomalous enhancement that is not in standard physics.
3. Therefore the only rational definition of "best geometry" is: **the configuration that places the largest number of candidate-anomaly sites (as located by the claims and theories with the most support) in front of the most sensitive, least artifact-prone detectors, with built-in controls.** Such a design also produces the most informative null if nothing is there.

This framing is not a hedge; it is the optimisation target. Designs that maximise "excess heat" alone are rejected because heat is ~10¹² times less sensitive than particle detection for standard D–D and is the most artifact-prone observable in the field's history.

## Two hypotheses the design must serve simultaneously

| | H1: enhanced conventional D–D | H2: anomalous "LENR" channel |
|---|---|---|
| Products | n (2.45 MeV) + ³He; p (3.02 MeV) + T | ⁴He + heat, no strong γ/n |
| Most sensitive signature | charged-particle spectroscopy, neutrons | heat correlated with ⁴He |
| Reactions/s for 5σ | ~10⁻³–10⁻² s⁻¹ | ~10¹⁰ s⁻¹ per 10 mW |
| Mainstream compatibility | Compatible in kind (screening/threshold resonance) | Requires new physics |

A single geometry that is blind to either hypothesis is inferior to one that tests both.

## Principles

1. **Quantify before choosing.** Every geometric choice cites a model in `docs/models/` or a source in `docs/research/`.
2. **Controls are geometry.** The H₂O/light-hydrogen twin, the Pt-cathode twin, and the detector blank are part of the configuration, not an afterthought.
3. **Energy spectra over totals.** A 3.02 MeV proton peak is worth more than any calorimetric excess.
4. **Tail, not mean.** Rates scale as ~Ue³⁰⁻⁵⁰; the rarest high-enhancement sites dominate. Geometry must multiply candidate sites and put them where products can escape to detectors.
5. **Safety is a hard constraint** (`docs/research/R6-safety-materials-build.md`).
6. **Pre-register** observables, analysis windows and stopping rules before data are taken.

## Deliverables

- `docs/research/R1…R7` — literature base with evidence grades.
- `docs/models/M0…M6` — quantitative models (code in `sim/`).
- `docs/decisions/ADR-*.md` — decision records.
- `docs/design/iteration-1.md` — the first full geometric configuration: drawings (ASCII/figures), dimensions, materials, detectors, controls, protocol, predicted sensitivity, and go/no-go criteria for iteration 2.
- `docs/design/red-team-*.md` — adversarial reviews and responses.
