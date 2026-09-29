# 01 — Master Plan

## Phases

| Phase | Output | Executed by | Status |
|---|---|---|---|
| 0. Charter + first-principles budget | `00-charter.md`, `models/M0-rate-budget.md` | lead | done |
| 1. Literature base | `research/R1…R7` | 7 parallel research agents | running |
| 2. Quantitative models | `models/M1…M6`, `sim/` | parallel cloud sessions (one per model) | queued |
| 3. Synthesis | `decisions/ADR-*`, `design/iteration-1.md` | lead | — |
| 4. Red team | `design/red-team-*.md`, revised iteration-1 | independent adversarial agents | — |

## Research streams (Phase 1)

| ID | Topic |
|---|---|
| R1 | Electrochemical Pd–D (Fleischmann–Pons, SRI, ENEA, NRL, SPAWAR, Google 2019) |
| R2 | Gas-phase loading, nanostructures, permeation (Arata, Kitamura/Takahashi, Iwamura/Clean Planet, Mizuno, Celani) |
| R3 | Low-energy nuclear physics in metals (screening, threshold resonance, ARPA-E, NASA LCF) |
| R4 | Theories and their geometric/dimensional/frequency predictions |
| R5 | Detection, backgrounds, statistics, artifact catalogue |
| R6 | Safety, materials, costs, build recipes |
| R7 | State of the art 2020–2026 |

## Candidate configurations under evaluation

These are the families the models and synthesis must rank. The lineage column shows the prior work each one descends from.

| ID | Configuration | Lineage | H1 sensitivity | H2 testability |
|---|---|---|---|---|
| C1 | Closed coaxial electrolytic cell: thin Pd wire cathode, concentric Pt anode, internal recombiner, calorimeter, ⁴He sampling, external neutron bank | Fleischmann–Pons, SRI/McKubre, ENEA | low (charged products absorbed in electrolyte) | high |
| C2 | Pd/D co-deposition on wire/mesh with CR-39 in contact | SPAWAR | medium (CR-39 only, artifact-prone) | low |
| C3 | **Detector-facing membrane (DFM):** thin Pd membrane loaded electrochemically from the back (Devanathan–Stachurski geometry); the front face is interface-engineered (PdO / CaO / multilayer / nanostructure) and faces Si charged-particle spectrometers in vacuum; current modulation drives flux and α/β cycling | Iwamura permeation, Lipson/Kasagi charged-particle work, Devanathan–Stachurski | **high** (spectroscopy, near-zero background) | low–medium |
| C4 | Gas-phase Pd/ZrO₂ or CuNi/ZrO₂ nanocomposite at 200–300 °C, calorimetry | Kitamura/Takahashi/Technova | low | medium |
| C5 | Ni/Cu multilayer thin film on Ni, D₂ gas, heated | Iwamura/Clean Planet | medium if detectors face film | medium |
| C6 | Temperature/pressure cycling of D-loaded Ti or Pd chips with neutron bank | Menlove/Jones, fracto-fusion | medium (neutrons only) | none |

Working hypothesis from M0, still to be tested by models and research: **a hybrid is optimal.** Candidates are C3 as the H1 instrument with C1-style calorimetry/⁴He on the same electrochemical loading circuit, or a single cell that serves both. The synthesis phase decides.

## Evaluation matrix (Phase 3)

Each candidate is scored 0–10 on each criterion. The weights reflect the charter (credible measurability first).

| Criterion | Weight | Measured by |
|---|---|---|
| H1 sensitivity: reactions/s for 5σ, and number of candidate sites visible to detectors | 0.25 | M0, M5 |
| H2 testability: heat resolution, ⁴He capture | 0.15 | M4, R5 |
| Candidate-site multiplicity × loading × non-equilibrium drive | 0.20 | M1, M3, M6, R1–R4 |
| Credibility: artifact exposure, control quality | 0.20 | R5, M5 |
| Evidence base: prior positive reports and their grades | 0.10 | R1, R2, R7 |
| Safety, cost, buildability | 0.10 | R6 |

## Modelling workstreams (Phase 2)

Briefs are in `docs/models/briefs/`. Each cloud session writes `docs/models/Mx-*.md`, `sim/mx_*.py` and `docs/models/figs/mx_*`, and pushes to branch `claude/lucid-davinci-gel4vu-mx`. The lead then merges.

| ID | Model |
|---|---|
| M1 | Site physics: D–D separations/enhancement at defect classes, crack-tip fields and fracto-emission energies, per-site-class yields |
| M2 | Electrochemistry: current distribution, loading vs overpotential, electrode geometry for uniform high loading |
| M3 | Hydrogen transport and mechanics: diffusion/permeation with α/β transition, loading time vs thickness, stress and cracking, cycling protocol |
| M4 | Calorimetry and thermal design |
| M5 | Detection: charged-particle transport and escape, Si/CR-39 geometry, neutron efficiency, backgrounds, run-time/power analysis |
| M6 | Surface micro-geometry: plasmon/roughness PSD, nanostructure field enhancement, mechanical/phonon mode engineering, laser coupling |
