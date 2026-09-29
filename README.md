# nuclear-fusion — Cold-Fusion Geometry Optimisation

This repository documents a design effort to find the geometric configuration of a **cold** (condensed-matter, room-temperature-class) nuclear experiment with the best chance of producing a **measurable, credible** signal. All reasoning is written down so others can audit it.

## Start here
1. [`docs/00-charter.md`](docs/00-charter.md): mission, scope, honest prior, the two hypotheses.
2. [`docs/01-master-plan.md`](docs/01-master-plan.md): phases, candidate configurations, evaluation matrix.
3. [`docs/models/M0-rate-budget.md`](docs/models/M0-rate-budget.md): first-principles numbers that set the design objective.

## Layout
| Path | Contents |
|---|---|
| `docs/research/` | Literature digests R1–R7 with evidence grades |
| `docs/models/` | Quantitative models M0–M6 (briefs in `docs/models/briefs/`) |
| `docs/decisions/` | Architecture/decision records (ADR) |
| `docs/design/` | Iteration designs and red-team reviews |
| `sim/` | Model code (`python3 sim/<file>.py`; needs numpy, scipy, matplotlib) |
