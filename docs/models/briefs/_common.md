# Common instructions for every modelling workstream (M1–M6)

You are one of several parallel modelling sessions working for the lead designer of a **cold-fusion (LENR) geometry-optimisation** project. The lead makes every final decision; your job is rigorous, quantitative models that the decision can rest on.

## Read first
1. `docs/00-charter.md` — mission, scope ("cold" only), honesty rules, the H1/H2 hypotheses.
2. `docs/01-master-plan.md` — candidate configurations C1–C6 and the evaluation matrix.
3. `docs/models/M0-rate-budget.md` — first-principles rate budget. Its key conclusions are that the rate is dominated by rare extreme sites, that charged-particle and neutron detection beats heat by ~10¹² for standard D+D, and that geometry must put many candidate sites within the escape depth facing low-background detectors.
4. Anything already present in `docs/research/` (literature digests; some may still be in progress).
5. Your own brief in `docs/models/briefs/`.

## Environment
- `pip install numpy scipy matplotlib sympy` (nothing is preinstalled). 4 CPUs and 15 GB RAM are available. Do not install heavy FEM suites unless essential. Finite-difference/finite-volume with `scipy.sparse` is preferred.
- Web access is available for looking up material constants and cross-sections. Cite each source (URL) next to every parameter value you use.

## Deliverables (only inside your namespace)
- `sim/mX_*.py`: runnable, deterministic scripts (`python3 sim/mX_....py` regenerates every number and figure).
- `docs/models/figs/mX_*.png` plus raw text outputs `docs/models/figs/mX_*.txt`.
- `docs/models/MX-<short-name>.md` with these sections:
  1. Question(s) answered
  2. Model and equations (with units)
  3. Parameters table (value, units, source URL, uncertainty)
  4. Verification (analytic limits, convergence checks)
  5. Results (tables and figures)
  6. **Design recommendations for the lead**: numbered, concrete, with numbers and tolerances (for example "cathode foil thickness 50 ± 10 µm because…")
  7. Sensitivities and uncertainties: which assumptions could flip a recommendation
  8. Open questions / hand-offs to other workstreams
- Do not edit files outside `sim/mX_*`, `docs/models/MX-*`, `docs/models/figs/mX_*`.

## Standards
- Quantitative over qualitative. Every claim is either computed or cited.
- Be honest: if a lever does nothing, say so. If a literature claim is weakly supported, say so.
- Keep the charter's definition of "cold" in mind: no plasma or beam energy sources in the device.

## Git
- Work on the branch named in your session prompt. Commit with clear messages and push with `git push -u origin <branch>`. On a network failure, retry up to 4 times with backoff (2 s, 4 s, 8 s, 16 s).
- Do NOT open pull requests. Do NOT push to any other branch.
- Finish by writing a ≤300-word summary at the top of your MX doc under "Summary for lead", then push.
