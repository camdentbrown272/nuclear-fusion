# ADR-003 — Objective function v2: claim-conditioned, multi-sample, ⁴He co-primary

**Status:** accepted · **Supersedes:** ADR-001 · **Date:** 2026-09-29 · **Basis:** red-team-0 §2 (findings F3, F4, F5, F10), revised M0

## Why ADR-001 was wrong
ADR-001 scored a geometry by the minimal localised enhancement it could detect. The revised M0 shows that this lever is weak:
- **Sensitivity is compressed by the exponential physics.** d ln λ / d ln U = 35–50, so 100× better sensitivity buys only 10–14 % in reachable U_eff.
- **Its value under an honest prior is small.** With the anomaly strength log-uniform over 40 decades, 100× better sensitivity adds ~5 % detection probability.
- **The order-unity levers were ignored:**
  - whether the configuration reproduces the conditions under which anomalies are claimed;
  - how many independent samples are run (with a per-sample success rate p = 0.2, 1 sample gives 0.20 and 8 give 0.83);
  - whether H2 (lattice-coupled ⁴He, no radiation) is visible at all. ADR-001 left H2 to calorimetry, which is ~10⁶× worse than ⁴He collection.

## Decision
A design d is scored by

  Score(d) = Σ_c w_c · P(cond_c | d) · [1 − (1 − p_c)^{N_samples}] · min(1, P_det(m_c/100 | d)) · Cred(d)  +  0.1 · S_generic(d)

where:
- c ranges over the claims (table below), each with an evidence weight w_c (B = 3, C = 2, D = 1);
- P(cond_c | d) is the probability that design d reproduces claim c's conditions **at a detector-visible location**;
- p_c is the claimants' per-sample success rate (default 0.2 where not stated);
- P_det is capped at 1 once the design reaches 1 % of the claimed magnitude m_c, because more sensitivity beyond that buys nothing for testing claims;
- Cred(d) is the credibility factor (artifact exposure, control quality);
- S_generic is the ADR-001 localised-anomaly sensitivity. It is kept at small weight for anomalies no one has claimed.

### Hypothesis set (every design must state what it sees for each)
| ID | Hypothesis | Signature | Primary detector |
|---|---|---|---|
| H1a | Enhanced screening, standard branching | p 3.02, t 1.01, ³He 0.82, n 2.45 MeV | Si ΔE–E telescopes; neutron bank |
| H1b | Crack- or transient-driven conventional fusion (fracto-emission) | Same products, **correlated with cracking** | Acoustic-emission sensor as a tag. This is a *confound* for H1a |
| H1c | e⁺e⁻ threshold-resonance channel (Czerski) | 511–511 keV coincidences, e± up to ~23 MeV | Back-to-back LaBr₃/NaI pair; high-energy γ window |
| H2 | Lattice-coupled D+D→⁴He | ⁴He at 2.6×10¹¹ J⁻¹; heat; possibly ≤ 20 keV particles | **⁴He**: static accumulation plus foil extraction. Calorimetry second |
| H3 | Transmutation (Iwamura), low-Z "fission daughters" (MIT/ARPA-E focus) | Isotopically anomalous surface species at 10¹²–10¹⁴ cm⁻² | Pre/post ToF-SIMS + HR-ICP-MS, with **isotope-tagged targets** (⁸⁶Sr → ⁹⁴Mo) |
| H4 | Reactions involving light hydrogen (p+d, Ni–H) | ³He; heat in H₂ | ³He mass spectrometry. The H₂O twin is a control **only under the stated assumption** that H is inert in Pd–D |

### Claim register (the c in the score)
| Claim | Grade | Conditions that must be reproduced at a visible location | Claimed magnitude m_c |
|---|---|---|---|
| SRI / ENEA electrolytic heat + ⁴He | B−/C | Bulk x ≥ 0.9 (resistance), i ≥ 0.1–0.4 A cm⁻², weeks, D flux (\|dx/dt\|), ENEA-processed Pd | 0.1–1 W → 10¹⁰–10¹¹ He s⁻¹ |
| Iwamura Pd/CaO permeation transmutation | B (contested) | D₂ gas, 1 atm, ~70 °C, multilayer on the **entry** face, J = 3–7×10¹⁷ D cm⁻² s⁻¹ | 10¹²–10¹⁴ product atoms cm⁻² |
| Kitamura/Takahashi CNZ/PNZ nanocomposites | B (in-network) | Ni-based ZrO₂ composites, 200–300 °C, H₂/D₂ | 5–25 W per 100–200 g |
| Tohoku/Clean Planet Ni/Cu multilayer | C | H₂ loaded at ~250 °C, then desorbed while heating to 500–900 °C | 1–6 W from 6 cm² |
| SPAWAR co-deposition | C/D | Pd from PdCl₂/LiCl co-deposited with D on Au/Ag | ~10³–10⁴ tracks cm⁻² over days |
| Lipson Pd/PdO:D desorption | C/D | Electrolytic loading, then desorption facing Si/CR-39 | ~10⁻²–1 s⁻¹ charged particles |
| Czerski beam-free 511 keV | D | D-loaded Pd/Zr, no beam, low background | not stated; take ~10⁻³ s⁻¹ |
| NTT out-diffusion bursts | D | Au on one face, out-diffusion into vacuum | bursts |

### Evaluation weights (master-plan matrix, revised)
| Criterion | Weight |
|---|---|
| Claim-condition reproduction at a visible location | 0.30 |
| H2 via ⁴He (and calorimetry) | 0.20 |
| Credibility (artifacts, controls) | 0.20 |
| H1 sensitivity (saturating at 1 % of the claimed magnitude) | 0.10 |
| Independent samples / statistics | 0.10 |
| Safety, cost, buildability | 0.10 |

## Consequences
- ⁴He is **co-primary** in every configuration. Foil extraction is pre-registered, with an H₂O-twin foil as the blank.
- Every configuration runs **≥ 4 independent samples**, with pre-registered per-sample denominators.
- The C3 platform is re-specified (ADR-002 rev. 2) to put claimed conditions in view of its detectors.
- Parallel low-cost arms are added where they raise Σ_c materially. See the iteration-1 design.
