# ADR-002 — Platform: electrochemically loaded, detector-facing permeation membrane with a segmented active skin

**Status:** accepted, **revision 2** (dimensions pending M3/M5/M6/M8) · **Date:** 2026-09-29 · **Basis:** R1–R7, M0, ADR-003, red-team-0

> **Revision 2 (after red-team-0 findings F1, F2, F6, F7, F10, F11).** The original text is kept below for the audit trail; these changes override it.
>
> 1. **Exit-face physics (F1).** A clean exit face in vacuum drains the membrane: detailed balance pins it at x ~ 10⁻⁵–10⁻². High loading, high flux and a clean face cannot coexist. Each exit-face sector therefore sets its *own* (x, J) locally, because lateral coupling reaches only about one membrane thickness. That is a feature, provided each sector's (x, J) is modelled and measured rather than assumed shared. Sectors are now defined by **regime**:
>    - **L (loading):** 20–50 nm Au cap. The foil reaches entry-face x ≥ 0.9 and J ≈ 0.
>    - **M (SRI regime):** tuned barrier (thin oxide or a patterned Au grid) giving x ≳ 0.9 *and* J ~ 10¹⁵–10¹⁶ D cm⁻² s⁻¹, i.e. McKubre's (x−x₀)²·flux product.
>    - **F (flux):** open Pd, PdO or Au/Pd/PdO (Lipson) skin. Maximum J at lower x.
>
>    The Au "null" of revision 1 was inverted: Au caps give the *most* highly loaded Pd. The null role now belongs to the H₂O twin and a flux-blocked twin.
> 2. **Thin membrane (F2).** 8–12 µm Pd, so the Si telescopes see protons from the **electrolyte-facing (entry) face**, which is where SRI/ENEA conditions exist. Protons emerge at 2.3–2.55 MeV from the entry face and 3.0 MeV from the exit face, which gives depth tomography. Consequence: Rn-progeny α from the electrolyte penetrate. Mitigate with Rn-free, sparged, sealed electrolyte and ΔE–E particle ID.
> 3. **Pressure-balanced front (F6).** The front (detector) chamber holds **~1 atm D₂**, not vacuum:
>    - Δp ≈ 0, so an unsupported thin membrane survives. Vacuum would give ~250 MPa against a 35–70 MPa annealed yield.
>    - The open exit sectors are clamped at x ≈ 0.65 (β) instead of ~10⁻⁴.
>    - The front becomes a static ⁴He accumulator (aliquots through a getter to the MS).
>    - Charged products lose ≤ 0.3 MeV over ≤ 5 mm of D₂.
>
>    Fast Δp interlock and a Si-bias trip are required.
> 4. **Electrical hygiene (F7):**
>    - the membrane is the detector-system hard ground and the anode is driven;
>    - a grounded ≤ 1 µm Al screen or Ni mesh sits between membrane and telescopes;
>    - a dummy Si detector behind 200 µm Al runs on the same electronics;
>    - waveform digitisation, linear supplies, and gauges/RGA off during counting windows;
>    - a **flux-blocked twin** (identical current waveform, no permeation);
>    - an **orthogonal flux modulation at constant current** (front D₂ pressure steps).
> 5. **Samples (F10).** An **array of ≥ 4 active membrane cells** plus twins, sharing one shielded detector house (neutron bank, 511 keV pair, muon veto). Each cell has its own sealed front volume, so ⁴He can be attributed per sample.
> 6. **Entry-side options (F2, F11):**
>    - an Iwamura-type gas-entry cell (**C3-G**: D₂ 1 atm, 70 °C, Pd/CaO multilayer and a ⁸⁶Sr-tagged target on the *entry* face, detectors facing the entry face through D₂);
>    - Pd/D co-deposition moves to a separate open-cell arm on a PSD scintillator (Cl₂ from a chloride bath would poison the closed-cell recombiner).
> 7. **Tags:** an acoustic-emission sensor on each membrane (the H1b confound); ¹⁰B/⁶LiF depth markers (ADR-004).

## Context: what the literature converges on
Seven independent research digests were compiled with limited primary-source access; tags in each file mark what was checked. They point to the same physical motif from different directions:

| Source | Finding | Grade |
|---|---|---|
| R1 (SRI) | Excess power ∝ (x−x₀)²(i−i₀)·**\|dx/dt\|**. Heat coincides with *changing* loading, i.e. deuterium flux | C |
| R1 (ENEA) | Thin (50 µm), annealed, etched foils reach D/Pd > 0.95 in > 90 % of samples. Surface PSD correlates with heat | C |
| R2 | 5 of 6 flux-type experiments see signals during out-flux or through-flux, not during static loading. Every flux success has a **thick reservoir + ~100 nm nanostructured active skin** | C→B |
| R2 (Iwamura/Toyota) | Pd/CaO multilayer permeation transmutation; Pr yield scales with D flux | B (contested) |
| R3 | Screening is metal- and surface-dependent: PdO ≈ 2× Pd; Au/Pd/PdO heterostructures ×10–50; vacancies raise U_e | B (existence) / C (values) |
| R4 | Theories agree on a **≤ 1 µm, vacancy-rich, highly loaded near-surface layer, driven by non-equilibrium flux** | theory |
| R5 | For charged-particle detection with electrolysis: "make the cathode a 10–25 µm Pd foil window" facing Si detectors | engineering |
| R7 | The 2020–26 credible results (UBC *Nature* 2025, UC Davis/LBNL *Nat. Commun.* 2026) all use a **dual-chamber foil**: electrochemical loading on the back, detectors on the front. The new variable they identify is **D flux** | B/C |

M0 adds the physics. The rate is carried by the rarest high-enhancement sites. Detection beats calorimetry by ~10¹² for conventional products. Active sites must sit within the charged-product escape depth (31 µm for 3.02 MeV p in Pd; 4.6 µm for 1.01 MeV t; R5).

## Decision
The iteration-1 platform is a **detector-facing permeation membrane (DFM)**:

1. **Membrane.** A thin Pd foil forms the wall between two chambers:
   - The **entry side** is a closed electrolytic cell (D₂O/LiOD, internal recombiner). It holds the entry face at the highest achievable deuterium chemical potential (x ≳ 0.9).
   - The **exit side** is a low-background detector chamber.

   Foil thickness is chosen so the whole foil is thinner than the 3 MeV proton range. Then the proton energy measured at the detector encodes the depth at which the reaction happened: entry face, bulk, or exit face. Target: 10–25 µm, set by ADR-004 after M3/M5.
2. **Segmented active skin** on the exit face. There are 4–6 sectors with different nanostructured skins, each viewed by its own collimated Si ΔE–E telescope. Candidate skins:
   - annealed bare Pd (baseline)
   - PdO
   - Au/Pd/PdO heterostructure
   - Pd/CaO multilayer
   - Ni/Cu multilayer
   - Au overlayer (null)

   All sectors share the same flux, temperature and electrical environment. **A sector difference is a built-in control** that common-mode artifacts cannot produce.
3. **Flux is the primary independent variable**, set by cell current and modulated on/off with a period matched to the foil diffusion time. D₂O ↔ H₂O swaps and Pd ↔ Au membranes serve as isotope and host controls.
4. **Detection covers all three hypothesis classes:**
   - **H1 (enhanced conventional D+D).** Si telescopes (p, t, ³He) with depth tomography; a neutron bank; e⁺e⁻/511 keV coincidence and a 3–25 MeV γ window (for the Czerski channel).
   - **H2 (heat + ⁴He).** Operate the exit chamber in a static, NEG-pumped mode. A non-evaporable getter pumps D₂ but not He, so any ⁴He leaving the exit face accumulates and can be counted by RGA/static mass spectrometry. Sensitivity is ~10⁸ atoms, compared with ~10¹⁵ atoms per mJ for calorimetry at 24 MeV/He. The entry-side cell sits in a calorimeter. Post-run, the membrane is melted for retained He.
   - **H3 (transmutation).** Isotope-tagged target deposits (e.g. ⁸⁶Sr → claimed ⁹⁴Mo) on selected sectors, with blinded external SIMS/ICP-MS.

## Rejected as the primary platform
| Config | Reason |
|---|---|
| C1 Fleischmann–Pons/SRI coaxial rod or wire | Charged products cannot cross 0.15 mm of electrolyte. It is blind to H1 and has no internal site-type control. Its calorimetric and ⁴He features are absorbed into the DFM entry cell |
| C2 co-deposition + CR-39 | CR-39 in electrolyte contact is the field's most artifact-prone measurement (R5). No spectroscopy |
| C4 nanopowder beds | Particles are buried in the bed with no line of sight to detectors. H ≈ D behaviour means no D–D ash (R2). Chemical confounders are large |
| C5 Ni/Cu film heat (Tohoku protocol) | Kept as **one skin sector**. The stand-alone 600–900 °C radiometric-calorimetry version is emissivity-limited (R2: Δε = 0.01 ↔ 0.9 W) |
| C6 thermal cycling of Ti chips | Neutrons only, bursty, EMI-prone. Its cracking/fracto-emission element is available in the DFM as a cycling protocol (M3) |
| KeV-beam or plasma-driven foil (UBC/LBNL) | Out of project scope: the charter requires "cold". R3/R4/R7 recommend it because it gives a guaranteed calibrated signal. We note the same membrane platform would accept such a module, but it is not part of this design |

## Consequences / open items
- ADR-003 (exit-chamber environment): vacuum with a perforated support grid, or pressure-balanced gas (D₂ or He at ~1 atm, so ΔP ≈ 0 and an unsupported foil is possible). Trade-offs: foil stress (M3), charged-particle energy loss (M5), ⁴He accumulation (He fill excluded), and hydrogen safety (R6). Decide after M3/M5.
- ADR-004: foil thickness, diameter and support. Decide after M2/M3/M5.
- ADR-005: skin set and sector layout. Decide after M6/M1.
- ADR-006: detector layout and shielding. Decide after M5/M7.
