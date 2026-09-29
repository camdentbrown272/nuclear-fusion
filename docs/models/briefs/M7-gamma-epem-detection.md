# Brief M7 — γ / e⁺e⁻ detection channel, and flux-modulation protocol

Branch: `claude/lucid-davinci-gel4vu-m7`. Read `_common.md` first, then `docs/research/R7-state-of-art-2020-2026.md` in full (especially §4.2 and §6).

## Why this workstream exists
R7 reports a claim by Czerski et al. (Phys. Rev. X 15, 041004, 2025; grade C): below ~5 keV, D–D fusion proceeds mostly through a 0⁺ threshold resonance in ⁴He that decays by **e⁺e⁻ pair emission** (total energy up to ~23.8 MeV). The same group reports 511 keV above background from Pd–D and Zr–D samples with **no beam** (D diffusion only; ICCF-26, grade D).

Our device is cold, so the relative energy is thermal, ~10⁻⁵ of 5 keV. If this channel exists, then in our regime the observable would be **e⁺e⁻ pairs, not protons or neutrons**. Two consequences:
- Annihilation photons escape from the entire active volume, so bulk sites become visible, not only those within the escape depth.
- The detector suite must be designed for this channel explicitly.

## Questions
1. What e⁺e⁻ energy sharing and angular correlation does E0 internal pair creation of a 0⁺ → 0⁺ (g.s.) transition at ~23.8 MeV in ⁴He produce? Give formulas and cite them (Rose/Oppenheimer–Schwinger IPC theory, or E0 pair-emission literature). Then determine, for positrons and electrons born in a 25–300 µm Pd foil with electrolyte on one side and vacuum on the other:
   - (a) where the positrons annihilate: in the foil, in the electrolyte, or in the chamber walls;
   - (b) the bremsstrahlung spectrum;
   - (c) the energy deposited in nearby detectors.
   A simplified MC is fine (continuous slowing down plus a multiple-scattering approximation). Use the ESTAR stopping powers and radiation lengths, and cite them.
2. **511 keV back-to-back coincidence.** Choose between LaBr₃(Ce), NaI(Tl), BGO and HPGe. Specify detector sizes, distance, and Pb/Cu/HDPE shielding, plus a plastic muon veto. Give the coincidence efficiency for annihilation at the foil, and the background coincidence rate at sea level. Background contributions:
   - cosmic-muon pair production in the shielding
   - ⁴⁰K and U/Th chains
   - radon
   - β⁺ emitters from activation
   Cite measured values from low-background labs. Compute the minimum detectable e⁺ production rate (5σ, 30 days) with and without an underground site.
3. **High-energy channel (5–25 MeV)**: a large BGO or NaI "calorimeter" for the pair energy/bremsstrahlung sum, with a cosmic veto. What are the efficiency, the background above 5 MeV, and the MDA?
4. **Flux modulation protocol.** The dual-chamber permeation foil lets us modulate the deuterium flux through the foil (electrolysis current on/off, or D₂O ↔ H₂O). Design the lock-in / on-off analysis:
   - optimal period given the foil diffusion time (τ ≈ L²/D, with L the foil thickness and D the diffusivity of D in Pd; compute it for 25–300 µm) and detector backgrounds
   - how to separate flux-correlated counts from temperature/EMI-correlated artifacts
   - how to combine the channels (p, n, 511 keV, high-E γ) with look-elsewhere control
5. **Calibration sources** that avoid NRC specific licensing where possible (exempt-quantity ²²Na, ¹³⁷Cs, ⁶⁰Co, ²⁴¹Am check sources), and what each calibrates.

## Outputs
- `docs/models/M7-gamma-epem-detection.md`, `sim/m7_*.py`, `docs/models/figs/m7_*`.
- A recommended γ/e± detector layout around a 10–25 mm diameter foil in a small UHV chamber. It must coexist with Si charged-particle telescopes facing the foil at 1–3 cm, and with a neutron moderator/³He bank or EJ-309 cells. Give a concrete arrangement with dimensions and a bill of materials with approximate prices.
- MDA table per channel (sea level vs shallow underground).
- The flux-modulation protocol specification.
