# Brief M6 — Surface micro-geometry, plasmonics, and mode engineering

Branch: `claude/lucid-davinci-gel4vu-m6`. Read `_common.md` first. Read these literature digests in full:
- `docs/research/R1-electrochemical-PdD.md`: especially the ENEA/Violante surface-PSD correlation and the Letts dual-laser THz results
- `docs/research/R4-theories-geometric-predictions.md`: the geometry-prediction matrix
- `docs/research/R2-gasphase-nano-permeation.md`

## Questions
1. **Surface roughness PSD.** ENEA reported that excess heat correlated with a specific spatial-frequency band in the cathode surface roughness PSD.
   - Reconstruct the claim quantitatively.
   - Compute surface-plasmon-polariton (SPP) dispersion for Pd, PdH/PdD (use the cited dielectric functions), Au on Pd, and the electrolyte/vacuum interface.
   - Determine which roughness wavevectors couple which photon energies. Assess whether the claimed band is physically consistent with plasmon coupling, and state which light sources (if any) would be required.
2. **Design a deterministic surface.** Specify a grating, a 2D lattice of nanoholes/nanopillars, or a random surface with a target PSD that reproduces the claimed band (or the plasmon-optimal band) on the vacuum-facing active face of the detector-facing membrane (C3) and on the C1 cathode. Give feature sizes, depths, periods and tolerances, plus fabrication routes and their cost: electrochemical etching, anodic-alumina templating, focused-ion-beam, interference lithography, laser texturing, sputtered multilayers. Remember the active layer must stay within the charged-particle escape depth; get the escape-depth numbers from M5 on branch `claude/lucid-davinci-gel4vu-m5` if available.
3. **Local field enhancement** at tips and nanogaps under laser illumination, as a factor |E_loc/E_0|. Compute realistic values. Evaluate the honest relevance to (a) nuclear barriers (expected negligible; quantify) and (b) surface chemistry and loading (possibly relevant).
4. **Phonon and mechanical mode engineering.**
   - PdD optical-phonon energies versus loading (cite neutron-scattering data), and comparison with the Letts beat frequencies (~8, 15, 20 THz).
   - Nanoparticle breathing modes versus diameter.
   - Membrane and wire flexural/acoustic modes for the candidate dimensions, and their Q.
   - Whether ultrasonic or megasonic drive, or dual-laser beat stimulation, can be integrated into C3 and C1 without compromising the detectors.
   - Hagelstein-type proposals: which frequencies, and what power density.
5. **Multilayer interfaces.** Iwamura-type Pd/CaO and Ni/Cu multilayers: layer thickness, number of periods, and interface density per cm². Compatibility with charged-particle escape.

## Outputs
- An active-surface specification sheet: roughness PSD target, pattern, multilayer stack, optional stimulation (laser wavelengths/powers/angles, or ultrasound frequency/power), each with its expected value justified by the literature and physics.
- A ranked list of which micro-geometric features are worth including in iteration 1, versus which are unsupported and should be left out.
