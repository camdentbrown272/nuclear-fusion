# M6: Surface micro-geometry, plasmonics, mode engineering and the segmented active face

Code: `sim/m6_*.py` (`python3 sim/m6_run_all.py` regenerates everything; `m6_gratings.py` takes about 40 min on 4 CPUs). Raw outputs: `docs/models/figs/m6_*.txt`. Optical data: `sim/m6_data/` (refractiveindex.info, CC0).

## Summary for lead

1. **The exit skin controls loading, and a bare exit face is unloaded (most important finding).** Any exit face with D₂ sticking s ≳ 10⁻⁴ (Pd, Pd/CaO, Ni/Cu) holds the vacuum face at x ≈ 10⁻⁵–10⁻³ at every achievable flux. The 50 µm membrane drains unless ~1 A cm⁻² of D is absorbed. As drawn, C3 puts α-phase Pd inside the 7 µm proton escape depth. Two fixes: a cap with s ≲ 10⁻⁷ (x ≥ 0.6), or 0.1 bar D₂ backfill on the detector side (x ≈ 0.61, ≤25 keV triton loss).
2. **Segmented face: geometry works, but sectors do not share flux.** High-s skins drain; Au and PdO sit at x ≈ 0.6. Sector comparisons therefore confound chemistry with loading and flux. Ø6 mm patches, ≥0.5 mm gaps and Ø6 mm detectors at 10 mm (knife-edge ≤1 mm) give ≤1 % cross-talk and 2 % efficiency, about 15× below one large detector.
3. **Plasmonics.** ENEA's band contains the SPP-coupling frequencies. But Pd SPPs are broad (Q 5–30), so any roughness couples, and dark cells have no SPPs (n < 10⁻¹⁶). A PSD correlation cannot be plasmonic. Optimal gratings (785 nm: Λ = 562/85 nm in D₂O, 761/101 nm in vacuum) absorb ~100 %. Barrier effect ~10⁻⁵. The chemistry effect is a real but modest photothermal flux modulation.
4. **THz.** A beat cannot drive the bulk Γ mode (Raman-forbidden) or L-point modes (momentum), and 30 mW gives ~10⁻⁹ of thermal. 15.3 THz lies above the PdD band.
5. **Iteration 1.** Include the loading valve, skins (a), (f), (b) and (d), AFM PSD as a covariate, and a support grid (span ≤1.8 mm). The laser serves only as a flux modulator. Exclude THz beats, breathing modes, hot spots and ultrasound.

---

## 1. Questions answered

| # | Question (brief) | Where answered |
|---|---|---|
| Q1 | ENEA roughness PSD; SPP dispersion for Pd, PdD, Au-on-Pd, electrolyte and vacuum; which roughness wavevectors couple which photons; consistency; light sources | §5.1, §5.2 |
| Q2 | A deterministic surface (grating, lattice or random PSD) with dimensions, tolerances and fabrication routes and costs | §5.3, §6 spec sheet |
| Q3 | Local field enhancement at tips and gaps; relevance to (a) barrier and (b) chemistry/loading | §5.4 |
| Q4 | PdD optical phonons vs Letts lines; breathing modes; membrane and wire modes and Q; ultrasound or beat integration; Hagelstein frequencies and power | §5.5 |
| Q5 | Iwamura Pd/CaO and Ni/Cu stacks: thickness, periods, interface density, escape compatibility | §5.6 |
| Q6 | Segmented active face: per-skin thickness and fabrication, k_r ordering, energy loss, cross-contamination, sector size, gap | §5.7, §5.8 |

## 2. Model and equations (units SI unless stated; e^{−iωt}, Im ε > 0 = loss)

### 2.1 SPP dispersion (`m6_spp.py`)
- Single interface: k_spp = k₀ √(ε_m ε_d/(ε_m + ε_d)).
- Propagation length: L = 1/(2 Im k_spp).
- Field decay into the dielectric: δ_d = 1/Im k_z,d, with k_z,d = √(ε_d k₀² − k_spp²).
- Three-layer mode (D₂O / Au film t / PdD): zero of (p₁+p₂)(p₂+p₃) + (p₁−p₂)(p₂−p₃) e^{2ik_{z2}t} = 0, with p_i = k_{z,i}/ε_i. Solved by Newton continuation in t.
- First-order grating coupling: G = 2π f = Re k_spp − n_d k₀ sinθ. At some angle this is possible for G ∈ [Re k_spp − n_d k₀, Re k_spp + n_d k₀].
- Thermal SPP occupation: n = 1/(e^{E/kT} − 1).

### 2.2 PdD dielectric function
- ε_PdD(λ) = F·ε_Pd(λ), with F = 0.75 (range 0.6–0.9) **[BK]**. This is a reading of von Rottkay et al. (1999) and Palm et al. (2018), who report a 20–40 % drop of |ε| on hydriding. The retrievable Palm file covers Pd only.
- H and D share ε at equal x, because electronic structure is isotope-independent to <10⁻³.

### 2.3 RCWA for 1-D gratings, TM (`m6_rcwa.py`)
- With u = η₀H_y and z′ = k₀z: dU/dz′ = iA⁻¹S_x and dS_x/dz′ = i(I − K E⁻¹ K)U.
- Hence d²U/dz′² = A⁻¹(K E⁻¹ K − I)U, where E = Toeplitz[ε], A = Toeplitz[1/ε] (Li's inverse rule) and K = diag(k_{x,n}/k₀).
- A top-down admittance recursion (S_x = Y U) handles the staircase sinusoid. Absorptance = 1 − Σ R_n (opaque substrate).
- Near-field intensity |E|²/|E_inc|² is evaluated on the plane touching the crests. This is a lower bound; the crest-plane field carries Gibbs ripple of about ±30 %.

### 2.4 Field enhancement (`m6_fields.py`, quasi-static, feature ≪ λ/2πn)
- Spheroid asperity (a hemispheroid on a plane equals a full spheroid by image symmetry): E_tip/E₀ = (ε/ε_d)/(1 + L(ε/ε_d − 1)), where L is the long-axis depolarisation factor.
- Sphere dimer: exact axisymmetric multipole solution. Response coefficients are B_l = C_l R^{2l+1} l(ε_d − ε)/(lε + (l+1)ε_d), coupled through the axial translation theorem P_l(cosθ_B)/r_B^{l+1} = Σ_n (∓1)^{l or n} C(l+n, n) r_A^n P_n(cosθ_A)/D^{l+n+1}.
- Barrier relevance, from M0: d ln λ/dU = ½√(E_G/U)/U, with E_G = 985.8 keV.
- Ponderomotive energy: U_p = e²E²/(4mω²).

### 2.5 Phonons and mechanics (`m6_modes.py`)
- Rocksalt Born–von Kármán model with central nearest-neighbour Pd–Pd, Pd–D and D–D force constants, fitted to three [BK] inelastic-neutron-scattering anchors. DOS on a 24³ Brillouin-zone grid.
- Lamb l = 0 breathing mode: ξ cot ξ = 1 − ξ²/(4β²), with β = c_T/c_L and f = ξc_L/(πd).
- Clamped plate: f₀₁ = (10.2158/2πa²)√(D/ρh), D = Eh³/12(1−ν²).
- Fluid loading on one side: f/√(1+β_f), with β_f = 0.6689 ρ_f a/(ρh).
- Thickness mode: c_L/2h.
- Cylinder radial mode: x J₀(x)/J₁(x) = 2c_T²/c_L².
- Large-deflection membrane under Δp: w₀ = 0.662a(pa/Eh)^{1/3} and σ = 0.423(Ep²a²/h²)^{1/3}.
- Acoustic intensity of a travelling wave: I = ½ρc_L³ε².
- Coherent optical-phonon cost: P = U ω/Q, with U = ½N_D m_D ω² u².

### 2.6 Stopping and escape (`m6_skins.py`)
- Elemental stopping powers come from pycatima (ATIMA with SRIM-type low-energy stopping), combined by Bragg additivity. CSDA ranges are integrated from these.
- Escape into a 2π hemisphere from depth z with E > E_thr: P(z) = ½(1 − z/R_thr).
- Escape-equivalent thickness: ∫P dz = R_thr/4.

### 2.7 Exit-face loading (detailed balance)
- The desorption flux into vacuum equals the absorption flux that would hold equilibrium with D₂ at fugacity f(x): J = 2 s Z₁ f(x).
- Z₁ = 7.58×10²³ D₂ cm⁻² s⁻¹ bar⁻¹ at 300 K.
- s is the effective sticking (= recombination) probability of the skin.
- f(x) is a piecewise-log isotherm for PdD at 300 K **[BK]**:
  - plateau 0.04 bar;
  - x = 0.65 at 1 bar, 0.70 at 10 bar, 0.80 at 3×10² bar, 0.90 at 10⁴ bar, 0.95 at 10⁵ bar.

### 2.8 Sector flux partition (`m6_sectors.py`)
- Steady ∇·(D n_Pd ∇x) = 0 in (lateral, depth).
- Entry boundary: J_in = J₀(1 − x/0.95), with J₀ = j_abs/e.
- Exit boundary: J_out = 2 s_i Z₁ f(x).
- Newton on a sparse finite-volume grid.
- Linear Fickian diffusion, with no α/β interface dynamics.

### 2.9 Collimation
- Point-source acceptance through a coaxial knife-edge aperture (height z_ap, radius r_ap) onto a detector disc (radius r_d, distance D), by polar quadrature.
- Cross-talk: counts from 6 hexagonal neighbour patches divided by own-patch counts, for uniform emission.

## 3. Parameters

| Parameter | Value | Units | Source | Uncertainty |
|---|---|---|---|---|
| ε(Pd) 0.19–1.94 µm | tabulated | – | Johnson & Christy, PRB 9, 5056 (1974); file https://github.com/polyanskiy/refractiveindex.info-database (data/main/Pd/nk/Johnson.yml), fetched this session | ±5–10 % |
| ε(Pd) cross-check | tabulated | – | Palm et al., ACS Photonics 5, 4677 (2018), https://doi.org/10.1021/acsphotonics.8b01243 (same repository) | Palm's ε″ is ~20 % higher than J&C |
| ε(Pd) > 1.9 µm | Lorentz–Drude | – | Rakić et al., Appl. Opt. 37, 5271 (1998), https://doi.org/10.1364/AO.37.005271 | ±20 % |
| ε(PdD) | 0.75 × ε(Pd) | – | **[BK]** von Rottkay, Rubin & Duine, JAP 85, 408 (1999); Palm 2018 | F = 0.6–0.9 |
| ε(Au) | tabulated | – | Johnson & Christy, PRB 6, 4370 (1972); Olmon et al., PRB 86, 235147 (2012) | ±5 % |
| ε(Ni), ε(Cu), ε(Ag) | tabulated | – | Johnson & Christy (1972, 1974) | ±10 % |
| n(D₂O) 0.5–1.6 µm | Sellmeier + k | – | Kedenburg et al., Opt. Mater. Express 2, 1588 (2012), https://doi.org/10.1364/OME.2.001588 | ±0.001 |
| n(H₂O) (used for D₂O outside 0.5–1.6 µm, −0.005) | tabulated | – | Hale & Querry, Appl. Opt. 12, 555 (1973) | ±0.005 (step visible at 0.78 eV in the light-line plot) |
| Stopping powers | ATIMA/SRIM-type | MeV cm² g⁻¹ | pycatima 1.98, https://github.com/hrosiak/catima | ±3–5 % |
| PdD lattice constant | 4.03 | Å | R4 §4 table **[BK]** | ±0.01 |
| Pd LA(X) | 6.7 | THz | **[BK]** Miiller & Brockhouse, Can. J. Phys. 49, 704 (1971) | ±0.2 |
| PdD₀.₆₃ optical band | 8.0 (Γ) – 11.5 (top) | THz | **[BK]** Rowe et al., PRL 29, 1250 (1972), https://doi.org/10.1103/PhysRevLett.29.1250; Ross et al., JPCM 10, 3219 (1998) | ±0.7 each |
| PdH/PdD optical-peak ratio | 1.5 (anharmonic; harmonic would be 1.41) | – | **[BK]** Ross 1998 | 1.41–1.55 |
| PdHₓ optical peak vs x | ~68 meV (α) → 56 meV (β), ~flat to x ≈ 1 | meV | **[BK]** Ross 1998 and earlier INS | ±3 meV |
| Letts beat lines | 8.3 (0.70), 15.3 (0.44), 20.4 (0.68) | THz (width) | R1, R4 † | as reported |
| Pd E, ν, ρ | 121 GPa, 0.39, 12.02 g cm⁻³ | – | CRC Handbook **[BK]** | ±5 % |
| PdD E, ρ | 109 GPa (−10 %), 11.07 g cm⁻³ | – | **[BK]**; ρ from lattice expansion | ±10 % |
| Au, Ni E, ν, ρ | 79/200 GPa; 0.44/0.31; 19.3/8.91 | – | CRC **[BK]** | ±5 % |
| Annealed Pd yield | ~50 | MPa | **[BK]** | 35–70 |
| D diffusivity in PdD | 5×10⁻⁷ | cm² s⁻¹ | R1 §3.2 | ×2 |
| D₂ sticking s: clean Pd | 0.1 | – | **[BK]** Conrad, Ertl & Latta, Surf. Sci. 41, 435 (1974) | 0.01–0.5 |
| s: contaminated Pd | 10⁻³ | – | **[BK]** assumption | 10⁻⁴–10⁻² |
| s: Ni | 10⁻² | – | **[BK]** | 10⁻³–0.1 |
| s: Au 20 nm on Pd | 10⁻⁷ | – | **[BK]** (bulk Au ≪ 10⁻¹⁵; thin films leak through pinholes and intermixing) | 10⁻¹⁰–10⁻⁵ |
| s: intact PdO | 10⁻⁶ | – | **[BK]** assumption | 10⁻⁹–10⁻⁴ |
| PdD isotherm anchors | see §2.7 | bar | **[BK]** high-pressure PdH/PdD data (Baranowski), electrochemical fugacity | ×10 per anchor |
| D₂ desorption activation energy (Pd) | 0.5 | eV | **[BK]** Behm, Christmann & Ertl, Surf. Sci. 99, 320 (1980) | 0.4–0.9 |
| Hot-carrier quantum yield | 10⁻⁴–10⁻² | – | **[BK]** Mukherjee et al., Nano Lett. 13, 240 (2013); Zhou et al., Science 362, 69 (2018) | – |
| Tunnelling limit of nanogaps | ~0.3–0.5 | nm | Savage et al., Nature 491, 574 (2012), https://doi.org/10.1038/nature11653 **[BK]** | – |
| Absorbed entry current j_abs (C3) | 0.05 | A cm⁻² | assumption (M2/M3 to supply) | 0.01–0.5 |

## 4. Verification

| Check | Result |
|---|---|
| RCWA, flat limit vs Fresnel (PdD/D₂O, 633 nm, 0/30/60°) | |r_p|² agrees to 6×10⁻⁷ |
| RCWA, lossless metal (ε = −15): energy conservation | Σ R = 1 to 10⁻¹⁴ for N = 10–40 |
| RCWA harmonic convergence (PdD/D₂O, 785 nm) | A = 0.668 (N = 10), 0.696 (30), 0.701 (40), 0.701 (60): 0.7 % at N = 30 |
| RCWA convergence (Au 20 nm/PdD) | 0.83–0.86 over N = 20–60: ±3 % (TM-metal staircase; quoted optima carry ±3 % in A and ±5 nm in period) |
| RCWA slices 8/16/32 | 0.694/0.696/0.691 |
| Au/air grating benchmark (Λ = 600 nm, h = 50 nm) | absorption peak at 630 nm, as expected from Λ·Re√(ε/(1+ε)) ≈ 615–630 nm |
| Three-layer SPP solver: 400 nm Au on PdD → Au/D₂O | identical to 6 digits |
| Three-layer SPP solver: film = substrate → PdD/D₂O | identical to 6 digits |
| Dimer translation theorem at random points | agrees to 8 digits |
| Dimer, gap = 100 R | 4.381 = single-sphere |3ε/(ε+2)| |
| Dimer multipole convergence (Au, g = 1 nm) | l_max = 50…300 gives 69.60 (converged) |
| Collimation quadrature vs analytic on-axis solid angle | 0.005532 vs 0.005532 |
| 2-D sector solver vs 1-D analytic (s = 0.1 and 10⁻⁷) | x_out equal to 4–5 digits |
| pycatima: Bragg-combined range vs direct Pd range (3.02 MeV p) | 31.20 vs 31.29 µm |
| PSD synthesis: realised rms vs target, Parseval | exact (30.00/36.00/13.00 nm) |

## 5. Results

### 5.1 SPP dispersion (Q1): `figs/m6_spp_dispersion.png`, `figs/m6_eps.png`, `figs/m6_spp.txt`

| Interface, 785 nm | n_eff | L (µm) | δ_d (nm) | Λ₀ (nm) | Q |
|---|---|---|---|---|---|
| Pd / vacuum | 1.0124 | 5.1 | 722 | 775 | 42 |
| PdD / vacuum | 1.0166 | 3.8 | 623 | 772 | 31 |
| Pd / D₂O | 1.3531 | 2.1 | 409 | 580 | 23 |
| PdD / D₂O | 1.3627 | 1.6 | 353 | 576 | 17 |
| Au 20 nm on PdD / D₂O | 1.3767 | 5.4 | 330 | 570 | 59 |
| Au (thick) / D₂O | 1.3779 | 17.6 | 328 | 570 | 194 |
| Ni / vacuum | 1.0116 | 4.5 | 724 | 776 | 37 |

- At 633 nm, PdD/D₂O has Λ₀ = 459 nm (f₀ = 2.18×10⁶ m⁻¹) and Q = 12. At 1064 nm it has Λ₀ = 790 nm and Q = 27.
- **Pd and PdD SPPs sit within 1–4 % of the light line and are strongly damped (Q ≈ 5–30).** A 20 nm Au overlayer triples Q.

### 5.2 ENEA band: reconstruction and consistency (Q1)

The digests give the band only as "micron to sub-micron, ~10⁵–10⁷ m⁻¹, not verified" (R1 G4). Both readings were tested:

- **Reading A, f = 1/Λ = 10⁵–10⁷ m⁻¹ (Λ = 0.1–10 µm).**
  - First-order coupling at normal incidence covers every photon from 0.66 to 3.4 eV, i.e. f = 0.7–4×10⁶ m⁻¹.
  - Above f ≈ 4×10⁶ m⁻¹ (Λ < 250 nm) no visible or NIR photon couples. Above ~5 eV no bound Pd SPP exists. The top half-decade of the band is therefore not plasmonic.
- **Reading B, q = 10⁵–10⁷ rad m⁻¹.** Normal-incidence coupling only works below 1.45 eV (λ > 855 nm) in D₂O.
- **Selectivity.** At grazing incidence, f as low as 2–12×10⁴ m⁻¹ couples, because (n_eff − n_d)/λ is tiny. The SPP linewidth in f is 0.5–3×10⁵ m⁻¹. So essentially any roughness between 10⁵ and 3×10⁶ m⁻¹ couples some visible or NIR photon at some angle. **A plasmonic mechanism cannot select a narrow band on Pd.**
- **Dark cells.** The thermal SPP occupation is 1.6×10⁻¹⁷ at 1 eV and 10⁻³³ at 1.96 eV. Without external light no SPP exists to couple. A PSD–heat correlation in a dark electrolytic cell, if real, is a proxy for metallurgy: grain size, texture, etch pits and dislocation outcrops, which ENEA's own recipe also controls.
- **Light sources required for any plasmonic role.**
  - An external p-polarised CW laser at 633–1064 nm, with the grating vector in the plane of incidence, through a window. The electrolyte path must stay ≤ a few mm at λ > 1.3 µm because of D₂O absorption.
  - Or, for the C3 vacuum face, light from the vacuum side.
  - Dual-laser pairs for the Letts beats: 670/682.7, 670/693.7 and 670/702.0 nm (8.3/15.3/20.4 THz); 785/802.4, 785/817.8 and 785/829.3 nm.

### 5.3 Deterministic gratings (Q2): `figs/m6_grating_maps.png`, `figs/m6_grating_spectra.png`, `figs/m6_gratings.txt`

Sinusoidal 1-D gratings, TM, normal incidence, optimised for maximum absorptance:

| Surface | λ (nm) | A_flat | Λ_opt (nm) | h_opt p-v (nm) | A_max | Λ window (half gain) | h window (A ≥ 0.9 A_max) | |E|²/|E_inc|² at crests |
|---|---|---|---|---|---|---|---|---|
| PdD / D₂O (C1) | 633 | 0.37 | 439 | 76 | 1.00 | 312–474 | 53–110 | 8.5 |
| PdD / D₂O (C1) | 785 | 0.32 | 562 | 85 | 1.00 | 480–589 | 60–124 | 11 |
| PdD / D₂O (C1) | 1064 | 0.28 | 779 | 106 | 1.00 | 715–800 | 74–148 | 16 |
| PdD / vacuum (C3) | 785 | 0.26 | 761 | 101 | 1.00 | 706–781 | 71–142 | 17 |
| PdD / vacuum (C3) | 1064 | 0.22 | 1045 | 125 | 0.995 | 997–1061 | 87–175 | 22 |
| Au 20 nm / PdD, D₂O | 785 | 0.10 | 566 | 38 | 1.00 | 559–570 | 30–49 | 90 |
| Au thick / D₂O | 785 | 0.03 | 568 | 27 | 0.995 | < ±4 | 22–37 | 185 |

- Angular acceptance (half gain) at 785 nm: ±4° for PdD and ±0.5° for Au/PdD.
- A PdD grating is a forgiving near-perfect absorber: ±10 % period and ±40 % depth. An Au overlayer gives about 8× more near-field intensity but needs ±1 % period and ±0.5° alignment.
- A random surface with an isotropic "ring" PSD (`figs/m6_psd.png`) centred at f₀ = 1/Λ_opt, with σ_f ≈ 5–8 %, couples at every azimuth.
  - A ring rms of about h/2 (30 nm for C1, 36 nm for C3) is equivalent to the optimal sinusoid.
  - A generic etched surface with 50 nm rms carries only 4.8 % of its variance in the 785 nm ring, i.e. an effective 11 nm rms.

### 5.4 Local fields (Q3): `figs/m6_fields.png`, `figs/m6_fields.txt`

| Structure (785 nm unless noted) | |E_loc/E₀| |
|---|---|
| PdD spheroid asperity, aspect 1/3/5/10, in D₂O | 3.4 / 13 / 19 / 15 |
| PdD asperity in vacuum, aspect 10 | 30 |
| Au asperity in D₂O, aspect 5 | 59 |
| Au asperity in vacuum, 1064 nm, aspect 10 | 660 (quasi-static resonance; retardation limits to ~100) |
| PdD dimer, R = 20 nm, gap 0.5/1/2/5 nm (D₂O) | 83 / 40 / 20 / 8 |
| Au dimer, R = 20 nm, gap 1 nm (D₂O, 633 / 785 nm) | 363 / 61 |
| Grating crests (from §5.3), |E| | PdD 3–5; Au 9–14 |

Realistic ceilings (tunnelling below ~0.5 nm) are |E/E₀| ≈ 10–30 on Pd/PdD and 100–200 on Au.

**(a) Nuclear barrier**
- For 30 mW on 1 mm², E₀ = 4.75×10³ V/m, so E_loc ≤ 5×10⁵ V/m. That is 8×10⁻¹⁵ of the deuteron field at 5 fm.
- The energy across 0.74 Å is 3.5×10⁻⁵ eV. The relative (gradient) shift on a D₂ pair at a 10 nm tip is 2.6×10⁻⁷ eV. U_p(e) is 1.7×10⁻⁹ eV.
- With M0's slope this multiplies the rate by exp(≤2.5×10⁻⁵).
- Even at the ps ablation threshold (out of the "cold" scope) the pair shift is ≤0.05 eV.
- **Verdict: irrelevant**, as R4 predicted.

**(b) Chemistry and loading**
- 30–300 mW over 0.5 cm² is 2.4×10¹⁷–2.4×10¹⁸ photons cm⁻² s⁻¹, comparable to the electron flux of 0.1 A cm⁻².
- With a hot-carrier quantum yield of 10⁻⁴–10⁻², photo-driven D events are 0.2 %–240 % of a 10¹⁶ cm⁻² s⁻¹ permeation flux.
- Photothermal heating is 0.07–5 K. Recombination sensitivity is 6.4–11.6 %/K.
- **Verdict: a laser is a credible, contact-free modulator of exit-face recombination and flux (percent to tens of percent), mostly thermal. It is not a nuclear lever.**

### 5.5 Phonons and mechanics (Q4): `figs/m6_phonons.png`, `figs/m6_breathing.png`, `figs/m6_modes.txt`

**Optical phonons vs Letts lines**
- The fitted model gives a PdD optical band of 8.0–11.5 THz (33–48 meV): Γ 8.0; X 9.9/11.5; L 8.9/11.5 THz. Shifting the anchors by ±0.7 THz moves the band to 7.3–12.2 THz, and a +1.5 THz top gives at most 13.0 THz.
- The PdH band (with the anharmonic factor) is 12.0–17.2 THz. The PdD two-phonon band is 16–23 THz.
- The Letts lines fall as follows:
  - **8.3 THz is inside the PdD band, near Γ.**
  - **15.3 THz is above the PdD band** under every anchor variation (it would need about +35 % hardening). It is inside the PdH band.
  - **20.4 THz is above PdH and inside the PdD two-phonon band.**
- Frequency matches are therefore non-unique and do not discriminate between assignments.
- Loading dependence [BK]: the optical peak softens from α to β (69 → 56 meV in PdH) and is roughly flat to x ≈ 1. Loading cannot move PdD to 15.3 THz.

**Selection rules for a two-laser beat**
- The rocksalt Γ optical mode is T₁u: IR-active and **first-order Raman-inactive**. A two-laser difference-frequency drive cannot couple to it in the bulk, only at surfaces and defects.
- |q_L| = 1.35×10¹⁰ m⁻¹ is ≥500× the maximum beat wavevector (2.5×10⁷ m⁻¹). L-point driving is kinematically excluded except through atomic-scale disorder.

**Cost of a coherent optical-phonon amplitude**
- A coherent 0.01 Å amplitude in a 1 mm² × 20 nm skin costs 5–26 kW (Q = 50–10).
- 30 mW gives ≤1.6×10⁶ coherent phonons against 1.2×10¹⁵ thermal ones.

**Breathing modes**
- The Pd l = 0 breathing mode follows f = 4.1 THz·nm/d: 820 GHz at 5 nm and 82 GHz at 50 nm.
- Reaching 8.3 THz needs d = 0.49 nm, which is not a nanoparticle. **Breathing modes cannot be tuned to the optical band.**

**Membrane and wire modes**
- C3 membrane f₀₁:
  - a = 1 mm, h = 50 µm: 80 kHz in vacuum, 52 kHz with D₂O on one side.
  - a = 10 mm: 0.8 / 0.21 kHz.
  - Thickness mode c_L/2h: 44 MHz at 50 µm.
- C1 wire, d = 0.5–1 mm, L = 30 mm:
  - longitudinal 52 kHz;
  - flexural 1.6–3.1 kHz;
  - radial 6.3–3.2 MHz.
- Q:
  - internal 10²–10³ **[BK]**;
  - thickness mode radiating into D₂O ≈ 50;
  - flexural in liquid 10–50.

**Pressure load (a mechanical constraint found here)**
- 1 bar across an unsupported 50 µm Pd membrane of radius 5 mm gives σ = 94 MPa, and 149 MPa at 10 mm. Both exceed the ~50 MPa yield.
- The maximum unsupported span for σ ≤ 30 MPa is 0.9 / 1.8 / 3.6 mm for h = 25 / 50 / 100 µm. **C3 needs a support grid.**

**Ultrasound and Hagelstein-type drive**
- A strain of 10⁻⁵ needs 4.8 W cm⁻² (stress 1.1 MPa).
- The MHz Hagelstein/Metzler proposal (1–3 MHz, 2.6×10¹⁵ quanta per 23.85 MeV) is classical at 300 K (occupation 3×10⁶).
- 1 W cm⁻² into a 44 MHz membrane resonance at Q = 20 gives a strain of 1.6×10⁻⁵, i.e. 1.6×10⁻⁹ eV per D.
- Si detectors and preamps are microphonic. **Not recommended in C3 iteration 1.** A Co-57 phonon–nuclear test belongs in a separate rig (R4 EV-4).

### 5.6 Multilayers and escape (Q5): `figs/m6_skins.txt`

**Ranges and escape (the M5 branch was not available, so they were computed here)**

| Particle | R(PdD₀.₈) | Escape-equivalent thickness (R_thr/4) |
|---|---|---|
| 3.02 MeV p | 33.4 µm | 6.8 µm (E > 1 MeV) |
| 1.01 MeV t | 4.7 µm | 0.93 µm (E > 0.2 MeV) |
| 0.82 MeV ³He | 1.35 µm | 0.22 µm |

**The two stacks**
- Pd/CaO: 40 + 5×(2+18) = 140 nm, with 10 buried interfaces.
- Ni/Cu: 6×(2+14) = 96 nm, with 11 interfaces.
- That is 1.5–1.7×10¹⁶ interface sites cm⁻² at one monolayer per interface.
- Both stacks sit well inside every escape depth.
- Interface sites are only **4.6×10⁻⁴** of the D sites within the proton escape volume (3.3×10¹⁹ cm⁻²). A stack wins only if an interface site is more than 10⁴–10⁵ times as active as a bulk site.

**Placement and operating regime**
- In Iwamura's experiments the stack is on the **upstream (high-pressure) side**. Placing it on the C3 exit face reverses that, and the stack then sits in the drained, low-x region (§5.7).
- Ni/Cu claims come from 250–900 °C gas cycling, which is not reachable in an electrolyte-backed C3.

### 5.7 Segmented face: skins, energy loss, exit loading (Q6): `figs/m6_skins.png`, `figs/m6_sector_flux.png`

| Skin | Thickness / route | ΔE p₃.₀₂ at 0° / 60° (keV) | ΔE t₁.₀₁ at 0° / 60° (keV) | s (exit k_r proxy) | x_out at J = 10¹⁶ | Notes |
|---|---|---|---|---|---|---|
| (a) bare Pd | anneal 850–900 °C in vacuum, etch; ~1.5 nm C/O | 0.03 / 0.06 | 0.14 / 0.27 | 0.1 (dirty: 10⁻³) | 2×10⁻⁵ (α) | baseline; drains |
| (b) PdO | 10 ± 3 nm thermal oxide (e.g. 450–550 °C in O₂ **[BK]**), masked | 0.48 / 0.96 | 1.6 / 3.3 | 10⁻⁶ (intact) | 0.006–0.61 | **Transient.** At J = 10¹⁶ its lifetime is 0.2–23 h (ε_red = 10⁻²–10⁻⁴). Monitor D₂O (m/z 20) |
| (c) Au/Pd/PdO | PdO 10 nm on exit; Au 20 nm on entry face of the sector | 0.48 / 0.96 | 1.6 / 3.3 | 10⁻⁶ | 0.01–0.08 (entry starved) | Back-face Au cuts the sector's supply (×0.1 assumed), so the sector drains laterally. Not a clean variant |
| (d) Pd/CaO | Pd 40 / [CaO 2 / Pd 18]×5, RF/DC sputter ±0.2 nm, shadow mask | 8.3 / 16.6 | 27 / 54 | 0.1 (top Pd) | 2×10⁻⁵ | Reversed vs Iwamura; drains |
| (e) Ni/Cu | 6×[Cu 2 / Ni 14] sputter, Ni on top **[assumed]** | 5.6 / 11.2 | 18 / 35 | 10⁻² | 6×10⁻⁵ | Outside its claimed T regime |
| (f) Au | 20 nm (or 50 nm) e-beam/sputter, no adhesion layer | 1.4 / 2.9 (50 nm: 3.6 / 7.2) | 3.9 / 7.8 (9.8 / 19.5) | 10⁻⁷ | 0.61 (β edge) | **Not a clean null.** The Pd under Au is the *most* loaded region. Au tests surface vs bulk origin, not "no D" |

- **k_r (exit recombination) ordering, for M3:** clean Pd ≈ Pd/CaO top > Ni/Cu > contaminated Pd ≫ intact PdO ≳ Au(20 nm) > Au(≥50 nm, pinhole-free).
- **Required s for a loaded exit face:**
  - x_out ≥ 0.6 needs s ≤ 1.5×10⁻⁷ at J = 10¹⁶ (1.5×10⁻⁶ at 10¹⁷);
  - x_out ≥ 0.8 needs s ≤ 2.2×10⁻¹¹ at 10¹⁶;
  - x_out ≥ 0.9 needs s ≤ 7×10⁻¹³.
- **The escape volume sits at x_out:** the loading drop across 7 µm is 2×10⁻⁴.
- **Draining without a cap:** keeping x_in = 0.9 behind a bare exit face at h = 50 µm needs 6×10¹⁸ D cm⁻² s⁻¹, i.e. about 1 A cm⁻² absorbed.
- **D₂ backfill alternative:**
  - 0.1 bar: x_out = 0.61, with losses over 2 cm of 4.6 keV (p), 25 keV (t) and 105 keV (³He).
  - 1 bar: x_out = 0.65, with losses of 46 keV (p) and 277 keV (t); ³He is stopped.
  - Si detectors at ≤100 V bias sit below the D₂ Paschen minimum (~270 V).

**Sector flux partition** (six 3 mm strips, 0.5 mm gaps, h = 50 µm, j_abs = 0.05 A cm⁻²)
- High-s sectors carry 3.0–3.4×10¹⁷ D cm⁻² s⁻¹ at x_out ≈ 0 (entry face x ≈ 0.04).
- PdO and Au sectors carry 0.7–1.1×10¹⁷ at x ≈ 0.55–0.65.
- Bare gaps take 27 % of the total flux; Au-capped gaps take 3 %.
- The edge influence width is ≤0.33 mm, so sectors ≥2 mm are independent of each other. They are not equivalent in (x, J).

**Cross-contamination**
- Lateral surface diffusion of Au, Ni and Cu on Pd at ≤350 K is negligible over µm in months **[BK]**.
- The real risks are process-induced:
  - Thermal oxidation for (b) and (c) at 450–550 °C oxidises every exposed sector unless it is done **first** under a SiO₂ or Au mask, or done locally.
  - Shadow-mask penumbra for sputtered layers is 50–100 µm at 0.1–0.2 mm mask gap. Keep gaps ≥0.5 mm.
  - Any anneal after Au deposition (>200 °C) interdiffuses Au–Pd over a few nm.
  - Cu and Ca from (d) and (e) can migrate during D₂O exposure of the back face. The back face must stay common and bare.
- **Sequence:** anneal and etch → PdO growth (masked) → sputter (d) and (e) through shadow masks → Au (f) last at room temperature.

### 5.8 Collimated detectors per sector (Q6): `figs/m6_sector_collimation.png`, `figs/m6_sectors.txt`

**Layout.** Hexagonal Ø6 mm patches. Each has a coaxial Ø6 mm (28 mm²) detector and a knife-edge aperture of equal radius.

| D (mm) | Aperture height (mm) | Efficiency Ω/4π (own patch) | Minimum gap for ≤1 % cross-talk (6 neighbours) |
|---|---|---|---|
| 10 | 0.5 | 2.0 % | ≤0.25 mm |
| 10 | 1.0 | 2.0 % | 0.5 mm |
| 10 | 2.0 | 2.0 % | 1.0 mm |
| 20 | 0.5–1.0 | 0.55 % | ≤0.25 mm |
| 30 | any ≤2 | 0.25 % | ≤0.25 mm |

- Reference: one uncollimated Ø25 mm detector at 5 mm sees 31 %.
- **Segmentation costs roughly ×15 in per-area efficiency** (about 2.7 e-folds of the M0 budget) in exchange for internal controls.
- Minimum sector: Ø6 mm (0.28 cm²) at D = 10 mm. Ø4 mm drops efficiency to 0.96 %.
- Six patches with 1 mm webs fit a Ø22 mm active face. The support-grid ribs (§5.5) can run along the webs.

## 6. Design recommendations for the lead

### Active-surface specification sheet

| Item | Specification | Expected value / justification |
|---|---|---|
| **Exit-face loading valve (C3)** | Either (i) a sector cap with s ≤ 10⁻⁷ (Au ≥50 nm, pinhole-free), or (ii) D₂ backfill of 100 ± 50 mbar on the detector side, with x_out verified by 4-wire R or ΔR of a witness strip | Without it the escape volume is α-phase (x ~10⁻⁴) and M0's "high loading within the escape depth" is not met. Backfill gives x ≈ 0.61 at ≤25 keV triton loss |
| **Support grid (C3)** | Perforated backing on the electrolyte side, open span ≤1.5 mm for h = 50 µm (≤0.8 mm for 25 µm); ribs under inter-sector webs | σ ≤ 30 MPa at 1 bar (§5.5) |
| **Sector layout** | 6 × Ø6.0 ± 0.1 mm patches, webs ≥1.0 mm (≥0.5 mm minimum), Au-capped webs | ≤1 % cross-talk, 2.0 % efficiency each at D = 10 mm |
| **Detectors** | 28 mm² Si, D = 10 ± 0.5 mm, knife-edge aperture Ø6.0 mm at 0.5–1.0 mm above the membrane, light-tight (≥0.2 µm Al on the aperture/window if a laser is used) | §5.8 |
| **Skins (iteration 1)** | (a) bare Pd; (f) Au 50 ± 5 nm; (b) PdO 10 ± 3 nm; (d) Pd 40/[CaO 2/Pd 18]×5 ± 0.2 nm; duplicate (a) and (f) positions | Energy loss ≤17 keV (p, 60°); (d) needs the t line-shape correction (27–54 keV) |
| **Roughness PSD (all cathodes)** | *Covariate, not a target*: AFM 5, 20 and 80 µm scans; report the radial PSD and the band rms in 10⁵–3×10⁶ m⁻¹; accept Ra ≤ 50 nm | ENEA band not plasmonic in the dark (§5.2) |
| **Optional plasmonic patch (only if a laser port is built)** | C3 vacuum face: 1-D sinusoidal or ring-PSD grating Λ = 761 ± 30 nm, depth 100 ± 30 nm (785 nm, p-polarised, ±4°); C1 in D₂O: Λ = 562 ± 40 nm, h = 85 ± 25 nm; Au/PdD only if ±1 % / ±0.5° can be held | A ≈ 1.0 vs 0.26–0.32 flat (×3–4 absorbed power) |
| **Optional stimulation** | 785 nm, 0.3 W CW, chopped 0.1–10 Hz, 1 mm spot on one sector, used as a thermal flux modulator (ΔT 0.7–5 K, Δk_r 5–50 %) | Lock-in test of flux-correlated emission. No nuclear-barrier effect expected (≤10⁻⁵) |
| **Stimulation to exclude** | Dual-laser THz beats, MHz or megasonic drive on C3, nanoparticle breathing-mode tuning | §5.5 |

### Ranked list for iteration 1

**Include, in order of value:**
1. Exit-loading valve (D₂ backfill or Au ≥50 nm cap). Without it, the other surface features sit on unloaded Pd.
2. Support grid with ribs under the sector webs. This is mechanically mandatory at 1 bar.
3. Segmented face with (a) and (f) duplicated. This is the cheapest internal control and splits surface-origin from bulk-origin emission.
4. PdO sector (b), with its reduction lifetime measured (D₂O at m/z 20). Kasagi's accelerator-regime U_e ≈ 2× Pd is the only peer-reviewed surface-chemistry lever.
5. Pd/CaO sector (d). Cheap if sputtering is available; it tests the Iwamura interface idea, but in a reversed and drained geometry.
6. AFM PSD and EBSD as recorded covariates on every membrane.

**Optional, low cost:**
7. One laser-port sector with a Λ = 761 nm grating and a chopped 785 nm laser, used as a flux/recombination modulator.

**Leave out (unsupported or counter-productive):**
8. Ni/Cu stack at room temperature.
9. Lipson-type heterostructure (c) with back-face Au.
10. Engineered ENEA-band PSD in dark operation.
11. Dual-laser THz beats.
12. Nanoparticle breathing-mode tuning.
13. Nanogap and hot-spot arrays as barrier modifiers.
14. MHz or megasonic drive on C3.

### Numbered recommendations

1. **Give C3 an exit-loading valve.** D₂ backfill of 0.1 ± 0.05 bar is preferred, because it is skin-independent and equalises sectors at x_out ≈ 0.61. Alternatively use an Au ≥50 nm cap on the active sectors. Log x_out per sector. Reason: every bare skin gives x_out < 10⁻³ (§5.7).
2. **Keep membrane thickness at h = 50 ± 10 µm on a ≤1.5 mm-span support grid.** Thinner membranes need proportionally smaller spans. Thicker ones drain more slowly but waste escape volume: only the top ~7 µm is seen.
3. **Six Ø6 mm sectors, ≥1 mm Au-capped webs, D = 10 mm, knife-edge at ≤1 mm.** Accept 2 % efficiency per sector, or reduce to four larger sectors if M5's sensitivity budget cannot absorb the ×15 cost.
4. **Skins for iteration 1:**
   - (a) and (f) twice each, as a position control;
   - (b) with an RGA m/z 20 log;
   - (d) as the only multilayer.
   - Drop (e) Ni/Cu: its claims are at 250–900 °C.
   - Drop (c) as specified: back-face Au starves the sector. Use (b) instead unless the back face can be patterned.
5. **Treat Au (f) as a surface-vs-bulk discriminator, not a null.** A true null is the H₂O twin (charter) or an unloaded Pd sector behind a D-impermeable back mask.
6. **Fabrication sequence:** anneal and etch → masked PdO → shadow-mask sputter (d) → Au (f) at room temperature. No anneals after Au. Web and penumbra ≥0.5 mm.
7. **Roughness:** record the PSD before and after; do not engineer a PSD band unless a laser port exists.
8. **If a laser port exists:** build the Λ = 761 nm grating on a single sector, specified as a sinusoid or ring PSD with a 36 nm ring rms. Use 785 nm at 0.3 W, chopped. The expected effect is modulation of flux and recombination only.
9. **Fabrication routes for the grating (cost [BK], per 25 mm membrane):**
   - fs-laser LIPSS: period ≈ 0.7–0.95 λ_fs, self-tuned to the SPP; about $0.5–2k at a job shop.
   - Interference lithography (Lloyd mirror, 405 nm, θ = 21° for 562 nm) plus Pd electroplating through resist: about $1–3k.
   - Nanoimprint: $5–10k for the master.
   - Anodic alumina: pitch ≤500 nm; suits nanopillars, not SPP gratings.
   - FIB: ≤0.01 mm² per hour; test coupons only.
   - Electrochemical etching: random PSD, band not controllable.

## 7. Sensitivities and uncertainties

- **Sticking and recombination values s** are order-of-magnitude [BK]. Recommendation 1 fails only if bare Pd in UHV had s ≲ 10⁻⁷, which is 5–6 orders of magnitude below literature values. It is robust.
- **The isotherm at x > 0.7** (×10 per anchor) moves the s required for x = 0.8 by ×10 but does not rescue bare Pd.
- **Linear-diffusion sector model:** α/β interfaces, stress and cracking (M3) will sharpen fronts, but they cannot equalise sector fluxes.
- **ε(PdD) scaling F:** varying F from 0.6 to 0.9 shifts Λ_opt by <2 % and h_opt by about 10 %. Conclusions are unchanged.
- **Phonon anchors ([BK])** are ±0.7 THz. The 15.3 THz exclusion survives up to +1.5 THz. If new INS data show PdD optical modes reaching 15 THz at x ≈ 0.95, the Letts assignment question reopens (the frequency argument only, not the power argument).
- **Collimation** assumes straight tracks. Edge scattering at knife edges adds ≲0.1 % tails. M5 should check with its transport code.
- **Energy losses:** stopping powers are ±5 %, and skin thickness tolerances dominate.
- **The ENEA band itself was not retrieved.** The analysis covers both unit conventions.
- **Assumed j_abs = 0.05 A cm⁻².** Raising it ×10 still leaves bare sectors drained (1 A cm⁻² needed).

## 8. Open questions and hand-offs

- **M3:** adopt the k_r ordering and s ranges (§5.7). Model the α/β front under each skin, and the support-grid-induced entry pattern. Stress at 1 bar vs span. PdO reduction kinetics.
- **M5:**
  - verify the escape numbers above (33.4 µm p, 4.7 µm t);
  - fold in the ×15 segmentation cost;
  - include line shapes for 27–54 keV triton shifts ((d) sector) and backfill losses;
  - check Si-detector operation at 0.1 bar D₂.
- **M2:** achievable *absorbed* current fraction at the C3 entry face, which sets J₀.
- **M1:** whether interface or oxide sites need to be >10⁴× more active than bulk sites for the multilayer sector to matter.
- **Lead and literature:** retrieve ENEA's actual PSD figure (Violante et al., ICCF-12 to 15) and a tabulated ε(PdHₓ) (von Rottkay 1999) to replace the [BK] scaling. Obtain INS data for PdD at x ≥ 0.9.
