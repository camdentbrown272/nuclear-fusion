# R4: LENR theories and what they predict about device geometry

**Scope.** Proposed LENR ("cold fusion") mechanisms, every geometric, dimensional, structural or frequency prediction they make that could shape a device, the strongest published critique of each, and our credence. Hot (thermonuclear or plasma) fusion is out of scope.

**Source caveat.** The egress proxy blocked direct fetches of lenr-canr.org, jcmns.org, arXiv and publisher PDFs. Numbers come from search-result abstracts, the established literature and our own calculations (Section 4). Items marked **†** come from memory or secondary summaries; verify them before using them as hard specifications.

**Credence** = our probability that the mechanism, as its authors state it, operates strongly enough to matter for a detectable nuclear signal in a benchtop device. It is not the probability that the physics ingredient exists: enhanced screening is almost certainly real at keV beam energies, but whether it matters at thermal energies is a separate question.

---

## 0. Executive summary

1. **Mainstream bounds** (Leggett–Baym equilibrium bound; Koonin–Nauenberg ~10⁻⁶⁴ s⁻¹ per D₂; Huizenga's "three miracles"; missing 23.8 MeV γ-rays and neutrons) constrain every mechanism. Each theory below is a proposal for evading one of them. Status: **mainstream**.
2. **Electron screening.** Measured U_e ≈ 300–800 eV in Pd, Zr, Ta, PdO versus ~25–100 eV theory. Czerski adds a 0⁺ threshold resonance and an e⁺e⁻ channel (PRC 2022/2024, PRX 2025); Berlinguette et al. (Nature 2025) found +15 % dd rate from electrochemical loading. Status: **mainstream-compatible in the keV regime**, *speculative* if extrapolated to thermal energies. The only family with peer-reviewed geometry/material dependencies (vacancies, oxide, loading, temperature).
3. **Hagelstein phonon–nuclear coupling ("fractionation").** Predicts D₂ in monovacancies at x ≳ 0.9–0.95 and a single coherent optical-phonon mode (PdD Γ ≈ 8.3 THz, L ≈ 15.3 THz); proposes MHz–THz vibration experiments. Explicit mathematics, unproven coupling strength, no independent replication. Status: **speculative**.
4. **Takahashi TSC.** 4d + 4e⁻ at a tetrahedral site collapse to fm scale in ~1.4 fs (⁸Be* → 2α); points to 2–10 nm Pd/Ni in ZrO₂ at 200–300 °C. Localizing the electrons costs ~37 MeV; the predicted 23.8 MeV alphas are absent. Status: **strongly disputed**.
5. **Widom–Larsen** (e⁻ + p → n + ν via "heavy" electrons). Needs β ≥ 2.53 mₑ, i.e. ~4×10¹¹ V/m at phonon frequencies (~10¹³ V/m optical); the screened lattice field is ~10⁶ V/m. Refuted by Ciuchi–Maiani–Polosa et al. (2012), Tennfors (2013), Hagelstein–Chaudhary (2008). Status: **strongly disputed**.
6. **Storms hydroton / NAE.** ~1 nm (<10 nm) nanocracks host linear H chains. Qualitative, no rate equation. Status: **speculative**.
7. **Kim BEC fusion.** nm grains, low temperature. Localized deuterons give T_c ≪ 1 K, and a BEC would not remove the short-range Coulomb hole. Status: **strongly disputed**.
8. **Coherence/band models** (Preparata, Chubb, Schwinger): qualitative needs (x → 1, ordered finite crystals, ~0.1 µm coherence domains) that contradict observed D localization. Status: **speculative to strongly disputed**. Swartz's "optimal operating point" is phenomenology.
9. **High fields (nanotips, plasmonic hot spots)** cannot touch the barrier: 10¹⁰ V/m is 10⁻¹⁰ of the Coulomb field at 5 fm, and a uniform field exerts no relative force on two deuterons (identical q/m). Tips only act as micro-accelerators (hot fusion; useful for detector calibration). Status: **ruled out by mainstream physics**.
10. **Robust targets** (several theories agree): vacancy-rich, high-loading (local D/Pd ≥ 0.9) near-surface Pd ≤1 µm thick; 2–10 nm features with ~1 nm gaps; driven D flux; diagnostics for dd-branch particles, 511 keV and ⁴He. **Highest expected value:** *low-energy (0.3–5 keV) D⁺ implantation into a loaded, defect-engineered target*, which gives a guaranteed, calibrated nuclear signal whose anomalous part can be measured against geometry (Section 5).

---

## 1. Framing: the bar every mechanism must clear

Our calculation (Section 4.1) uses the shifted-energy WKB model with the dd reaction constant A = S(0)/(π μ c α) = 1.5×10⁻¹⁶ cm³ s⁻¹. It gives these requirements for 1 cm³ of PdD in which *every* deuteron pair is active:

- **1 neutron/s** (a credible signal with a good detector) needs a *sustained static* screening energy of U_e ≈ 200 eV (range 110–440 eV once the prefactor uncertainty of ±10 orders of magnitude is included).
- **1 W** through ⁴He at 23.85 MeV needs U_e ≈ 500 eV (220–2100 eV).
- If only 10⁻⁶ of sites are active (an "NAE" scenario), 1 W needs U_e ≈ 1 keV.

The measured beam-regime screening energies (300–800 eV) sit *uncomfortably close* to these numbers, the most interesting fact in this review. The catch: a static U_e of 200 eV means the D–D potential stays flat down to a 7 pm turning point, which would collapse the observed equilibrium D–D spacing (0.74 Å in D₂, 2.85 Å between PdD O-sites). That is the physical content of Leggett–Baym.

A viable mechanism therefore needs (a) non-equilibrium, localized ~0.2–1 keV per pair without keV beams; (b) a new coupling that bypasses tunnelling (Hagelstein); or (c) a different nuclear process (Widom–Larsen weak interaction; Czerski resonance/e⁺e⁻ channel).

---

## 2. Theories

### 2.1 Hagelstein (MIT): phonon–nuclear coupling, lossy spin-boson model, fractionation

**Core claim.** In Fleischmann–Pons (FP) excess heat, D₂ → ⁴He (Q = 23.85 MeV) occurs *coherently*. The large quantum is fractionated into ~10⁸–10⁹ optical phonons instead of being emitted as a γ-ray or as energetic particles. The model system is a lossy spin-boson Hamiltonian:

H = ΔE·S_z/ħ + ħω₀ a†a + V(a + a†)·2S_x/ħ − iħΓ(E)/2

The loss term Γ(E) (decay channels open at high energy) destroys the destructive interference that normally blocks exchange between one large quantum and many small ones. Exchange is fast when g ≈ V√n/ΔE ~ 1 (n = phonon occupation), helped by a Dicke factor √N from many D₂ sites sharing one mode. Since ~2012 the coupling is attributed to a relativistic term from boosting a many-particle Dirac Hamiltonian with the lattice motion, H_int = **a·cP** (P = phonon-induced momentum); electron-mediated couplings were found too weak. The framework also predicts vibration-induced **excitation transfer** and **up-conversion**, which explains Karabut's collimated ~1.5 keV X-rays as coherent excitation of the 1565 eV level of ²⁰¹Hg (the lowest excited state of any stable nucleus), emitted as a phased array normal to the cathode.

**Numbers (our calculations).**
- 23.85 MeV = **6.95×10⁸** quanta at 8.3 THz (34.3 meV), 2.6×10¹⁵ at 2.21 MHz. 1565 eV = 4.6×10⁴ quanta at 8.3 THz, 1.7×10¹¹ at 2.21 MHz.
- Thermal occupation at 8.3 THz, 300 K: n ≈ 0.36; the model needs n ≫ 1 in one coherent mode.
- Sustaining a coherent 0.001 Å, 8.3 THz displacement over 1 cm³ PdD costs ~1.6×10⁹ W at Q = 100. Such a mode must be self-excited or confined to tiny volumes; tens-of-mW lasers cannot drive it.

**Explicit geometric, dimensional and frequency predictions**

| Parameter | Predicted optimum or range | Units | How to realize in hardware |
|---|---|---|---|
| Bulk D/Pd loading x | ≥ 0.90; ~0.95 for spontaneous superabundant vacancies | atom ratio | Thin films (≤1 µm), codeposition; monitor R/R₀ (4-wire) |
| Pd vacancy fraction | Up to ~25 % (δ-phase, Pd₃VacD₄-like) at x ≳ 0.95; monovacancies host molecular D₂ | site fraction | Pd/D **codeposition** (Szpak-type), high-pressure H anneal (Fukai superabundant vacancies), irradiation |
| Phonon mode (FP heat) | Compressional **optical** mode, zero group velocity: Γ ≈ 8.3 THz, L ≈ 15.3 THz (PdD); ≈ 20 THz attributed to PdH or D in Au | THz | Dual-laser beat on Au-coated PdD; THz difference-frequency or FEL sources; ion/current-pulse phonon injection |
| Phonon mode (Karabut X-ray, host-nucleus transitions) | Acoustic modes, MHz–GHz ("frequency as high as possible") | Hz | Piezo at 1–3 MHz (Metzler/Tanzella used ≈2.21 MHz); pulsed-discharge excitation; heavy holder as acoustic resonator |
| Mode Q / coherence | A single, highly excited, coherent mode shared by all active sites; zero group velocity to keep energy local. No numeric Q is given. **Our estimate:** PdD optical-phonon Q ≈ 10–50 (ps lifetimes), which works against the model | — | Choose Γ/L-point modes; minimize inhomogeneous broadening (uniform loading) |
| Dicke number N | As large as possible within one coherent mode volume | — | Maximize the density of vacancy–D₂ sites inside one phonon coherence volume |
| Electrochemical current density | Threshold behaviour (McKubre: P_xs ∝ (x − x₀)²(i − i₀)·\|i_D\|, x₀ ≈ 0.875, i₀ ~ 0.1–0.25 A cm⁻² †) | A cm⁻² | Galvanostatic drive above threshold, together with deuterium flux across the surface |
| X-ray beam geometry (Karabut) | Emission collimated **normal to the cathode surface**; divergence ~λ/D (0.83 nm/1 mm ≈ 10⁻⁶ rad) | rad | Flat, polished cathodes; X-ray detectors on the surface normal; trace Hg (²⁰¹Hg) impurity |

**Testable predictions.** Heat–⁴He correlation near 24 MeV/atom with no energetic particles (alphas ≲20 keV); heat responds to THz stimulation at Γ/L frequencies; MHz vibration alters Co-57→Fe-57 emission (excitation transfer, delocalization, 14.4 keV anisotropy); collimated keV X-rays from vibrated surfaces with low-lying nuclear levels.

**Strongest critique.**
- The a·cP coupling has never been shown to be large enough; Hagelstein's matrix-element calculations (JCMNS 2013–2025) are unfinished.
- Dicke enhancement assumes coherence across 10⁸–10⁹ phonons and many sites, which ~ps phonon lifetimes make implausible.
- The experimental anchors are unreplicated or ambiguous. The Letts heat was **not reproduced** (0 events in 231 laser-triggered trials, 9 cathodes). The Metzler–Hagelstein–Lu Co-57 run (JCMNS 27, 2018) did not see the vibration-induced effect it sought and reported "non-exponential" count histories instead, which are prone to artefacts. The Hagelstein–Tanzella vibrating-copper run was also negative. Karabut is a single group.
- No hostile peer review in mainstream journals.

**Credence: 4 %.** Internally consistent, with falsifiable frequency and vacancy predictions (rare here), but nothing independent supports the huge coupling. Still the most *design-relevant* LENR-specific theory.

---

### 2.2 Takahashi: tetrahedral symmetric condensate (TSC) and multibody fusion

**Core claim.** Four deuterons on the O-sites around an empty tetrahedral (T) site, plus four electrons on the interpenetrating tetrahedron, form a neutral pseudo-particle. Semi-classical Langevin calculations (EQPET) have it condense from ~74–100 pm to ~10–20 fm in **≈1.4 fs**, where 4D fusion gives ⁸Be* (47.6 MeV). Originally ⁸Be* → 2α at 23.8 MeV each; later a "burst of low-energy photons" (BOLEP) cascade to the ⁸Be ground state gives two ~46 keV alphas. A 4H/TSC version covers Ni–H. Reported rate: 5.5×10⁻⁸ f s⁻¹ per cluster (time-averaged); projected "several to a few 100 W/cm³".

**Predictions**

| Parameter | Predicted optimum or range | Units | Realization |
|---|---|---|---|
| Local site | T-site with 4 filled neighbouring O-sites (local D/M ≈ 1); D–D spacing before collapse 2.85 Å (PdD), O→T distance 1.75 Å | Å | High local loading at the surface and subsurface |
| Particle size | ~2–10 nm Pd, Ni or Pd–Ni "nano-cores" in an oxide matrix; "sub-nano-holes" (SNH, ~0.3–1 nm surface defects) and a "global mesoscopic potential well" on the particle | nm | Melt-spun, oxidized Pd₁Ni₁₀/ZrO₂ (PNZ) and Cu₁Ni₇/ZrO₂ (CNZ) composites (Kitamura/Takahashi; NEDO-MHE project) |
| Temperature | ~200–350 °C for gas loading | °C | Heated D₂/H₂ gas cells |
| Symmetry | Tetrahedral (4D) or octahedral (8D) | — | fcc hosts (Pd, Ni) |
| Drive | Phonon excitation or D desorption bursts (non-equilibrium) | — | Temperature/pressure cycling, net desorption mode |
| Reported heat | 80–400 W/kg sustained for weeks at ~300 °C (PNZ with D₂) | W/kg | — |

**Testable predictions.** ⁴He is the main ash, with neutrons and tritium at very low levels. The original version predicts 23.8 MeV alphas. H and D gas both work, with different ash.

**Strongest critique.**
- Localizing four electrons within 20 fm costs ≈9.4 MeV *each* (≈37 MeV; Section 4 (other checks)), which chemistry cannot supply. The "quadruplet electron" state is unproven, and the quantum many-body problem is treated semi-classically.
- A symmetric four-body collapse is what Leggett–Baym-type bounds exclude in equilibrium.
- 23.8 MeV alphas would cause copious secondary radiation and are not seen; the ≲20 keV alpha bound rules out the original version. BOLEP is ad hoc (no known process damps 47.6 MeV by a "black-body" cascade).
- The gas-loading heat has no correlated nuclear ash.

**Credence: 1 %.** The 200–300 °C nanocomposite recipe is still cheap to test.

---

### 2.3 Widom–Larsen: ultra-low-momentum neutrons (ULMN)

**Core claim.** Collective oscillations of surface H/D "patches" on hydride surfaces create strong local fields that raise the mass of surface-plasmon-polariton (SPP) electrons, enabling ẽ⁻ + p → n + ν_e. Born in a patch-scale collective mode, the neutrons have ultra-low momentum and huge capture cross-sections, so they are captured locally (transmutation). Heavy electrons supposedly absorb the resulting γ-rays.

**Equations and numbers.**
- Threshold m̃c² ≥ (m_n − m_p)c² ⇒ **β = m̃/mₑ ≥ 2.531**. With β = √(1 + a₀²), a₀ = eE/(mₑcω), this needs a₀ ≥ 2.33.
- Required field: **3.6×10¹¹ V/m** at ħω = 60 meV (PdH optical phonon); 9.3×10¹² V/m at 1.55 eV, i.e. ~10¹⁹ W/cm² laser intensity. That is the relativistic laser-plasma regime, which yields MeV electrons, not ULMNs.
- Realistic field: displacing all H in PdH coherently by 0.1 Å gives 1.2×10¹⁰ V/m unscreened; conduction-electron screening by (ω/ω_p)² ≈ (0.06/7)² leaves ~10⁶ V/m, so β − 1 ≈ 10⁻¹¹.

**Predictions**

| Parameter | Predicted optimum | Units | Realization |
|---|---|---|---|
| Collective-proton "patch" size | nm to ~µm (not fixed quantitatively); the ULMN wavelength is of order the patch size | m | Nanostructured surfaces, particles, cracks |
| Local field | ≳10¹¹ V/m claimed; we calculate ≥3.6×10¹¹ V/m | V/m | Sharp features, plasmonic hot spots, very high H flux |
| Surface roughness / nanoparticles | Features ~10–100 nm to concentrate SPP fields | nm | Roughened Pd/Ni, Au/Ag nanoparticle decoration |
| Non-equilibrium drive | High H/D flux through the surface, current or laser excitation | — | Electrolysis, glow discharge, laser |

**Testable predictions.** Neutron-rich transmutations and isotope shifts; no free neutrons or γ ("heavy-electron γ shield"); easily detected β-emitters such as ¹⁰⁹Pd (13.7 h).

**Strongest critique.**
- Ciuchi, Maiani, Polosa, Riquer, Ruocco & Vignati (EPJC 72, 2193, 2012), with proper collective wavefunctions: "little room for such a remarkable effect". Maiani–Polosa–Riquer (EPJC 74, 2843, 2014): rates in plasmas are lower still.
- Tennfors (EPJ Plus 128, 15, 2013): neutron-scattering data on H in Pd bound the mass increase to **<1 %**. Hagelstein–Chaudhary (J. Phys. B 41, 125001, 2008): far smaller shift in Coulomb gauge.
- γ-absorption by heavy electrons contradicts QED. Each neutron needs 0.78 MeV focused onto one electron from ~10⁷ phonons of 0.06 eV: fractionation in reverse, with no mechanism.

**Credence: 1 %.**

---

### 2.4 Storms: the nuclear active environment (NAE) and the "hydroton"

**Core claim.** LENR happens not in the bulk lattice but in a rare separate structure, the NAE: nanocracks formed by stress relief during loading. In a gap of the right width, H nuclei and electrons form a linear chain (the **hydroton**) with metallic-hydrogen-like bonding and strong screening. Axial resonance brings nuclei periodically close, and mass-energy leaves as *many coherent low-energy photons* ("slow fusion"), not one γ-ray.

**Predictions**

| Parameter | Predicted optimum or range | Units | Realization |
|---|---|---|---|
| Gap width | ~1 nm, <10 nm: narrow enough that H₂ (kinetic diameter 0.29 nm) cannot form ("gaps larger than a few atomic diameters" do), wide enough for single hydrons. Our reading: ~0.3–1 nm | nm | Loading/deloading cycles (stress-relief cracking), codeposition, dealloyed nanoporous metals, electromigrated break-junction arrays, ALD-spaced nanogaps |
| Location | Near-surface layer (µm-scale), chemically independent of the host lattice | µm | Thin films; surface treatments |
| Chain | Linear chain of H, D or T nuclei with reduced spacing; length unspecified | — | — |
| Resonance frequency | Axial chain vibration; **not quantified**. Our guess: 10–100 THz range for H chains | THz | — |
| Temperature | Rate rises with temperature (Arrhenius-like) | K | Heated cells |
| Isotope → product | D chains → ⁴He; H chains → D (weak process); mixed H+D → T | — | — |
| D/Pd loading | Not fundamental. High loading matters mainly because it produces cracks | — | — |

**Testable predictions.** Heat/ash scale with density of correctly sized cracks, not bulk loading; soft X-rays from active sites; tritium from H+D mixtures.

**Strongest critique.**
- No Hamiltonian or rate equation, so it cannot be falsified quantitatively.
- A tighter chain lowers the barrier only through screening; at ppm site density that needs ≈0.3–1 keV sustained (Section 4.1), 10–30× what metallic electrons give.
- Emitting 23.8 MeV as small photons faces the fractionation objection, with less mathematics than Hagelstein.
- Cracks are ubiquitous; the effect is rare.

**Credence: 1.5 %.** Useful idea: vary *engineered ~1 nm gaps* systematically; cheap, and relevant to screening too.

---

### 2.5 Kim (Purdue): Bose–Einstein condensation nuclear fusion (BECNF)

**Core claim.** Deuterons trapped in micro/nano-scale grains form a BEC in a harmonic trap (equivalent linear two-body method, optical-theorem formulation; PRC 55, 801, 1997; Naturwissenschaften 96, 803, 2009). Fusion is collective, energy is shared by the condensate (no γ), and ⁴He is the main product.

**Predictions**

| Parameter | Predicted | Units | Realization |
|---|---|---|---|
| Grain or particle size | Nano-scale. Smaller grains favour condensate formation. Kim cites the Arata–Zhang ZrO₂–Pd nanoparticles (~5–10 nm †) | nm | Pd in ZrO₂ nanocomposites |
| Temperature | Lower is better for the BEC, but D mobility needs finite temperature. Kim proposed low-temperature tests † | K | Cryogenic gas loading |
| D density | High | m⁻³ | High loading |

**Strongest critique (our calculations, Section 4 (other checks)).**
- *Free* deuterons at PdD density: T_c ≈ 13 K; at 300 K nλ³ = 0.024 versus the 2.61 needed.
- D in Pd is localized (band mass ≥10²–10⁴ m_D), so T_c ≈ 0.1–0.001 K.
- A BEC of charged bosons keeps its short-range Coulomb hole; fusion depends on g(r → fm), which condensation does not change (Leggett–Baym again).

**Credence: <1 %.**

---

### 2.6 Coherence and band theories: Preparata, Chubb & Chubb, Schwinger (plus Swartz's OOP)

**Preparata** (Bressani, Del Giudice & Preparata, Nuovo Cimento A 101, 845, 1989; *QED Coherence in Matter*, 1995). Matter forms QED "coherence domains" (CDs) oscillating in phase with a trapped EM mode. Above a loading threshold (x ≳ 0.7–0.85 †) deuterons form a coherent plasma screened by the d-electron plasma. Geometry: CD size ≈ mode wavelength (~0.1 µm for a ~10 eV mode); a loading threshold; longitudinal currents/fields along Pd wires to drive loading †. **Critique:** the CD ground-state instability is not accepted in mainstream QED, and no quantitative prediction has been independently confirmed (see Chechin et al., nucl-th/0303057). **Credence: 1 %.**

**Chubb & Chubb ion band states** (Fusion Technol. 20, 93, 1991). Near stoichiometry, D⁺ (and ⁴He⁺⁺) occupy Bloch-like states across a finite, perfectly periodic crystal; overlapping delocalized charge "fuses" with energy shared by the lattice (Mössbauer-like). Predictions: ordered, defect-free **finite nanocrystals** ("particular nanoscale crystals turn on faster" †), x → 1, and a slow Zener-like onset of ion conduction under sustained field. **Critique:** H tunnelling bandwidths in Pd are ≲meV and destroyed by phonons at 300 K; delocalizing single-particle densities does not reduce the *two-particle* short-range correlation; contradicts all defect-centred theories. **Credence: <1 %.**

**Schwinger** (Z. Naturforsch. 45a, 756, 1990): p + d → ³He in HD impurities, 5.5 MeV passed to lattice phonons. Geometry: only lattice periodicity; predicts ³He scaling with H contamination. **Critique:** no MeV→phonon coupling mechanism (the problem Hagelstein later took up). **Credence: <1 %.**

**Swartz "optimal operating point" (OOP)** (JCMNS 6, 149, 2012): excess power (or ⁴He) versus input power peaks in a narrow window, for PHUSOR (spiral Pd, high-impedance D₂O) and NANOR (preloaded ZrO₂–PdNiD nanocomposite) devices. Phenomenology, not mechanism; peaked responses also arise from artefacts (recombination, calorimeter nonlinearity). **Implication:** sweep drive power finely. **Credence it reflects nuclear physics: 3 %.**

---

### 2.7 Electron screening and barrier-lowering models

**Core claim (mainstream observation).** Low-energy dd reactions in deuterated metals show enhanced yields, parametrized as σ(E) → σ(E + U_e), i.e. enhancement f = exp(πη·U_e/E).

| Host | U_e (measured) | Source |
|---|---|---|
| D₂ gas / atomic (theory ≈ 25–28 eV) | ~25 eV | adiabatic limit |
| Ta | 309 ± 12 eV | Raiola et al. (EPJ A 2002) |
| Pd | 310 ± 20 ± 50 eV (Lipson); up to ~800 eV (Raiola †) | Lipson et al.; Raiola et al. |
| Au/Pd/PdO heterostructure | 601 ± 23 eV; PdO 600 ± 20 ± 75 eV | Lipson et al. |
| Zr | ~115 eV to 497 ± 7 eV depending on surface/UHV condition; 340 eV (ZrD₂, 2024) | Huke/Czerski et al. PRC 78, 015803 (2008); J. Phys. G 35, 014012; arXiv 2409.02112 |
| Cu, Ag, Au (noble) | small (Au ≈ 70 eV) | Raiola; Lipson |

**Models.**
- **Debye plasma model** (Rolfs group): U_e ∝ (n_eff/T)^½, so U_e falls with temperature (reported); physically questionable, since degenerate metallic electrons should follow Thomas–Fermi screening.
- **Czerski:** screening plus a **0⁺ threshold resonance** in ⁴He near the d + d threshold (EPL 113, 22001, 2016; PRC 106, L011601, 2022); **electron emission** from it (PRC 109, L021601, 2024); **e⁺e⁻ branching** measured down to 5 keV (PRX 15, 041004, 2025); U_e rises with **lattice defects** (positron-annihilation data on deuterated Zr, Materials 16, 6255, 2023); a thermal "plateau" below 2.5 keV attributed to phonon-heated ion tracks (arXiv 2409.02112).
- **Lipson/Miley:** PdO/Pd heterostructures boost U_e; low-level emissions during exothermic D desorption.
- **Ichimaru** (RMP 65, 255, 1993): strongly coupled plasma screening, important in dense stars, only tens of eV at metallic densities. **Frisone:** enhancement at microcracks/dislocations with a temperature optimum † (qualitative). **Electron clusters / non-equilibrium electrons:** no quantitative model; hot electrons *reduce* static screening.
- **Berlinguette et al.** (Nature 2025, "Thunderbird"): in-situ electrochemical loading raised the dd rate by **15(2) %** under plasma-immersion D⁺ implantation, the first high-profile demonstration that *loading modulates a real nuclear rate*.

**Predictions**

| Parameter | Predicted optimum | Units | Realization |
|---|---|---|---|
| Host metal | Pd, Zr, Ta, PdO-capped Pd (high U_e); avoid Cu/Ag/Au | — | Thin films or foils |
| Vacancy / defect density | Higher is better (Czerski, positron-annihilation data) | site fraction | Ion pre-irradiation, cold work, codeposition |
| Surface oxide | Conflicting: PdO raises U_e (Lipson), while uncontrolled oxidation and contamination is the dominant *uncertainty* in UHV data (Czerski) | — | Deliberately compare clean-UHV and PdO-capped surfaces |
| Temperature | Lower (if the Debye model holds) | K | Cooled target stage |
| Projectile energy | 0.5–10 keV. Below ~5 keV, U_e of 300 eV gives f ≈ 1.5–100 (Section 4.2) | keV | Low-energy UHV ion beam or plasma-immersion implantation |
| Target thickness | > ion range (~10–100 nm at keV); ≤1 µm to allow fast loading | nm | Sputtered or evaporated films on a cooled substrate |
| D loading of target | Higher is better (+15 % observed) | — | In-situ electrochemical loading from the back (Thunderbird geometry) |

**Strongest critique.** Measured U_e is 3–10× theory and *unexplained* ("screening puzzle"); surface contamination and stopping-power systematics are serious. The enhancement is seen only with keV projectiles; extrapolating to thermal energies needs a static potential incompatible with observed D–D spacings (Section 1).

**Credence.**
- That enhanced screening of hundreds of eV is real at keV energies: **85 %**.
- That Czerski's threshold resonance / e⁺e⁻ channel is real: **35 %**.
- That any screening effect produces measurable fusion *without* keV projectiles: **<1 %**.

---

### 2.8 Heavy-electron and high-field ideas: nanotips and plasmonic hot spots

Claim: intense fields at tips, cracks or plasmonic hot spots modify the barrier. Quantitatively (Section 4 (other checks)):

1. **Scale.** The deuteron's Coulomb field is 5.8×10¹⁹ V/m at 5 fm and 5.8×10¹¹ V/m at 0.5 Å. The best tip field (~10¹⁰ V/m, near field evaporation) is **10⁻¹⁰** of the nuclear-scale field and ~2 % of the field at atomic spacing.
2. **No relative force.** For identical q/m, a uniform field moves only the centre of mass; tunnelling sees only the *gradient*. At a 10 nm tip (~10¹⁸ V/m²) that is 7×10⁷ V/m across 0.74 Å, a ~5 meV shift.
3. Even if the full 0.74 eV acted on the pair (U_eff = 47 eV), the rate changes ×3, on a rate 50–60 orders too small.
4. **Plasmonic hot spots** (10–100× enhancement) at CW intensities of 10⁴–10⁸ W/m² reach only ~10⁴–10⁷ V/m.
5. **What fields do achieve is ion acceleration:** 10 keV per µm at 10¹⁰ V/m. Pyroelectric tip fusion (Naranjo, Gimzewski & Putterman, Nature 434, 1115, 2005) is *hot* beam–target fusion, out of scope, but an excellent **neutron-detector calibration source**.

**Credence that high fields modify the barrier in a useful way: <0.5 %.**

---

### 2.9 Phonon and lattice resonance: Letts–Cravens and PdD dispersion

**Observation.** Two weak (~tens of mW) tunable lasers, p-polarized (surface-plasmon coupling), beat on an Au-coated, electrolytically loaded PdD cathode. Across 170 runs (2007–2008), excess heat responded at beat frequencies **8.3 THz (width 0.70), 15.3 THz (0.44) and 20.4 THz (0.68)** (elsewhere 8.2, 15.1, 20.8 THz). Hagelstein had predicted the optical-band edges, where compressional modes have zero group velocity, and assigns 8.3 THz to Γ and 15.3 THz to L. The 20 THz line was attributed to H contamination or D in Au vacancies.

**PdD phonon facts.**
- Pd acoustic branches reach ≲6–7 THz (Θ_D ≈ 274 K ⇒ 5.7 THz).
- The PdH optical fundamental is ≈55–60 meV (≈13.5–14.5 THz), from inelastic neutron scattering †.
- The PdD optical band lies at ≈35–40 meV at low wavevector (8.5–9.6 THz) †, dispersing upward at high loading.
- Anharmonicity makes the H/D frequency ratio differ from √2.

**Our check:** 20.8/15.1 = 1.377 and 20.4/15.3 = 1.333, compared with √2 = 1.414. Assigning the 20 THz line to PdH at the "L-point" is *roughly* consistent, but not exact.

**Predictions for hardware.** Excite at 8.3 ± 0.35 THz and 15.3 ± 0.22 THz for PdD (and ≈20.4 THz for PdH). Use p-polarized light on a thin Au overlayer (tens of nm) on highly loaded PdD. This can be combined with the Hagelstein vacancy recipe.

**Strongest critique.** The only replication found no excess heat (of the ~100 mW scale Letts reported) in 231 triggered trials. Tens of mW cannot coherently drive a macroscopic THz mode (Section 4 (other checks)). The 15.3 THz (63 meV) assignment lies above most PdD optical-band estimates †.

**Credence: 2 %.**

---

### 2.10 Mainstream critiques that bound all mechanisms

- **Leggett & Baym** (PRL 63, 191, 1989): *in equilibrium*, tunnelling to small separations is rigorously bounded by the Born–Oppenheimer potential; for D in metals the bound is "far too small". **Loophole:** non-equilibrium (flux, phonon pumping, ion tracks), which every serious theory invokes.
- **Koonin & Nauenberg** (Nature 339, 690, 1989): D₂ rate ~3×10⁻⁶⁴ s⁻¹; FP rates would need an effective electron-mass enhancement of ≳5–10×.
- **Huizenga's three "miracles"** (1992):
  1. a barrier-penetration enhancement of ~40–50 orders of magnitude;
  2. a change in dd branching from the ~50:50 n + ³He / p + t to ~100 % ⁴He;
  3. ⁴He produced without its 23.8 MeV γ-ray (the normal branching ratio is ~10⁻⁷).

  Our numbers: 1 W from conventional dd would emit **8.6×10¹¹ neutrons/s**, which is lethal and easily detected. 1 W from ⁴He at 23.85 MeV means 2.6×10¹¹ He/s. With normal branching it would also produce ~1.7×10⁵ γ-rays/s at 23.8 MeV.
- **Energetic-particle constraint** (Hagelstein, Naturwissenschaften 97, 345, 2010): the absence of secondary signals limits any FP-reaction alphas to ≲20 keV. This rules out TSC's 23.8 MeV alphas and requires any mechanism to hand ~24 MeV to the lattice.
- **Momentum/energy conservation:** d + d → ⁴He needs a third body or a γ. "The lattice takes it" needs MeV→meV coupling unknown to established physics; Czerski's e⁺e⁻/electron emission is the only mainstream-published third-body candidate.
- **Reproducibility:** DOE reviews (1989, 2004) unpersuaded; Berlinguette et al. (Nature 570, 45, 2019) found no anomalies and showed how hard D/Pd > 0.9 is; the Letts replication failed.
- **Metzler, Hunt, Hagelstein & Galvanetto** (New J. Phys. 26, 101202, 2024): observable solid-state dd fusion needs >40 orders of enhancement; known mechanisms give up to ~30 orders each, so cascades are *not excluded in principle*, but none is demonstrated.

**Credence that these bounds correctly exclude LENR at the claimed watt level in near-equilibrium systems: ~97 %.**

---

## 3. Consolidated geometry-prediction matrix

Abbreviations: **Hag** = Hagelstein; **TSC** = Takahashi; **WL** = Widom–Larsen; **Sto** = Storms; **Kim** = Kim BEC; **Coh** = Preparata/Chubb; **Scr** = screening family (Czerski, Raiola, Lipson, Berlinguette); **HF** = high-field; **Let** = Letts/THz. Credence is shown in brackets in the header.

| Geometric feature | Hag [4%] | TSC [1%] | WL [1%] | Sto [1.5%] | Kim [<1%] | Coh [<1%] | Scr [85% keV / <1% thermal] | HF [<0.5%] | Let [2%] |
|---|---|---|---|---|---|---|---|---|---|
| Nanoparticle diameter | Small particles/codeposits OK (vacancy-rich) | **2–10 nm** Pd/Ni in ZrO₂ | 10–100 nm (plasmonic) | n/a (interparticle gaps) | **few nm** (smaller better) | finite ordered nanocrystals (size unspecified) | n/a | n/a | n/a |
| Crack/gap width | n/a | SNH **0.3–1 nm** | nm–µm patches | **~1 nm (<10 nm)** | n/a | n/a (defects harmful) | defects raise U_e | n/a | n/a |
| Surface roughness wavelength | n/a | n/a | **10–100 nm** | n/a | n/a | n/a | smooth, clean surfaces for systematics | tip radius 1–100 nm | Au overlayer for surface plasmons (optical λ 0.6–0.8 µm) |
| Grain size | n/a (vacancies dominate) | nm | n/a | small grains → more cracks | nm | **large, perfect** (conflict) | more grain boundaries/defects → higher U_e (plausible) | n/a | n/a |
| Film thickness | ≤1 µm (fast loading to x ≥ 0.9) | n/a | surface layer | µm surface layer | n/a | n/a | **> ion range (10–100 nm), ≤1 µm** | n/a | Au tens of nm on PdD |
| Tip radius | n/a | n/a | sharp features help | n/a | n/a | n/a | n/a | 1–100 nm (only accelerates ions) | n/a |
| Vacancy concentration | **High; up to ~25 % at x ≥ 0.95** | local defects (SNH) | n/a | voids yes, vacancies insufficient | n/a | harmful | **higher → higher U_e** | n/a | vacancies in Au (20 THz) |
| Local D/M loading | **≥ 0.90–0.95** | **≈1 locally** | high flux more than loading | proxy for cracking | high | **→1** | higher → +15 % | n/a | high |
| Mechanical/phonon frequency | **8.3 & 15.3 THz** optical (PdD); MHz–GHz acoustic (host nuclei, Karabut) | phonon-driven O→T motion (optical, ~8 THz implied) | proton oscillation ~10¹⁴ s⁻¹ (~60 meV) | chain axial resonance (unquantified; our guess 10–100 THz) | n/a | coherent plasma modes | ion-track phonons (Czerski) | n/a | **8.3, 15.3, 20.4 THz** |
| Cavity/plasmon resonance | n/a | n/a | **surface-plasmon resonance (1.5–3.5 eV)** | n/a | n/a | CD ~0.1 µm (≈10 eV mode) | n/a | plasmonic hot spots (useless) | surface-plasmon coupling via p-polarization |
| Temperature | neutral (loading falls with T) | **200–350 °C** | n/a | higher | **lower** | low (band coherence) | **lower (Debye model)** | n/a | room temperature |
| Drive / flux | current > threshold, D flux | desorption bursts | high H flux | loading cycles | n/a | sustained field (Chubb) | keV D⁺ + in-situ loading | ion acceleration | two lasers |
| H vs D | H poisons D₂→⁴He | both (4H/TSC) | both | both, different products | D | D | d+d, p+d | — | D (H at 20 THz) |

**Agreement (robust targets):** (1) near-surface (≤1 µm), vacancy/defect-rich Pd, Pd–Ni or Zr (Hag, Sto, TSC, Scr, Let); (2) local D/M ≳ 0.9 (Hag, TSC, Coh, Scr; Sto indirectly); (3) 2–10 nm features with ~0.3–1 nm gaps/voids (TSC, Kim, Sto; WL loosely); (4) non-equilibrium drive (every theory escaping Leggett–Baym); (5) 8–15 THz PdD lattice modes as the stimulus band (Hag, Let; TSC implicitly).

**Conflict (each is a discriminating experiment):** perfect periodicity (Chubb) vs defects (all others); low temperature (Kim, Debye) vs 200–350 °C (Takahashi, Storms); H as poison (Hagelstein) vs fuel (WL, TSC, Storms); PdO surface helps (Lipson) vs oxidation as systematic error (Czerski).

---

## 4. Quantitative sanity checks

Computed with `python3`/numpy (outputs copied below; core of 4.1–4.2 reproduced here). Model: shifted-energy Gamow/WKB, λ = (A/V_conf)·exp(−2πη(U_e)), A = S(0)/(πμcα), S(0) = 110 keV·b (both dd branches), V_conf = (4/3)π(0.74 Å)³. The prefactor is uncertain by many orders, so required U_e is quoted for prefactors ×10⁻¹⁰ to ×10¹⁰.

```python
import numpy as np
hbar,c,e,amu,alpha=1.0546e-34,2.998e8,1.602e-19,1.6605e-27,1/137.036
mu=2.01410*amu/2; S=110e3*e*1e-28; A=S/(np.pi*mu*c*alpha)       # 1.5e-22 m^3/s
V=4/3*np.pi*(0.74e-10)**3
G=lambda E_eV: 2*np.pi*alpha*c/np.sqrt(2*E_eV*e/mu)              # 2*pi*eta
lam=lambda Ue: A/V*np.exp(-G(Ue))                                 # per pair, s^-1
Ue_needed=lambda lt,pref=A/V: 0.5*mu*(2*np.pi*alpha*c/np.log(pref/lt))**2/e
f_beam=lambda E_eV,Ue: np.exp(0.5*G(E_eV)*Ue/E_eV)               # beam enhancement
print(lam(300), Ue_needed(2/6.8e22), f_beam(1e3,300))           # 1e-17, ~200 eV, ~111
```

**4.1 Screening needed at thermal energy (1 cm³ PdD, 6.8×10²² D).**

| U_e (eV) | 2πη | λ per pair (s⁻¹) |
|---|---|---|
| 27.7 | 188.7 | 1×10⁻⁷⁴ |
| 100 | 99.3 | 7×10⁻³⁶ |
| 200 | 70.2 | 3×10⁻²³ |
| 300 | 57.3 | 1×10⁻¹⁷ |
| 500 | 44.4 | 5×10⁻¹² |
| 1000 | 31.4 | 2×10⁻⁶ |

| Target | All sites active: U_e needed | 10⁻⁶ of sites active: U_e needed |
|---|---|---|
| 1 n/s | 200 eV (113–443) | 310 eV (156–886) |
| 1 mW via ⁴He | 372 eV (177–1215) | 694 eV (267–4583) |
| 1 W via ⁴He | 496 eV (216–2118) | 1040 eV (341–16 376) |

The turning point for a static U_e is 52 pm at 27.7 eV and 7 pm at 200 eV. A static potential of this kind is incompatible with measured D–D spacings. *Conclusion:* screening near the measured values would suffice *if* it were available to thermal pairs. That caveat is the whole debate.

**4.2 Beam regime f = exp(πη·U_e/E).**

| E (keV) | U_e = 25 eV | U_e = 300 eV | U_e = 800 eV |
|---|---|---|---|
| 1 | 1.48 | 111 | 2.9×10⁵ |
| 2 | 1.15 | 5.3 | 85 |
| 5 | 1.04 | 1.52 | 3.1 |
| 10 | 1.01 | 1.16 | 1.49 |

The anomalous signal is therefore largest for 1–3 keV projectiles. This sets the design energy of the recommended experiment.

**4.3 Branch bookkeeping per watt.** Conventional dd gives 1.7×10¹² reactions/s, i.e. **8.6×10¹¹ n/s**. ⁴He gives 2.6×10¹¹ /s. At a branching ratio of 10⁻⁷, that is 1.7×10⁵ γ-rays/s at 23.8 MeV. A heat-only claim of 1 W with <10 n/s implies branching distortion ≥10¹¹.

**4.4–4.10 Other checks** (used in Section 2):

| Check | Result |
|---|---|
| Coulomb field of d | 5.8×10¹⁹ V/m (5 fm), 1.4×10¹⁵ (1 pm), 5.8×10¹¹ (0.5 Å); atomic unit 5.1×10¹¹ V/m |
| 10¹⁰ V/m nanotip | 0.74 eV across 0.74 Å; *differential* field 7×10⁷ V/m ⇒ ~5 meV pair shift; ×3.1 rate even if the full 0.74 eV acted at U_eff = 47 eV; D⁺ gains 10 keV/µm (hot fusion) |
| Widom–Larsen | β = 2.531 ⇒ a₀ = 2.33 ⇒ 3.6×10¹¹ V/m at 60 meV (≈1.7×10¹⁶ W/cm² equiv.), 9.3×10¹² V/m at 1.55 eV (1.2×10¹⁹ W/cm²); screened lattice field ~9×10⁵ V/m ⇒ β − 1 ≈ 10⁻¹¹ |
| TSC localization | K_e ≈ 9.4 MeV at 20 fm (×4 = 37 MeV); 37 keV at 1 pm; 0.4 keV at 10 pm |
| Fractionation | quanta per 23.85 MeV: 6.95×10⁸ (8.3 THz), 3.8×10⁸ (15.3 THz), 2.6×10¹⁵ (2.21 MHz); n_thermal(300 K) = 0.36 (8.3 THz) |
| Coherent THz drive | 0.01 Å at 8.3 THz in 1 cm³ stores 0.31 J (5.6×10¹⁹ phonons) and costs 1.6×10¹¹ W (Q = 100), 1.6×10⁹ W (Q = 10⁴); 0.001 Å costs 1.6×10⁹ W (Q = 100). External drive can pump only nm–µm volumes |
| BEC | T_c = 13.3 K (m* = m_D), 0.13 K (100 m_D), 1.3 mK (10⁴ m_D); nλ³ = 0.024 (300 K), 0.19 (77 K) versus 2.61 |
| Letts ratios | 20.8/15.1 = 1.377, 20.4/15.3 = 1.333 versus √2 = 1.414; 8.3/15.3/20.4 THz = 34.3/63.3/84.4 meV |
| Length scales | PdD a ≈ 4.03 Å, O–O 2.85 Å, O→T 1.75 Å; H₂ bond 0.74 Å, kinetic diameter 2.89 Å; "few atomic diameters" ≈ 8–11 Å (Storms' ~1 nm); surface-plasmon resonance 1.5–3.5 eV = 354–827 nm; 10 eV coherence domain ≈ 124 nm |

---

## 5. Design implications, ranked by expected value

Expected value (EV) ≈ P(real) × P(geometry matters | real) × P(credible detection | real) + information value. "Credible" = a nuclear product (neutrons, charged particles, 511 keV, ⁴He/³He, tritium) measured with blanks and calibrations; heat alone does not qualify.

1. **Low-energy D⁺ implantation (0.3–5 keV) into a D-loaded, defect-engineered thin-film target. Highest EV.**
   - *Geometry:* 100 nm–1 µm Pd, Zr/ZrD₂ or PdO-capped Pd film on a membrane loaded from the back (electrochemically or by gas; Thunderbird-style). Vary vacancy density (10⁻⁴–10⁻² site fraction via pre-irradiation or codeposition) and grain size (nm vs µm); UHV-clean and PdO-capped regions side by side; stage at 77–400 K.
   - *Diagnostics:* Si detectors (3.02 MeV p, 1.01 MeV t, 0.82 MeV ³He); ³He/BF₃ or EJ-309 neutron counter (2.45 MeV); NaI/HPGe for 511 keV and energetic electrons (Czerski e⁺e⁻ channel).
   - *Why:* the nuclear signal is guaranteed, and the *anomalous* enhancement over 25–100 eV theory is measured against exactly the variables Hagelstein, Storms, Takahashi and Czerski name. P(real anomaly) ≈ 0.85; P(informative geometry dependence) ≈ 0.3.
   - *Caveat:* this is "lattice-assisted" beam–target fusion, not thermal cold fusion; the project lead should confirm scope. We recommend it as the calibrated backbone against which any LENR claim is benchmarked.
2. **Vacancy-rich, high-loading Pd surface layer driven far from equilibrium, with nuclear diagnostics.** This is the consensus target of Hag, Sto, TSC and Scr.
   - *Geometry:* Pd/D codeposited film (≤1–5 µm) on Au or Cu, or thin foil loaded to D/Pd ≥ 0.9 (resistance-verified); current density swept 0.05–1 A/cm² in fine steps (Swartz OOP); loading/deloading cycles for ~1 nm cracks.
   - *Diagnostics:* in-situ Si detectors (not CR-39 alone), neutron counters, tritium assay, ⁴He mass spectrometry that resolves D₂ from ⁴He. P(real) ≈ 0.04; low cost.
3. **2–10 nm Pd or Pd–Ni in ZrO₂ nanocomposites under D₂ gas at 200–350 °C.** TSC, Kim and Storms (gaps between particles) agree here.
   - *Geometry:* re-calcined Pd₁Ni₁₀/ZrO₂ or Cu₁Ni₇/ZrO₂ with particle size (2, 5, 10 nm) as a variable; H₂ vs D₂ as a control (theories disagree). *Diagnostics:* ⁴He in gas, neutrons, γ; calorimetry only secondary. P(real) ≈ 0.02.
4. **Phonon–nuclear coupling test with radioactive tracers (Hagelstein–Metzler), independent of fusion.**
   - *Geometry:* Co-57 bonded to a steel/Fe plate or foil; piezo drive at 1–3 MHz on a high-Q acoustic resonance; matched *undriven* twin; detectors at several angles for the predicted anisotropy/delocalization of the 14.4 keV line. Clean counting statistics make a replicated positive the most credible LENR-relevant evidence possible; guard against piezo heating, microphonics, pile-up. P(real) ≈ 0.03, high information per dollar.
5. **THz stimulation of loaded PdD at 8.3 ± 0.35 and 15.3 ± 0.22 THz** (≈20.4 THz as the PdH control).
   - Add-on to item 2: Au overlayer ~20–50 nm, p-polarized dual lasers or a difference-frequency THz source, scanned with off-resonance controls. Failed replication once; P ≈ 0.02.
6. **Engineered ~0.3–1 nm nanogaps** (Storms): electromigrated break-junction arrays, dealloyed nanoporous Pd, ALD-spaced gaps, with gap width as the swept variable. P ≈ 0.015. The same samples double as high-defect targets for item 1.
7. **Plasmonic roughness and nanotips** (Widom–Larsen, high-field): P < 0.01 as barrier modifiers; use a tip or pyroelectric crystal only as a cheap **hot-fusion neutron calibration source**.

**Cross-cutting rules implied by the theory review:**
- Detect dd-branch products (p, t, n) *and* ⁴He/511 keV: screening mechanisms predict normal branching (neutrons are the easiest credible signal); Hagelstein/TSC predict ⁴He without energetic particles (needs heat–⁴He correlation).
- Use H/D substitution as a built-in control that splits the theories.
- Log loading (R/R₀), temperature, flux and defect density as *primary* variables; every theory predicts steep thresholds.

---

## 6. References

Primary sources could not be fetched in this environment; links are as found through search. † in the text marks values to verify.

**Hagelstein and co-workers**
- Hagelstein, "Current status of the theory and modeling effort based on fractionation", JCMNS 19, 98 (2016). https://dspace.mit.edu/bitstream/handle/1721.1/122637/ICCF19review.pdf ; https://jcmns.org/article/72378-current-status-of-the-theory-and-modeling-effort-based-on-fractionation
- Hagelstein & Chaudhary, "Energy exchange in the lossy spin-boson model", JCMNS 5, 52 (2011). https://jcmns.org/article/72149-energy-exchange-in-the-lossy-spin-boson-model/attachment/150135.pdf
- Hagelstein & Chaudhary, "Including nuclear degrees of freedom in a lattice Hamiltonian", arXiv:1201.4377. https://arxiv.org/abs/1201.4377
- Hagelstein, "Coupling between a deuteron and a lattice", arXiv:1204.2159. https://arxiv.org/pdf/1204.2159
- Hagelstein & Chaudhary, "Models for phonon–nuclear interactions and collimated X-ray emission in the Karabut experiment". https://www.semanticscholar.org/paper/86deac7b9c9e09bf4c9358ff2af474ec58380692 ; US patent 10660191. https://patents.google.com/patent/US10660191B1/en
- Hagelstein, "Theory and experiments in condensed matter nuclear science", JCMNS 35, 49 (2022). https://jcmns.org/article/72582-theory-and-experiments-in-condensed-matter-nuclear-science
- Hagelstein, "Relativistic phonon–nuclear coupling matrix element for the D₂/⁴He transition", JCMNS (2025). https://jcmns.org/article/134000-relativistic-phonon-nuclear-coupling-matrix-element-for-the-d2-4he-transition
- Hagelstein, "Arguments for dideuterium near monovacancies in PdD" (2009). https://newenergytimes.com/v2/library/2009/2009Hagelstein-ArgumentsForDideuterium.pdf
- Hagelstein, "Constraints on energetic particles in the Fleischmann–Pons experiment", Naturwissenschaften 97, 345 (2010). https://www.researchgate.net/publication/41417374
- Metzler, Hagelstein & Lu, "Observation of non-exponential decay in X-ray and γ emission lines from Co-57", JCMNS 27, 46 (2018). https://jcmns.org/article/72479
- Metzler & Hagelstein, "Developing phonon–nuclear coupling experiments with vibrating plates and radiation detectors". https://jcmns.org/article/72439
- "Phonon-mediated nuclear excitation transfer" (ICCF21). https://dspace.mit.edu/bitstream/handle/1721.1/122616/ICCF21ExTrans.pdf
- Patent WO2019236455A1. https://patents.google.com/patent/WO2019236455A1/en
- Hagelstein & Tanzella, vibrating-copper experiment (ICCF19 summary). https://coldfusionnow.org/hagelstein-and-tanzellas-vibrating-copper-experiment/
- Metzler, Hunt, Hagelstein & Galvanetto, "Known mechanisms that increase nuclear fusion rates in the solid state", New J. Phys. 26, 101202 (2024). https://iopscience.iop.org/article/10.1088/1367-2630/ad091c ; https://arxiv.org/abs/2208.07245

**Takahashi and co-workers**
- Takahashi, "Physics of cold fusion by TSC theory", JCMNS. https://jcmns.org/api/v1/articles/72286-physics-of-cold-fusion-by-tsc-theory.pdf
- Takahashi, "Nuclear products of cold fusion by TSC theory". https://jcmns.org/article/72320.pdf
- Takahashi, "Time-dependent EQPET analysis of TSC". https://www.osti.gov/etdeweb/biblio/21061111
- Kitamura, Takahashi et al., "Anomalous heat effects induced by metal nano-composites". https://jcmns.org/article/72496.pdf ; "Brief summary report of MHE project". https://www.researchgate.net/publication/322160963

**Widom–Larsen and critiques**
- Widom & Larsen, "Ultra low momentum neutron catalyzed nuclear reactions on metallic hydride surfaces", Eur. Phys. J. C 46, 107 (2006). https://arxiv.org/abs/cond-mat/0505026
- Ciuchi, Maiani, Polosa, Riquer, Ruocco & Vignati, Eur. Phys. J. C 72, 2193 (2012). https://link.springer.com/article/10.1140/epjc/s10052-012-2193-9 ; https://arxiv.org/abs/1209.6501
- Maiani, Polosa & Riquer, Eur. Phys. J. C 74, 2843 (2014). https://link.springer.com/article/10.1140/epjc/s10052-014-2843-1
- Tennfors, Eur. Phys. J. Plus 128, 15 (2013). DOI 10.1140/epjp/i2013-13015-3
- Hagelstein & Chaudhary, "Electron mass shift in nonthermal systems", J. Phys. B 41, 125001 (2008). https://arxiv.org/abs/0801.3810 ; Widom–Larsen reply: https://arxiv.org/pdf/0802.0466
- Hagelstein, "Electron mass enhancement and the Widom–Larsen model", JCMNS. https://jcmns.org/article/72227-electron-mass-enhancement-and-the-widom-larsen-model

**Storms**
- Storms, "How basic behavior of LENR can guide a search for an explanation", JCMNS 20, 100 (2016). https://jcmns.org/article/72408
- Storms, "The role of voids as the location of LENR", JCMNS 11, 123 (2013). https://jcmns.org/article/72221-the-role-of-voids-as-the-location-of-lenr
- Storms, "Explaining cold fusion", JCMNS 15, 295 (2015). https://jcmns.org/article/72347-explaining-cold-fusion.pdf
- Storms, *The Explanation of Low Energy Nuclear Reaction* (Infinite Energy Press, 2014).

**Kim, Chubb, Preparata, Schwinger, Swartz**
- Kim, "Theory of Bose–Einstein condensation mechanism for deuteron-induced nuclear reactions in micro/nano-scale metal grains and particles", Naturwissenschaften 96, 803 (2009). https://link.springer.com/article/10.1007/s00114-009-0537-6
- Kim, "Bose–Einstein condensate theory of deuteron fusion in metal", JCMNS 4, 188 (2011). https://jcmns.org/article/72132
- Chubb & Chubb, "Cold fusion as an interaction between ion band states", Fusion Technol. 20, 93 (1991). https://www.tandfonline.com/doi/abs/10.13182/FST91-A29646 ; "Context for understanding why particular nanoscale crystals turn-on faster". https://www.osti.gov/etdeweb/biblio/21061109
- Bressani, Del Giudice & Preparata, Nuovo Cimento A 101, 845 (1989); Preparata, *QED Coherence in Matter* (World Scientific, 1995). https://link.springer.com/article/10.1007/BF02844877
- Chechin et al., "Critical review of theoretical models for anomalous effects in deuterated metals". https://arxiv.org/abs/nucl-th/0303057
- Schwinger, "Cold fusion: a hypothesis", Z. Naturforsch. 45a, 756 (1990).
- Swartz, "LANR nanostructures and metamaterials driven at their optimal operating point", JCMNS 6, 149 (2012). https://jcmns.org/article/72167

**Electron screening**
- Raiola et al., "Electron screening in d(d,p)t for deuterated metals and the periodic table", Phys. Lett. B (2002). https://www.sciencedirect.com/science/article/pii/S0370269302027740
- Raiola et al., "Enhanced electron screening in d(d,p)t for deuterated Ta", Eur. Phys. J. A (2002). https://link.springer.com/article/10.1007/s10050-002-8766-5
- Raiola et al., "Electron screening in d(d,p)t for deuterated metals: temperature effects", J. Phys. G 31 (2005). https://iopscience.iop.org/article/10.1088/0954-3899/31/11/002
- Huke, Czerski et al., "Enhancement of deuteron-fusion reactions in metals and experimental implications", Phys. Rev. C 78, 015803 (2008). https://link.aps.org/doi/10.1103/PhysRevC.78.015803
- Czerski et al., J. Phys. G 35, 014012 (2008). https://iopscience.iop.org/article/10.1088/0954-3899/35/1/014012
- Czerski et al., "Deuteron–deuteron nuclear reactions at extremely low energies", Phys. Rev. C 106, L011601 (2022). https://link.aps.org/doi/10.1103/PhysRevC.106.L011601
- Czerski et al., "Indications of electron emission from the deuteron–deuteron threshold resonance", Phys. Rev. C 109, L021601 (2024). https://link.aps.org/doi/10.1103/PhysRevC.109.L021601
- Dubey, Czerski et al., Phys. Rev. X 15, 041004 (2025). https://link.aps.org/doi/10.1103/chlp-b215 ; https://arxiv.org/abs/2408.07567
- Czerski et al., "Observation of thermal deuteron–deuteron fusion in ion tracks". https://arxiv.org/abs/2409.02112 ; "Thermal deuteron–deuteron fusion in metallic targets". https://arxiv.org/abs/2605.27438
- Threshold resonance, Acta Phys. Pol. B 48, 489 (2017). https://www.actaphys.uj.edu.pl/R/48/3/489/pdf
- Lipson et al., screening in Pd and Au/Pd/PdO. https://www.osti.gov/etdeweb/biblio/20172333
- "Electron screening in laboratory nuclear reactions", Particles 7, 50 (2024). https://www.mdpi.com/2571-712X/7/3/50
- Berlinguette et al., "Electrochemical loading enhances deuterium fusion rates in a metal target", Nature (2025). https://www.nature.com/articles/s41586-025-09042-7
- Ichimaru, "Nuclear fusion in dense plasmas", Rev. Mod. Phys. 65, 255 (1993).

**Phonons, lasers, high fields**
- Letts, Cravens & Hagelstein, "Terahertz difference frequency response of PdD in two-laser experiments". https://www.researchgate.net/publication/265804225
- Letts, "Stimulation of optical phonons in deuterated palladium". https://www.lenr-canr.org/acrobat/LettsDstimulatio.pdf
- "Attempted replication of excess heat in the Letts dual-laser experiment", JCMNS. https://jcmns.org/article/72403.pdf
- Naranjo, Gimzewski & Putterman, "Observation of nuclear fusion driven by a pyroelectric crystal", Nature 434, 1115 (2005). DOI 10.1038/nature03575
- Fukai & Okuma, "Formation of superabundant vacancies in Pd hydride under high hydrogen pressures", Phys. Rev. Lett. 73, 1640 (1994).

**Mainstream critiques and reviews**
- Leggett & Baym, Phys. Rev. Lett. 63, 191 (1989). https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.63.191
- Koonin & Nauenberg, "Calculated fusion rates in isotopic hydrogen molecules", Nature 339, 690 (1989). DOI 10.1038/339690a0
- Huizenga, *Cold Fusion: The Scientific Fiasco of the Century* (Univ. Rochester Press, 1992).
- U.S. DOE, *Report of the Review of Low Energy Nuclear Reactions* (2004).
- Berlinguette et al., "Revisiting the cold case of cold fusion", Nature 570, 45 (2019). DOI 10.1038/s41586-019-1256-6
