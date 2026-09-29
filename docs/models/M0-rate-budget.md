# M0 — First-Principles Rate Budget

Code: [`sim/m0_rate_budget.py`](../../sim/m0_rate_budget.py) · Raw output: [`figs/m0_output.txt`](figs/m0_output.txt)

This model answers one question before any geometry is drawn: **how large an effect has to exist for any detector to see it, and which geometric levers actually move the rate?**

## 1. Framework

For two deuterons the fusion rate is set by Coulomb-barrier penetration.

- Gamow energy: $E_G = 2\mu c^2(\pi\alpha Z_1Z_2)^2 = 985.8$ keV for D+D.
- Cross-section: $\sigma(E) = \frac{S}{E}\,e^{-\sqrt{E_G/E}}$ with $S_{dd,n}\approx55$, $S_{dd,p}\approx57$ keV·b.
- Bound pair (molecule, lattice pair): $\lambda = A\,\rho_0\,P$, where $A = \dfrac{S\,c}{\pi\alpha\,\mu c^2} = 1.56\times10^{-16}$ cm³/s, $\rho_0$ is the pair probability density where tunnelling starts, and $P$ is the barrier penetration.
- We express every environment through one lumped number, the **effective screening energy** $U_{e,\mathrm{eff}}$, defined by $P \equiv \exp(-\sqrt{E_G/U_{e,\mathrm{eff}}})$. It absorbs electron screening *and* how close the lattice lets the pair sit.

**Calibration.** The D₂ molecule (Koonin & Nauenberg 1989: λ ≈ 3×10⁻⁶⁴ s⁻¹) corresponds to $U_{e,\mathrm{eff}} = 34$ eV with $\rho_0 = 1.5\times10^{26}$ cm⁻³ from the D₂ zero-point width.

**Caution: $U_{e,\mathrm{eff}}$ is not the accelerator $U_e$.** Accelerator experiments measure screening at keV energies (turning points ~100 fm). For a thermal pair, a Yukawa-shaped screening of accelerator strength $U_e$ gives a much weaker effective value (exact WKB):

| Accelerator-style $U_e$ (eV) | 100 | 300 | 500 | 800 | 1000 | 2000 |
|---|---|---|---|---|---|---|
| $U_{e,\mathrm{eff}}$ for a freely-approaching thermal pair (eV) | 41 | 123 | 206 | 331 | 414 | 841 |

## 2. What each detector needs (14-day run, 5σ)

| Channel | Assumptions | Reactions/s needed |
|---|---|---|
| Neutrons (2.45 MeV) | ³He bank, ε = 10 %, bkg 0.05 cps | **2×10⁻²** |
| Protons (3.02 MeV), Si in vacuum | 5 % solid angle × 30 % escape, bkg 2×10⁻⁵ cps in window | **3×10⁻³** |
| Heat, standard D+D branching | 10 mW resolution | **1.7×10¹⁰** |

- For standard D+D, **heat is ~10¹² times less sensitive than particle detection.** A heat-first design is the wrong optimisation for H1.
- 1 W of standard D+D means 8.4×10¹¹ n/s, or **~10 Sv/h at 1 m**. Watt-level heat claims without lethal neutrons are only possible through a non-standard channel (H2). That is why H2 needs its own signature: heat correlated with ⁴He.

## 3. Required enhancement

Required $U_{e,\mathrm{eff}}$ (eV) to reach the neutron threshold, as a function of the number $N$ of active pairs:

| N pairs | ρ₀=10²⁴ | ρ₀=10²⁵ | ρ₀=10²⁶ |
|---|---|---|---|
| 10⁹ | 521 | 470 | 426 |
| 10¹² | 388 | 355 | 326 |
| 10¹⁵ | 300 | 277 | 257 |
| 10¹⁸ | 239 | 223 | 208 |
| 10²¹ (≈ whole 0.5 g cathode) | 195 | 183 | 172 |

**Exponent budget.** Molecular D₂ has penetration exponent X = 170. A neutron-detectable rate needs X ≈ 53–73. The gap is **97–118 e-folds**. The engineering levers are small by comparison:

- 10⁹× more active sites buys 20.7 e-folds.
- 100× better detection buys 4.6 e-folds.

**Null-limit constraint.** Bulk PdD is known not to fuse at ≳10⁻²⁵ /pair/s (1989–90 null searches; to be verified in R5). In the free-gas picture this forces $U_{e,\mathrm{eff}}(\text{thermal}) \lesssim 147$ eV in bulk. So the large accelerator screening values do not carry over to thermal pairs in ordinary bulk metal, and bulk averages are not where to look.

## 4. Steepness and tail dominance: the key geometric insight

$\dfrac{d\ln\lambda}{d\ln U_{e,\mathrm{eff}}} = \tfrac12\sqrt{E_G/U_{e,\mathrm{eff}}}$, which is **35–50** in the relevant range. Locally the rate scales as $U^{35\text{–}50}$, so +10 % enhancement gives 28–113× the rate.

If enhancement varies from site to site (Gaussian, mean 100 eV):

| σ (eV) | share of all reactions from the top 0.1 % of sites |
|---|---|
| 5 | 16 % |
| 10 | 59 % |
| 20 | 90 % |

**Consequence for geometry.** The rate is set by the rarest, most extreme sites (defects, crack faces, vacancy clusters, interfaces, surface asperities, nanoparticles), not by bulk averages. The design objective is therefore:

1. **Multiply** candidate extreme sites: high interface/defect density and high-loading nanostructure.
2. **Place them within the escape depth** of charged products (tens of µm in Pd for 3 MeV protons) facing the detectors.
3. **Load them to the highest possible deuterium chemical potential**, with flux and cycling (non-equilibrium widens the site distribution).
4. **Keep everything else identical in a light-hydrogen control.**

## 5. What this rules in and out

- **Out:** designs whose only observable is bulk excess heat (insensitive under H1, artifact-prone under H2).
- **Out:** very thick bulk cathodes as the primary sensor (sites deep in the metal are invisible to charged-particle detectors, and loading is slow).
- **In:** thin, interface-rich, high-loading active layers directly facing low-background charged-particle spectrometers, surrounded by neutron detection; plus a calorimetric/⁴He channel to test H2.

## 6. Limitations

- ρ₀ and the potential shape are uncertain by orders of magnitude. The exponent dominates, and conclusions change by less than ±15 % in required $U_{e,\mathrm{eff}}$ across ρ₀ = 10²⁴–10²⁶.
- The model ignores resonances (for example Czerski's proposed near-threshold 0⁺ state in ⁴He). A resonance would multiply $S$, and so λ, linearly by the enhancement factor. Even 10⁴× equals only 9 e-folds, so it does not change the conclusion but does raise the value of charged-particle and e⁺e⁻ (511 keV) detection. R3 will quantify this.
- Non-thermal local energies (fracto-emission fields at crack tips, desorption transients) could raise the relative energy $E$ from 0.03 eV to eV–keV. This is the one "cold" lever that acts directly on the exponent. It is modelled in M1.

![rate vs Ue](figs/m0_rate_vs_Ue.png)
