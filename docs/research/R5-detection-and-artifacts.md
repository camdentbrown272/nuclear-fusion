# R5 — Measurement methods, detector geometry, backgrounds and the artifact catalog for LENR experiments

**Scope.** How to measure a claimed LENR (cold fusion) effect so that the result holds up under expert review: calorimetry, neutrons, charged particles, tritium, ⁴He, gammas and X-rays, and transmutation. For each: sensitivity, backgrounds, geometry, cost and known artifacts. Also the statistical design and a recommended integrated detector layout for a single cell.

**Where the numbers come from (read this first).** Only limited web checking was possible in this session: most publisher and preprint domains were blocked, and the search quota ran out after a few queries. Values are tagged as follows:
- **[C]**: computed here and reproducible. Ion ranges and residual energies come from the ATIMA/`pycatima` stopping-power code (SRIM-type low-energy stopping). Poisson, Currie and Feldman–Cousins numbers were computed with numpy/scipy. Isotope masses, reaction Q-values and simple transport estimates are from first principles.
- **[V]**: literature or vendor values quoted from domain knowledge. Check each against the cited primary source before relying on it for a design decision.
- Untagged values are standard textbook constants.

---

## 0. Executive summary

1. **Use the per-watt yardstick.** One watt of d+d fusion by the known branches would give about 8.6×10¹¹ n/s, roughly **10 Sv/h at 1 m** [C]. By the claimed d+d→⁴He channel it would give **2.6×10¹¹ ⁴He/s** [C]. Every LENR heat claim therefore implies either a lethal neutron field (never observed) or a helium production rate that has to be measured quantitatively.
2. **Tritium and neutrons are extremely sensitive but, historically, the wrong scale.** 1 mW·day of the t+p branch in 50 mL of electrolyte gives about 4.8×10³ Bq/mL [C], against a liquid-scintillation (LSC) detection limit of about 0.15 Bq/mL. Neutron upper limits of about 0.01 n/s at 50 mW mean that neutrons are suppressed by ≥10¹³ relative to heat [C]. Use these channels to set upper limits and branching ratios, not as the main evidence for heat.
3. **Helium-4 is the only claimed product expected at a level comparable to the heat, and it has the worst background: air, at 5.24 ppm.** One W·day of ⁴He equals the helium in **0.16 cm³ of air** [C]. The key geometric lever is a **sealed all-metal cell with a small headspace (≤30 cm³)**. In that volume, 100 mW for 2–4 days pushes the headspace ⁴He **above the ambient 5.24 ppm**, which in-leaking air can never do [C].
4. **Calorimetry must be position-independent and done in a closed cell.** Use a closed cell with an internal recombiner and a pressure transducer. Use a Seebeck-envelope or mass-flow calorimeter whose calibration changes by less than 0.2% when a heater is moved between the cathode, recombiner and headspace positions. Aim for σ ≈ max(5 mW, 0.3% of P_in). Keep input power low (≤2–5 W), because most errors scale with P_in. This design addresses Shanahan's calibration-constant-shift critique and the recombination (Jones/Hansen) critique directly.
5. **A heat claim needs a size criterion as well as a significance criterion.** 100 eV per Pd atom in a 1 g cathode is 25 Wh, i.e. about 10 days at 100 mW [C]. That is roughly 600× the enthalpy of loading the Pd with deuterium.
6. **Neutrons: use a ³He well counter, not single tubes or bubble detectors.** Specifically: a 4π HDPE well counter with 12 or more ³He tubes behind about 4 cm HDPE-equivalent moderator, with a Cd liner, borated-PE shield and plastic-scintillator muon veto, in a basement. Expected efficiency is 15–35% [V]. The detection limit is about 0.01 n/s per day at 0.05 cps background [C]. **Systematic stability of the background, not counting time, sets the floor**: a 2% background systematic at 0.05 cps limits discovery to ≥0.017 n/s regardless of run time [C]. That is why a twin control station and cell swaps are mandatory.
7. **CR-39 is cheap but has the weakest credibility of the nuclear channels.** Use a **differential Mylar filter stack** on each chip (0/6/25/60/150 µm). The windows separate p, t, ³He and α [C]. The 150 µm window is an internal null: it stops 3 MeV protons, so it should show only background. Note that the 6 µm Mylar used in past SPAWAR-type work stops ³He but passes tritons and all radon alphas [C]. Count every chip blind.
8. **Gamma spectra: HPGe resolution matters.** In NaI or LaBr₃, the ²¹⁴Bi line at **2204 keV** (radon daughter) cannot be separated from the **2223 keV** n-capture-on-hydrogen line. LaBr₃ also has an internal ²²⁷Ac alpha background at 1.6–2.8 MeV electron-equivalent (MeVee). The window above 3 MeV (up to 25 MeV) is nearly background-free once muons are vetoed. That window is the best place to test the 23.8 MeV and e⁺e⁻/bremsstrahlung channels.
9. **Transmutation evidence based on elements appearing with natural isotope ratios is contamination until proven otherwise.** SIMS and ICP-MS isotope "anomalies" are dominated by hydride and oxide interferences (e.g. ¹⁰⁴PdD⁺ at mass 106), matrix effects and ordinary chemical fractionation, such as Li isotopes during electrolysis.
10. **Integrated design, in order of credibility per dollar:**
    - closed-cell dual-method calorimetry;
    - sealed-headspace ⁴He sampling with getter purification and high-resolution MS;
    - distilled-sample LSC tritium with a full inventory;
    - ³He well counter with muon veto and a "dummy" (neutron-blind) tube;
    - blind CR-39 filter stacks;
    - a LaBr₃/CeBr₃ back-to-back pair or an HPGe detector.

    Run everything twice: identical **twin stations** (D₂O/Pd test cell vs. H₂O/Pd or D₂O/Pt control), a pre-registered analysis, blinded data, and at least one cell swap between stations.

---

## 1. The yardstick: what a real signal would look like

| Quantity (per watt unless noted) | Value | Note |
|---|---|---|
| d+d→⁴He (23.85 MeV) | 2.62×10¹¹ ⁴He/s; 2.26×10¹⁶ He per W·day = 0.84 µL(STP) | [C] |
| d+d→t+p (4.03 MeV) | 1.55×10¹² t/s; 1 W·day of tritium = 2.4×10⁸ Bq | [C] |
| d+d→³He+n (3.27 MeV) | 1.91×10¹² n/s if that branch alone | [C] |
| Conventional d+d (≈50/50 branches) | 1.71×10¹² fusions/s → 8.6×10¹¹ n/s; H*(10) ≈ 10 Sv/h at 1 m | [C], using ~400 pSv·cm² for 2.5 MeV neutrons [V: ICRP-74] |
| He in air | 5.24 ppm → 1.41×10¹⁴ He atoms/cm³ of air | [C] |
| Air volume equal to 1 W·day of ⁴He | 0.16 cm³ | [C] |
| Energy to reach 100 eV/Pd atom in 1 g Pd | 9.1×10⁴ J = 25 Wh (10.5 days at 100 mW) | [C]; loading enthalpy is only ~150 J [C] |

**Design consequence.** Any observed neutron or tritium rate should be reported as a *ratio to the simultaneously measured excess power*. A 50 mW claim with a 0.01 n/s upper limit means n/heat ≤ 2×10⁻¹³ of the conventional value [C].

---

## 2. Calorimetry

### 2.1 Methods compared

| Type | Principle | Realistic accuracy | Time constant | Position dependence | Cost |
|---|---|---|---|---|---|
| Isoperibolic Dewar (Fleischmann–Pons) | Cell in a partially silvered Dewar inside a ±0.01 K bath. Heat leaves mostly by radiation through the unsilvered section, k_eff ≈ 4k_R·T³ ≈ 0.09 W/K [C] | F&P claimed ~0.1% [V]. Independent analyses (Miles; GE/Wilson) put it at ≈ ±1% of input or ±10–20 mW [V] | τ = C/k ≈ 400 J/K ÷ 0.086 W/K ≈ **1.3 h** [C] | **High**: depends on electrolyte level, which changes k, and on stratification at low current | $3–10k |
| Isoperibolic, conduction-dominated (Miles-type) | Well-defined conductive thermal resistance | ±20 mW or ±1% [V: Miles 1994] | 0.5–2 h | Moderate | $3–10k |
| Mass-flow | P = ṁ·c_p·ΔT, with coolant jacket fully enclosing the cell | ±0.3–0.5% of input [V: SRI/McKubre] | 5–20 min | Low if jacket and insulation are complete | $10–25k |
| Seebeck envelope (heat-flux) | Thermopile / thermoelectric modules covering >90% of an enclosing box, held in a ±0.01 K bath | 0.1–0.5% of input; a few mW floor [V: Storms; Thermonetics-class commercial units] | 5–20 min | **Lowest** when coverage is complete | $5–15k DIY; $30–60k commercial [V] |
| Twin / differential | Two identical cells in one bath, one of them a control | ~0.1% of input on the difference [V] | same as base type | Cancels bath drift, not position effects | +50% |

**Flow-calorimeter numbers [C].** At ṁ = 0.5 g/s and P = 5 W, ΔT = 2.39 K, so a 1 mK sensor mismatch is 0.04%. At 1 W the same 1 mK is 0.21%. The rules: choose ṁ so that ΔT ≥ 1 K, measure flow gravimetrically to 0.1%, and use a matched thermistor pair calibrated to ≤1 mK.

### 2.2 Recombination and cell type
- **Open cell.** The calculation assumes all evolved D₂ and O₂ leave the cell, so P_net = I(V − E_tn), with E_tn = 1.527 V for D₂O. If a fraction *f* recombines inside the cell, the apparent excess is f·I·E_tn. At I = 0.5 A and f = 5% this gives **38 mW**; at 1 A and 10% it gives **153 mW** [C]. Jones, Hansen et al. (1995) measured Faradaic efficiencies below 100% in open cells at low current density and showed this could account for many reported excesses [V].
- **Closed cell** with an internal recombiner (Pt or Pd on alumina, or Pd/C): P_in = V·I exactly, with no Faradaic assumption. Add a **pressure transducer**: any pressure rise means incomplete recombination, which invalidates the energy balance. The recombiner releases the gas enthalpy at the recombiner, not at the cathode, so the calorimeter must be position-independent (see 2.3).

### 2.3 The major critiques and how the design answers each

| Critique | Content | Design answer |
|---|---|---|
| **Shanahan CCS / ATER** (Thermochim. Acta 2002, 2005) | A 1–3% shift of the calibration constant, caused by heat moving between the electrode region and the recombiner or headspace (at-the-electrode recombination), can explain hundreds of mW of apparent excess in Storms' Seebeck data [V]. Storms (2006) replied, and Marwan, Krivit et al. (2010) responded in J. Environ. Monit. | Map the calibration constant with a Joule heater at ≥5 positions (cathode centre, cathode end, recombiner, headspace, bottom), **with electrolysis running at the same total power**. Require spread < 0.2%, and report it as a systematic. |
| **Wilson et al. (GE), J. Electroanal. Chem. 332 (1992)** | Criticized how F&P determined the heat-transfer coefficient and handled enthalpy and evaporation terms. Concluded the claimed excess was substantially overstated; F&P rebutted in the same issue [V] | Use a calorimeter with a directly measured calibration, not a modelled heat-transfer coefficient. Account for vapor enthalpy only in closed cells. |
| **Miles et al. error analysis** (J. Phys. Chem. 98, 1948, 1994) | Error budget for isoperibolic cells: ~±20 mW, dominated by electrolyte level, bath stability and sensor calibration [V] | Hold the electrolyte level constant with a closed cell. Bath at ±0.01 K. Calibrate sensors in situ. |
| **Boil-off / spray** | Entrained liquid (foam, spray) counted as evaporated water inflates the vaporization enthalpy | Avoid boiling regimes. Condense and weigh everything that leaves. |
| **Input-power error** | ⟨V⟩⟨I⟩ ≠ ⟨VI⟩ with ripple, pulses or galvanostat oscillation | Digitize V and I at ≥100 kHz and compute ⟨VI⟩, or use a true power analyzer. Check the galvanostat for oscillation with a scope. |
| **Gas-loading / high-temperature artifacts** | Steam quality (wet steam counted as dry), infrared-emissivity errors in thermography (e.g., the 2014 Lugano E-Cat report), thermocouple pickup from AC heaters [V] | Do not use thermography-based calorimetry. Use enclosed flow or Seebeck calorimetry only. |
| **Heat after death / transients** | Heat-capacity and heat-transfer changes during dry-out or current steps | Integrate energy only across states with a valid calibration. Report integrated energy, not peak power. |

### 2.4 Recommended calorimeter (numbers)
- **Dual method in series.** Heat flows from the cell through a Seebeck envelope of thermoelectric modules on an Al or Cu can, into a water-cooled jacket that forms a mass-flow calorimeter. The **same heat is measured twice with independent systematics**. The two readings must agree within their errors. The water jacket (≈2 cm) also serves as the first neutron moderator layer (§11).
- Input 1–5 W. Target σ = max(5 mW, 0.3% P_in). Bath or reservoir at ±0.01 K. Outer enclosure at ±0.1 K.
- Calibration heater: a 4-wire DC Joule heater inside a dummy-cathode geometry, plus a heater at the recombiner. Calibrate before, during (as superimposed pulses) and after each run, over the full power range.
- **Salting.** A third party injects hidden heater pulses and hidden offsets into the data stream. The analysis is unblinded only after the pipeline has been frozen.
- Record cell pressure, electrolyte level and room temperature and humidity. Record D/Pd loading via the 4-wire resistance ratio R/R₀.

---

## 3. Neutron detection

### 3.1 Detector options

| Detector | Physics | Intrinsic efficiency for 2.45 MeV | Background / issues | Cost [V] |
|---|---|---|---|---|
| **³He proportional tubes in HDPE** | ³He(n,p)t, Q = 764 keV, σ_th = 5330 b. Tubes 2.54 cm × 30–100 cm at 4 atm | Thermal ~70–90% per tube. System efficiency is set by geometry and moderation (see 3.3) | Excellent gamma rejection with a threshold (~10⁻⁶–10⁻⁸). Microphonics and EMI; HV discharge in humid air | $1.5–3k per tube; $20–100k for a system |
| BF₃ tubes | ¹⁰B(n,α), Q = 2.79 MeV, 3840 b | ~1/2–1/3 of ³He per tube | Toxic gas. More gamma-tolerant. Same microphonic and EMI issues | $1–2k per tube |
| Liquid scintillator with pulse-shape discrimination (PSD): EJ-301/EJ-309; stilbene; EJ-276 plastic | n–p elastic recoil; PSD separates neutrons from gammas | 15–30% for a 5 cm cell at a PSD threshold of ~100–200 keVee (≈0.8–1.2 MeV proton recoil) [V] | **Gives energy information**: the recoil edge at 2.45 MeV distinguishes d–d neutrons from the cosmic spectrum. Higher background from PSD leakage and cosmics; ns timing | $5–10k per detector plus a $10k digitizer |
| Capture-gated spectrometer (organic scintillator + ⁶Li glass or ¹⁰B) | Recoil signal followed by delayed capture | ~1% (the BYU 1989 system) [V] | Good background rejection | custom |
| Bubble detectors (BTI BD-PND) | Superheated droplets; threshold ~100–200 keV | ~1 bubble per 1–2.5×10³ n/cm² [V] | 1 n/s at 5 cm gives only **~0.1 bubbles/day**, while the cosmic background gives ~0.2–0.4 bubbles/day [C]. Temperature-, shock- and age-sensitive. Useful only for ≥10² n/s | $0.3–0.5k each |
| Activation foils | Indium: ¹¹⁵In(n,n′)¹¹⁵ᵐIn (4.5 h, 336 keV) or ¹¹⁵In(n,γ)¹¹⁶ᵐIn after moderation (54 min, 162 b). Gold: ¹⁹⁷Au(n,γ) (2.7 d, 412 keV) | 1 g In at φ = 1 n/cm²/s: A_sat = 1.5×10⁻³ Bq (fast) or 0.8 Bq (thermalized) [C] | Forensic only: needs φ ≳ 10–100 n/cm²/s | <$1k plus gamma counting |

### 3.2 Backgrounds
- **Cosmic-ray neutrons, sea level, New York City reference (Gordon et al. 2004):** integral flux above 10 MeV = 3.6×10⁻³ cm⁻²s⁻¹. The total over all energies is ~1.3×10⁻² cm⁻²s⁻¹ [V]. The spectrum has a thermal peak, a 1/E region, an evaporation peak at 1–2 MeV that overlaps 2.45 MeV, and a cascade peak near 100 MeV.
- **Variation:**
  - barometric coefficient ≈ −0.7%/hPa [V], so typical weather swings of ±20 hPa change the rate by **±14%**;
  - solar-cycle modulation ~5–15% at mid-latitudes, and Forbush decreases of up to −10–20% within a day [V];
  - altitude: ×3–4 at 1.6 km and ×10 at 3 km [V];
  - overburden: the hadronic attenuation length is ~150 g/cm², so **2–3 m of concrete or soil (~500–700 g/cm²) reduces the rate ~×20–100** [C]. A basement is the cheapest background reduction available.
- **Muons:** ~1 cm⁻²min⁻¹ at sea level (barometric coefficient ≈ −0.15%/hPa) [V]. Muon capture and spallation in lead or steel near the detector produce **multi-neutron bursts**, which is exactly what "neutron burst" claims look like. Keep high-Z material out of the moderator volume, or veto it.
- **Accidental coincidences:** at B = 0.05 cps with a 100 µs gate, about 8 accidental doublets per year; at B = 1 cps, about 3×10³ per year [C].
- **Other sources:** people standing nearby (hydrogenous bodies reflect and moderate neutrons, changing the rate by a few %), water level changes in the adjacent flow calorimeter, neutron sources used elsewhere in the building, and radon (via gamma pile-up and PSD leakage in scintillators).

### 3.3 Geometry recommendations (tabletop)
- **Efficiency ladder** [V; confirm with Monte Carlo and a ²⁵²Cf source at the cathode position]:

  | Configuration | Typical efficiency |
  |---|---|
  | Single moderated tube at 10–20 cm | 0.1–1% |
  | Slab of 4–8 tubes facing the cell at 5 cm | 2–8% |
  | Single-ring 4π well counter (8–18 tubes) | 15–35% |
  | 2–3 rings (epithermal-multiplicity-counter class) | 50–65% |

- **Moderator:** 3–6 cm HDPE-equivalent between the cell and the tube centres, where water counts as ≈0.85× HDPE. Add ≥5 cm HDPE reflector behind the tubes, 1 mm Cd outside that, then 5–10 cm of 5% borated PE. Die-away time is ~30–60 µs, so use a coincidence gate of ~64 µs and a veto window of ~150–200 µs.
- **Muon veto:** 5 cm plastic-scintillator panels (e.g. EJ-200) on the top and four sides. At ~15 Hz per 30×30 cm panel the dead time is ≲0.5% [C].
- **Controls built into the counter:** split the tubes into two interleaved groups, A and B, with separate electronics; a real signal must appear in both. Add **one neutron-blind "dummy" tube** (⁴He-filled or unfilled, with the same HV and preamp) mounted on the same frame to record EMI and microphonic events. Digitize full waveforms and reject events with abnormal rise times.
- **Spectroscopic confirmation:** two EJ-309 or stilbene detectors in radial ports at ~15 cm. A d–d signal has to show a recoil edge at 2.45 MeV.

### 3.4 Minimum detectable emission (Currie L_D, paired blank) [C]

| Background B (cps) | 1 day, ε = 10% | 1 day, ε = 30% | 1 week, ε = 30% |
|---|---|---|---|
| 0.005 (underground / deep basement) | 0.012 n/s | 0.0038 n/s | 0.0014 n/s |
| 0.05 (basement + veto + borated PE) | 0.036 | 0.012 | 0.0045 |
| 0.5 (surface, shielded) | 0.11 | 0.037 | 0.014 |
| 2.0 (surface, unshielded) | 0.22 | 0.075 | 0.028 |

### 3.5 History of LENR neutron claims
- **Fleischmann & Pons (1989):** claimed ~4×10⁴ n/s based on a 2.2 MeV gamma peak. Petrasso et al. (MIT) showed the peak's width and energy, and the absent Compton edge, were inconsistent with capture gammas. The claim was withdrawn [V].
- **Jones et al., BYU (Nature 1989):** a capture-gated spectrometer with ~1% efficiency recorded a ~10⁻³ cps excess at 2.5 MeV, a marginal ~3σ result implying order 0.1–1 n/s. It was never reproduced at higher significance [V].
- **De Ninno/Scaramuzzi, ENEA Frascati (1989):** Ti chips in pressurized D₂ with thermal cycling; BF₃ counters recorded bursts. The effect was intermittent and not confirmed with better detectors [V]. Thermal cycling and pressure changes are exactly the conditions that cause microphonics.
- **Menlove et al., LANL/BYU (1990):** high-efficiency ³He counters recorded bursts from Ti/D₂. Cosmic spallation multiplets, which a veto or underground site reduces, are a known burst mimic. No reproducible signal was established [V].
- **Harwell, Caltech, MIT, Utah (Salamon) (1989–90):** null results with upper limits [V].
- **Berlinguette et al. (Google/UBC/MIT/LBNL, Nature 2019):** found "no evidence" for neutrons or 2.223 MeV gammas above background. *Pull their exact detection limits from the paper's supplementary information before citing them.*
- **Conventional but anomalous-looking sources of real neutrons:** fracto-emission ("fractofusion") in cracking TiDₓ or PdDₓ, and pyroelectric acceleration (Naranjo et al., Nature 2005). These produce genuine *hot* d–d neutrons. A real neutron signal correlated with cracking or thermal cycling is therefore not evidence of LENR.

---

## 4. Charged-particle detection

### 4.1 Ranges (µm) [C: ATIMA/pycatima]

| Particle | Mylar | CR-39 | Al | Pd | PdD₀.₈ | D₂O | Air (cm) |
|---|---|---|---|---|---|---|---|
| p 3.02 MeV (d+d) | 114 | 118 | 82 | 31 | 33 | 148 | 14.5 |
| t 1.01 MeV (d+d) | 11.9 | 12.2 | 10.0 | 4.6 | 4.7 | 15.7 | 1.6 |
| ³He 0.82 MeV (d+d) | 3.4 | 3.5 | 2.9 | 1.4 | 1.4 | 4.6 | 0.48 |
| α 3.7 MeV (d+³He) | 18.3 | 18.7 | 14.3 | 6.1 | 6.2 | 23.6 | 2.4 |
| α 5.49 MeV (²²²Rn) | 32.9 | 33.7 | 24.7 | 10.1 | 10.5 | 42.5 | 4.2 |
| α 7.69 MeV (²¹⁴Po) | 55.4 | 56.9 | 40.7 | 16.2 | 16.9 | 71.7 | 7.1 |
| p 14.7 MeV (d+³He) | 1834 | 1895 | 1240 | 412 | 438 | 2383 | 229 |

**Consequences:**
- Only protons from the outer ~30 µm of a Pd cathode can escape.
- Tritons escape from the outer ~5 µm, and ³He from the outer ~1.4 µm.
- **Nothing escapes through electrolyte thicker than 0.15 mm.** CR-39 immersed in electrolyte, or separated from the cathode by a liquid gap, sees only particles born within microns of its surface.

### 4.2 Filter stack: residual energy (MeV) behind Mylar [C]

| Mylar (µm) | 6 | 12 | 25 | 40 | 60 | 100 | 150 |
|---|---|---|---|---|---|---|---|
| p 3.02 | 2.93 | 2.83 | 2.61 | 2.33 | 1.93 | 0.84 | stopped |
| t 1.01 | 0.54 | stopped | – | – | – | – | – |
| ³He 0.82 | stopped | – | – | – | – | – | – |
| α 5.49 (Rn) | 4.80 | 4.03 | 1.89 | stopped | – | – | – |
| α 7.69 (²¹⁴Po) | 7.15 | 6.59 | 5.21 | 3.24 | stopped | – | – |
| α 8.78 (²¹²Po) | 8.29 | 7.78 | 6.58 | 4.98 | 2.01 | stopped | – |

Equivalents in Al: 25 µm stops ≤5.5 MeV α and leaves a 3 MeV proton at 2.42 MeV; 45 µm Al stops ²¹⁴Po alphas and leaves the proton at 1.84 MeV [C].

**Recommended stack on each chip:** bare / 6 / 25 / 60 / 150 µm Mylar. How each window is read:
- **60 µm window: the proton window.** It stops all natural alphas up to 8.78 MeV and leaves 3 MeV protons at ~1.9 MeV, near the track-diameter optimum.
- **150 µm window: the null.** It stops 3.02 MeV protons, so it should show only intrinsic background and neutron recoils.
- **6–12 µm step:** distinguishes tritons.
- **Bare window:** the only one that sees ³He, which has a 3.5 µm range.

### 4.3 CR-39 protocol, sensitivity, background
- **Etch:** 6–6.5 M NaOH at 70 °C for 6–7 h, the SPAWAR/Lipson-type protocol [V]. The bulk etch rate is ~1.2–1.5 µm/h, so ~8–10 µm is removed per face [V]. Measure it for every batch from chip thickness or the fission-fragment method.
- **Track diameters (indicative [V]; calibrate in-house):** α 3–8 MeV → 8–15 µm. Protons 1–3 MeV → 3–8 µm, shrinking with energy; a 3 MeV proton barely develops unless degraded or etched longer. A 1 MeV triton behaves roughly like a 0.33 MeV proton. **0.82 MeV ³He has a 3.5 µm range, less than the ~9 µm etched away**, so it survives only as shallow round pits. Detecting ³He needs a short 1–2 h etch or sequential etching.
- **Calibration:** ²⁴¹Am with Mylar degraders; an accelerator proton or deuteron beam; a DD neutron generator for recoils. The MIT ICF group's calibration method (Séguin et al. 2003) is the reference.
- **Background:** intrinsic ~10–60 tracks/cm² for fresh sheet [V]. Radon plating on bare surfaces adds about 2–3 tracks/cm² per kBq·h/m³ [V]; 50 Bq/m³ for one week gives **~17–25 tracks/cm²** [C]. Store and expose chips in N₂-purged, radon-tight containers. Include travel blanks and in-room blanks.
- **Artifacts:**
  - chemical attack by electrolyte (Cl⁻, Li⁺, PdCl₂ plating bath, oxidizers);
  - gas bubbles and mechanical contact from wires or pressure;
  - heating: bulk etch rate and track fading change above ~50–60 °C;
  - UV or gamma exposure altering etch response;
  - Pd deposits catalysing local etching;
  - scanner or observer bias.
- **SPAWAR triple tracks** (Mosier-Boss et al. 2009) were attributed to ¹²C(n,n′)3α, which requires E_n ≥ 9.6 MeV. Rough estimate [C, σ ≈ 0.25 b, sensitive depth 2×10 µm]: cosmic neutrons above 10 MeV (≈4×10³ n/cm² in 2 weeks) give ~0.08 triple tracks/cm². Ten per cm² needs ~6×10⁵ n/cm² of ≥10 MeV neutrons.
  - If those neutrons came from secondary d–t reactions (d–t/d–d ≈ 10⁻³), the primary d–d rate would be ~10³ n/s. A ³He well counter would see that trivially, and the recoil track density would be enormous.
  - **Cross-modality consistency is the test.** Kowalski and others showed that pits resembling tracks form through chemical and mechanical damage [V]. SPAWAR responded with tracks behind 6 µm Mylar and on the back side of chips, but 6 µm Mylar does not stop radon alphas or tritons (table 4.2).

### 4.4 Silicon detectors (surface-barrier / PIPS)
- **Advantages:**
  - resolution ~20 keV for α in vacuum;
  - background < 0.05 counts/h above 3 MeV for 450 mm² [V: vendor alpha-PIPS spec];
  - a minimum-ionizing muon deposits only ~0.1 MeV in 300 µm, so a threshold at ≥1 MeV rejects gammas and muons.
- **Geometry:** a disk detector of radius R = 12 mm at d = 10 mm covers Ω/4π = 0.18; at 5 mm, 0.31 [C].
- **Use:** only in vacuum or He gas. For electrolysis, make the cathode a **10–25 µm Pd foil window**: 3 MeV protons emerge at 2.4 MeV through 10 µm and 1.0 MeV through 25 µm [C]. D₂ permeating into the detector chamber must be pumped away.
- **Prior work:** Lipson et al. reported ~3 MeV protons and 11–16 MeV alphas from Pd/PdO:Dₓ under desorption or electrolysis. Kasagi, and the Rolfs/Raiola and Czerski groups, measured enhanced electron screening in d+d on metal targets at 5–20 keV beam energies [V]. That accelerator result is real and credible, but the reaction channels are the conventional ones.
- **Artifacts:** EMI from the electrolysis supply producing noise bursts (digitize waveforms), radon daughters plated on the detector (²¹⁰Po builds up over years), and light leaks.

---

## 5. Tritium

- **Measurement:** LSC of 1–8 mL of **distilled** water. Sample three places: the electrolyte, the recombiner and condensate water, and the cathode at end of run (vacuum extraction at 500–900 °C with oxidation, or dissolution).
- **Detection limits [C]:**
  - standard LSC (1 mL, 60 min, B ≈ 20 cpm, ε ≈ 30%): L_D ≈ 9 dpm ≈ **0.15 Bq/mL**;
  - low-level counter (8 mL, 500 min, B ≈ 1.5 cpm): ≈ **2 Bq/L**;
  - ³He ingrowth mass spectrometry is more sensitive still (1 Bq yields 3.2×10⁷ ³He atoms per year) but needs resolving power ≥510 to separate ³He from HD (table 6).
- **Electrolytic enrichment [C].** The T/D separation factor α ≈ 2 (range 1.5–2.5 [V]). In an open cell topped up with fresh D₂O at constant volume, x/x₀ = α − (α−1)·e^(−n/α), where n is the number of cell volumes electrolyzed. For α = 2: 1.39 after 1 volume, 1.78 after 3, 1.99 as n → ∞. Batch electrolysis without top-up follows x/x₀ = (V/V₀)^(1/α−1); reducing the volume 10× gives ×3.2. **Any tritium increase smaller than ~3× the starting level, without a full inventory balance, is not evidence.** A closed cell with a recombiner removes enrichment by gas loss, but tritium still partitions into the Pd.
- **Contamination:**
  - D₂O lots vary widely, typically 10¹–10³ dpm/mL [V]; reactor-derived D₂O can be orders of magnitude higher;
  - some as-received Pd lots have been reported to contain tritium [V: Texas A&M/Wolf];
  - laboratory tritium sources and luminous items;
  - the spiking allegations around the 1990 Texas A&M claims.
- **LSC artifacts:** chemiluminescence and photoluminescence from alkaline LiOD samples, and quench from Pd colloids. **Always distill and neutralize before counting**, dark-adapt the vials, and repeat each count.
- **Scale check:** 1 mW·day of the t+p branch in 50 mL gives 4.8×10³ Bq/mL [C]. Historical LENR tritium claims correspond to ≪10⁻⁶ of the heat and cannot be tied to it energetically.

---

## 6. Helium-4

| Mass pair | ΔM (u) | Resolving power M/ΔM needed at FWHM [C] |
|---|---|---|
| ⁴He (4.002603) vs D₂ (4.028204) | 0.02560 | 156 |
| ⁴He vs HT (4.023874) | 0.02127 | 188 |
| ⁴He vs H₂D (4.029752) | 0.02715 | 147 |
| ³He (3.016029) vs HD (3.021927) | 0.00590 | 511 |
| ³He vs H₃ (3.023475) | 0.00745 | 405 |

- **Practical note:** separating a trace ⁴He peak on a D₂ peak 10⁶× larger needs R ≈ 500+ together with good abundance sensitivity, or else removing the D₂ first.
- **Recommended analysis chain:**
  1. an all-metal sample cylinder (VCR or bellows valves, baked, evacuated);
  2. Zr-V-Fe or Ti getter to remove D₂, O₂ and other reactive gases (factor ≥10⁶);
  3. activated charcoal at 77 K to remove Ar;
  4. static-mode magnetic-sector noble-gas MS or a high-resolution quadrupole, with spike calibration (known ⁴He added to a blank).

  Instrument blanks are ~10⁹ atoms [V], compared with 1 mW·day = 2.3×10¹³ ⁴He [C]. **Sensitivity is not the problem; the air blank is.**
- **Leak budget [C]:**
  - keeping in-leakage below 1% of a 1 W signal requires < 2.6×10⁹ He/s, i.e. an air leak < 2×10⁻⁵ atm·cc/s. A helium-leak-checked metal system at ≤10⁻⁹ mbar·L/s meets this easily;
  - permeation through elastomer O-rings from atmospheric helium is ~10⁶–10⁷ atoms/s [V], also negligible;
  - **glass sample flasks are the exception.** Pyrex (K ≈ 10⁻¹¹ cm³·mm·s⁻¹·cm⁻²·cmHg⁻¹ [V: Altemose 1961]; 100 cm², 2 mm wall) takes in ~5×10⁶ He/s, or **~10 ppb-equivalent per month in a 50 mL flask** [C]. Early heat–helium samples stored in glass were criticized on exactly this point.
- **Dissolved air helium:** 100 mL of electrolyte in equilibrium with air holds 1.2×10¹⁴ He atoms, equal to 8 minutes of 1 W [C]. Purge the electrolyte with He-free Ar or D₂ first and discard the early samples.
- **The leak-immunity criterion (key geometric design).** Air in-leakage can only drive the headspace toward 5.24 ppm, never above it. With a 30 cm³ headspace, exceeding ambient needs 1.6×10⁴ J at full release, or 3.2×10⁴ J at 50% release (1.9–3.7 days at 100 mW) [C]. With a 500 cm³ headspace it would take 31–62 days. **Keep the headspace ≤20–30 cm³.** Also: **ban helium use in the laboratory during runs** (no He leak-testing, no He cylinders or balloons nearby) and log room-air helium.
- **Retention.** ⁴He produced in Pd is partly retained. SRI recovered the missing fraction by deloading, anodic stripping or heating, reporting ~104 ± 10% of the value expected at 23.85 MeV per ⁴He in experiment M4 [V]. Miles (China Lake) reported ~10¹¹ ⁴He per J, of order the theoretical 2.6×10¹¹/J, with helium present in heat-producing cells and absent in controls [V]. Critics note that those concentrations were **below ambient** (ppb in effluent gas), so leaks could not be excluded; the small-headspace design is meant to remove that objection. At end of run, **anneal or dissolve the cathode to release retained helium**.

---

## 7. Gamma and X-ray detection

| Detector | Resolution | Strengths | Pitfalls | Cost [V] |
|---|---|---|---|---|
| NaI(Tl) 3″×3″ | ~7% at 662 keV | Cheap, large volume | **Cannot separate 2204 keV (²¹⁴Bi) from 2223 keV (n-H capture)**; PMT gain drifts about −0.3%/°C | $2–4k |
| LaBr₃(Ce) 2″×2″ | ~2.7–3% | ~ns timing (good for 511–511 coincidence) | ¹³⁸La lines at 789 and 1436 keV; **²²⁷Ac alphas appear at 1.6–2.8 MeVee, over the 2.2 MeV region**. Prefer CeBr₃ or apply PSD | $10–15k |
| BGO | ~10% | High-Z; efficient for 511 keV and >10 MeV | Poor resolution | $5–10k |
| HPGe (coaxial, 40%) | ~0.15% (≈2 keV at 1332 keV) | Unambiguous line identification. **Also works as a neutron monitor** via the 2223 keV line, ¹⁰B(n,αγ) at 478 keV, ¹¹³Cd(n,γ) at 558 keV, and the asymmetric ⁷²Ge(n,n′) 691 keV and ⁷⁴Ge 596 keV peaks from fast neutrons | Needs LN₂ or cryocooler, lead shield, radon purge | $50–120k |

**Lines to know:**
- ²¹⁴Pb: 295, 352 keV. ²¹⁴Bi: 609, 1120, 1764, 2204 keV. These vary with radon.
- ⁴⁰K 1461 keV; ²⁰⁸Tl 2614 keV (the natural-background endpoint).
- 511 keV annihilation: always present, from cosmic pair production in lead and from β⁺ emitters.
- d+d→⁴He+γ at 23.8 MeV (branching ~10⁻⁷, never seen in LENR); p+d→³He+γ at 5.5 MeV (relevant because D₂O contains H); 2223 keV n-H capture.

**Czerski e⁺e⁻ channel.** Czerski et al. propose a 0⁺ threshold resonance in ⁴He that decays by internal pair creation, with ≈22.8 MeV shared between e⁺ and e⁻ (EPL 2016 and later work). The signature is **back-to-back 511–511 keV coincidences** (two LaBr₃/CeBr₃ or BGO detectors at 180°, ±5 ns window) plus high-energy e± and their **bremsstrahlung continuum at 3–20 MeV**. Above 2.614 MeV, natural background is almost entirely cosmic muons, which a veto removes, so that window is quiet. Background 511–511 pairs come from β⁺ activation and pair production in nearby lead: keep lead outside the veto and far from the coincidence axis.

**Other artifacts:**
- tritium bremsstrahlung below 18.6 keV, which fakes "soft X-rays" (Pd K X-rays cannot be excited, since the K-edge is 24.35 keV);
- high-voltage discharge devices produce genuine keV X-rays;
- in autoradiography, hydrogen and H₂O₂ fog emulsions chemically (the Russell effect).

---

## 8. Transmutation analysis

| Artifact | Example |
|---|---|
| Hydride and deuteride molecular ions | ¹⁰⁴PdD⁺ and ¹⁰⁵PdH⁺ at mass 106; ¹⁰⁶PdD⁺ at 108 (mimics ¹⁰⁸Pd enrichment or Ag); ⁷LiH⁺ at 8 |
| Oxides, argides, dimers, doubly charged ions (SIMS and ICP-MS) | ⁴⁰Ar¹⁶O / ⁴⁰Ca¹⁶O at 56 (Fe); ⁴⁰Ar²³Na at 63 (Cu); ¹²C¹⁶O / ¹⁴N₂ / ²⁸Si at 28; ⁵⁶Fe²⁺ at m/z 28 |
| Matrix effects | SIMS ion yields vary by 10³–10⁴ between elements and with surface state, so "concentrations" without standards are meaningless |
| Contamination from cell parts | Pt from anode dissolution plating onto the cathode (the classic case). Si, B, Na, Ca, Al leached from borosilicate glass. Fe, Cr, Ni from stainless steel. Cu, Zn, Pb from solder and fittings. Au from contacts. K and Na from skin. Laboratory dust |
| Natural chemical fractionation | Electrochemical Li-isotope fractionation (⁶Li is preferentially incorporated), hence a ⁶Li/⁷Li shift without any nuclear process; sputtering mass bias; detector dead time |
| Surface segregation | Impurities in the bulk Pd diffuse to the surface during loading and appear "new" in surface-sensitive analyses (XPS, SIMS, EDX) |

**Rules:**
- Analyze every input material (cathode, anode, electrolyte, glass) by ICP-MS beforehand.
- Run blank cells with an H₂O or Pt swap.
- Measure the same spot before and after (fiducial marks).
- Use ≥2 independent techniques, e.g. high-mass-resolution SIMS plus ICP-MS with a collision cell, or neutron activation analysis.
- Spike a control with the suspected contaminant.
- **Elements that appear with natural isotope ratios count as contamination.** Only isotopically anomalous products that exceed interference limits, seen with multiple techniques, deserve attention (cf. the Iwamura Cs→Pr permeation experiments, which are still debated).

---

## 9. Master artifact catalog

| # | Artifact | Affects | Mechanism | Historical example | Mitigation | Design implication for geometry |
|---|---|---|---|---|---|---|
| 1 | In-cell recombination | Open-cell calorimetry | Faradaic efficiency <1 inflates computed excess by f·I·E_tn (38 mW at 0.5 A, 5%) | Jones/Hansen 1995 vs open-cell claims | Closed cell with recombiner and pressure monitoring | Recombiner inside the calorimetric boundary; small headspace |
| 2 | Calibration constant shift | Any calorimeter | Heat source moves (electrode ↔ recombiner ↔ gas) | Shanahan vs Storms 2002–2010 | Map heater positions under operating conditions; <0.2% spread | Seebeck envelope or full flow jacket; no single-point sensors |
| 3 | Electrolyte level / stratification | Isoperibolic | Heat-transfer coefficient depends on level; poor mixing | F&P Dewar corrections; GE analysis | Closed cell, stirring, constant level | Avoid designs whose heat path changes with level |
| 4 | Input-power mismeasurement | All calorimetry | Ripple, pulses, oscillation | Pulsed and plasma electrolysis claims | ⟨VI⟩ at ≥100 kHz; power analyzer | Short, shielded leads; 4-wire sensing at the cell |
| 5 | Steam quality / emissivity | High-temperature gas-loading calorimetry | Wet steam; wrong infrared emissivity | E-Cat demos 2011; Lugano 2014 | Enclosed flow calorimetry | No open-air or thermography calorimetry |
| 6 | Chemical energy storage | Heat | Stored D₂/O₂, PdD enthalpy, Li alloying | "Heat after death" debates | 100 eV/atom criterion; integrate over full cycles | Know the cathode mass; minimize gas inventory |
| 7 | Barometric and solar modulation | Neutron background | −0.7%/hPa; ±14% with weather | Many marginal neutron excesses | Barometer covariate; twin station | Identical twin station in the same room |
| 8 | Cosmic spallation multiplets | Neutron "bursts" | Muons and hadrons in Pb or steel produce multiple neutrons | Burst claims 1989–91 | Muon veto; no high-Z near counter; underground | Veto panels enclose everything; basement location |
| 9 | Microphonics / EMI / HV discharge | ³He, BF₃, Si counts | Vibration, switching supplies, humidity-driven discharge | ENEA and BYU-era bursts (suspected) | Waveform digitization; dummy tube; A/B groups | Mechanically isolate counters from pumps; separate power |
| 10 | Nearby hydrogenous mass | Neutron efficiency and background | People and water albedo | Common in low-level counting | Access log; fixed water volumes | Rigid enclosure; fixed flow-jacket water |
| 11 | Gamma mis-assignment | Gamma "neutron" evidence | 2204 keV ²¹⁴Bi vs 2223 keV; bad calibration | F&P 1989 gamma spectrum (Petrasso) | HPGe; calibration sources; radon monitor | N₂ purge around the gamma detector |
| 12 | Internal activity of scintillators | Gamma | ¹³⁸La, ²²⁷Ac in LaBr₃ | — | CeBr₃; measure without cell | Put the gamma detector in its own shield |
| 13 | 511 keV background | e⁺e⁻ searches | Pair production in lead, β⁺ | Every HPGe spectrum | Back-to-back coincidence; veto | Pb outside veto; opposed detector pair |
| 14 | Radon plate-out | CR-39, Si | α 5.5–7.7 MeV | Common CR-39 background | N₂ purge; filter stack; blanks | Sealed detector chambers |
| 15 | Chemical and mechanical pits | CR-39 | Electrolyte attack, bubbles, contact | Kowalski vs SPAWAR | No electrolyte contact; controls; blind count | Detector behind a window, not immersed |
| 16 | Filter too thin | CR-39 particle ID | 6 µm Mylar passes t and Rn alphas | SPAWAR interpretation debate | 0/6/25/60/150 µm stack | Stack built into every chip holder |
| 17 | Cosmic triple tracks and recoils | CR-39 "neutrons" | ¹²C(n,n′)3α from >10 MeV cosmic neutrons | — | Blanks exposed for the same time | Blank chips in the same shield |
| 18 | Electrolytic T enrichment | Tritium | α ≈ 2 separation factor | Early tritium claims | Inventory model; closed cell | Recombiner water sampling port |
| 19 | Tritium contamination | Tritium | D₂O lots, Pd lots, spiking | Texas A&M 1989–90 | Assay every lot; chain of custody | — |
| 20 | LSC chemiluminescence | Tritium | Alkaline samples, colloids | Early electrolyte counts | Distill and neutralize | — |
| 21 | Air in-leak / dissolved helium | ⁴He | 5.24 ppm air; 1.2×10¹⁴ He per 100 mL | Heat–He critiques (below-ambient data) | All-metal system; above-ambient criterion; purge | Headspace ≤30 cm³, metal seals |
| 22 | Glass permeation | ⁴He samples | ~10 ppb per month in glass flasks | Early heat–helium flasks | Metal sample cylinders | — |
| 23 | D₂ / HT at mass 4 | ⁴He MS | Unresolved peaks | Low-resolution QMS reports | Getter; R ≥ 500 | — |
| 24 | Laboratory helium | ⁴He | Leak testing, He cylinders | — | Ban helium; monitor room air | — |
| 25 | Retained helium | He/heat ratio | He trapped in Pd | SRI M4 stripping | Anneal or dissolve cathode | Cathode removable for extraction |
| 26 | Mass interferences | Transmutation | Hydrides, oxides, argides, 2⁺ ions | Many SIMS claims | High mass resolution; two methods | — |
| 27 | Material contamination | Transmutation | Anode, glass, steel, dust | Pt, Fe, Cu, Zn "products" | Pre-assay; blanks; spikes | Minimize dissimilar materials in the electrolyte |
| 28 | Li-isotope fractionation | Transmutation | Electrochemistry | ⁶Li/⁷Li "shifts" | Control cell; fractionation model | — |
| 29 | Real hot fusion | Neutrons, p, t | Fracto-emission, pyroelectric fields | Naranjo 2005 | Correlate with cracking and fields | Avoid dielectric crystals and high fields near the D-loaded metal |
| 30 | Look-elsewhere / optional stopping | All | Searching many bins or windows | Burst and heat-spike claims | Pre-registration; global p-values | — |
| 31 | Temperature drift | PMT, NaI, PSD, thresholds | Gain drifts ~0.3–0.5%/°C | — | ±0.5 °C detector enclosure; gain stabilization | Thermal isolation between cell and detectors |
| 32 | Observer bias | CR-39, heat analysis | Unblinded choices | Many | Blind coding; salting | — |

---

## 10. Statistics and experimental design

### 10.1 Formulas
- **Currie** (α = β = 0.05), with gross counts compared to a paired blank of equal time:
  - critical level L_C = 2.33√B, detection limit L_D = 2.71 + 4.65√B (B = expected blank counts);
  - for a well-known blank: L_C = 1.645√B, L_D = 2.71 + 3.29√B;
  - minimum detectable source rate: MDA = L_D / (ε·t).
- **Discovery significance** (Asimov approximation, Cowan et al. 2011), for known background: Z = √{2[(s+b)ln(1+s/b) − s]}. With an off-measurement of time ratio τ, use the Li–Ma likelihood ratio.
- **Background systematics.** Z = S / √[(S+B)/t + (δB)²], so Z_max = S/(δB). **Discovery requires S > 5δB no matter how long you count.**
- **Feldman–Cousins 90% CL intervals** [C, own construction; the published FC table gives 1.08 for n=0, b=3]:

  | Observed n | b = 0 | b = 1 | b = 3 |
  |---|---|---|---|
  | 0 | [0, 2.44] | [0, 1.60] | [0, 0.95–1.08] |
  | 1 | [0.11, 4.36] | [0, 3.35] | [0, 1.88] |
  | 3 | [1.10, 7.42] | [0.10, 6.42] | [0, 4.42] |
  | 5 | [1.84, 9.98] | [1.25, 8.98] | [0, 6.99] |

### 10.2 Worked examples [C]
1. **Source 1 n/s, ε = 10%, B = 0.05 cps (S = 0.1 cps).**

   | Assumption | Time to 5σ |
   |---|---|
   | Gaussian, known background (S·t/√(B·t) = 5) | 125 s |
   | Poisson Asimov, known background | **193 s** |
   | Signal variance included | 375 s |
   | Equal-time on/off, Li–Ma | 478 s on + 478 s off |

   Easy. The 5σ systematic floor is δ < 40%.
2. **Source 0.01 n/s (typical of claimed LENR levels):**

   | Setup | Time to 5σ (known B) | Systematic check |
   |---|---|---|
   | Single tube, ε = 1%, B = 0.01 cps | **290 days** | S/B = 1% ≈ the barometric systematic → **not achievable** |
   | Well counter, ε = 30%, B = 0.5 cps (surface) | 16 days | S/B = 0.6%, below realistic δ → **not achievable** |
   | Well counter, ε = 30%, B = 0.05 cps (basement, veto) | 1.6 days | S/B = 6%; needs δ ≲ 1% (twin station) |
   | Well counter, ε = 30%, B = 0.005 cps (underground) | 0.2 days | robust |
3. **Burst search.** 30 days of 1-minute bins = 43,200 trials. A local 3σ excess is expected **58 times** and a local 4σ excess ~1.4 times, by chance alone. With a mean of 3 counts per bin, a 12-count minute looks like 5.2σ in the Gaussian approximation but is really 3.8σ local (Poisson), and has global p ≈ 0.95. Requiring 5σ global means ≥16 counts in a single minute, which is 5.2σ local.
4. **Heat.** With P_in = 2 W and σ = max(5 mW, 0.3%) = 6 mW, 50 mW is an 8σ excess *per calibrated state*. But the credible σ is the day-to-day reproducibility of the calibration, including position mapping. Require that and the 100 eV/atom criterion (25 Wh for a 1 g cathode).
5. **Helium.** 50 mW for 3 days at 50% release gives 1.7×10¹⁵ ⁴He, i.e. 2.1 ppm in a 30 cm³ headspace [C]. A pre-purged blank is ≪0.05 ppm, so this is highly significant statistically. The leak-immune (above-ambient) threshold needs ~3.7 days at 100 mW.

### 10.3 Design protocol
- **Pre-register:**
  - primary endpoint: integrated excess energy > 5σ_sys *and* > 100 eV/Pd atom;
  - secondary endpoints: ⁴He/energy ratio vs 2.6×10¹¹ /J; upper limits on n/W, t/W and e⁺e⁻/W;
  - the analysis windows, cuts and stopping rule.
- **Blind:**
  - third-party salting of calorimetry data;
  - coded cells (disguise D₂O vs H₂O by hiding cell masses);
  - coded CR-39, LSC and helium samples.
- **Controls:**
  - H₂O/Pd and D₂O/Pt in an **identical twin station** running at the same time;
  - swap cells between stations mid-run (crossover design);
  - dummy calorimeter runs with a resistor only.
- **Modulation:** use controlled current on/off and loading/deloading steps. Test the lagged cross-correlation against pre-registered lags only.

---

## 11. Recommended integrated measurement suite for a single cell

### 11.1 Geometry (indicative radii from the cathode axis; check with an MCNP/Geant4/OpenMC model and ²⁵²Cf calibration)

```
 r (cm)  layer
 0–3.5   Sealed metal cell (Inconel/SS with PTFE liner; all-metal gas boundary), ~100 mL electrolyte,
         HEADSPACE ≤30 cm³, internal Pt/Al2O3 recombiner, pressure transducer, Joule heaters at cathode
         and recombiner, R/R0 loading leads. Bottom port: 10–25 µm Pd-foil cathode window ->
         small vacuum/N2 chamber with Si PIPS (d≈10 mm) OR CR-39 behind 0/6/25/60/150 µm Mylar stack.
 3.5–4.5 Seebeck envelope (thermoelectric modules on Cu can)             } dual-method calorimeter:
 4.5–6.5 Water flow jacket, 0.5–1 g/s, ±0.01 K reservoir                  } Seebeck + mass-flow in series
 6.5–7.5 Low-density insulation (foam/aerogel)
 7.5–9.5 HDPE (≈2 cm) → total ≈4 cm HDPE-equivalent to tube centres
 ~10.8   Ring of 12 ³He tubes (2.54 cm, 4 atm, ~50 cm active), groups A/B interleaved + 1 dummy tube
 12–17   HDPE reflector; axial HDPE plugs top/bottom (penetrated only by feedthroughs)
 17      1 mm Cd
 17–25   5% borated PE
 outside 5 cm plastic-scintillator muon veto on top + 4 sides; N2-purged, radon-tight enclosure
 ports   4 radial Ø6 cm ports at 90°: 2 × LaBr3/CeBr3 (or BGO) at 180° for 511–511 + 3–25 MeV window;
         2 × EJ-309/stilbene at ~15 cm for the 2.45 MeV recoil edge. (Port losses ≈10–20% of ε.)
 site    basement ≥2 m concrete/soil overburden; no He use in the lab; barometer, radon monitor,
         T/RH, accelerometer, mains-EMI probe logged with the same timestamps.
 twin    Identical Station B with the control cell; swap cells mid-run.
```

Expected performance [C/V]:
- neutron ε ≈ 15–30%, B ≈ 0.05–0.2 cps (basement, vetoed), 1-day MDA ≈ 0.01–0.03 n/s;
- calorimetry σ ≈ 5–15 mW at 1–5 W input;
- ⁴He above-ambient threshold after ≈2–4 days at 100 mW;
- tritium MDA ≈ 0.15 Bq/mL;
- CR-39 proton window sensitive to a few tracks/cm² above blank.

An HPGe detector can replace one LaBr₃ once the station is mature, for unambiguous 2204/2223 keV and Ge(n,n′) identification.

### 11.2 Ranking by credibility per dollar

| Rank | Measurement | Cost [V] | What it proves | Credibility per $ |
|---|---|---|---|---|
| 1 | ⁴He from a sealed ≤30 cm³ headspace; getter + high-resolution MS; cathode extraction at end | $10–30k (plus ~$0.5–2k per outsourced sample) | Quantitative product/energy ratio; leak-immune above 5.24 ppm | **Highest, if heat is present** |
| 2 | Closed-cell dual-method calorimetry with position mapping and twin control | $15–40k | Necessary gate. Not nuclear-specific by itself | High (indispensable) |
| 3 | Tritium by LSC with distillation and full inventory | $2–5k | Extreme sensitivity to the t+p branch; branching-ratio limits | High for limits; moderate for positives |
| 4 | ³He well counter + veto + dummy tube + A/B groups | $40–100k | Cleanest yes/no nuclear signature; n/W limit | Moderate–high |
| 5 | CR-39 differential filter stack, blind-counted | $3–10k | Screening; particle-ID hints | Moderate (cheap, weak evidence) |
| 6 | LaBr₃/CeBr₃ back-to-back pair (511 coincidence, 3–25 MeV) | $25–40k | Tests the Czerski e⁺e⁻ channel, 23.8 MeV and bremsstrahlung | Moderate |
| 7 | Si detectors behind a Pd-foil window | $10–20k | Charged-particle spectroscopy | Moderate (thin-film cells only) |
| 8 | HPGe | $60–120k | Line identification; Ge(n,n′) as a fast-neutron check | Lower per $ |
| 9 | Transmutation (ICP-MS, SIMS) | $5–30k | Rarely conclusive | Lowest |

**Bottom line.** No LENR report to date has met the combined bar of: position-independent closed-cell heat above 100 eV/atom; ⁴He above ambient in a sealed metal system, in the 23.85 MeV ratio; null controls in a twin station; and blind analysis. The suite above is designed so that a positive result on that combination could not plausibly be explained by any artifact in §9. Equally, a null result would give branching-ratio limits that are informative in their own right.

---

## 12. References
(DOI links are given where the identifier is certain. Otherwise a Scholar search link is given. [V] values should be checked against these sources.)

1. Fleischmann, Pons, Hawkins, J. Electroanal. Chem. 261, 301 (1989). https://doi.org/10.1016/0022-0728(89)80006-3
2. Jones et al., "Observation of cold nuclear fusion in condensed matter," Nature 338, 737 (1989). https://doi.org/10.1038/338737a0
3. Petrasso et al., "Problems with the γ-ray spectrum in the Fleischmann et al. experiments," Nature 339, 183 (1989). https://doi.org/10.1038/339183a0
4. Salamon et al., "Limits on the emission of neutrons, γ-rays, electrons and protons from Pons/Fleischmann electrolytic cells," Nature 344, 401 (1990). https://doi.org/10.1038/344401a0
5. Williams et al. (Harwell), Nature 342, 375 (1989). https://doi.org/10.1038/342375a0 ; Lewis et al. (Caltech), Nature 340, 525 (1989). https://doi.org/10.1038/340525a0
6. Berlinguette et al., "Revisiting the cold case of cold fusion," Nature 570, 45 (2019). https://doi.org/10.1038/s41586-019-1256-6 (accessible summary: https://www.researchgate.net/publication/333405089_Revisiting_the_cold_case_of_Cold_Fusion)
7. De Ninno et al., Europhys. Lett. 9, 221 (1989). https://scholar.google.com/scholar?q=De+Ninno+1989+Europhysics+Letters+neutron+titanium+deuterium
8. Menlove et al., J. Fusion Energy 9, 495 (1990). https://scholar.google.com/scholar?q=Menlove+1990+neutron+emission+Ti+Pd+pressurized+D2
9. Wilson et al., "Analysis of experiments on the calorimetry of LiOD–D2O electrochemical cells," J. Electroanal. Chem. 332, 1 (1992). https://scholar.google.com/scholar?q=Wilson+1992+calorimetry+LiOD-D2O+electrochemical+cells
10. Fleischmann & Pons, Phys. Lett. A 176, 118 (1993). https://scholar.google.com/scholar?q=Fleischmann+Pons+1993+calorimetry+Pd-D2O+simplicity+complications
11. Miles, Bush, Stilwell, "Calorimetric principles and problems in measurements of excess power during Pd-D2O electrolysis," J. Phys. Chem. 98, 1948 (1994). https://scholar.google.com/scholar?q=Miles+Bush+Stilwell+1994+calorimetric+principles+problems
12. Miles et al., "Correlation of excess power and helium production during D2O and H2O electrolysis using palladium cathodes," J. Electroanal. Chem. 346, 99 (1993). https://scholar.google.com/scholar?q=Miles+1993+correlation+excess+power+helium+production
13. Jones, Hansen et al., "Faradaic efficiencies less than 100% during electrolysis of water can account for reports of excess heat in 'cold fusion' cells," J. Phys. Chem. 99, 6973 (1995). https://scholar.google.com/scholar?q=Faradaic+efficiencies+less+than+100%25+cold+fusion+cells
14. Shanahan, Thermochim. Acta 382, 95 (2002); 428, 207 (2005). https://scholar.google.com/scholar?q=Shanahan+Thermochimica+Acta+calibration+constant+shift
15. Storms, "Comment on papers by K. Shanahan…," Thermochim. Acta (2006). https://www.osti.gov/etdeweb/biblio/20829521 ; Shanahan reply: https://www.osti.gov/etdeweb/biblio/20829522
16. Marwan, Krivit et al., "A new look at LENR research: a response to Shanahan," J. Environ. Monit. 12, 1765 (2010). https://www.lenr-canr.org/acrobat/MarwanJanewlookat.pdf
17. McKubre et al., SRI heat–helium (M4) and EPRI TR-107843 (1998). https://scholar.google.com/scholar?q=McKubre+helium+M4+heat+SRI
18. Hagelstein et al., "New physics effects in metal deuterides" (DOE review, 2004). https://scholar.google.com/scholar?q=Hagelstein+McKubre+Nagel+Chubb+Hekman+new+physics+effects+metal+deuterides
19. Mosier-Boss et al., "Triple tracks in CR-39 as the result of Pd–D co-deposition," Naturwissenschaften 96, 135 (2009). https://link.springer.com/article/10.1007/s00114-008-0449-x ; DT-neutron comparison: https://www.osti.gov/etdeweb/biblio/21416357
20. Kowalski, "Comments on 'Use of CR-39 in Pd/D co-deposition experiments'," Eur. Phys. J. Appl. Phys. 44, 287 (2008). https://scholar.google.com/scholar?q=Kowalski+comments+use+of+CR-39+Pd%2FD+co-deposition
21. Séguin et al., "Spectrometry of charged particles from ICF plasmas" (CR-39 calibration), Rev. Sci. Instrum. 74, 975 (2003). https://doi.org/10.1063/1.1518141
22. Lipson et al., CR-39/Si charged-particle emission from Pd/PdO:Dx. https://scholar.google.com/scholar?q=Lipson+charged+particle+emission+Pd+PdO+deuterium+CR-39
23. Kasagi et al., J. Phys. Soc. Jpn. 71, 2881 (2002); Raiola et al., Eur. Phys. J. A 19, 283 (2004). https://scholar.google.com/scholar?q=Kasagi+2002+energetic+protons+alpha+150+keV+deuteron+PdDx
24. Czerski et al., "Screening and resonance enhancements of the 2H(d,p)3H reaction yield in metallic environments," EPL 113, 22001 (2016). https://doi.org/10.1209/0295-5075/113/22001 ; later e⁺e⁻ work: https://scholar.google.com/scholar?q=Czerski+electron+positron+pair+deuteron+fusion+threshold+resonance
25. Naranjo, Gimzewski, Putterman, "Observation of nuclear fusion driven by a pyroelectric crystal," Nature 434, 1115 (2005). https://doi.org/10.1038/nature03575
26. Iwamura et al., Jpn. J. Appl. Phys. 41, 4642 (2002). https://doi.org/10.1143/JJAP.41.4642
27. Gordon et al., "Measurement of the flux and energy spectrum of cosmic-ray induced neutrons on the ground," IEEE Trans. Nucl. Sci. 51, 3427 (2004). https://www.researchgate.net/publication/3139171 ; JEDEC JESD89A; PNNL-30228 cosmic flux tool: https://www.pnnl.gov/main/publications/external/technical_reports/PNNL-30228.pdf
28. Currie, "Limits for qualitative detection and quantitative determination," Anal. Chem. 40, 586 (1968). https://doi.org/10.1021/ac60259a007
29. Feldman & Cousins, Phys. Rev. D 57, 3873 (1998). https://doi.org/10.1103/PhysRevD.57.3873
30. Cowan, Cranmer, Gross, Vitells, Eur. Phys. J. C 71, 1554 (2011). https://doi.org/10.1140/epjc/s10052-011-1554-0 ; Li & Ma, ApJ 272, 317 (1983). https://doi.org/10.1086/161295
31. Klein & Roodman, "Blind analysis in nuclear and particle physics," Annu. Rev. Nucl. Part. Sci. 55, 141 (2005). https://doi.org/10.1146/annurev.nucl.55.090704.151521
32. Altemose, "Helium diffusion through glass," J. Appl. Phys. 32, 1309 (1961). https://doi.org/10.1063/1.1736226
33. Ensslin et al., *Passive Nondestructive Assay of Nuclear Materials* (LANL "PANDA", NUREG/CR-5550) — ³He counter design, die-away, multiplicity. https://scholar.google.com/scholar?q=Passive+Nondestructive+Assay+of+Nuclear+Materials+NUREG%2FCR-5550
34. Knoll, *Radiation Detection and Measurement*, 4th ed. (Wiley, 2010).
35. ATIMA/catima stopping-power code (used for all [C] ranges): https://github.com/hrosiak/catima ; NIST PSTAR/ASTAR for cross-checks: https://physics.nist.gov/PhysRefData/Star/Text/intro.html
36. Eljen EJ-301/EJ-309 PSD scintillators: https://eljentechnology.com/products/liquid-scintillators/ej-301-ej-309
37. ICRP Publication 74 (1996), fluence-to-dose conversion coefficients.
38. Storms, *The Science of Low Energy Nuclear Reaction* (World Scientific, 2007); Taubes, *Bad Science* (1993) — history of the tritium and contamination controversies.
