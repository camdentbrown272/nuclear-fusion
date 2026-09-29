# R6 — Safety, Materials, Sourcing, Costs, and Build Constraints

*Research date: 2026-09-29. Scope: the practical and safety envelope that any LENR (Pd/D electrolysis, co-deposition, or gas-loading) geometry must fit inside.*

> **Source-quality note.** Web page fetching was blocked in this research session, and most vendor catalogs (Fisher, Sigma, Goodfellow) put prices behind a login. Items marked **(sourced)** come from a cited page or search-indexed snippet. Items marked **(est.)** are engineering estimates from metal spot value plus typical fabrication and catalog markups. Get quotes before committing money. The hazard calculations are first-principles and are shown in the text so they can be checked.

---

## 1. Executive summary

- **The main hazard is chemical: D₂/O₂ in the headspace.** It is not radiation. Every documented LENR accident was a gas-phase explosion or a stored-energy release: SRI 1992 (one death), the Fleischmann–Pons cube in 1985, Biberian 2004, and Mizuno 2005. One litre of stoichiometric 2D₂+O₂ at 1 bar holds about **6.7 kJ, roughly 1.6 g TNT-equivalent**.
- **Keep the headspace small and never let it become stoichiometric.** In open cells, sweep with an inert gas at **≥0.8 L/min per amp** (holds D₂ ≤1%, i.e. 25% of the LFL). Closed cells should be D₂-prefilled, fitted with a recombiner, and **tripped when ΔP ≥ 10% of the absolute fill pressure**, which keeps O₂ ≤3%, below the ~5% limiting oxygen concentration.
- **Closed-cell pressure rating:** if a detonable mixture cannot be ruled out, the design must survive about **20× the initial absolute pressure** (Chapman–Jouguet), or about 50× where reflected detonation is possible. Otherwise prevent the mixture by design and fit a burst disk.
- **Oxygen-clean everything.** The SRI forensic analysis (LLNL) found hydrocarbon machining oil inside the metal cell, in an O₂-enriched headspace. Use no organics except fluoropolymers.
- **Limit the stored-energy inventory.** Pd loaded to D/Pd ≈ 0.9 holds about **12–15 kJ of combustible D per cm³**. Keep cathodes **≤0.3 cm³**, avoid sealed hollow cathodes, and deload before opening a cell.
- **Keep the electrical system SELV:** ≤60 V DC at maximum current. This limits the electrode gap to roughly ≤5 mm for a plate at 1 A/cm² in 0.1 M LiOD.
- **Palladium is cheap; detectors and calorimetry are not.** Pd spot is **$1,198–1,265/ozt ($38.5–40.7/g)** (Sept 2026), so a 1 mm × 10 cm wire holds about $36 of metal. D₂O costs about **$0.3–0.65/g** in bulk and $1–3/g at retail (est.).
- **Add a neutron monitor with automatic shutdown.** A hypothetical 1 W of textbook D–D fusion would give about 9 Sv/h at 1 m. It has never been observed, but the alarm is cheap and the monitor doubles as the signal detector.
- **Regulatory load is light** with no licensed sources, no accelerator or X-ray generator, and tritium below the 1,000 µCi exempt quantity. D₂O needs no US purchase license; exports fall under 10 CFR 110.
- **Recommended tier: about $50k.** It buys Seebeck calorimetry (±20–60 mW), a closed-cell safety system, He-3, NaI and CR-39 detection, a vacuum anneal furnace, and outsourced SEM. $5k cannot produce credible calorimetry. Most of the $500k tier goes to staff, parallel cells, and ⁴He mass spectrometry.

---

## 2. Incident record and lessons

| Incident | What happened | Cause | Lesson for the geometry |
|---|---|---|---|
| **SRI International, 2 Jan 1992** (Andrew Riley killed, three injured) | Riley lifted the metal cell out of its water bath and was waiting for it to drain when it exploded. | D₂ + O₂ had accumulated **despite an internal recombiner**. LLNL forensics found **hydrocarbon (machining) oil** inside the cell and no explosives. Oil plus an O₂-enriched headspace can react explosively. | Oxygen-clean all parts. Trend pressure to catch recombiner failure. **Never handle, tilt or open a cell until it is deloaded and purged.** |
| **Fleischmann & Pons, Feb 1985** | A 1 cm Pd cube reportedly melted and partly vaporised, destroying the bench and pitting the concrete floor about 4 in deep. | Disputed. F&P offered an ignition interpretation and urged "extreme caution". | Keep Pd volume small and fit a thermal cutoff for unattended runs. Loaded Pd is a **fuel store**. |
| **Biberian, 2004** (glass tube, hollow Pd cathode) | The glass tube was shattered into many small pieces. The electrodes were unaffected. | Gas-phase explosion, probably initiated by **SWACER** from a reaction in the hollow cathode (Ruer & Biberian, JCMNS 26, 2018). | **No hollow cathodes.** Glass cells need secondary containment. |
| **Mizuno, 2005** (plasma electrolysis, glass cell) | Electrolyte went from 25 to 70 °C in about 10 s, then the cell exploded. Mizuno, about 1 m away, was cut by glass and deafened for a week. | H₂/O₂ in the headspace suspected. Not established. | High-voltage plasma electrolysis needs remote operation and a blast enclosure. |

A SPAWAR co-deposition cell also had a "catastrophic thermal event" after 3 days. Co-deposition cells need the same thermal cutoffs and enclosures.

---

## 3. Key hazard physics (design numbers)

**Faradaic gas generation** (D₂O, 25 °C, 1 bar), per amp of cell current:

- D₂: 5.18 µmol/s, or **7.6 mL/min**
- O₂: 2.59 µmol/s, or 3.8 mL/min
- Total: **11.4 mL/min (0.68 L/h)**
- D₂O consumed in an open cell: **9.0 g per amp-day**. A 60-day run at 1 A uses about 540 g.
- Recombination heat inside a closed cell: **1.53 W per amp** (thermoneutral voltage of D₂O ≈ 1.527 V)

**Flammability.** H₂ in air: LFL 4%, UFL 75%, detonable 18.3–59%, minimum ignition energy **0.017 mJ**, autoignition about 585 °C. H₂ in O₂: about 4–94%. Payman & Titman (Nature, 1936) measured D₂ limits in air and O₂ and found them similar. **Treat D₂ as H₂ for design.** The limiting O₂ concentration is about 5%, so a D₂-rich atmosphere is non-flammable only while O₂ stays ≤5%. Industrial recombiners start working at 0.2–4.4% H₂, depending on temperature, pressure and steam.

**Headspace energy** (stoichiometric 2D₂+O₂, 1 bar, reacting to D₂O vapour, 249 kJ/mol):

| Headspace volume | 10 mL | 50 mL | 100 mL | 500 mL |
|---|---|---|---|---|
| Energy | 67 J | 0.34 kJ | 0.67 kJ | 3.4 kJ |
| TNT-equivalent | 0.016 g | 0.08 g | 0.16 g | 0.8 g |

**Pressure rise from ignition in a closed vessel:**

- Adiabatic isochoric complete combustion: about **8–10× P₀**
- Chapman–Jouguet detonation of stoichiometric H₂/O₂: about **18–19× P₀**
- Reflected detonation: about **45–50× P₀**

**Closed cell with a failed recombiner:** dP/dt = 0.23 bar/min per amp in 50 mL of headspace (it scales as I/V). If the cell starts O₂-free at absolute pressure P_fill, the O₂ mole fraction after a fractional pressure rise r = ΔP/P_fill is x_O₂ = r / 3(1+r). **Keeping x_O₂ ≤3% requires r ≤ 0.10**, whatever the volume or current. Example: 2 bar D₂ fill, 50 mL, 1 A. The trip at ΔP = 0.2 bar is reached in about 105 s, with O₂ at 3%. Without a trip, O₂ reaches 5% in about 3 min.

**Stored chemical energy in loaded Pd:** 1 cm³ of Pd holds 0.113 mol Pd. At D/Pd = 0.9 that is 0.051 mol D₂, or **12.7–15 kJ/cm³** if burned. A 1 mm × 2 cm wire (0.016 cm³) holds about 0.2 kJ. The F&P 1 cm cube held about 13 kJ.

**Hypothetical radiation dose (a safety bound, not an expectation):** 1 W of textbook D–D fusion is about 8.5×10¹¹ n/s, or about 6.8×10⁶ n/cm²/s at 1 m, which is about **9 Sv/h**. Emission of this size has never been seen; reported neutron signals sit at or near background. A neutron monitor wired to a power trip still belongs in the design.

**Hydride expansion:** β-PdD expands by about 11% in volume (about 3.5% linear), so a 20 mm wire grows about 0.7 mm. Storms and a Pd-alloy patent recommend cathodes that expand ≤12% in volume.

---

## 4. Hazard register

Likelihood: **H** = expected at some point over a multi-month campaign without the mitigation; **M** = plausible; **L** = unlikely.

| # | Hazard | Cause | Consequence | Likelihood | Mitigation | Design implication (geometry) |
|---|---|---|---|---|---|---|
| 1 | D₂/O₂ deflagration or detonation, open cell | Off-gas builds up in the headspace, enclosure or room; vent blocked by condensate or a valve | Glass fragments, burns, hearing damage | M (3 known events) | Continuous N₂/Ar sweep; H₂ sensor alarm at ≤1% (NFPA 2); flame arrestor; shield | Headspace **≤20 mL**. Sweep **≥0.8 L/min per A** (≤1% D₂), ideally 2 L/min per A (≤0.4%). Vent bore ≥6 mm, no valves, falls toward the cell so condensate drains back. |
| 2 | Recombiner failure, closed cell | Catalyst flooded, wetted or poisoned (Cl⁻, S, Si, organics); undersized catalyst | Stoichiometric mixture builds, then ignites on a hot spot or when the catalyst dries out and reactivates | M–H | Hydrophobic Pt/Pd catalyst (PTFE-bound); pressure-trend interlock; thermocouple on the catalyst bed | Recombiner at the **top of the headspace, ≥20 mm above maximum liquid level, behind a splash baffle**. Capacity ≥2× the stoichiometric gas rate at I_max. Trip at **ΔP ≥ 0.1·P_fill**. |
| 3 | Recombiner hot-spot ignition | A burst of accumulated gas drives the catalyst above about 560–585 °C | Ignition inside the vessel | M | Limit catalyst mass; metal-mesh heat spreader; sintered-metal flame arrestor around the bed; keep the bed below 180 °C (hydrophobic coatings decompose above that) | Thermal path from the catalyst to the cell wall. Arrestor separating the catalyst from the bulk headspace. |
| 4 | Organic residue in O₂-enriched gas (the SRI lesson) | Machining oil, silicone grease, o-ring lubricant, fingerprints | Explosive oxidation; also contamination of the cathode surface | M | Oxygen-clean to ASTM G93; solvent plus ultrasonic clean; vacuum bake; lint-free gloves | Wetted and headspace materials limited to **PTFE, PFA, FEP, PCTFE (Kel-F), quartz, Pt, Pd, 316L**. No elastomers in the headspace; use metal or PTFE seals. |
| 5 | Stored-energy release from loaded Pd | Large cathode; hollow cathode; sudden exposure of loaded Pd to air or O₂ | Meltdown or explosion (F&P 1985, Biberian 2004) | L–M | Deload by anodic stripping or slow outgassing under inert gas before opening; no unattended operation without interlocks | Cathode volume **≤0.3 cm³** (≤4 kJ). **No sealed voids.** Supports must allow 3.5% linear growth. |
| 6 | Pressure-vessel rupture, closed cell | Over-pressure from gas generation, heating or ignition | Shrapnel | L–M | Relief valve plus burst disk vented to the hood; hydrotest at 1.5× MAWP; burst ≥4× MAWP | MAWP ≥10× P_fill (isochoric combustion) when the mixture is controlled, **≥20× P_fill** if detonation is credible. Headspace ≤50 mL. Metal body (316L) or thick PTFE inside a steel shell. |
| 7 | Electrical shock or arc | DC supplies of 10–300 V, especially in plasma electrolysis; wet benches | Electrocution, burns | M | Design for SELV **≤60 V DC**; current limit; fuses; GFCI on the AC side; enclosure interlock; E-stop | Electrode gap and electrolyte conductivity chosen so V_cell ≤ 30–40 V at I_max (see §9). Plasma electrolysis only inside an interlocked enclosure. |
| 8 | Caustic electrolyte; Li metal | 0.1–1 M LiOD (pH about 13); Li + D₂O is exothermic and releases D₂ | Eye or skin burns; fire | M | Goggles plus face shield, nitrile or neoprene gloves, eyewash; add Li in ≤100 mg pieces to chilled D₂O under Ar in a hood (69 mg Li per 100 mL at 0.1 M releases about 120 mL D₂) | Fill ports sized so filling and draining can be done without tipping the cell. |
| 9 | Aqua regia (etching) | NOCl and NO₂ fumes; gas pressure in a closed container; reacts violently with organics | Inhalation injury; bottle bursts | M | Mix fresh, small volumes, in a hood; never cap; quench into water; dedicated waste stream | Cathode fixtures that can be removed for etching and remounted without touching the surface. |
| 10 | Pyrophoric Pd black, nanopowder or co-deposit | High-surface Pd loaded with H or D dries out in air | Ignites solvents, ignites itself | M | Keep wet; passivate slowly with dilute air in N₂; handle under inert gas | Co-deposition cells must be drained and flushed under inert gas. |
| 11 | High-pressure D₂ and hydrogen embrittlement (Sieverts rigs, gas loading, cell prefill) | Cylinder at 2,000+ psi; regulator failure; high-strength steel, Ti or Ni-alloy fittings | Jet fire, over-pressure, cracked fittings | M | Lecture-bottle cylinder in a ventilated cabinet; regulator with relief; flow restrictor; He leak check; 316L with VCR fittings | All downstream parts 316L and rated above the relief set pressure. |
| 12 | Glass failure | Thermal shock; pressure; fragments thrown by a gas event | Lacerations | M | Polycarbonate enclosure ≥6 mm thick; prefer PTFE/PFA cell bodies | Secondary containment around every glass part. |
| 13 | Electrolyte dry-out or boiling | 9 g/(A·day) electrolysis loss plus evaporation (D₂O boils at 101.4 °C) | Electrodes exposed, arcing, ignition | H over long runs | Level sensor; automatic top-up; 90–95 °C trip | Reserve for ≥30 days at I_max, or automatic top-up. |
| 14 | Ionising radiation (hypothetical n/γ; X-rays from >5 kV discharges) | Nuclear process (unconfirmed); high-voltage discharge | Dose | L | Neutron plus gamma monitor with alarm and trip at 10× background; bubble dosimeters; no licensed sources | Room for a moderated He-3 tube or EJ-309 within 10–30 cm, and Pb around the NaI. |
| 15 | Tritium | Background in D₂O (61–2,500 Bq/L), enriched up to about 5× by electrolysis | Negligible dose; **false-positive** risk | M (for data) | Assay each D₂O lot; vent off-gas | Sampling port; closed cell keeps the tritium inventory accounted for. |
| 16 | Unattended multi-week runs | Power cut, PC crash, sensor failure | Any of the above, undetected | H | UPS; watchdog that fails safe (power off); remote alarms | Every trip wired in hardware, not only in software. |

Minor: PdCl₂ is a skin sensitiser, and the hot 6.25 M NaOH used to etch CR-39 is caustic. Use PPE and a hood.

---

## 5. Materials and price table

Metal spot prices: **Pd $1,198/ozt (Kitco, 2026-09-28); $1,265/ozt (Umicore, 2026-09-25)**, which is $38.5–40.7/g. **Pt $1,726/ozt (2026-09-28)**, which is $55.5/g. Pd density 12.02 g/cm³; Pt 21.45 g/cm³. Small fabricated pieces usually sell at 2–5× their melt value (est.). Scrap Pd can be sold back to refiners for most of its metal value.

| Item | Spec | Suppliers | Metal content → melt value | Expected purchase price | Source / date |
|---|---|---|---|---|---|
| Pd wire 0.5 mm | 99.95–99.99% | Surepure, Goodfellow, Thermo (Alfa) | 0.024 g/cm → 10 cm ≈ $9 | $40–100 per 10 cm (est.) | Goodfellow/Surepure listings (sourced availability) |
| Pd wire 1.0 mm | 99.9% (Thermo AA10280BS, 10 cm); 99.95% | Thermo/Fisher, Goodfellow | 0.094 g/cm → 10 cm ≈ $36 | $100–250 per 10 cm (est.) | Fisher listing (sourced availability; price behind login) |
| Pd wire 2.0 mm | 99.95% (Thermo AA14761BU, 25 cm) | Thermo/Fisher | 0.378 g/cm → 25 cm ≈ $364 | $600–1,200 (est.) | Fisher listing |
| Pd foil 25 µm, 25×25 mm | 99.95% "light tested" | Sigma/Goodfellow (GF03781863) | 0.19 g ≈ $7 | $60–150 (est.) | Sigma listing |
| Pd foil 50–125 µm, 25×25 mm | 99.95%, as rolled | Sigma/Goodfellow (GF86358077, GF68759061); Thermo 011515 (0.1 mm) | 0.38–0.94 g ≈ $15–36 | $100–300 (est.) | Sigma, Thermo listings |
| Pd rod | 4 mm × 20 mm | Goodfellow, ESPI, Surepure | 3.0 g ≈ $116 | $250–500 (est.) | est. |
| Pd sponge or powder | 99.95%, −60 mesh | Thermo, Sigma, ESPI | $/g ≈ spot | 1.3–2× spot (est.) | est. |
| Pd nanopowder | <25 nm | Sigma, US Research Nanomaterials | — | $50–150/g (est.); **pyrophoric when loaded** | est. |
| Pd–Ag 77/23 foil | Hydrogen-membrane alloy | Goodfellow, Johnson Matthey | — | 1–1.5× Pd foil (est.) | est. |
| Pd–Rh, Pd–B, Pd–Ce | Research alloys (NRL Pd–B; ENEA Pd–Ce) | **Custom melt only** (Goodfellow, ESPI, Ames Lab) | — | $1,500–5,000 per melt, with lead time (est.) | est. |
| PdCl₂ (co-deposition) | 99.9% | Sigma, Thermo | 0.13 g per 25 mL of 0.03 M | $10–30 per run (est.) | est. |
| Pt wire anode | 0.5 mm × 1 m, 99.95% | Thermo (AA10285, 1.0 mm), Goodfellow | 4.2 g ≈ $234 | $450–900 (est.) | Fisher listing |
| Pt gauze / mesh | 52 mesh, 25×50 mm | Thermo, Goodfellow | — | $300–700 (est.) | est. |
| **D₂O 99.9 atom % D** | Low-tritium grade preferred | Sigma 151882, Cambridge Isotope, United Nuclear, Heavy Water Board (India) | — | Bulk ex-factory global average **$646/kg (2024)**; research grade **$300–500/L (2025)**; retail 100 g about $1–3/g (est.); United Nuclear lists small lots at $29–39 | IndexBox; DataM Intelligence; unitednuclear.com |
| LiOD | 7.5 wt% in D₂O, ≥98 atom % D (Sigma 347450), **or** Li metal dissolved in D₂O | Sigma, SCBT | — | $150–400 per 25–50 g solution (est.); Li metal $2–5/g (est.) | Sigma listing |
| Au foil (co-deposition substrate) | 25–50 µm | Goodfellow, Sigma | — | $100–250 per 25×25 mm (est.) | est. |
| Ni mesh; ZrO₂, CaO powders; constantan (Cu55Ni44Mn1) wire 0.1–0.2 mm | Minor or alternative materials | McMaster, Goodfellow, Sigma, Omega | — | $30–100 each (est.) | est. |
| PTFE/PFA cell body and parts | Machined | Local shop | — | $100–600 (est.) | est. |
| Recombiner catalyst | 0.5% Pt or Pd on alumina, hydrophobised | Sigma, Alfa | — | $50–200 (est.) | est. |
| Burst disk, relief valve, 316L VCR fittings | Rated ≥ MAWP | Swagelok, Fike | — | $500–2,000 per closed cell (est.) | est. |
| D₂ gas | 99.8%, lecture bottle | Airgas, Linde | — | $300–700 (est.) | est. |

Avoid LiOH monohydrate as a LiOD substitute because it brings H into the electrolyte. Anhydrous LiOH adds about 0.1 at% H at 0.1 M.

---

## 6. Detector and instrumentation price table

| Item | Role | Price | Source / date |
|---|---|---|---|
| RadiaCode 103 (CsI(Tl), 8.4% FWHM) | Hobby gamma spectrometer and dose logger | **$319** | radiacode.com (sourced, 2026) |
| 2"×2" NaI(Tl) with digital MCA (Ortec digiBASE) | Gamma spectroscopy | **About $6,000 plus $2,000 for Maestro software; $9–10k complete** | Maximus Energy (2022) |
| 2"×2" NaI plus low-cost MCA (Spectrum Techniques UCS-30, used Bicron) | Gamma spectroscopy | $2,000–5,000 (est.) | est. |
| He-3 proportional tube | Thermal neutron counting in a moderator | Used $200–1,200 (e.g., SNM-16 $199; NEUTRON-LITE kit **$1,195**). New 4 atm tube $3–8k (est.). The He-3 gas shortage pushed prices from $40–85/L before 2009 to over $2,000/L | LabX, eBay, Chemistry World, GAO-11-753 |
| Ludlum 42-41L PRESCILA | Proton-recoil fast/thermal neutron dose probe | **$2,880** (2017 list) | Ludlum price list 2017 |
| Ludlum 12-4 rem ball | Neutron dose-equivalent survey meter | About $2,000+ (older forum figure); likely $4–6k new (est.) | fusor.net |
| BTI BD-PND bubble detectors | Personal neutron dosimeter; gamma-blind | About **$300 each** (2015 forum); likely $350–500 now (est.) | bubbletech.ca, lenr-forum |
| EJ-309 liquid scintillator, 2"×2" with PMT | n/γ pulse-shape discrimination, fast neutrons | $3–6k (est.); digitizer (e.g., CAEN desktop) $6–12k (est.) | Eljen (high flash point, 144 °C) |
| CR-39 (TASL Tastrak; Fukuvi) | Charged-particle track detector | $3–15 per chip; $100–300 per sheet (est.) | tasl.co.uk |
| CR-39 readout | Microscope (manual) or TASLImage (automated) | $300–2,000 (manual, est.); automated readers are far dearer | est. |
| H₂ detector (catalytic or thermal-conductivity, ≤1% alarm) | NFPA 2 compliance | $200–1,500 (est.) | est. |
| Gamry Reference 3000 | Potentiostat/galvanostat with EIS | **$17,995** (2015 list) | Scribd price list |
| Gamry Interface 1010E; BioLogic SP-200 | Potentiostat | $10–30k (est.) | est. |
| Keithley 2450/2460 SourceMeter | Precision constant current with 4-wire sense | $6–10k (est.) | est. |
| Bench PSU (Rigol DP832, Siglent SPD3303X) | Hobby constant current | $300–600 (est.) | est. |
| Keysight DAQ970A plus multiplexer cards | Precision DAQ (RTDs, V, I) | $3–5k (est.) | est. |
| Pt100 Class A 4-wire RTDs / calibrated thermistors | Calorimetry temperatures | $20–100 each (est.) | est. |
| TEC modules, 40×40 mm, for a Storms-style Seebeck calorimeter | Heat-flux sensing | $3–15 each; 20–60 needed (est.) | Storms, "How to make a cheap and effective Seebeck calorimeter" |
| Recirculating chiller (±0.02 °C) | Calorimeter heat sink | $3–8k (est.) | est. |
| Vacuum tube furnace, 1200 °C, with rotary pump | Pd annealing | $4–10k (est.) | est. |
| Glovebox / glove bag (Ar) | D₂O and Li handling | $100 (bag) to $40–60k (box) (est.) | est. |

**Outsourced services (est.):**

| Service | Price |
|---|---|
| SEM/EDS | $100–400/h |
| EBSD grain and texture maps | $200–500 per sample |
| AFM | $100–200/h |
| XRD texture | $100–300 per sample |
| ICP-MS impurity analysis | $50–200 per sample |
| GDMS of bulk Pd | $500–1,000 per sample |
| Vacuum heat-treat | $100–300 per batch |
| Pd sputtered films | $500–2,000 per run |
| Tritium by liquid scintillation counting | $50–300 per sample |
| **⁴He by high-resolution mass spectrometry** (must resolve D₂ at 4.0282 u from ⁴He at 4.0026 u; few labs offer it) | $500–2,000 per sample |

---

## 7. Budget tiers

### (a) About $5k: hobbyist or garage

| Item | $ |
|---|---|
| Pd: 2× 1 mm × 5–10 cm wire plus 1 foil | 500 |
| Pt anode wire | 500 |
| D₂O, 250–500 g | 500 |
| LiOD solution or Li metal; PdCl₂, LiCl, Au substrate (co-deposition) | 400 |
| PTFE cell bodies, fittings, vent, bubbler | 400 |
| Bench PSU with 4-wire resistance measurement (for D/Pd) | 500 |
| RTDs plus 24-bit ADC logger; insulated isoperibolic jacket | 500 |
| CR-39 chips, NaOH etch setup, microscope | 600 |
| RadiaCode 103 plus 2 BD-PND bubble detectors | 1,000 |
| H₂ sensor, polycarbonate enclosure, PPE, eyewash | 600 |
| **Total** | **about 5,500** |

**Can do:** Galileo co-deposition with CR-39, a loading demonstration by resistance ratio, and isoperibolic calorimetry at about ±50–150 mW (2–5%). **Cannot do:** a credible heat claim, neutron counting above background, a safe closed cell, or cathode annealing. Open, swept, vented cells only, never left unattended.

### (b) About $50k: small lab (recommended)

| Item | $ |
|---|---|
| Pd from 2–3 vetted lots (99.95–99.99%), wires and foils; one Pd–Ag or alloy sample | 4,000 |
| D₂O, 2 kg (low-tritium lot, pre-assayed) plus H₂O controls | 2,000 |
| Pt anodes and mesh | 1,500 |
| Two closed PTFE-lined 316L cells: recombiner, burst disk, relief valve, pressure transducer, VCR fittings | 5,000 |
| Storms-style Seebeck calorimeter, ×2 (active and H₂O control): TECs, Al box, insulation | 3,000 |
| Recirculating chiller | 5,000 |
| Keithley 2460 or Gamry 1010E (constant current plus EIS) | 10,000 |
| Keysight DAQ970A with RTDs, pressure logging, and hardware interlock relay board | 5,000 |
| Vacuum tube furnace and pump | 6,000 |
| Used He-3 tube in HDPE moderator plus counting electronics | 2,000 |
| 2"×2" NaI with MCA | 3,000 |
| CR-39 program plus bubble detectors | 1,000 |
| Outsourced SEM/EBSD/ICP-MS | 2,500 |
| Tritium LSC on D₂O lots and runs (about 10 samples) | 1,000 |
| Safety: H₂ detection, enclosure, UPS, watchdog, gas cabinet | 2,000 |
| **Total** | **about 53,000** |

**Can do:** Seebeck calorimetry at about ±20–60 mW (0.2–0.7%, per Storms) with paired D₂O and H₂O cells. Also in-house annealing and etching, loading tracked by resistance plus EIS, closed cells with the ΔP interlock, neutron and gamma trips, and tritium assays. **Missing:** ⁴He measurement (outsourcing adds $5–15k), parallel statistics, and staff.

### (c) About $500k: professional

| Item | $ |
|---|---|
| 1.5 FTE for 2 years (electrochemist plus technician) | 250,000 |
| 6–8 parallel cells with dual calorimetry (Seebeck plus mass-flow), multichannel potentiostat (BioLogic VMP-3e class) | 70,000 |
| Pd program: custom melts (Pd–B, Pd–Rh, Pd–Ce), 5 lots, sputtered films | 25,000 |
| ⁴He: high-resolution QMS or sector-field MS with getter purification, **or** outsourced contract | 60,000 |
| Neutrons: EJ-309 array (4×) plus digitizers; He-3 long counter | 40,000 |
| Gammas: HPGe (used or entry-level) | 30,000 |
| Glovebox, tube furnaces, electropolishing station | 25,000 |
| Outsourced characterisation (EBSD, GDMS, TEM), tritium LSC | 15,000 |
| Safety engineering: gas cabinet, ventilation to NFPA 2, pressure-vessel review | 10,000 |
| Contingency | about 30,000 |

**What it can do:** statistically powered runs (n ≥ 20 cathodes with controls), heat measured against ⁴He, blind analysis, and external review. Most of the spending is staff. The instruments are only about a third of the total.

---

## 8. Cathode preparation recipes

### 8.1 ENEA (Violante), as used by Energetics/SKINR and Letts

Reported excess-heat reproducibility above 60%. Energetics runs using Violante foils showed excess heat in about 73% of runs.

1. **Material:** high-purity Pd, 99.95–99.99%, from a *vetted lot*. Pd varies greatly from batch to batch, and the patent notes there is no standard protocol.
2. **Cold roll to foil**, rotating the sample **120° about the surface normal between passes**. The target is about **50 µm** thick. Intermediate anneals as needed.
3. **Ultrasonic cleaning** in solvent and then DI water.
4. **Vacuum anneal:** ribbons at **450 °C for 30 min, then 850 °C for 90 min** in vacuum. The final 850 °C anneal enlarges grains and relieves stress.
5. **Target microstructure:** surface grains mostly **15–40 µm, ideally about 35 µm** (acceptable range 5–100 µm), with strong **<100> texture**. Verify by EBSD or XRD.
6. **Etch before every run:** **50% aqua regia for 2 min** at room temperature, or warmer if the oxide resists (maximum about 100 °C). Etching roughens the grains, deepens grain boundaries, and raises double-layer capacitance and effective surface area, which gives faster and higher loading.
7. **Characterise:** SEM/AFM for surface morphology, XPS for surface contamination.

### 8.2 Letts–Cravens (ASTI 2004; patents US 10,559,831 and 11,563,217)

1. Cold-work to about **50% thickness reduction** (for example 0.5 mm → 0.25 mm) to create dislocations.
2. Polish with a rotary brush and Al₂O₃ abrasive.
3. Etch in aqua regia for about 2 min at up to 100 °C.
4. Anneal at **800–900 °C for ≥30 min** (typically 850 °C), adjusting per lot to reach about 35 µm grains.
5. Use small cathodes (about 0.15 g) and load to roughly **10⁷ C per mole of Pd** before expecting heat (Cravens & Letts, "Enabling Criteria", a review of 167 papers from 1989–2007).
6. Prefer alloys and treatments that give **≤12% volume expansion** on loading.

### 8.3 SRI (McKubre)

1. Twice vacuum-melted Pd, machined to a wire or rod.
2. **Vacuum anneal at 800 °C for 3 h**, or 850 °C for 2.5 h, then **backfill with D₂ at temperature**.
3. Run criteria: **D/Pd ≥ 0.875–0.9** held for an incubation period. Current density above threshold: excess power scales with *i* above about **100 mA/cm²**. The threshold is clearer and higher with a wire cathode and concentric anode than with parallel plates.
4. Measure D/Pd in situ by 4-wire resistance ratio R/R₀, which requires a wire at least about 2 cm long.

### 8.4 Fleischmann–Pons (baseline geometry)

- Cathode: Pd rod, **2 mm diameter × 12.5 mm**, spot-welded to a Pt lead.
- Anode: Pt wire, concentric.
- Electrolyte: **0.1 M LiOD** in D₂O.
- Cell: Pyrex Dewar, silvered in the upper part; electrodes held by a Kel-F plug.
- Operation: constant current, stepped from **10 up to 1,000 mA/cm²**, with D₂O topped up.

### 8.5 Storms' guidelines ("How to produce the Pons–Fleischmann effect")

- Reject cathodes that crack or expand grossly on a trial load–deload cycle.
- Load near room temperature at moderate current density, then raise the temperature step by step toward boiling. Limit high-current pre-electrolysis "activation" to a couple of cycles.
- Surface deposits (Li, Pt from the anode, Si from glass) shape the outcome. Control them.

### 8.6 SPAWAR Pd/D co-deposition (Szpak/Mosier-Boss; Galileo protocol)

1. Electrolyte: **0.03 M PdCl₂ + 0.3 M LiCl in D₂O**. Substrate: Au foil or wire (Ag, Ni and Pt were also used). Anode: Pt.
2. Galvanostatic steps: **1.0 mA/cm² for 8 h → 3 mA/cm² for 8 h → 5.0 mA/cm² until all Pd²⁺ is reduced** (the solution turns colourless). Then **30–50 mA/cm² for 2–3 h** to keep gas visibly evolving, and afterwards step up the current. An external magnetic or electric field is optional.
3. CR-39 goes against the cathode. Standard etch is 6.25 M NaOH at 70 °C for about 6 h. Run H₂O and blank (no Pd) controls to separate tracks from chemical or mechanical pitting (see "Interpreting CR-39 Detectors used in Pd/D Co-deposition", JCMNS).
4. **Safety:** the co-deposit is pyrophoric when dry, and one catastrophic thermal event was reported.

### 8.7 Other programs

- **Google/UBC/MIT:** Benck et al. (Chem. Mater. 2019) reached **H/Pd near 1** by electrochemical insertion. Composition is a dynamic insertion/evolution balance, so it has to be measured operando. The Nature 2019 perspective found no excess heat. The 2025 UBC follow-up found that electrochemical loading **raised plasma-driven D–D fusion rates by about 15%**, but that route requires a registered radiation-producing machine (§9).
- **MFMP:** open, live replications of Celani's constantan wires reported 12.5% and then 5.3% apparent excess (quantumheat.org). Useful as a model for open data, not as confirmation of an effect.

### Common anneal and etch parameters (synthesis)

| Parameter | Value |
|---|---|
| Atmosphere | ≤10⁻⁵ mbar vacuum or flowing Ar/H₂ (use D₂ backfill if it is available) |
| Temperature | **850 °C** (range 800–900) |
| Time | **1.5–3 h** |
| Cooling | Furnace cool |
| Handling | Fluoropolymer or quartz boats only (no graphite, no Ni); gloved |
| Etch | 50% aqua regia, 2 min, just before mounting; rinse in DI water, then D₂O |

---

## 9. Regulatory notes (US focus)

- **D₂O:** no NRC license is needed to buy or possess it domestically. Large acquisitions may be reported by vendors. Institutional vendors (Sigma, Fisher, Cambridge Isotope) generally require a verified institutional account. **Export** is controlled under 10 CFR 110: a general license covers ≤10 kg deuterium (about 50 kg D₂O) per shipment and ≤200 kg per year per country to non-embargoed destinations.
- **D₂ gas:** shipped as a DOT Class 2.1 flammable gas (UN1957). Storage and use fall under NFPA 55, NFPA 2 and OSHA 29 CFR 1910.103.
- **Tritium:** the exempt quantity for H-3 is **1,000 µCi (37 MBq)** in 10 CFR 30.71 Schedule B (about 2×10¹⁶ atoms). Claimed LENR tritium yields (10¹¹–10¹⁵ atoms) fall well below it. Commercial D₂O already carries **61–2,500 Bq/L** of tritium, rising with purity, and electrolysis enriches it up to about 5×. **Assay every lot.**
- **Radiation sources:** use only exempt check sources, or borrow calibration sources at a licensed facility. **Avoid AmBe, Cf-252 and other licensed neutron sources.**
- **Radiation-producing machines:** anything that accelerates electrons or ions to >5–10 kV (glow discharge, fusor, ion beam) can produce X-rays or neutrons. Most states require it to be registered with the state radiation control program. That keeps accelerator-based D–D measurements out of the low tiers.
- **Codes and standards:** NFPA 2 (H₂ detection at ≤1% by volume; ventilation ≥1 scf/min/ft² of floor area where H₂ is stored or used); NFPA 45 (labs); NFPA 55; NFPA 70/70E (electrical); OSHA 1910.1450 (chemical hygiene plan); ASME BPVC VIII (vessels over 15 psig; very small-diameter vessels are exempt, but follow its design margins anyway); ASTM G93 (oxygen cleaning).
- **Waste:** recover Pd through a refiner. Neutralise LiOD. Aqua regia is hazardous waste (corrosive and oxidiser). Tritiated water below exempt levels needs no special handling, but document it.

---

## 10. Design constraints the geometry must satisfy

1. **Headspace volume:** ≤20 mL for open cells and ≤50 mL for closed cells. This caps the worst-case stoichiometric energy at about 0.13 kJ (open) or 0.34 kJ (closed) at 1 bar.
2. **Open-cell vent path:** inert sweep ≥0.8 L/min per amp at I_max (≥2 L/min per A preferred). Vent bore ≥6 mm, with no valves or restrictions. Condensate must drain back into the cell. Include a bubbler or check valve and a flame arrestor, and exhaust to a hood or outdoors.
3. **Closed-cell pressure rating:** MAWP ≥10× the absolute fill pressure, or ≥20× if detonation cannot be excluded. Hydrotest at 1.5× MAWP. Relief valve set ≤ MAWP and a burst disk at 1.1× relief, both vented to the exhaust.
4. **Closed-cell interlock:** trip cell current at **ΔP ≥ 0.10 × P_fill** (absolute), with a response time under 60 s. Start the cell O₂-free (evacuate, then backfill with D₂).
5. **Recombiner placement:** top of the headspace, ≥20 mm above maximum liquid level, behind a splash baffle, inside a sintered-metal flame arrestor, and thermally sunk to the wall. Thermocouple-monitored, with a trip at 150 °C. Hydrophobic catalyst sized for ≥2× the stoichiometric gas rate at I_max.
6. **Cathode volume:** ≤0.3 cm³, with no enclosed voids or hollow cathodes. Mounts must allow **3.5% linear growth** (≈0.7 mm per 20 mm).
7. **Wire cathodes for loading measurement:** ≥20 mm active length, with separate current and voltage leads (4-wire) for R/R₀.
8. **Current-density uniformity:** concentric anode (cylinder or helix) around a wire cathode, with the gap uniform to ±10%. Anode area ≥3× cathode area. Parallel plates only with guard geometry.
9. **SELV:** V_cell ≤60 V DC at I_max. For 0.1 M LiOD (ρ ≈ 60 Ω·cm in D₂O; bubbles can double it), a 1 mm × 20 mm wire with a 10 mm-radius concentric anode is about 14 Ω, or about 8 V at 0.31 A (500 mA/cm²). A 1 cm² plate needs a gap of **≤5 mm** at 1 A.
10. **Wetted and headspace materials:** PTFE, PFA, FEP, PCTFE, quartz, Pt, Pd, 316L only, oxygen-cleaned. No elastomers or greases in the headspace. Avoid borosilicate in prolonged contact with hot LiOD if Si contamination is to be controlled.
11. **Electrolyte inventory:** a reserve of ≥30 days × 9 g/(A·day) for open cells, or automatic top-up. Level sensor interlocked to the current.
12. **Thermal limits:** electrolyte ≤95 °C (D₂O boils at 101.4 °C), with a hardware thermal cutoff independent of the DAQ.
13. **Secondary containment:** polycarbonate enclosure ≥6 mm thick around every glass part and around closed cells. Keep ≥1 m of operator standoff during current steps and pressure changes. **No opening, tilting or lifting of any cell until it has been deloaded, purged and cooled.**
14. **Detector access:** leave room for a moderated neutron detector within 10–30 cm of the cathode, a CR-39 position adjacent to the cathode, and a lead-shielded NaI. Avoid large metal masses between cathode and detectors. Neutron or gamma alarm at 10× background trips the power.
15. **Loading gas systems** (if gas-loaded or prefilled): 316L with metal-gasket fittings, a regulator with relief, a flow restrictor, and a lecture-bottle-scale cylinder in a ventilated cabinet. Ceiling-level H₂ sensor with a ≤1% alarm.
16. **Fail-safe control:** every trip (pressure, temperature, level, H₂, radiation, E-stop, watchdog) removes cell power in hardware. UPS on the monitoring system.
17. **Calorimetric boundary:** open cells must correct for the enthalpy carried off by gas (I × 1.527 V). Closed cells keep the recombination heat inside the boundary, at the price of constraints 3–5.

---

## 11. References

- Forensic analyses of explosion debris from the 2 Jan 1992 Pd/D₂O incident at SRI (LLNL): https://www.osti.gov/biblio/212468 ; https://digital.library.unt.edu/ark:/67531/metadc665463/ ; https://inis.iaea.org/records/pjf88-kav54
- The January 2, 1992, explosion in a deuterium/palladium electrolytic system at SRI: https://inis.iaea.org/records/jw89d-fvx92
- Deseret News, "Equipment failure blamed..." (1992): https://www.deseret.com/1992/1/4/18960119/equipment-failure-blamed-in-lab-explosion-that-killed-cold-fusion-researcher-br/
- Michael McKubre (incident account): https://en.wikipedia.org/wiki/Michael_McKubre
- Ruer & Biberian, "Reanalysis of an Explosion in a LENR Experiment," JCMNS 26 (2018): https://jcmns.org/article/72472-reanalysis-of-an-explosion-in-a-lenr-experiment
- Mizuno 2005 explosion: https://sciencespot.co.uk/2005-3-10-18423-3423.html ; https://en.wikipedia.org/wiki/Tadahiko_Mizuno
- F&P 1985 meltdown: https://matheasy.substack.com/p/martin-fleischmann-recounts-the-1985 ; https://encyclopedia.pub/entry/31286
- Hydrogen safety data: https://en.wikipedia.org/wiki/Hydrogen_safety ; https://www.aiche.org/sites/default/files/docs/pages/the_elemental_-_hydrogen_flammability.pdf ; https://www.osti.gov/servlets/purl/1721461
- Payman & Titman, "Limits of Inflammability of Hydrogen and Deuterium in Oxygen and in Air," Nature 137 (1936): https://www.nature.com/articles/137190a0
- Recombiner catalyst behaviour: https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/4350610 ; https://glc.ans.org/nureth-16/data/papers/14071.pdf
- NFPA 2 guidance: https://www.labmanager.com/lab-gas-generator-safety-managing-hydrogen-flammability-pressure-risks-and-ventilation-requirements-35520 ; https://h2tools.org/bestpractices/laboratory-safety/laboratory-design/ventilation ; https://altfuelgarage.frontierenergy.com/wp-content/uploads/2018/05/Hydrogen-Requirements-Best-Practices-NFPA-2.pdf
- Pd price: https://www.kitco.com/charts/palladium ; https://pmm.umicore.com/en/prices/palladium/
- Pt price: https://tradingeconomics.com/commodity/platinum
- Pd/Pt catalog listings: https://www.fishersci.com/shop/products/palladium-wire-1-0mm-0-04in-dia-99-9-metals-basis-thermo-scientific/AA10280BS ; https://www.fishersci.com/shop/products/palladium-wire-2-0mm-0-08-in-dia-99-95-metals-basis-thermo-scientific-1/AA14761BU ; https://www.sigmaaldrich.com/US/en/product/aldrich/gf68759061 ; https://www.thermofisher.com/order/catalog/product/011515.FF ; https://www.goodfellow.com/usa/palladium-spooled-wire-group ; https://www.surepure.com/Palladium-Wire-Rod/a/8,1
- D₂O pricing and market: https://www.indexbox.io/search/price-for-heavy-water-deuterium-oxide-the-united-states/ ; https://www.datamintelligence.com/research-report/heavy-water-market ; https://unitednuclear.com/chemicals-metals-c-69/deuterium-oxide-p-135.html ; https://www.sigmaaldrich.com/US/en/product/aldrich/151882
- LiOD: https://www.sigmaaldrich.com/catalog/product/aldrich/347450
- 10 CFR 110 export: https://www.ecfr.gov/current/title-10/chapter-I/part-110 ; https://biennialsandeducation.org/can-you-buy-heavy-water
- 10 CFR 30 exemptions: https://www.ecfr.gov/current/title-10/chapter-I/part-30 ; https://www.morganlewis.com/blogs/upandatom/2026/06/material-changes-nrc-proposes-overhaul-of-byproduct-materials-regulations
- Tritium in commercial D₂O: https://link.springer.com/article/10.1007/s44211-024-00615-6 ; https://files.ncas.org/erab/apx3b.htm
- He-3 shortage: https://www.chemistryworld.com/news/shortages-spur-race-for-helium-3-alternatives-/3003620.article ; https://www.gao.gov/assets/a585515.html ; https://www.labx.com/item/neutron-lite-neutron-detector-spectrometer-mca-he3-helium-3/scp-221792-f5fbc219-815d-4f6a-aed7-89f9568ee49e
- Detectors: https://radiacode.com/products/radiacode-103?lang=en ; https://maximus.energy/index.php/2022/07/02/ortec-digibase-e-gamma-spectrometer/ ; https://device.report/m/0f634454662db1c2ad91fcaf781fd0774eb9940cb1b3e566877286b340879f87.pdf ; https://bubbletech.ca/product/bd-pnd-personal-neutron-dosimeter/ ; https://www.lenr-forum.com/forum/thread/2230-bd-pnd-neutron-detectors/ ; https://eljentechnology.com/products/liquid-scintillators/ej-301-ej-309 ; https://www.tasl.co.uk/tastrak-padc.php
- Gamry pricing (2015): https://www.scribd.com/document/309879595/Ceni-IntlSalesPriceList-August2015
- Seebeck calorimetry (Storms): https://lenr-canr.org/acrobat/StormsEhowtomakea.pdf ; https://lenr-canr.org/acrobat/StormsEdescriptioa.pdf ; https://www.osti.gov/etdeweb/biblio/21066223
- Storms, "How to Produce the Pons–Fleischmann Effect": https://www.lenr-canr.org/acrobat/StormsEhowtoprodu.pdf
- Letts & Cravens, cathode fabrication: https://www.lenr-canr.org/acrobat/LettsDcathodefab.pdf ; https://lenr-canr.org/acrobat/LettsDlaserstimu.pdf
- Cravens & Letts, "Enabling Criteria": https://www.lenr-canr.org/acrobat/CravensDtheenablin.pdf ; https://www.lenr-canr.org/acrobat/CravensDfactorsaff.pdf
- Exothermically responsive cathodes patents: https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/10559831 ; https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/11563217
- ENEA/SKINR: https://www.researchgate.net/publication/283861408_Sidney_Kimmel_Institute_for_Nuclear_Renaissance ; https://mospace.umsystem.edu/xmlui/handle/10355/37370 ; https://www.researchgate.net/publication/237284370_Material_Science_on_Pd-D_System_to_Study_the_Occurrence_of_Excess_Power ; https://jcmns.org/article/72458.pdf
- McKubre calorimetry and loading: https://lenr-canr.org/acrobat/McKubreMCHcalorimetr.pdf ; https://lenr-canr.org/acrobat/KunimatsuKdeuteriuml.pdf
- F&P calorimetry and cell: https://lenr-canr.org/acrobat/Fleischmancalorimetr.pdf ; https://www.lenr-canr.org/acrobat/LonchamptGreproducti.pdf
- Co-deposition: https://www.lenr-canr.org/acrobat/SzpakScalorimetra.pdf ; https://www.lenr-canr.org/acrobat/MosierBosscharacteri.pdf ; https://jcmns.org/article/72567.pdf ; https://www.researchgate.net/publication/283569283_Condensed_Matter_Nuclear_Science_Using_PdD_Co-Deposition
- Google/UBC/MIT: https://www.nature.com/articles/s41586-019-1256-6 ; https://pubs.acs.org/doi/abs/10.1021/acs.chemmater.9b01243 ; https://phys.org/news/2025-08-room-temperature-reactor-electrochemistry-boost.html
- MFMP: http://www.quantumheat.org/ ; https://www.researchgate.net/publication/286466066_Martin_Fleischmann_Memorial_Project_status_review
- 316L hydrogen embrittlement: https://www.sciencedirect.com/science/article/abs/pii/S0921509399003196
