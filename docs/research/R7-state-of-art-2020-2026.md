# R7 — State of the Art in LENR / Condensed-Matter Nuclear Science, 2020–2026

*Research note for the cold-fusion/LENR geometry-design project. Compiled 2026-09-29.*

## 0. Provenance and how to read this document

- **Evidence grades** (as requested):
  **A** = independently replicated, mainstream-accepted;
  **B** = replicated by several groups but still debated;
  **C** = single group (possibly high-quality venue);
  **D** = anecdotal / disputed / commercial claim without independent data;
  **N** = credible null result.
- **Source tags.** During this session, direct page fetching (WebFetch) was blocked by the network egress proxy for every domain tried (nature.com, arxiv.org, aps.org mirrors, pmc, lenr-canr.org, newenergytimes, arpa-e.energy.gov, physicsworld). The session-wide web-search budget was also exhausted partway through. Therefore:
  - **[V]** = the fact was confirmed in this session from search-engine result snippets (title, venue, date, key numbers). The primary PDF was not read.
  - **[BK]** = from the author's background knowledge (training data through mid-2026). It was **not** re-verified this session. Downstream agents should check [BK] items before relying on exact numbers.
- **Bottom line on honesty.** From 2020 to 2026, the credible, reproducible physics in this field was almost entirely **beam- or plasma-driven deuterium fusion in metal hydrides at keV and sub-keV energies**. Those are conventional nuclear products at rates enhanced by the solid-state environment, and they are far too small for energy. The "excess heat" branch still has **no independently replicated result in a mainstream venue**.

---

## 1. Executive summary

1. **[A] Metal lattices enhance keV-energy D–D fusion through electron screening.** Many independent accelerator groups have measured this since ~2000 (Kasagi, Rolfs/Raiola, Czerski/Huke, and others). Screening potentials of hundreds of eV in metals are well above the ~25 eV seen for gas targets. Work in 2020–2026 refined this effect and did not overturn it. The size of the effect is still not explained quantitatively by theory.
2. **[B] A sub-keV "yield plateau" was reported by two independent groups, 2024–2026.**
   - UC Davis + LBNL (Karahadian, Colborne, Persaud, Schenkel, Munday), funded by ARPA-E. Their dual-chamber platform combines electrochemical loading with a low-energy ion beam. They report that D–D yield in Pd and Ti foils stops falling exponentially below ~2–2.5 keV. At the lowest energies they quote an enhancement of **>10¹⁸ over bare-nucleus expectations**. Published in *Nature Communications* 17:8845, 18 July 2026 [V].
   - Szczecin (Czerski et al.) independently sees a 2H(d,p)3H plateau down to ~1 keV in ZrD₂, Ti and Pd. They interpret it as thermal fusion in ion tracks [V].
   - The two groups interpret the effect differently, and absolute rates are tiny.
3. **[C] Electrochemical loading measurably increases a nuclear reaction rate.** UBC's "Thunderbird" (Berlinguette group, *Nature*, Aug 2025) used plasma-immersion D implantation into one face of a 300 µm Pd foil, with electrochemical D loading from the other face. Electrolysis raised the neutron rate by **15 ± 2 %** (peak 188 ± 2 n/s). Output was ~10⁻⁹ W of fusion power against ~15 W of input [V]. The result is real and reproducible within the group, but "modest" (Chemistry World), and it is conventional D–D fusion.
4. **[C] A claimed new D–D channel emits e⁺e⁻ pairs.** Dubey, Czerski et al. (*Phys. Rev. X* 15, 041004, 7 Oct 2025) report annihilation and bremsstrahlung radiation consistent with e⁺e⁻ pairs up to ~23 MeV from a 0⁺ threshold resonance in ⁴He. They claim this channel dominates D–D fusion below ~5 keV [V]. If true, low-energy D–D fusion would be largely "silent" in neutrons and protons. It has not been independently replicated.
5. **[C] The most consistent excess-heat claims come from Tohoku University (Iwamura, Itoh) and Clean Planet, using Ni–Cu nano-multilayers with H₂.** Reported figures [V]:
   - 4–6 W excess and 460 ± 120 kJ over 80 h, i.e. ≥410 ± 108 keV per H atom (photon-radiation calorimetry, arXiv 2311.18347 / JCMNS).
   - More than 10 keV per absorbed H, with no neutrons or gammas (*Jpn. J. Appl. Phys.* 63, 2024).
   - Infrared imaging of micron-scale hot spots (ICCF-26, 2025).

   No peer-reviewed independent replication outside the Tohoku / Clean Planet consortium was found.
6. **[D] Commercial claims have no independent public calorimetry.** These include:
   - Clean Planet's "QHe": joint boiler development with Miura Co., a 120 cm "IKAROS" module targeting 24 kW, a ¥500 M raise in Jan 2026, and mass-production targets for 2026–27 [V, press/secondary].
   - Brillouin, Rossi's E-Cat, and Aureon [BK].
7. **[N, provisional] ARPA-E's LENR Exploratory Topic produced no published excess heat or anomalous products.** The program ran ~$10 M across 8 projects from 2023 to ~2025/26. We found no publication reporting excess heat, transmutation, or nuclear products beyond conventional beam-target fusion. Its one high-profile peer-reviewed output (LBNL/UC Davis) is the plateau result in item 2. We found no public program-level final conclusion from ARPA-E.
8. **[N] The 2019 Google-funded *Nature* study remains the reference null.** It found no excess heat and no nuclear products, and electrochemistry could not reliably push loading above D/Pd ≈ 0.95. The consortium members (Berlinguette, Schenkel, Munday, Chiang) then turned to beam-plus-electrochemistry hybrids, which produced items 2 and 3.
9. **[C/D] Several single-lab nuclear-emission claims from 2020–2026 remain unreplicated:**
   - NASA Glenn "Lattice Confinement Fusion": bremsstrahlung-driven fusion in ErD₃/TiD₂ (PRC 2020 [BK]).
   - UIUC CR-39 tracks read as ~138 keV alphas from a Pd–D₂ glow discharge [V].
   - Maximus Energy: neutrons from cavitation of deuterated Ti powder (*Sci. Rep.* 2024, single author) [V].
   - Claims of ³He or D production in hydrogen–metal systems [V].
10. **The field changed methods.** Its center of gravity moved from "excess heat in Pd electrolysis" to "nuclear-rate enhancement in metal hydrides." The new work is detection-first, hypothesis-driven, publishes in *Nature*, PRX and *Nature Communications*, and runs mostly on ARPA-E, EU H2020 (CleanHME, HERMES) and philanthropic money.
11. **The leading hypotheses are:**
    - enhanced screening together with defects and locally concentrated D ("materials-driven fusion");
    - the ⁴He 0⁺ threshold resonance with e⁺e⁻ decay;
    - thermal-spike (ion-track) fusion;
    - phonon–nuclear energy transfer (Hagelstein/Metzler "nuclear Dicke" models);
    - reactions triggered by hydrogen diffusion flux (Iwamura; Czerski's 511 keV from D diffusion).
12. **Design implication.** The configuration most likely to give a *credible, measurable* signal is a **dual-chamber thin-foil geometry**:
    - a Pd foil (with Ti or Au as controls);
    - electrochemical D loading on the back face;
    - low-energy (0.3–20 keV), mass-analysed D⁺ ions or plasma on the front face;
    - close-geometry charged-particle, neutron and 511 keV/bremsstrahlung detection;
    - interleaved on/off and H/D controls.

    A heat-only design (Ni–Cu multilayer) is a lower-credibility secondary option.

---

## 2. Timeline, 2020–2026

| Date | Group | Result | Venue | Grade | URL |
|---|---|---|---|---|---|
| 2020-04 | NASA Glenn (Pines, Steinetz, Benyo, Forsley, Mosier-Boss et al.) | 2.9 MeV electron-linac bremsstrahlung on ErD₃/TiD₂ produces D–D neutrons; screened fusion and Oppenheimer–Phillips stripping proposed ("Lattice Confinement Fusion") [BK] | Phys. Rev. C 101, 044609 & 044610 | C | https://journals.aps.org/prc/abstract/10.1103/PhysRevC.101.044610 |
| 2021 | Co-deposition group (Mosier-Boss/Forsley lineage) | Neutrons from Pd/D co-deposition measured with bubble detectors | J. Electroanal. Chem. (2021) | C/D | https://www.sciencedirect.com/science/article/abs/pii/S1572665721000503 |
| 2021 (Jun) | ICCF-23 (virtual, hosted from Xiamen) [BK] | Conference | — | — | — |
| 2021-09 | US House Armed Services Cttee | FY2022 NDAA committee report asks DoD (USD R&E) for a briefing on LENR's state and military utility [BK; verify exact text] | H. Rept. 117-118 | — | https://www.congress.gov/congressional-report/117th-congress/house-report/118 |
| 2021-10 | ARPA-E | LENR workshop; NASA shows LCF gas-cycling experiments | ARPA-E workshop | — | https://arpa-e.energy.gov/sites/default/files/migrated/2021LENR_workshop_Benyo.pdf |
| 2022-07 | Czerski (Szczecin) | Theory plus data: 0⁺ threshold resonance in ⁴He with large e⁺e⁻ width explains anomalous low-energy D–D rates | Phys. Rev. C 106, L011601 | C | https://link.aps.org/doi/10.1103/PhysRevC.106.L011601 |
| 2022-07 | ICCF-24, Mountain View CA (Anthropocene Institute / Solid State Energy Summit) [BK] | Conference | — | — | https://en.wikipedia.org/wiki/International_Conference_on_Cold_Fusion |
| 2022-09 | ARPA-E | LENR Exploratory Topic call (~$10 M) | FOA | — | https://arpa-e-foa.energy.gov/Default.aspx?Search=lenr&SearchType= |
| 2023-02-17 | ARPA-E | 8 projects selected (7 organisations; U. Michigan holds two) | DOE press | — | https://www.greencarcongress.com/2023/02/20230218-lenr.html |
| 2023-04 | MIT (Metzler, Hunt, Messinger, Galvanetto, Hagelstein) | Review of claimed neutrons and "fission-daughter" low-Z elements in gas-loaded, laser-irradiated metal hydrides; motivates the ARPA-E program | SSRN preprint | — | https://dx.doi.org/10.2139/ssrn.4411160 |
| 2023-08-27/31 | ICCF-25, Szczecin | Conference (CleanHME host) | — | — | https://newenergytimes.com/v2/conferences/2023/ICCF25/ICCF-25-Book-of-Abstracts-2023.07.04.pdf |
| 2023-09 | Szczecin | Positron annihilation spectroscopy (PAS) and XRD of lattice defects in deuterated Zr with O/C impurities; defect-screening link | Materials 16, 6255 | C | https://doi.org/10.3390/ma16186255 |
| 2023-11 | Tohoku / Clean Planet (Kasagi, Itoh, Iwamura) | Photon-radiation calorimetry of Ni–Cu multilayer: 4–6 W excess; 460 ± 120 kJ in 80 h; ≥410 keV/H | arXiv 2311.18347; JCMNS | C | https://arxiv.org/abs/2311.18347 |
| 2024-01/02 | UIUC (Ziehm, Miley) | CR-39 tracks consistent with 138 ± 21 keV alphas from Pd electrodes in 10 Torr D₂ discharge at about −500 V; ~100× H₂/He controls | arXiv 2402.05117 | C/D | https://arxiv.org/abs/2402.05117 |
| 2024-02 | Szczecin | "Indications of electron emission from the D–D threshold resonance" | Phys. Rev. C 109, L021601 | C | https://link.aps.org/doi/10.1103/PhysRevC.109.L021601 |
| 2024-03-04 | Tohoku / Clean Planet (Iwamura, Itoh, Yamauchi, Takahashi) | Heat bursts in H-preloaded Ni-based multilayers on rapid heating; >10 keV/H; no n/γ | Jpn. J. Appl. Phys. 63 (3) | C | https://iopscience.iop.org/article/10.35848/1347-4065/ad2622 |
| 2024-05-21 | Maximus Energy (Fomitchev-Zamilov) | Neutrons (peak >6500 CPM, >10⁴× background) during acoustic cavitation of TiD powder in oil; author correction July 2024 | Sci. Rep. | D | https://www.nature.com/articles/s41598-024-62055-6 |
| 2024 | Unnamed group | D/H ratio ~280× natural after TiHₓ thermal cycling (mass spectrometry) | Symmetry 16, 1542 (MDPI) | D | https://www.mdpi.com/2073-8994/16/11/1542 |
| 2024-08/09 | Szczecin | Preprints: new e⁺e⁻ channel (2408.07567); thermal D–D fusion in ion tracks, ZrD₂ plateau to ~1 keV (2409.02112) | arXiv | C | https://arxiv.org/abs/2409.02112 |
| 2024-09 | Tohoku group et al. | ³He detected in Ni–Cu/ZrO₂ nanocomposites exposed to H₂ at high temperature | arXiv 2409.05382 | D | https://arxiv.org/pdf/2409.05382 |
| 2024-09-05 | CleanHME | Workshop at European Parliament, Strasbourg | — | — | https://cleanhme.eu/ |
| 2024-10 | MIT (Metzler, Hunt, Hagelstein, Galvanetto) | Known solid-state mechanisms each give up to ~10³⁰ enhancement; cascading them could reach the >10⁴⁰ needed for observable D–D rates | New J. Phys. 26 | theory | https://iopscience.iop.org/article/10.1088/1367-2630/ad091c |
| 2025-01 | MIT/Cambridge/Zurich (Hagelstein et al.) | 234-page "Models for nuclear fusion in the solid state" (nuclear Dicke model, fusion–fission energy transfer) | arXiv 2501.08338 | theory | https://arxiv.org/abs/2501.08338 |
| 2025-01-23/24 | CleanHME | Final event, Szczecin (H2020 project ends) | — | — | https://cordis.europa.eu/project/id/951974/results |
| 2025-03 | Szczecin | Electron screening in D–D on Zr with O/C contamination | Materials 18, 1331 | C | https://doi.org/10.3390/ma18061331 |
| 2025 | Biberian et al. (CleanHME) | "Excess heat in nanoparticles of nickel alloys in hydrogen" | JCMNS 38 | D | https://cleanhme.eu/?page_id=27 |
| 2025-05-26/30 | ICCF-26, Morioka, Japan | 141 attendees; 81 abstracts (44 oral, 35 posters). Highlights: Iwamura IR hot spots; ³He from H+D; 511 keV above background from Pd-D and Zr-D underground | Conference | C/D | https://iccf26.org/ ; https://www.sciengine.com/JMCC/doi/10.16084/j.issn1001-3555.2025.06.009 |
| 2025-08 | UBC (Berlinguette group) + collaborators | Thunderbird: electrochemistry raises plasma-implantation D–D neutron rate by 15 ± 2 %; max 188 ± 2 n/s | Nature | C | https://www.nature.com/articles/s41586-025-09042-7 |
| 2025-08 | ARPA-E | FY2023 Annual Report to Congress (lists LENR topic) | Report | — | https://arpa-e.energy.gov/sites/default/files/2025-09/ARPA-E%20FY%202023%20Annual%20Report.pdf |
| 2025-09 | Tohoku (Kasagi, Itoh, Shibasaki, Iwamura) | Radiant spectrum during excess heat in NiCu and Ni films | arXiv 2509.13847 | C | https://arxiv.org/pdf/2509.13847 |
| 2025-10-07 | Szczecin (Dubey, Czerski, Das H., et al.) | e⁺e⁻ channel of D–D below ~5 keV; ZrD target; 5–20 kV D⁺/D₂⁺ at 40–60 µA; NaI(Tl) + HPGe | Phys. Rev. X 15, 041004 | C | https://link.aps.org/doi/10.1103/chlp-b215 |
| 2025-12 | UC Davis + LBNL | Preprint of sub-keV plateau | arXiv 2512.06212 | B (with Szczecin) | https://arxiv.org/abs/2512.06212 |
| 2026-01 | Messinger, Metzler, Price | "Gatekeeping: a partial history of cold fusion" (history/philosophy of science) | arXiv 2601.09996 | — | https://arxiv.org/abs/2601.09996 |
| 2026-01 | Clean Planet | ~¥500 M strategic equity raise tied to QHe commercialization | Press (secondary) | D | https://newfireenergy.substack.com/p/clean-planet-moves-deeper-into-commercialization |
| 2026-05-22 | Szczecin (Czerski, Dubey, Das, Thulichery, Kowalska, Targosz-Ślęczka, Valat) | Thermal D–D fusion in metallic targets: plateau in Zr, Ti, Pd; thermal-spike + screening + resonance model | arXiv 2605.27438 | C (B with UCD/LBNL) | https://arxiv.org/abs/2605.27438 |
| 2026-07-18 | UC Davis + LBNL (ARPA-E-funded) | >10¹⁸ enhancement over bare nuclei at lowest energies; plateau below ~2 keV in Pd and Ti | Nat. Commun. 17:8845 | B | https://www.nature.com/articles/s41467-026-74421-1 |
| 2026 (mid) | Independent analyst (ResearchGate, not peer-reviewed) | "Beyond scalar screening": identifiability critique of the sub-keV plateau | Preprint | — | https://www.researchgate.net/publication/412147776 |

---

## 3. The ARPA-E LENR Exploratory Topic (2022–2025/26)

**Framing.**
- ARPA-E held a two-day LENR workshop in October 2021 [V]. It issued the Exploratory Topic call in September 2022 and announced ~$10 M for 8 projects on 17 February 2023 [V].
- The stated aim was to "break the stalemate": to establish whether on-demand, repeatable LENR with nuclear diagnostics could be achieved, or rule it out decisively [V].
- The program deliberately favoured **hypothesis-driven** tests over open-ended replication, and **peer-reviewed, top-journal publication** [V].
- It included two "capability teams" (radiation diagnostics at U. Michigan; materials characterisation at Texas Tech) to serve the experimental teams [V].

| # | Team (location) | Award | Approach | Published outcome found (as of 2026-09) |
|---|---|---|---|---|
| 1 | **Lawrence Berkeley National Lab** (T. Schenkel; with UC Davis, J. Munday) | ~$1.5 M / 2.5 yr [V] | Electrochemical loading of Pd combined with ion beams (<500 eV d⁺) for loading and defect engineering; optional laser-driven plasmon excitation. Diagnostics: neutrons, MeV ions, mass spectrometry of low-energy ions, ex-situ tritium and transmutation analysis [V] | **Yes.** Karahadian et al., *Nat. Commun.* 2026: sub-keV yield plateau in Pd and Ti, >10¹⁸ over bare-nucleus rates [V]. No excess heat or transmutation reported (as far as found). |
| 2 | **MIT** (F. Metzler, P. Hagelstein; J. Messinger, N. Galvanetto) | ~$1.5–1.8 M (not verified; inferred from total) | Platform to "thoroughly and reproducibly test claims of nuclear anomalies in gas-loaded metal-hydrogen systems", focused on unambiguous indicators such as neutrons; targets earlier claims of neutrons and low-Z "fission daughters" in gas-loaded / laser-irradiated hydrides [V] | Theory and review papers (NJP 2024; arXiv 2501.08338; SSRN 2023) [V]. **No experimental positive or null result located.** |
| 3 | **Stanford University** | $1.5 M [V] | "Nuclear product detection from deuterated nanoparticles under phonon stimulation": tests whether LENR-active sites in metal nanoparticles form on D₂ exposure [V] | **None located.** |
| 4 | **Energetics Technology Center** (Indian Head, MD) — "CATHODE" | $1.5 M [V] | Pd/D co-deposition (PdCl₂ + LiCl in D₂O) onto a metal film conformed onto a **plastic scintillator**, so the detector sits microns from the active cathode [V] | **None located.** |
| 5 | **University of Michigan** — gas cycling | ~$1.15 M [V] | D₂ gas cycling through a chamber of nanocrystalline Pd; variables: temperature, crystallite size, laser wavelength [V] | **None located.** |
| 6 | **University of Michigan** — diagnostics capability team | ~$0.9–1.1 M [V; sources conflict] | Neutron, gamma and ion detection support for all teams [V] | **None located** (likely internal support). |
| 7 | **Texas Tech University** — materials capability team | $1.15 M [V] | Materials fabrication, characterisation and nuclear-product detection as a shared resource [V]. TTU's Center for Emerging Energy Sciences (R. Duncan) [BK] | **None located.** |
| 8 | **Amphionic LLC** (Dexter, MI) | $295,924 [V] | Nanostructured Pd–aramid-nanofiber (ANF) composites for "controlled LENR" | **None located.** |

**ARPA-E's own conclusions.**
- We found no public program-level final report or summary of conclusions. The FY2023 Annual Report (Aug 2025) lists the topic [V]. We could not read its text.
- With a ~24–30-month performance period, the program should have wrapped up around 2025–26.

**Our interpretation (clearly an inference).** The program was explicitly designed around publication in top journals, and it had dedicated diagnostic teams. Its only prominent paper reports an enhanced-screening-type phenomenon: conventional D–D products at enhanced rates. Nothing published claims heat or new nuclear products. That is weak but real evidence (**grade N, provisional**) that the most-cited "excess heat / transmutation" conditions did not appear on demand under ARPA-E-grade diagnostics. Unpublished null results would be consistent with this, but we cannot confirm them.

---

## 4. Most-credible positive results (detailed)

### 4.1 Sub-keV D–D yield plateau in metal hydrides — UC Davis/LBNL (Nat. Commun. 2026) — Grade B

- **Who / funding.** M. E. Karahadian, M. Colborne, A. Persaud, T. Schenkel, J. N. Munday (UC Davis ECE and LBNL ATAP). Funded by the ARPA-E LENR topic [V]. Munday and Schenkel were both in the 2019 Google consortium.
- **Geometry.** A "dual-chamber platform". A metal foil (Pd or Ti) is the wall between an **electrochemical cell**, which loads D from the back, and a **vacuum chamber**, where a **low-energy D ion beam** strikes the front face [V].
  - We could not retrieve beam current, spot size, foil thickness or detector geometry. [BK: LBNL low-energy ion-beam systems typically deliver µA-level mass-analysed beams, and fusion protons are detected with silicon detectors on the beam side.]
- **Result.** The D–D yield falls with energy as expected (Gamow suppression moderated by screening) down to ~2–2.5 keV. Below that it **plateaus** instead of continuing to fall exponentially. At the lowest energies the yield is **>10¹⁸ above the bare-nucleus expectation** [V].
  - Ti, normally regarded as weakly screening with lower D diffusivity, also shows the plateau [V].
  - The authors frame this as "materials-driven fusion": electrons, defects and locally concentrated deuterium create reaction environments that do not exist in plasmas [V].
- **Scale check (our calculation).** The bare D–D Gamow factor is exp(−2πη) with 2πη ≈ 31.3·√(1.007/E_cm[keV]):
  - ~2×10⁻¹⁴ at E_cm = 1 keV;
  - ~5×10⁻²⁰ at 0.5 keV.

  So a flat yield through this range is equivalent to enhancement factors of 10¹⁰ to 10¹⁸ or more. The absolute rate, however, stays at "countable-with-silicon-detectors" levels, **nowhere near** any energy relevance. Physics World and the press releases stress this [V].
- **Caveats.**
  - A small beam-energy tail (charge exchange in the column, neutral components, or molecular-ion fragments with more energy per nucleon than intended) could mimic a plateau.
  - So could a slowly accumulating D inventory plus a small fraction of full-energy ions.
  - An independent (non-peer-reviewed) "identifiability" critique exists [V].
  - Szczecin's independent observation of a similar plateau in Zr/Ti/Pd (4.2) is the main reason this rates **B** rather than C.

### 4.2 Szczecin programme (Czerski et al.): threshold resonance, e⁺e⁻ channel, thermal-spike fusion — Grade C (plateau B)

- **Facility.** The eLBRUS ultra-high-vacuum accelerator at the University of Szczecin, the coordinator of the EU H2020 CleanHME project (2020–Jan 2025) [V].
- **Beams and detectors.** D⁺ and D₂⁺ at **5–20 kV, 40–60 µA** onto Zr (ZrD₂), Ti and Pd targets. Detectors: large NaI(Tl) and HPGe for gammas, plus silicon charged-particle detectors [V].
- **Threshold resonance (PRC 2022).** Czerski proposed a single-particle **0⁺ resonance in ⁴He at the d+d threshold**. Its e⁺e⁻ (internal pair) width would exceed its proton width, and its interplay with screening explains why metal-target D–D rates at a few keV exceed extrapolations [V].
- **Electron emission (PRC 2024).** The measured electron energy spectrum and the electron/proton branching ratio agree with e⁺e⁻ decay of the 0⁺ state to the ground state, supported by Monte Carlo simulation [V].
- **PRX 2025 (Dubey et al.).** Both **511 keV annihilation** and **bremsstrahlung** were seen, consistent with electrons and positrons of **up to ~23 MeV**. The paper claims that **below ~5 keV the e⁺e⁻ channel dominates D–D fusion** [V].
- **Thermal-spike fusion (arXiv 2409.02112; 2605.27438, May 2026).**
  - A 2H(d,p)3H yield plateau appears in ZrD₂ at beam energies as low as ~1 keV. The protons come from a near-resting centre of mass, which the authors interpret as thermal fusion in beam-induced ion tracks [V].
  - Ti and Pd, which differ in thermal properties and screening, fit the thermal-spike + screening + resonance model [V].
- **Materials papers.** Oxygen and carbon contamination and lattice defects (measured by PAS/XRD) change screening energies (Materials 2023, 2025) [V]. **Surface chemistry and defects are experimental variables, not nuisances.**
- **ICCF-26 claim.** The conference summary reports 511 keV gamma counts **above background in Pd-D and Zr-D samples measured in an underground laboratory, with no beam**, attributed to D diffusion alone [V, from conference summary; D-grade until published].
- **Caveats.**
  - A 23.8 MeV-scale e⁺e⁻ channel dominating at low energy is extraordinary.
  - It has not been seen by gas-target D–D measurements (as far as we know [BK]).
  - It has not been independently confirmed.
  - 511 keV lines are ubiquitous in background: cosmic-ray pair production, and β⁺ activation such as ¹³N/¹⁷F from (d,n) on C/O contaminants, although this is small at keV energies.

### 4.3 Electrochemistry boosts plasma-implantation fusion — UBC "Thunderbird" (Nature 2025) — Grade C

- **Geometry** [V]. Three parts:
  1. a **plasma thruster / plasma-immersion ion implantation (PIII)** source firing D⁺ at the front of a **300 µm Pd foil** in a vacuum chamber;
  2. the vacuum chamber;
  3. an **electrochemical cell on the back face** loading D from heavy water.

  The same foil is both electrochemical cathode and fusion target.
- **Results** [V].
  - D–D fusion neutrons at ~130–140 counts/s typical, peak **188 ± 2 n/s**.
  - Switching electrolysis on raised the fusion rate by **15 ± 2 %** on average.
  - The authors equate ~1 V of electrochemical driving force to ~800 atm of D₂ pressure.
  - Described as the first PIII-driven D–D fusion and the first reproducible demonstration that electrochemical loading increases a nuclear rate.
- **Energetics.** ~10⁻⁹ W of fusion against ~15 W input [V].
- **Why it matters for us.** It is the cleanest template for **"electrochemical knob → nuclear observable"** with an internal control: electrolysis toggled on and off while the beam stays constant. The effect is explained by higher near-surface D density (more targets) rather than new physics.

### 4.4 Tohoku University / Clean Planet Ni–Cu multilayer heat — Grade C (commercial: D)

- **Configuration** [V for the concept; BK for dimensions].
  - Nano-multilayer Ni/Cu (and earlier Pd/Ni with D₂) sputtered onto a Ni substrate. [BK: ~25 × 25 × 0.1 mm, ~6 bilayers of ~2 nm Cu / ~14 nm Ni, often coated on both faces.]
  - The film is **preloaded with H₂ gas**, the chamber is evacuated, and the sample is **heated rapidly**. Heat appears as hydrogen diffuses out through the layered interfaces.
  - Heat bursts can be induced deliberately [V].
- **Numbers** [V].
  - JJAP 2024: maximum energy released per absorbed H **>10 keV**; **no gammas or neutrons** observed.
  - Photon-radiation calorimetry (three detectors, 0.3–5.5 µm): **4–6 W** excess; **460 ± 120 kJ** over 80 h; **≥410 ± 108 keV per H atom**.
- **Hot spots.** ICCF-26 IR imaging reportedly shows excess heat arising in **micron-scale regions** where local temperature exceeds the melting point of Ni [V, conference summary].
- **Related claims.** ³He detected in Ni–Cu/ZrO₂ nanocomposites (arXiv 2409.05382), attributed at ICCF-26 to H+D reactions involving trace D [V]. Grade D.
- **Replication status.** Earlier NEDO-funded "mutual verification" (~2015–2018) among Japanese groups (Tohoku, Kyushu, Nagoya, Kobe, Technova, Nissan) reported AHE in some runs, with limited reproducibility and partners collaborating rather than independent [BK]. **No independent, peer-reviewed replication outside the consortium was found for 2020–2026.**
- **Why it is still the best heat candidate.** The claimed energy per H is **far above chemical (~eV)**. Radiometric and thermocouple calorimetry agree. It is gas-phase, with no electrolyte recombination ambiguity. It is published in a mainstream applied-physics journal. It shows an internal dependence on multilayer structure.
- **Why it is not better than C.**
  - Single ecosystem with a commercial interest.
  - Calorimetry at 500–900 °C in vacuum depends heavily on emissivity and heat-flow modelling. A changing surface emissivity during desorption is an obvious systematic, which the photon-calorimetry paper is designed to address.
  - No nuclear ash has been correlated quantitatively with the heat.
- **Commercial layer** [V, secondary/press].
  - Joint development agreement with **Miura Co.** for an industrial boiler [BK: announced ~2021].
  - "QHe IKAROS" rod module, 120 cm tall, **target** 24 kW, stackable to MW.
  - ~¥500 M raise in January 2026; mass production aimed at 2026–27.
  - No third-party calorimetry of these modules has been published. **Grade D** until it is.

### 4.5 NASA Glenn "Lattice Confinement Fusion" — Grade C

- **Setup** [BK]. ErD₃ and TiD₂ samples irradiated with bremsstrahlung from a ~2.9 MeV electron linac.
- **Observations** [BK]. Neutrons consistent with D–D fusion and some higher-energy neutrons.
- **Proposed mechanism** [BK]. Photodisintegration of D makes energetic neutrons, which "heat" deuterons, which then fuse in a screened lattice; Oppenheimer–Phillips stripping is also proposed.
- **Publications.** PRC 101, 044609/044610 (2020) [BK]. A 2021 ARPA-E workshop presentation covered follow-on **gas-cycling** LCF experiments (no beam) [V].
- **Assessment.** This is essentially beam-driven fusion with a condensed-matter enhancement: plausible, single-group, and not a heat source. We found no peer-reviewed report of a positive gas-cycling-only result.

### 4.6 Lower-credibility nuclear-emission claims (C/D)

- **UIUC glow discharge** (Ziehm & Miley, 2024) [V].
  - Conditions: 10 Torr D₂, 5–7 mm gap, 20–40 mA/cm², cathode about −500 V.
  - CR-39 tracks consistent with **138 ± 21 keV alphas** from the Pd electrodes, ~100× H₂/He controls.
  - Weaknesses: CR-39 track interpretation (sputtered-ion damage, chemical etching artefacts), and alphas at 138 keV have no obvious nuclear source. Grade C/D.
- **Cavitation of TiD powder** (Fomitchev-Zamilov, Sci. Rep. 2024) [V].
  - Neutron counts up to >6500 CPM (>10⁴× background), sustained for hours, only during violent secondary acoustic pressure peaks.
  - The author notes that spallation and other processes still need to be ruled out.
  - Single author and company; acoustic microphonics and EMI in gas-proportional neutron counters are well-known false-positive sources. Grade D.
- **Pd/D co-deposition neutrons** (bubble detectors, J. Electroanal. Chem. 2021) [V]. Bubble detectors are temperature- and shock-sensitive. Grade C/D. This is the lineage behind ETC's ARPA-E CATHODE project.

---

## 5. Most-credible null results and skeptical analyses

### 5.1 Google-funded consortium, *Nature* 2019 (context for 2020–2026) — Grade N

- **Who and effort.** Berlinguette, Chiang, Munday, Schenkel, Fork, Koningstein, Trevithick, et al. About 30 researchers over ~4 years and ~$10 M [BK].
- **Three tracks** [BK]:
  1. Pd electrochemistry at high loading;
  2. Pd powder under hydrogen at elevated temperature/pressure;
  3. low-energy D beams/plasmas on metal hydrides.
- **Result** [BK]. No excess heat and no anomalous nuclear products. Most electrochemical cathodes could not exceed **D/Pd ≈ 0.95**, which is near or below the ~0.9–0.95 thresholds historically claimed necessary for Fleischmann–Pons-type heat.
- **Why it matters.** It is a credible null for "casual" reproduction of Fleischmann–Pons. The consortium's later path (4.1, 4.3) shows where the credible signal actually lives: **enhanced conventional fusion under ion bombardment**.

### 5.2 ARPA-E program (2023–2026) — Grade N (provisional)

See Section 3. Eight projects, including two dedicated diagnostics and materials capability teams, published no positive heat, transmutation or anomalous-radiation results that we could locate. Treat this as an absence of evidence, not a formal null publication.

### 5.3 Mainstream commentary on the positive results

- *Chemistry World* called the Thunderbird enhancement "modest" [V].
- Physics World, phys.org and the UC Davis/LBNL press releases on the 2026 plateau stress that rates remain **far too small for energy** [V].
- The broader mainstream reading is that these results are **solid-state nuclear physics (screening, defects, local density, possibly thermal spikes)**, not vindication of Fleischmann–Pons excess heat.
- Popular Science ("Cold fusion is making a scientific comeback") and Atomic Insights ("How hot is cold fusion?") document renewed but cautious interest [V].
- **Sociology.** Messinger, Metzler and Price's "Gatekeeping: a partial history of cold fusion" (2026) and "Risk and scientific reputation: lessons from cold fusion" (arXiv 2201.03776, 2022) argue that the reputational cost of working on LENR suppressed replication attempts [V]. That is relevant to why nulls are under-published.

### 5.4 Known systematic failure modes that null or down-grade heat claims [BK]

These matter for any heat claim, including Iwamura's:

- **Calibration-constant shift** (Shanahan critique of Storms/McKubre-type calorimetry).
- **D₂/O₂ recombination in open cells** (NRL-era analyses, e.g. Kidwell et al., of gas-loaded Pd–ZrO₂ heat bursts).
- **Emissivity changes** in radiometric/high-temperature calorimetry.
- **Heater-power measurement errors** with pulsed or AC waveforms. This is central to the E-Cat, Parkhomov-type and Mizuno-type disputes.
- **Replications of the Mizuno & Rothwell (2019) Pd-rubbed Ni-mesh claims (hundreds of W)** by independent amateur and small-lab groups in 2019–2022 did not report robust excess heat [BK]. Grade N/D; the setups were not high-quality, so the nulls are weak too.

### 5.5 Status of named programmes and companies (independent verification?)

| Actor | 2020–2026 claim / activity | Independent verification? | Grade |
|---|---|---|---|
| SKINR (U. Missouri; Kimmel gift 2012) | Pd-D electrolysis incl. Energetics "SuperWave" replication; materials work [BK] | No robust published replication; profile declined by early 2020s [BK] | D/N |
| Texas Tech (R. Duncan, CEES) | ARPA-E capability team; Pd-D electrochemistry and materials [V/BK] | No positive peer-reviewed result located | — |
| ENEA (Frascati; Violante lineage) | Pd-cathode metallurgy (grain texture, impurities) controlling reproducibility; collaborations with SRI/Energetics historically [BK] | Historical partial cross-lab agreement (pre-2020); little new 2020–26 output found | C/D |
| Italian groups (Celani Constantan wires; CleanHME/HERMES partners) | Few-W excess in Cu–Ni–Mn wires under H₂/D₂ [BK]; Ni-alloy nanoparticle heat (Biberian, JCMNS 2025) [V] | No independent mainstream replication | D |
| Parkhomov (Russia) | Ni–LiAlH₄ Rossi-type reactors, long-run excess-heat claims [BK] | Calorimetry disputed; some partial repeats by enthusiasts | D |
| Brillouin Energy | "Controlled Electron Capture Reaction", pulsed Ni–H; SRI (Tanzella) measurements reported COP>1 in some runs (pre-2020) [BK] | Only commissioned tests; no peer-reviewed publication | D |
| Industrial Heat | Investor in LENR companies after the Rossi dispute (lawsuit 2016–17) [BK] | Tested Rossi's device and found no excess heat (per its litigation position) [BK] | (N for Rossi) |
| E-Cat (Rossi) | SK/SKL/"SKLep" demonstrations, electricity-output claims [BK] | None credible | D |
| Aureon | Commercial LENR claims (details not verified this session) | None located | D |
| Clean Planet | See 4.4 | Consortium-internal only | C (science) / D (product) |
| NASA GRC | See 4.5 | No | C |
| US Navy / DoD | FY2022 NDAA House report asked for an LENR briefing [BK]. Navy-adjacent work continues via ETC (Indian Head) co-deposition and scintillator cathodes [V]. NRL's earlier gas-loading work leaned toward recombination explanations [BK]. DTRA 2020–26 activity not verified | No positive DoD publication found | — |

---

## 6. Dominant hypotheses (2020–2026) and what each predicts for detection

1. **Enhanced screening + defects + local D density ("materials-driven fusion")** — UC Davis/LBNL, Szczecin, NASA.
   - Predicts conventional D–D products (p+t, n+³He in roughly equal branches) at rates that depend on the material, its loading, vacancies, surface contamination and temperature.
   - Detectable with Si detectors and neutron counters.
   - **Mainstream-compatible.**
2. **⁴He 0⁺ threshold resonance with e⁺e⁻ decay** — Czerski.
   - Predicts that at E ≲ 5 keV most reactions emit e⁺e⁻ pairs (up to ~23–24 MeV total), giving 511 keV coincidences and high-energy bremsstrahlung, with a rising e/p ratio at lower energy.
   - **Need back-to-back 511 keV coincidence and ≥20 MeV-range calorimetric gamma/electron detection.**
3. **Thermal-spike / ion-track fusion** — Czerski 2024–26.
   - Predicts yield plateaus at ~1 keV, proton spectra from near-rest centre of mass, and dependence on the lattice's thermal conductivity and D diffusivity.
4. **Phonon–nuclear coupling / nuclear Dicke model** — Hagelstein, Metzler et al.
   - Predicts D–D energy transferred to the lattice or to host nuclei without energetic particles, possibly fusion–fission of Pd into low-Z products, and a need for strong, coherent lattice excitation.
   - **Hard to test. The best test is a search for low-Z elements with anomalous isotope ratios together with neutron/charged-particle absence.**
5. **Diffusion-flux-triggered reactions** — Iwamura multilayers; Czerski's D-diffusion 511 keV claim.
   - Predicts that heat or radiation correlates with hydrogen **flux** across interfaces (desorption/permeation), not static loading.
   - This aligns with "permeation" geometries: foil membranes with pressure or electrochemical gradients.
6. **Legacy models** (Widom–Larsen ULM neutrons, Takahashi TSC, Storms' "Hydroton" in cracks, Staker's superabundant-vacancy δ-phase) remain in the LENR literature. No 2020–2026 mainstream support was found. Widom–Larsen in particular is widely rejected [BK].

---

## 7. What the 2020–2026 literature implies for the best configuration to try now

**Guiding principle.** Every credible signal from 2020 to 2026 came from **detection-first experiments that measure a nuclear observable directly, with an internal on/off control, in a thin-foil metal hydride under ion/plasma bombardment, with electrochemical or gas loading as an enhancement knob.** Heat-only claims still sit at C/D. A project optimising for *credible and measurable* should therefore build a **nuclear-product experiment first**. Heat is at most a secondary channel.

### 7.1 Recommended primary configuration: "dual-chamber permeation-foil target"

**Geometry.**
- A circular **Pd foil**, ~25–300 µm thick and ~10–25 mm diameter, forms the pressure and vacuum wall between two chambers.
- **Back chamber:** electrochemical cell (e.g. 0.1 M LiOD in D₂O, Pt anode; current density a few to tens of mA/cm²). Alternatively a D₂ gas side (1–10 bar) for a gas-permeation variant.
  - Thinner foils give faster back-to-front D flux, and flux may matter (hypothesis 5).
  - Thunderbird used 300 µm.
- **Front chamber:** UHV (≤10⁻⁷ Torr). Szczecin shows that O/C surface contamination changes screening, so bake-out and in-situ sputter-cleaning with Ar⁺ are recommended.
  - Beam: a **mass-analysed** low-energy D⁺ beam, **0.3–20 keV**, 10–100 µA on ~0.2–1 cm². Alternatively a PIII / glow-discharge plasma as in Thunderbird (simpler, but its energy distribution is less well defined).
- **Material matrix** (same geometry each time):
  - Pd (strong screening, high D mobility);
  - **Ti** (weak screening control; UCD/LBNL still see a plateau);
  - **Zr** (Szczecin reference);
  - **Au or Cu** (non-hydriding null control).

**Detectors** (arranged around the front face):
1. **Charged particles.** Si PIPS/surface-barrier detectors, or ΔE–E telescopes, at backward angles (~120–150°) a few cm from the spot.
   - Thin Al or Mylar windows (~1–10 µm, sized per product) stop scattered D and light. Windows must still pass 3.02 MeV p, 1.01 MeV t and 0.82 MeV ³He.
   - Reaction depth at keV energies is ≲100 nm, so charged products escape the beam face. **Rear-face detection through the foil is impossible**: a 3 MeV proton has a range of tens of µm in Pd.
2. **Neutrons.** Moderated ³He counters plus EJ-309/stilbene with pulse-shape discrimination, placed close (≤10 cm) for solid angle, with timestamps. Thunderbird counted ~10² n/s. Plateau-regime rates will be far lower.
3. **511 keV and bremsstrahlung.**
   - Two LaBr₃/NaI or HPGe detectors in **180° coincidence** for positron annihilation.
   - A large NaI/BGO for the 5–25 MeV range (e⁺e⁻ channel test).
   - Pb shielding plus a cosmic veto. An underground site if available.
   - Calibrate with ²²Na and check β⁺ activation of C/O contaminants.
4. **Loading diagnostics.** Foil resistance ratio (R/R₀ versus D/Pd), electrochemical charge counting, and gas-side pressure.
5. **Post-run analysis.** Tritium (liquid scintillation of electrolyte and foil), ³He/⁴He mass spectrometry, and SIMS/ICP-MS for low-Z elements and isotope ratios. This tests hypothesis 4 cheaply.

**Protocol.**
- **Interleave** the variables on minute-to-hour cadences with the beam held constant: electrolysis on/off (the Thunderbird method), D versus H (D₂O versus H₂O electrolyte; H⁺ beam), beam energy steps, and foil temperature (20–200 °C).
- Measure a **yield-versus-energy curve from 20 keV down to ≤0.5 keV** for loaded versus unloaded and Pd versus Ti.
- Target observables:
  - (a) screening energy U_e;
  - (b) presence and height of a plateau;
  - (c) e⁺e⁻/proton ratio;
  - (d) dependence on **back-side D flux**, which is the genuinely new variable that could separate "LENR-like" diffusion effects from plain screening.
- **Pre-register analyses** and blind the on/off labels. Run two independent detector chains per product.

**Why this geometry maximises credibility.**
- It reproduces two results that have peer-reviewed support in *Nature*, *Nature Communications* and PRX (B/C grade), which calibrates the apparatus.
- It adds variables (loading flux, material, surface state) where the dominant hypotheses make **different, falsifiable predictions**.
- Any anomaly, such as a rate increase correlated with permeation flux beyond what D density explains, or an e⁺e⁻ excess, is measured in counts with controls, not inferred from calorimetry.

### 7.2 Secondary (lower-credibility) configuration: gas-phase Ni–Cu multilayer heat

- **When to use it.** If the project also wants a heat channel.
- **Build.** Replicate the Tohoku/Clean Planet approach: sputtered Ni/Cu nano-multilayers on a thin Ni substrate [BK: ~25 × 25 × 0.1 mm, ~6 × (2 nm Cu / 14 nm Ni)], loaded with H₂ (and separately D₂) at moderate temperature, then evacuated and heated rapidly to several hundred °C.
- **Measure** with **two independent calorimetry methods**: radiometric plus thermocouple/heat-flow. Use an identical blank Ni substrate and an Ar/He-loaded dummy, and measure emissivity before and after.
- **Instrument for nuclear products at the same time** (neutron, gamma, ³He/⁴He sampling) so that any heat can be tied to ash.
- **Credibility threshold.** Reproduce ≥ several W and ≥10 keV/H with blanks at zero.
- **Expected grade.** C at best until it is independently replicated.

### 7.3 Things the 2020–2026 literature says *not* to prioritise

- **Classic open-cell Pd–D₂O electrolysis for heat alone.** Loading is hard to reach and control (Google 2019 null), and recombination and calibration systematics dominate. Electrolysis is valuable **as a loading knob for a beam target** (7.1), not as a stand-alone heat source.
- **CR-39-only, bubble-detector-only or single-counter neutron claims.** These carry the highest false-positive history. Always use at least two independent detection physics.
- **Nickel–LiAlH₄ "Rossi-type" reactors and commercial black boxes.** No independent verification exists.
- **Cavitation or sonofusion neutron claims.** Single-author claims have a history of artefacts; only pursue them with coincidence-gated, pulse-shape-discriminated neutron detection.

### 7.4 Quantitative expectations to set before building

- **Beam-target regime, 5–20 keV, µA–mA currents.** D–D rates of ~10–10³ events/s are achievable (Thunderbird ~10² n/s). This is enough to measure a 15 % modulation in minutes to hours.
- **Sub-keV plateau regime.** Rates are much lower, so detector efficiency and background (underground or shielded, with coincidence) dominate the design.
- **Energy.** Every credible nuclear result so far gives ≲10⁻⁸ of input power. **No credible configuration from 2020–2026 approaches net energy.** Treat "measurable nuclear signal" and "useful heat" as different goals.

---

## 8. References (with URLs)

**Beam/plasma + solid-state enhanced fusion (most credible)**
1. Karahadian M.E., Colborne M., Persaud A., Schenkel T., Munday J.N., "Enhanced nuclear fusion in the sub-keV energy regime", *Nat. Commun.* 17, 8845 (2026). https://www.nature.com/articles/s41467-026-74421-1 ; arXiv: https://arxiv.org/abs/2512.06212
2. LBNL news, "When it comes to fusion, materials matter" (23 Jul 2026). https://newscenter.lbl.gov/2026/07/23/when-it-comes-to-fusion-materials-matter/ ; UC Davis: https://www.ucdavis.edu/blog/when-it-comes-fusion-materials-matter ; phys.org: https://phys.org/news/2026-07-materials-fusion-reaction.html ; Physics World: https://physicsworld.com/a/nuclear-fusion-persists-at-ultralow-energies-inside-metal-foils/
3. "Electrochemical loading enhances deuterium fusion rates in a metal target", *Nature* (Aug 2025), Berlinguette group, UBC. https://www.nature.com/articles/s41586-025-09042-7 ; eScholarship PDF: https://escholarship.org/content/qt9r65z0pt/qt9r65z0pt.pdf ; group page: https://groups.chem.ubc.ca/cberling/thunderbird-reactor/
4. Chemistry World, "Electrochemistry offers 'modest' boost to deuterium fusion reaction". https://www.chemistryworld.com/news/electrochemistry-offers-modest-boost-to-deuterium-fusion-reaction/4022070.article ; phys.org: https://phys.org/news/2025-08-room-temperature-reactor-electrochemistry-boost.html
5. Dubey R., Czerski K., Das H. G., et al., "Experimental signatures of a new channel of the deuteron–deuteron reaction at very low energy", *Phys. Rev. X* 15, 041004 (2025). https://link.aps.org/doi/10.1103/chlp-b215 ; arXiv: https://arxiv.org/abs/2408.07567
6. Czerski K., "Deuteron-deuteron nuclear reactions at extremely low energies", *Phys. Rev. C* 106, L011601 (2022). https://link.aps.org/doi/10.1103/PhysRevC.106.L011601
7. Czerski K. et al., "Indications of electron emission from the deuteron-deuteron threshold resonance", *Phys. Rev. C* 109, L021601 (2024). https://link.aps.org/doi/10.1103/PhysRevC.109.L021601 ; arXiv: https://arxiv.org/abs/2305.17101
8. Czerski K., Dubey R., Kowalska A., Das H. G., Kaczmarski M., Targosz-Ślęczka N., Valat M., "Observation of thermal deuteron-deuteron fusion in ion tracks", arXiv:2409.02112. https://arxiv.org/abs/2409.02112
9. Czerski K., Dubey R., Das H. G., Thulichery S., Kowalska A., Targosz-Ślęczka N., Valat M., "Thermal deuteron-deuteron fusion in metallic targets", arXiv:2605.27438 (May 2026). https://arxiv.org/abs/2605.27438
10. Szczecin group, "Electron screening in deuteron–deuteron reactions on a Zr target with oxygen and carbon contamination", *Materials* 18, 1331 (2025). https://doi.org/10.3390/ma18061331
11. Szczecin group, "Crystal lattice defects in deuterated Zr in presence of O and C impurities studied by PAS and XRD for electron screening effect", *Materials* 16, 6255 (2023). https://doi.org/10.3390/ma16186255
12. "Beyond scalar screening: transfer-function closure and identifiability of the sub-keV deuterium-fusion plateau" (ResearchGate preprint, 2026). https://www.researchgate.net/publication/412147776

**Heat claims (Japan) and commercial**
13. Iwamura Y., Itoh T., Yamauchi S., Takahashi T., "Anomalous heat generation that cannot be explained by known chemical reactions produced by nano-structured multilayer metal composites and hydrogen gas", *Jpn. J. Appl. Phys.* 63 (2024). https://iopscience.iop.org/article/10.35848/1347-4065/ad2622
14. Kasagi J., Itoh T., Shibasaki Y., Iwamura Y., "Photon radiation calorimetry for anomalous heat generation in NiCu multilayer thin film during hydrogen gas desorption", arXiv:2311.18347 / JCMNS. https://arxiv.org/abs/2311.18347 ; https://jcmns.org/article/134004-photon-radiation-calorimetry-for-anomalous-heat-generation-in-nicu-multilayer-thin-film-during-hydrogen-gas-desorption
15. Kasagi J. et al., "Measurement of radiant spectrum for excess heat generation in NiCu and Ni thin film during hydrogen gas desorption", arXiv:2509.13847. https://arxiv.org/pdf/2509.13847
16. "Detections of He-3 in Ni-based binary metal nanocomposites with Cu in zirconia exposed to hydrogen gas at elevated temperatures", arXiv:2409.05382. https://arxiv.org/pdf/2409.05382
17. Frontiers in Materials (2024), "Changes in the structure and composition of nano-sized Ni-Cu multilayer films with increasing temperature in an atmosphere with and without hydrogen". https://www.frontiersin.org/journals/materials/articles/10.3389/fmats.2024.1407810/full
18. Clean Planet technology page. https://www.cleanplanet.co.jp/en/technology/ ; Miura JDA: https://fuelcellsworks.com/subscribers/japan-miura-co-ltd-clean-planet-conclude-joint-development-agreement-for-industrial-boiler-using-quantum-hydrogen-energy ; commercialization coverage (secondary): https://newfireenergy.substack.com/p/clean-planet-moves-deeper-into-commercialization ; https://www.lenrbuyerguide.com/posts/2026-05-14-clean-planet-s-qhe-lenr-nears-commercial-reality-in-2026

**ARPA-E program**
19. ARPA-E press release (Feb 2023). https://arpa-e.energy.gov/news-and-events/news-and-insights/us-department-energy-announces-10-million-funding-projects-studying-low-energy-nuclear-reactions
20. Green Car Congress summary of the 8 projects. https://www.greencarcongress.com/2023/02/20230218-lenr.html ; ANS Nuclear Newswire: https://www.ans.org/news/article-4769/arpae-picks-eight-teams-to-proveor-debunklowenergy-nuclear-reactions/ ; ExecutiveGov: https://executivegov.com/2023/02/arpa-e-funds-8-projects-to-break-impasse-in-low-energy-nuclear-reactions-field/ ; New Energy Times: https://news.newenergytimes.net/2023/03/01/u-s-department-of-energys-arpa-e-funds-low-energy-nuclear-reaction-research/
21. ARPA-E project pages: CATHODE (ETC) https://arpa-e.energy.gov/programs-and-initiatives/search-all-projects/cathode-cathode-scintillator-detector-electrochemistry ; Stanford https://arpa-e.energy.gov/programs-and-initiatives/search-all-projects/nuclear-product-detection-deuterated-nanoparticles-under-phonon-stimulation ; Texas Tech https://arpa-e.energy.gov/programs-and-initiatives/search-all-projects/advanced-materials-characterization-and-nuclear-product-detection-lenr ; Amphionic https://arpa-e.energy.gov/programs-and-initiatives/search-all-projects/nanostructured-pd-anf-composites-controlled-lenr-exploitation
22. LBNL ATAP, "Berkeley Lab to lead ARPA-E low energy nuclear reactions project". https://atap.lbl.gov/news/berkeley-lab-to-lead-arpa-e-low-energy-nuclear-reactions-project/ ; LBNL IBT project page: https://www-ibt.lbl.gov/projects/low-energy-nuclear-reaction
23. ARPA-E FY2023 Annual Report to Congress (Aug 2025). https://arpa-e.energy.gov/sites/default/files/2025-09/ARPA-E%20FY%202023%20Annual%20Report.pdf
24. Benyo T. et al., "Lattice confinement fusion gas cycling experiments", ARPA-E LENR workshop (2021). https://arpa-e.energy.gov/sites/default/files/migrated/2021LENR_workshop_Benyo.pdf

**Theory, reviews, history**
25. Metzler F., Hunt C., Hagelstein P.L., Galvanetto N., "Known mechanisms that increase nuclear fusion rates in the solid state", *New J. Phys.* 26 (2024). https://iopscience.iop.org/article/10.1088/1367-2630/ad091c
26. Hagelstein P.L., Metzler F., Lilley M.K., Messinger J.F., Galvanetto N., "Models for nuclear fusion in the solid state", arXiv:2501.08338. https://arxiv.org/abs/2501.08338
27. Metzler F., Hunt C., Messinger J., Galvanetto N., Hagelstein P., "Probing neutrons and purported fission daughter products from gas-loaded, laser-irradiated metal-hydrogen targets", SSRN (2023). https://dx.doi.org/10.2139/ssrn.4411160
28. Messinger J.F., Metzler F., Price H., "Gatekeeping: a partial history of cold fusion", arXiv:2601.09996 (2026). https://arxiv.org/abs/2601.09996
29. "Risk and scientific reputation: lessons from cold fusion", arXiv:2201.03776. https://arxiv.org/pdf/2201.03776
30. "A review of experiments reporting non-conventional phenomena in nuclear matter…", arXiv:2308.13533. https://arxiv.org/pdf/2308.13533

**Single-lab nuclear-emission claims**
31. Ziehm E.P., Miley G.H., "On the detection of alpha emission from a low-voltage DC deuterium discharge with palladium electrodes", arXiv:2402.05117. https://arxiv.org/abs/2402.05117
32. Fomitchev-Zamilov M., "Observation of neutron emission during acoustic cavitation of deuterated titanium powder", *Sci. Rep.* (2024); correction https://www.nature.com/articles/s41598-024-67868-z . https://www.nature.com/articles/s41598-024-62055-6
33. "An experimental study on deuterium production from titanium hydride powders subjected to thermal cycles", *Symmetry* 16, 1542 (2024). https://www.mdpi.com/2073-8994/16/11/1542
34. "Electrolytic co-deposition neutron production measured by bubble detectors", *J. Electroanal. Chem.* (2021). https://www.sciencedirect.com/science/article/abs/pii/S1572665721000503
35. Pines V. et al.; Steinetz B. et al., NASA GRC lattice confinement fusion, *Phys. Rev. C* 101, 044609/044610 (2020) [BK]. https://journals.aps.org/prc/abstract/10.1103/PhysRevC.101.044610

**Conferences, EU projects, press**
36. ICCF-26 (Morioka, 26–30 May 2025). https://iccf26.org/ ; abstracts: https://iccf26.org/img_common/abstract/iccf26_abstract.pdf ; summary: https://www.sciengine.com/JMCC/doi/10.16084/j.issn1001-3555.2025.06.009 ; photos: https://solidstatefusion.org/2025/06/event-gallery-iccf-26/
37. ICCF-25 book of abstracts (Szczecin, 2023). https://newenergytimes.com/v2/conferences/2023/ICCF25/ICCF-25-Book-of-Abstracts-2023.07.04.pdf
38. CleanHME (H2020 951974). https://cordis.europa.eu/project/id/951974/results ; https://cleanhme.eu/ ; publications: https://cleanhme.eu/?page_id=27 ; HERMES (H2020 952184): https://cordis.europa.eu/project/id/952184
39. Popular Science, "Cold fusion is making a scientific comeback". https://www.popsci.com/science/cold-fusion-low-energy-nuclear-reaction/ ; Atomic Insights, "How hot is cold fusion?": https://atomicinsights.com/how-hot-is-cold-fusion/
40. Conscious Energy, "Mapping global LENR funding" (Feb 2025). https://conscious.energy/2025/02/09/mapping-global-lenr-funding/
41. Berlinguette C.P. et al., "Revisiting the cold case of cold fusion", *Nature* 570, 45–51 (2019) [BK]. https://www.nature.com/articles/s41586-019-1256-6

---

*Items to verify first (highest decision value):*
1. Beam energies, currents, detector geometry and absolute rates in Karahadian et al. 2026.
2. The Thunderbird PIII voltage and neutron-detector type.
3. The Czerski e/p branching-ratio values and whether any independent group has tested the 0⁺ resonance.
4. Any ARPA-E final program summary or team publications from late 2025–2026.
5. The exact dimensions of the Tohoku/Clean Planet multilayer samples.
6. The exact FY2022 NDAA report language on LENR.
