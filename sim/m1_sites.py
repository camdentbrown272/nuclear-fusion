"""
M1 Part A - site catalogue and per-site D-D rates.

For every site class: minimum D-D separation, local electron density, vibrational
quantum, binding energy, achievable number density, and U_e,eff from
  (1) Thomas-Fermi (Yukawa) and full Lindhard (RPA) static screening at the
      local electron density, solved by exact WKB at the site's relative energy;
  (2) the empirical model: accelerator-measured U_e of the host (R3), scaled to the
      site by (n_local/n_bulk)^1/2 (the Debye-type scaling used to fit the
      accelerator data), times an optional defect factor, and mapped to thermal
      energy by the M0 Yukawa-WKB method (m0.wkb_yukawa).
Also: the bulk-null constraint on each class.

Run: python3 sim/m1_sites.py  -> docs/models/figs/m1_sites.txt, m1_sites.png
"""
import os
import warnings
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import m1_physics as P

warnings.filterwarnings("ignore")
kT300 = P.kB * 300

# Accelerator screening energies used by the empirical model (eV).
# "lo" = lowest credible lab value, "hi" = highest (R3 Table 2.1 / R4 Table 2.7).
U_ACC = {
    "Pd":  (310, 800),    # Tohoku (Kasagi 2002) / Bochum (Raiola 2004)          [R3, m]
    "PdO": (600, 600),    # Tohoku PdO; Lipson Au/Pd/PdO 601+-23                  [R3, R4]
    "Ti":  (30, 300),     # Bochum hydride at 20C (<~30) / JINR-Tomsk TiD2 100-300 [R3, m]
    "Zr":  (105, 340),    # Szczecin with resonance / ZrD2 thick-target 2024     [R3, v]
    "Ni":  (150, 380),    # lo = Pd-like lab ratio 310/800 applied [guess]; hi Bochum [R3, m]
    "Cu":  (180, 470),    # same scaling [guess] / Bochum                         [R3, m]
    "LiD": (30, 60),      # insulators "small effects" (Bochum)                  [R3, v]
    "oxide_other": (310, 600),  # CaO/Pd, ZrO2/Pd: no data; bracket Pd..PdO      [guess]
}
N_BULK = {"Pd": 0.020, "PdO": 0.010, "Ti": 0.030, "Zr": 0.028, "Ni": 0.030, "Cu": 0.025,
          "LiD": 0.001, "oxide_other": 0.010}

# ---------------------------------------------------------------------------
# Site catalogue.  Units: d_min (A); n_loc (e/bohr^3, all valence electrons at the
# D position); hw (eV, D vibrational quantum); E_b (text, eV); dens_typ / dens_hi
# (close D-D pairs per cm^3 of host material); acc = key into U_ACC; fdef = defect
# multiplier on U_acc for the "defects raise U_e" claim (Szczecin, grade C);
# mol = molecular D2 pair (Koonin-Nauenberg calibration applies).
# ---------------------------------------------------------------------------
SITES = [
    dict(k="bulkO", name="Pd bulk octahedral (beta-PdD, x~0.9)", host="Pd", d_min=2.85,
         n_loc=0.020, hw=0.037, E_b="-0.20 per D vs 1/2 D2 (beta formation enthalpy)",
         dens_typ=6 * 6.1e22, dens_hi=6 * 6.8e22, acc="Pd", fdef=1.0, mol=False,
         ctrl="loading x (electrochemical overpotential / gas pressure)",
         src="a=4.03-4.07 A; hw 35-40 meV INS (Rowe et al. PRL 29, 1250, 1972) [BK]; "
             "n_loc DFT/EMT interstitial density ~0.01-0.03 a.u. (Puska & Nieminen PRB 29, 5382, 1984) [BK]; "
             "dH ~ -19 kJ/mol H (Fukai, The Metal-Hydrogen System, 2005) [BK]; 6 NN pairs per D"),
    dict(k="bulkT", name="Pd tetrahedral (T-site D next to O-site D)", host="Pd", d_min=1.75,
         n_loc=0.025, hw=0.070, E_b="+0.1 to +0.2 above O-site (DFT)",
         dens_typ=4 * 6e19, dens_hi=4 * 6e20, acc="Pd", fdef=1.0, mol=False,
         ctrl="thermal occupancy exp(-0.1 eV/kT): temperature, x->1",
         src="O-T distance a*sqrt(3)/4 = 1.75 A [calc]; T-O energy 0.1-0.2 eV (DFT, e.g. Kamakoti & Sholl "
             "J. Membr. Sci. 225, 145, 2003) [BK]; occupancy 1e-3..1e-2 of D [calc from Boltzmann]"),
    dict(k="vacD6", name="Vacancy-D_n cluster (SAV, VD6)", host="Pd", d_min=2.55,
         n_loc=0.012, hw=0.050, E_b="0.23 per D trap energy (H-vacancy, Pd)",
         dens_typ=12 * 6.8e19, dens_hi=12 * 6.8e21, acc="Pd", fdef=1.3, mol=False,
         ctrl="cathodic charging (1e-3 near surface), codeposition (1e-2), GPa anneal (0.1-0.25)",
         src="D-D in VH6 ~2.5 A (R2; DFT Nazarov et al. PRB 89, 144108, 2014; Isaeva et al. IJHE 36, 1254, 2011) [BK]; "
             "trap energy 0.23 eV (Besenbacher et al. J. Less-Common Met. 130, 475, 1987) [BK]; "
             "SAV fractions Fukai & Okuma PRL 73, 1640 (1994) [R2]; 12 NN D-D pairs per VD6"),
    dict(k="vacD2", name="D2 molecule inside monovacancy (Hagelstein conjecture; DFT-disfavoured)",
         host="Pd", d_min=0.80, n_loc=0.004, hw=0.371, E_b="unbound vs 2 isolated trapped D in DFT",
         dens_typ=6.8e19, dens_hi=6.8e21, acc=None, fdef=1.0, mol=True,
         ctrl="same as SAV; exists only if the conjecture holds",
         src="Hagelstein 2009 (R4 ref); DFT finds H2 not stable in Pd monovacancy (Isaeva 2011) [BK]"),
    dict(k="void", name="Divacancy/void holding D2 fluid (~1 GPa)", host="Pd", d_min=0.74,
         n_loc=0.0, hw=0.371, E_b="D2 bond 4.56 eV; void pressure = D fugacity",
         dens_typ=4.3e22 * 1e-5, dens_hi=4.3e22 * 1e-3, acc=None, fdef=1.0, mol=True,
         ctrl="void volume fraction 1e-5 (cycled) .. 1e-3 (blistered); fugacity via overpotential",
         src="D2 fluid at 1 GPa: molar volume ~14 cm3/mol -> 4.3e22 D2/cm3 (H2 EOS, Loubeyre 1996) [BK]; "
             "bond length unchanged <1% at 1 GPa [BK]; molecular rate: Koonin & Nauenberg 1989 via M0"),
    dict(k="disl", name="Dislocation core (edge, tensile side)", host="Pd", d_min=2.80,
         n_loc=0.015, hw=0.035, E_b="0.1-0.6 trap energy",
         dens_typ=6 * 1.1e19, dens_hi=6 * 1.1e20, acc="Pd", fdef=1.3, mol=False,
         ctrl="rho = 1e11 (hydride-cycled) .. 1e12 cm^-2 (heavily cycled/cold-worked)",
         src="trap energies Kirchheim, Acta Metall. 29, 835 (1981) [BK]; rho after alpha/beta cycling "
             "1e11-1e12 cm^-2 (Flanagan & Oates, Annu. Rev. Mater. Sci. 21, 269, 1991) [BK]; "
             "~3 core sites per b (2.75 A) -> 1.1e8 sites/cm of line"),
    dict(k="gb", name="Grain boundary", host="Pd", d_min=2.70,
         n_loc=0.013, hw=0.035, E_b="0.1-0.3 segregation energy",
         dens_typ=6 * 3.0e18, dens_hi=6 * 3.0e21, acc="Pd", fdef=1.3, mol=False,
         ctrl="grain size: 20 um (annealed foil) .. 20 nm (nanocrystalline/codeposit); S_v=3/d",
         src="GB sites 2 layers x 1.5e15 cm^-2 x S_v [calc]; segregation (Mütschele & Kirchheim, "
             "Scripta Metall. 21, 135, 1987) [BK]"),
    dict(k="surf", name="Free surface (fcc hollow, adsorbed D)", host="Pd", d_min=2.75,
         n_loc=0.008, hw=0.070, E_b="0.45-0.5 per D vs 1/2 D2 (chemisorption)",
         dens_typ=3 * 1.5e15, dens_hi=3 * 1.5e15, acc="Pd", fdef=1.0, mol=False, areal=True,
         ctrl="real area (roughness factor 1..100); coverage theta ~1 at high fugacity",
         src="Pd(111) 1.53e15 cm^-2; adsorption Christmann, Surf. Sci. Rep. 9, 1 (1988) [BK]; "
             "3 NN pairs per adsorbed D; hw(perp) ~ 70-100 meV for D [BK]"),
    dict(k="subsurf", name="Subsurface octahedral layer", host="Pd", d_min=2.85,
         n_loc=0.018, hw=0.040, E_b="~ -0.1 to -0.2",
         dens_typ=6 * 1.5e15, dens_hi=6 * 1.5e15, acc="Pd", fdef=1.0, mol=False, areal=True,
         ctrl="as surface", src="one (111) layer below surface [calc]"),
    dict(k="crack", name="Crack face / nanogap, w = 0.3-1 nm (Storms NAE)", host="Pd", d_min=1.2,
         n_loc=0.005, hw=0.070, E_b="as surface",
         dens_typ=1.5e15 * 2 * 10 * 0.1, dens_hi=1.5e15 * 2 * 1e5 * 0.1, acc="Pd", fdef=1.0, mol=False,
         ctrl="crack area S_v: 10 cm^-1 (cycled foil) .. 1e5 cm^-1 (codeposit/dealloyed); 10% of "
              "gap area in the 0.3-1 nm window",
         src="cross-gap D-D = w - 2x0.9 A adsorption height -> 1.2 A at w=0.3 nm [calc]; "
             "S_v and the 10% window fraction are [guess]; gap electron density from spill-out overlap [guess]"),
    dict(k="PdO", name="PdO/Pd interface (D-rich Pd side)", host="Pd", d_min=2.85,
         n_loc=0.010, hw=0.040, E_b="interface accumulation (Lipson)",
         dens_typ=6 * 1.3e15, dens_hi=6 * 1.3e15, acc="PdO", fdef=1.0, mol=False, areal=True,
         ctrl="thermal/anodic oxide 2-20 nm; area = active face",
         src="PdO U_e ~600 eV at keV (Kasagi 2002; Lipson 601+-23) [R3,R4]; areal site density [calc]"),
    dict(k="oxint", name="CaO/Pd or ZrO2/Pd interface", host="Pd", d_min=2.85,
         n_loc=0.010, hw=0.040, E_b="unknown",
         dens_typ=6 * 1.3e15, dens_hi=6 * 1.3e15, acc="oxide_other", fdef=1.0, mol=False, areal=True,
         ctrl="multilayer count (Iwamura: 5 CaO/Pd bilayers); NP/ZrO2 contact area",
         src="no screening data; bracketed Pd..PdO [guess]"),
    dict(k="NP", name="Pd nanoparticle 2-10 nm (core)", host="Pd", d_min=2.85,
         n_loc=0.018, hw=0.038, E_b="weaker than bulk; x=0.3-0.5 at 1 bar",
         dens_typ=6 * 2.7e22, dens_hi=6 * 3.4e22, acc="Pd", fdef=1.0, mol=False,
         ctrl="particle size; support (ZrO2) interface counted separately",
         src="size-dependent capacity (R2 item 1; Pundt & Kirchheim Annu. Rev. Mater. Res. 36, 555, 2006) [BK]"),
    dict(k="TiD2", name="TiD2 tetrahedral (fluorite)", host="Ti", d_min=2.22,
         n_loc=0.030, hw=0.100, E_b="-0.6 per D (TiH2 formation)",
         dens_typ=3 * 1.1e23, dens_hi=3 * 1.1e23, acc="Ti", fdef=1.0, mol=False,
         ctrl="gas loading 400-600 C; cracks on hydriding",
         src="a(TiH2)=4.45 A -> D-D a/2 = 2.22 A; H optical ~140 meV -> D ~100 meV [BK]"),
    dict(k="ZrD2", name="ZrD2 tetrahedral", host="Zr", d_min=2.40,
         n_loc=0.028, hw=0.095, E_b="-0.8 per D",
         dens_typ=3 * 8.6e22, dens_hi=3 * 8.6e22, acc="Zr", fdef=1.0, mol=False,
         ctrl="as Ti", src="epsilon-ZrH2 lattice [BK]; U_e Szczecin (R3)"),
    dict(k="NiD", name="Ni(D) octahedral / Ni-Cu interface trap", host="Ni", d_min=2.64,
         n_loc=0.030, hw=0.062, E_b="+0.1 per D (endothermic solution)",
         dens_typ=6 * 1e18, dens_hi=6 * 1e20, acc="Ni", fdef=1.0, mol=False,
         ctrl="solubility ~1e-4 at 1 bar/300 C; interfaces trap more",
         src="NiH a=3.73 A -> 2.64 A; H optical 88 meV -> D 62 meV [BK]"),
    dict(k="LiD", name="LiD (ionic insulator, D- anion)", host="LiD", d_min=2.88,
         n_loc=0.001, hw=0.060, E_b="-0.95 per D (formation)",
         dens_typ=6 * 5.5e22, dens_hi=6 * 5.5e22, acc="LiD", fdef=1.0, mol=False,
         ctrl="fracture only (fracto-emission host)",
         src="rock-salt a=4.07 A [BK]; TF n/a (insulator) - n_loc is a nominal placeholder"),
]
SITE = {s["k"]: s for s in SITES}

KN_UEFF, KN_RHO0 = 34.0, 1.5e26          # M0 D2 calibration (Koonin & Nauenberg 1989)
NULL_PER_D = 1e-25                       # fusions per D per s, bulk Pd nulls [M0; BK Gai et al. Nature 340, 29 (1989)]

_lind_cache = {}


def E_rel(hw, T=300):
    """Mean energy of the relative-coordinate oscillator mode, eV."""
    return 0.5 * hw / np.tanh(hw / (2 * P.kB * T))


def site_models(s):
    """Return dict of U_eff (eV) and rates per pair (s^-1) for all models."""
    out = {}
    if s["mol"]:
        for m in ("TF", "Lind", "Elo", "Ehi"):
            out["U_" + m] = KN_UEFF
        rho0 = KN_RHO0
    else:
        E = E_rel(s["hw"])
        n = s["n_loc"]
        out["U_TF"] = P.Ueff_yukawa(P.U_TF(n), E)
        if n not in _lind_cache:
            _lind_cache[n] = P.LindhardPotential(n)
        L = _lind_cache[n]
        out["U_Lind"] = P.Ueff_from_P(P.wkb_penetration(L.V, E)[0])
        lo, hi = U_ACC[s["acc"]]
        nb = N_BULK[s["acc"]] if s["acc"] in ("PdO", "oxide_other") else N_BULK[s["host"]]
        scale = np.sqrt(n / nb) if s["acc"] not in ("PdO", "oxide_other") else 1.0
        out["Uacc_lo"] = lo * scale
        out["Uacc_hi"] = hi * scale * s["fdef"]
        out["U_Elo"] = P.Ueff_yukawa(out["Uacc_lo"], E)
        out["U_Ehi"] = P.Ueff_yukawa(out["Uacc_hi"], E)
        rho0 = P.rho0_zero_point(s["hw"])
        # informational: harmonic confinement penalty to reach the TF turning point
        _, rtp = P.m0.wkb_yukawa(P.U_TF(n), E=E)
        s1 = P.hbarc / np.sqrt(2 * P.m_d_c2 * s["hw"])
        out["conf_TF"] = np.exp(-max(s["d_min"] - rtp, 0) ** 2 / (2 * 2 * s1 ** 2))
    out["rho0"] = rho0
    for m in ("TF", "Lind", "Elo", "Ehi"):
        out["lam_" + m] = P.rate_pair(out["U_" + m], rho0)
    return out


def null_limit_per_pair(s, ref_pairs_per_D=None):
    """
    Upper limit on the per-pair rate of class s implied by 1989-90 bulk Pd nulls
    (<= 1e-25 fusions/D/s), if class s was present at its typical density in the
    null samples (cold-worked electrolytic Pd, ~0.5 cm3, ~10 cm2 surface per cm3).
    Returns None for classes not present in those samples.
    """
    present = {"bulkO", "bulkT", "disl", "gb", "surf", "subsurf", "void", "crack"}
    if s["k"] not in present:
        return None
    nD = 6.1e22
    dens = s["dens_typ"] * (10.0 if s.get("areal") else 1.0)   # 10 cm2 of surface per cm3
    if s["k"] == "vacD6":
        dens = 12 * 6.8e17    # thermal + cold-work vacancies ~1e-5
    return NULL_PER_D * nD / dens


def main():
    lines = []
    Pr = lambda s="": (print(s), lines.append(s))
    Pr("M1 Part A - site catalogue (Pd-D unless stated)")
    Pr("Columns: d_min (A) | n_loc (a.u.) | hw (meV) | pairs/cm3 typ..hi (per cm2 if areal) |"
       " U_eff (eV) TF / Lindhard / Emp-lo / Emp-hi | log10 lambda per pair (s^-1) TF / Lind / Elo / Ehi")
    rows = []
    for s in SITES:
        r = site_models(s)
        rows.append((s, r))
        unit = "cm^-2" if s.get("areal") else "cm^-3"
        Pr(f"\n[{s['k']}] {s['name']}")
        Pr(f"   d_min={s['d_min']:.2f} A  n_loc={s['n_loc']:.3f} a.u. ({P.n_au_to_cm3(s['n_loc']):.1e} cm^-3)"
           f"  hw={1e3*s['hw']:.0f} meV  E_b: {s['E_b']}")
        Pr(f"   pairs: {s['dens_typ']:.1e} .. {s['dens_hi']:.1e} {unit}; control: {s['ctrl']}")
        if not s["mol"]:
            Pr(f"   U(r->0): TF {P.U_TF(s['n_loc']):.1f} eV, Lindhard {_lind_cache[s['n_loc']].U0:.1f} eV;"
               f" accelerator-scaled lo/hi {r['Uacc_lo']:.0f}/{r['Uacc_hi']:.0f} eV;"
               f" E_rel={1e3*E_rel(s['hw']):.0f} meV; harmonic confinement factor (not applied) {r['conf_TF']:.1e}")
        Pr(f"   U_eff:  TF {r['U_TF']:6.1f}  Lind {r['U_Lind']:6.1f}  Emp-lo {r['U_Elo']:6.1f}  Emp-hi {r['U_Ehi']:6.1f}"
           f"   rho0={r['rho0']:.1e}")
        Pr("   log10 lambda/pair: " + "  ".join(f"{m} {np.log10(r['lam_'+m]):7.1f}" for m in ("TF", "Lind", "Elo", "Ehi")))
        Pr("   log10 rate per " + ("cm2" if s.get("areal") else "cm3") + " (hi density): " +
           "  ".join(f"{m} {np.log10(r['lam_'+m]*s['dens_hi']):7.1f}" for m in ("TF", "Lind", "Elo", "Ehi")))
        lim = null_limit_per_pair(s)
        if lim is not None:
            Uc = P.Ue_required(lim, r["rho0"])
            flag = "EXCLUDED by bulk nulls" if r["lam_Ehi"] > lim else "allowed"
            Pr(f"   bulk-null limit: lambda <= {lim:.1e} /pair/s  (U_eff <= {Uc:.0f} eV); Emp-hi is {flag}")
        else:
            Pr("   bulk-null limit: none (class absent or negligible in 1989-90 null samples)")

    Pr("\nRequired U_eff for a 5-sigma/30-day proton signal (C3-like: 2.7e-3 fusions/s from M0 scaled) "
       "vs number of pairs, rho0 = 5e24:")
    for N in (1e15, 1e18, 1e21):
        Pr(f"   N={N:.0e}: U_eff >= {P.Ue_required(2.7e-3 / N, 5e24):.0f} eV")
    Pr("TF cannot exceed ~35 eV (U_eff ~13 eV) for any metallic density; reaching U(r->0)=300 eV "
       f"needs n = {((14.4/300/ P.a0_A)**-2 * np.pi / 4)**3 / (3*np.pi**2):.1e} a.u. "
       "(>1e4 x any metal): the static-linear-screening family is closed.")

    # -------- figure: U_eff by site and model
    fig, ax = plt.subplots(figsize=(9, 5.2))
    names = [s["k"] for s, _ in rows]
    y = np.arange(len(rows))
    cols = {"TF": "#2a78d6", "Lind": "#1baf7a", "Elo": "#eda100", "Ehi": "#eb6834"}
    lab = {"TF": "Thomas-Fermi (Yukawa)", "Lind": "Lindhard (RPA)",
           "Elo": "accelerator-scaled, low lab value", "Ehi": "accelerator-scaled, high lab value"}
    for i, m in enumerate(("TF", "Lind", "Elo", "Ehi")):
        ax.scatter([r["U_" + m] for _, r in rows], y + (i - 1.5) * 0.15, s=36, color=cols[m],
                   label=lab[m], zorder=3, edgecolor="white", linewidth=0.8)
    ax.axvline(KN_UEFF, color="grey", lw=1, ls=":")
    ax.text(KN_UEFF * 1.03, len(rows) - 0.4, "D$_2$ molecule (34 eV)", fontsize=7, color="grey")
    ax.axvspan(147, 165, color="grey", alpha=0.12)
    ax.text(150, -0.9, "bulk-null ceiling\n(147-165 eV)", fontsize=7, color="#52514e")
    ax.axvspan(208, 355, color="#eb6834", alpha=0.07)
    ax.text(212, len(rows) - 0.6, "needed for 5$\\sigma$/30 d\n(10$^{18}$-10$^{12}$ pairs)", fontsize=7,
            color="#52514e")
    ax.set_xscale("log")
    ax.set_yticks(y, names, fontsize=8)
    ax.set_xlabel("U$_{e,eff}$ seen by a thermal pair (eV)")
    ax.set_title("M1-A: effective screening by site class and model", fontsize=10)
    ax.grid(axis="x", color="#e6e6e6", lw=0.6)
    ax.set_axisbelow(True)
    ax.legend(fontsize=7, loc="lower right", frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(P.OUT, "m1_sites.png"), dpi=140)
    with open(os.path.join(P.OUT, "m1_sites.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    return rows


if __name__ == "__main__":
    main()
