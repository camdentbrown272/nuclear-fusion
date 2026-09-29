"""
M1 Part C - yield ladder per configuration C1-C6 and the ADR-001 sensitivity matrix
[configuration x site class]: the minimal localised enhancement at site class k that
configuration C detects at 5 sigma in 30 days.

M3 (cracks/vacancies) and M5 (detection) were not available on their branches when this
was written; site inventories are estimated here (documented per entry) and detection
follows M0/R5 assumptions.

Run: python3 sim/m1_ladder.py -> docs/models/figs/m1_ladder.txt, m1_matrix.png, m1_ladder.png
"""
import os
import warnings
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

import m1_physics as P
import m1_sites as S
import m1_nonthermal as NT

warnings.filterwarnings("ignore")
T = P.T30
CLASSES = ["bulkO", "bulkT", "vacD6", "vacD2", "void", "disl", "gb", "surf", "subsurf",
           "crack", "PdO", "oxint", "NP", "TiD2", "NiD"]
lines = []


def Pr(s=""):
    print(s)
    lines.append(s)


# ------------------------------------------------------------------ escape efficiency
_R = {}


def useful_range(product, host="PdD0.9"):
    """Path length (um) over which a product keeps more than its analysis threshold."""
    key = (product, host)
    if key not in _R:
        Z1, M1, E0, Ethr = {"p": (1, 1.007, 3020, 1500), "t": (1, 3.016, 1010, 300),
                            "he3": (2, 3.016, 820, 400)}[product]
        _R[key] = P.ion_range_um(E0, host, Z1, M1) - P.ion_range_um(Ethr, host, Z1, M1)
    return _R[key]


def eps_charged(z_um, f_omega, host="PdD0.9", water_um=0.0):
    """
    Detected charged products per fusion, site at depth z (um) below a face viewed by a
    detector covering solid-angle fraction f_omega (cone about the normal).
    Isotropic emission; a product counts if theta < theta_det and z/cos(theta) < useful range:
      fraction = (1 - max(cos theta_det, z/R_u)) / 2.
    p-branch (1/2): p or t (back-to-back, mutually exclusive); n-branch (1/2): 3He.
    water_um: electrolyte film between face and detector (CR-39 case) reduces the useful range.
    """
    cos_d = 1 - 2 * f_omega
    frac = {}
    for prod in ("p", "t", "he3"):
        Ru = useful_range(prod, host)
        if water_um > 0:
            Ru = Ru * max(1 - water_um / useful_range(prod, "D2O"), 0.0)
        frac[prod] = 0.0 if (Ru <= 0 or z_um >= Ru) else (1 - max(cos_d, z_um / Ru)) / 2
    return 0.5 * (frac["p"] + frac["t"]) + 0.5 * frac["he3"]


def eps_profile(z_lo, z_hi, f_omega, host="PdD0.9", water_um=0.0, n=60):
    zs = np.linspace(z_lo, z_hi, n) if z_hi > z_lo else np.array([z_lo])
    return float(np.mean([eps_charged(z, f_omega, host, water_um) for z in zs]))


# ------------------------------------------------------------------ configurations
def dens(k, which="typ"):
    return S.SITE[k]["dens_" + which]


def make_configs(b5):
    """
    Each config: inventory {class: (pairs, (z_lo, z_hi) depth window in um or None)},
    charged-channel params, neutron/511 channel params, D inventory for cosmic floor,
    D-specific non-fusion neutron rate (B5), and the drives.
    All numbers are [guess]-grade engineering estimates unless tagged; see M1 write-up sec. 3.
    """
    Dn_PdD = b5["PdD"]["n2n"] + b5["PdD"]["gn"]            # per 6.1e22 D
    Dn_D2O = b5["D2O"]["n2n"] + b5["D2O"]["gn"]            # per 50 mL
    ko_PdD, ko_D2O = b5["PdD"]["knock"], b5["D2O"]["knock"]
    C = {}

    # ---- C1 coaxial electrolytic cell: Pd wire 1 mm x 5 cm, 50 mL D2O, neutrons + 511 only
    V, A = np.pi * 0.05 ** 2 * 5, np.pi * 0.1 * 5
    C["C1"] = dict(
        desc="Pd wire 1 mm x 5 cm (0.039 cm3, 1.6 cm2), x=0.9, annealed 20 um grains, 50 mL D2O",
        inv={"bulkO": (dens("bulkO") * V, None), "bulkT": (dens("bulkT") * V, None),
             "vacD6": (12 * 6.8e19 * A * 1e-4, None), "vacD2": (6.8e19 * A * 1e-4, None),
             "void": (dens("void") * V, None), "disl": (dens("disl") * V, None),
             "gb": (dens("gb") * V, None), "surf": (dens("surf") * A, None),
             "subsurf": (dens("subsurf") * A, None), "crack": (dens("crack") * V, None)},
        f_omega=0.0, cp_bkg=0.0, cp_sys=0.0, host="PdD0.9", water_um=0.0,
        eps_n=0.10, n_bkg=0.05, n_sys=0.01, eps_511=0.01, b511=2e-3,
        Dn=Dn_PdD * V + Dn_D2O, knock=ko_PdD * V + ko_D2O,
        drives="galvanostatic loading; ~1 alpha/beta cycle per week; no exit face (no net flux)")

    # ---- C2 Pd/D codeposition, CR-39 in contact
    V, A = 5e-4 * 2 * 0.5, 2.0              # 5 um x 2 cm2, 50 % porosity
    C["C2"] = dict(
        desc="Pd/D codeposit 5 um on 2 cm2 mesh facing CR-39 through ~20 um electrolyte; 20 nm grains",
        inv={"bulkO": (dens("bulkO") * V, (0, 5)), "bulkT": (dens("bulkT") * V, (0, 5)),
             "vacD6": (dens("vacD6", "hi") * 0.1 * V, (0, 5)), "vacD2": (dens("vacD2", "hi") * 0.1 * V, (0, 5)),
             "void": (dens("void", "hi") * V, (0, 5)), "disl": (dens("disl", "hi") * V, (0, 5)),
             "gb": (dens("gb", "hi") * V, (0, 5)), "surf": (dens("surf") * A * 20, (0, 5)),
             "subsurf": (dens("subsurf") * A * 20, (0, 5)), "crack": (1.5e15 * 2 * 1e4 * 0.1 * V, (0, 5))},
        f_omega=0.5, cp_bkg=60 / T, cp_sys=0.30, host="PdD0.9", water_um=20.0, cp_name="CR-39",
        eps_n=0.01, n_bkg=0.01, n_sys=0.05, eps_511=0.0, b511=1.0,
        Dn=Dn_D2O * 0.6, knock=ko_D2O * 0.6,
        drives="continuous deposition at high D fugacity; no controlled cycling")

    # ---- C3 detector-facing membrane (DFM)
    A, t = 2.0, 25.0
    V = A * t * 1e-4
    Vf = A * 1e-4                           # front 1 um nanostructured/codeposited layer
    C["C3"] = dict(
        desc="Pd membrane 25 um x 2 cm2, back-loaded (20 mL D2O), front: 1 um nanocrystalline "
             "codeposit + PdO 5 nm (or CaO/Pd x5) + 5 nm Pd NP decoration, Si in vacuum",
        inv={"bulkO": (dens("bulkO") * V, (0, t)), "bulkT": (dens("bulkT") * V, (0, t)),
             "vacD6": (dens("vacD6", "hi") * 0.1 * Vf, (0, 1)),
             "vacD2": (dens("vacD2", "hi") * 0.1 * Vf, (0, 1)),
             "void": (dens("void", "hi") * V, (0, t)), "disl": (dens("disl", "hi") * V, (0, t)),
             "gb": (dens("gb", "hi") * Vf, (0, 1)), "surf": (dens("surf") * A * 10, (0, 0)),
             "subsurf": (dens("subsurf") * A * 10, (0, 0.001)),
             "crack": (1.5e15 * 2 * 100 * 0.1 * V + 1.5e15 * 2 * 1e4 * 0.1 * Vf, (0, t)),
             "PdO": (dens("PdO") * A, (0.005, 0.005)), "oxint": (5 * dens("oxint") * A, (0.01, 0.1)),
             "NP": (dens("NP") * A * 5e-7 * 0.5, (0, 0.005))},
        f_omega=0.18, cp_bkg=6e-5, cp_sys=0.05, host="PdD0.9", water_um=0.0, cp_name="Si",
        eps_n=0.15, n_bkg=0.05, n_sys=0.01, eps_511=0.02, b511=2e-3,
        Dn=Dn_PdD * V + Dn_D2O * 0.4, knock=ko_PdD * V + ko_D2O * 0.4,
        drives="current modulation: permeation flux up to J_max (Part D), 24 alpha/beta cycles/day")

    # ---- C4 gas-phase Pd/ZrO2 nanocomposite, 10 g Pd as 5 nm particles, 250 C
    Vpd, d = 10 / 12.0, 5e-7
    Asurf = 6 / d * Vpd
    C["C4"] = dict(
        desc="100 g Pd/ZrO2 composite (10 g Pd, 5 nm particles), D2 gas 250 C, calorimeter + outside n bank",
        inv={"NP": (dens("NP", "typ") * 0.55 * Vpd, None), "bulkT": (dens("bulkT") * 20 * 0.5 * Vpd, None),
             "surf": (dens("surf") * Asurf * 0.5, None), "subsurf": (dens("subsurf") * Asurf, None),
             "oxint": (dens("oxint") * Asurf * 0.5, None), "vacD6": (12 * 6.8e19 * Vpd, None),
             "crack": (1.5e15 * 2 * 1e5 * 0.1 * Vpd, None)},
        f_omega=0.0, cp_bkg=0.0, cp_sys=0.0, host="PdD0.9", water_um=0.0,
        eps_n=0.05, n_bkg=0.05, n_sys=0.02, eps_511=0.005, b511=2e-3,
        Dn=Dn_PdD * Vpd * 0.55, knock=ko_PdD * Vpd * 0.55,
        drives="temperature cycling 200-300 C, absorption/desorption bursts")

    # ---- C5 Ni/Cu multilayer on Ni, D2 gas, heated; Si facing one face
    A = 6.25
    int_pairs = 12 * A * 1.3e15 * 0.1 * 6        # 12 interfaces per face, theta=0.1 [guess]
    C["C5"] = dict(
        desc="25x25x0.1 mm Ni with 6x(2 nm Cu/14 nm Ni) on both faces, D2-loaded then heated; Si at one face",
        inv={"NiD": (2 * int_pairs, (0, 0.1)), "surf": (2 * 4.5e15 * A * 0.1, (0, 0))},
        vis_frac={"NiD": 0.5, "surf": 0.5},
        f_omega=0.02, cp_bkg=1e-4, cp_sys=0.05, host="PdD0.9", water_um=0.0, cp_name="Si",
        eps_n=0.05, n_bkg=0.05, n_sys=0.02, eps_511=0.005, b511=2e-3,
        Dn=Dn_PdD * 5e18 / 6.1e22, knock=0.0,
        drives="heating ramps to 500-900 C: desorption flux 1e14-1e15 D/cm2/s through interfaces")

    # ---- C6 Ti chips, T/P cycling, 3He well counter
    V = 100 / 4.5
    C["C6"] = dict(
        desc="100 g Ti chips (22 cm3) as TiD1.5, 12 LN2<->RT cycles/day, 3He well counter eps 0.3",
        inv={"TiD2": (dens("TiD2") * 0.75 * V, None), "surf": (dens("surf") * 1e4, None),
             "crack": (1.5e15 * 2 * 1e3 * 0.1 * V, None), "disl": (dens("disl") * V, None),
             "gb": (dens("gb") * 2 * V, None)},
        f_omega=0.0, cp_bkg=0.0, cp_sys=0.0, host="TiD2", water_um=0.0,
        eps_n=0.30, n_bkg=0.05, n_sys=0.01, eps_511=0.0, b511=1.0,
        Dn=Dn_PdD * 0.75 * V * 1.1e23 / 6.1e22, knock=ko_PdD * 0.75 * V * 1.1e23 / 6.1e22,
        drives="thermal cycling; fracto events in metallic TiD2 are energetically dead (B1)",
        null_host="Ti")
    return C


# ------------------------------------------------------------------ detection model
def channel_thresholds(c, f_ee=0.0):
    """
    Minimal detected-signal counts in 30 d for each channel of config c and the per-fusion
    detection efficiency factor that multiplies a pair's escape efficiency.
    f_ee: fraction of fusions going to the e+e- channel (Czerski scenario).
    """
    ch = {}
    if c["f_omega"] > 0:
        b = c["cp_bkg"] * T
        ch["cp"] = dict(smin=P.s_min_5sigma(b, c["cp_sys"]), scale=1 - f_ee)
    b = c["n_bkg"] * T + c["eps_n"] * c["Dn"] * T
    sysn = np.hypot(c["n_sys"] * c["n_bkg"] * T, 0.3 * c["eps_n"] * c["Dn"] * T) / max(b, 1e-30)
    ch["n"] = dict(smin=P.s_min_5sigma(b, sysn), scale=c["eps_n"] * 0.5 * (1 - f_ee))
    if c["eps_511"] > 0 and f_ee > 0:
        ch["ee"] = dict(smin=P.s_min_5sigma(c["b511"] * T, 0.02), scale=c["eps_511"] * f_ee)
    return ch


def eff_pairs(c, k, ch_name):
    """Detection-weighted pair count N_k * eps_k for channel ch_name."""
    N, win = c["inv"][k]
    if ch_name == "cp":
        if win is None:
            return 0.0
        vis = c.get("vis_frac", {}).get(k, 1.0)
        return N * vis * eps_profile(win[0], win[1], c["f_omega"], c["host"], c["water_um"])
    return N


def lam_min(c, k, f_ee=0.0):
    """Minimal per-pair rate (s^-1) at class k detectable at 5 sigma in 30 d (best channel)."""
    best = np.inf
    for name, ch in channel_thresholds(c, f_ee).items():
        Ne = eff_pairs(c, k, "cp" if name == "cp" else "x") * ch["scale"]
        if Ne > 0:
            best = min(best, ch["smin"] / (T * Ne))
    return best


def events_per_day(c, lam_by_class, f_ee=0.0):
    """Detected events/day per channel for given per-class rates (plus cosmic floor in n)."""
    out = {}
    for name, ch in channel_thresholds(c, f_ee).items():
        tot = sum(eff_pairs(c, k, "cp" if name == "cp" else "x") * ch["scale"] * lam_by_class.get(k, 0.0)
                  for k in c["inv"])
        if name == "n":
            tot += c["knock"] * c["eps_n"] * 0.5
        out[name] = (tot * P.DAY, ch["smin"] / 30.0)
    return out


def main():
    rates = {s["k"]: S.site_models(s) for s in S.SITES}
    NT.lines.clear()
    b5 = NT.part_B5()
    C = make_configs(b5)

    Pr("M1 Part C - yield ladder and sensitivity matrix (5 sigma, 30 days)")
    Pr("Useful ranges in PdD0.9 (um): " + ", ".join(f"{p}: {useful_range(p):.1f}" for p in ("p", "t", "he3")))
    for cn, c in C.items():
        Pr(f"\n[{cn}] {c['desc']}")
        Pr(f"   drives: {c['drives']}")
        for name, ch in channel_thresholds(c).items():
            Pr(f"   channel {name}: s_min(30 d) = {ch['smin']:.0f} counts, efficiency factor {ch['scale']:.3g}")
        Pr(f"   D-specific non-fusion neutrons {c['Dn']:.1e} n/s; cosmic knock-on fusions {c['knock']:.1e} /s")
        Pr("   class      pairs      eps_cp     log10 lam_min  log10 n(TF)  U_req(eV)")
        for k, (N, win) in c["inv"].items():
            e = eff_pairs(c, k, "cp") / N if c["f_omega"] > 0 else 0.0
            lm = lam_min(c, k)
            n_TF = np.log10(lm / rates[k]["lam_TF"])
            Ur = P.Ue_required(lm, rates[k]["rho0"])
            Pr(f"   {k:8s} {N:10.2e}  {e:9.2e}   {np.log10(lm):9.1f}     {n_TF:8.1f}   {Ur:7.0f}")

    # ---------------- ladder: events/day under scenarios
    Pr("\nYield ladder: detected events per day [best channel]  (threshold events/day for 5 sigma in 30 d)")
    Pr("   config | (i) standard (TF) + cosmic floor | (ii-lo) accel. low | (ii-hi) accel. high | threshold")
    ladder = {}
    for cn, c in C.items():
        row = {}
        for scen, key in (("i", "lam_TF"), ("iilo", "lam_Elo"), ("iihi", "lam_Ehi")):
            lam = {k: rates[k][key] for k in c["inv"]}
            ev = events_per_day(c, lam)
            best = max(ev.items(), key=lambda kv: kv[1][0] / kv[1][1])
            row[scen] = (best[0], best[1][0], best[1][1])
        ladder[cn] = row
        Pr(f"   {cn}     | {row['i'][1]:9.2e} ({row['i'][0]})               | {row['iilo'][1]:9.2e} ({row['iilo'][0]})"
           f"   | {row['iihi'][1]:9.2e} ({row['iihi'][0]})    | {row['iihi'][2]:.2g}")
    Pr("   (ii-hi) is already EXCLUDED by prior nulls for most classes (m1_sites); rerun with every class capped"
       " at its prior-null limit (the most optimistic rate still allowed by data):")
    for cn, c in C.items():
        lam = {}
        for k in c["inv"]:
            lim = S.null_limit_per_pair(S.SITE[k], c.get("null_host", "Pd"))
            lam[k] = rates[k]["lam_Ehi"] if lim is None else min(rates[k]["lam_Ehi"], lim)
        ev = events_per_day(c, lam)
        top = sorted(((eff_pairs(c, k, 'cp' if 'cp' in ev else 'x') * lam[k], k) for k in c["inv"]), reverse=True)[:2]
        Pr(f"   {cn}: " + ", ".join(f"{n}: {v[0]:.2e}/d (thr {v[1]:.2g})" for n, v in ev.items())
           + f"   dominant classes: {[t[1] for t in top]}")

    # ---------------- the matrix
    Pr("\nSENSITIVITY MATRIX (standard D-D branching). Entry = log10 of the minimal enhancement over the"
       " Thomas-Fermi standard rate | required U_eff (eV) | '*' = required rate already excluded by prior nulls")
    hdr = "   class    " + "".join(f"{cn:>16s}" for cn in C)
    Pr(hdr)
    M_U = np.full((len(CLASSES), len(C)), np.nan)
    M_n = np.full_like(M_U, np.nan)
    M_null = np.zeros_like(M_U, dtype=bool)
    for i, k in enumerate(CLASSES):
        cells = []
        for j, (cn, c) in enumerate(C.items()):
            if k not in c["inv"]:
                cells.append(f"{'-':>16s}")
                continue
            lm = lam_min(c, k)
            if not np.isfinite(lm):
                cells.append(f"{'blind':>16s}")
                continue
            n = np.log10(lm / rates[k]["lam_TF"])
            Ur = P.Ue_required(lm, rates[k]["rho0"])
            lim = S.null_limit_per_pair(S.SITE[k], c.get("null_host", "Pd"))
            star = lim is not None and lm > lim
            M_U[i, j], M_n[i, j], M_null[i, j] = Ur, n, star
            cells.append(f"{n:7.0f} |{Ur:5.0f}{'*' if star else ' '}  ")
        Pr(f"   {k:8s}" + "".join(cells))
    Pr("\nSame matrix, entry = log10 of minimal enhancement over the (ii-lo) 'Tohoku-static' rate (negative ="
       " (ii-lo) would already be detected):")
    Pr(hdr)
    for k in CLASSES:
        cells = []
        for cn, c in C.items():
            if k not in c["inv"] or not np.isfinite(lam_min(c, k)):
                cells.append(f"{'-':>16s}")
                continue
            cells.append(f"{np.log10(lam_min(c, k) / rates[k]['lam_Elo']):16.1f}")
        Pr(f"   {k:8s}" + "".join(cells))

    # ---------------- e+e- scenario
    Pr("\nCzerski e+e- scenario (f_ee = 0.91): minimal per-pair rate, log10 (standard branching -> e+e-):")
    for cn, c in C.items():
        ks = [k for k in ("vacD6", "PdO", "NP", "TiD2", "surf") if k in c["inv"]][:3]
        Pr(f"   {cn}: " + ", ".join(f"{k}: {np.log10(lam_min(c, k)):.1f} -> {np.log10(lam_min(c, k, 0.91)):.1f}"
                                   for k in ks))

    # ---------------- sensitivities
    Pr("\nSensitivities (log10 lam_min at the named class; U_req in eV):")
    c3 = C["C3"]
    base = {k: lam_min(c3, k) for k in ("PdO", "vacD6", "NP")}
    for lab, mod in (("C3 baseline", {}), ("C3 close geometry f_omega=0.31", {"f_omega": 0.31}),
                     ("C3 Si background x10", {"cp_bkg": 6e-4}),
                     ("C3 Si background systematic 20%", {"cp_sys": 0.20})):
        cc = dict(c3, **mod)
        Pr(f"   {lab:36s}: " + ", ".join(f"{k} {np.log10(lam_min(cc, k)):.2f} (U {P.Ue_required(lam_min(cc, k), rates[k]['rho0']):.0f})"
                                         for k in base))
    cc = dict(c3, inv={k: (v[0] * 5, v[1]) for k, v in c3["inv"].items()})
    Pr(f"   {'C3 area 2 -> 10 cm2 (5 detectors)':36s}: " + ", ".join(
        f"{k} {np.log10(lam_min(cc, k)):.2f} (U {P.Ue_required(lam_min(cc, k), rates[k]['rho0']):.0f})" for k in base))
    for cn in ("C1", "C6"):
        for sy in (0.01, 0.001):
            cc = dict(C[cn], n_sys=sy)
            k = "bulkO" if cn == "C1" else "TiD2"
            Pr(f"   {cn} neutron background systematic {sy:.1%}: s_min {channel_thresholds(cc)['n']['smin']:.0f},"
               f" {k} log10 lam_min {np.log10(lam_min(cc, k)):.2f}")
    for f in (0.1, 10):
        Pr(f"   rho0 x{f}: C3 PdO U_req {P.Ue_required(lam_min(c3, 'PdO'), rates['PdO']['rho0'] * f):.0f} eV"
           f" (baseline {P.Ue_required(lam_min(c3, 'PdO'), rates['PdO']['rho0']):.0f})")
    lam_pdo = rates["PdO"]["lam_Elo"]
    for ER, lab in ((0.04, "Yukawa-mapped at 0.04 eV (baseline)"), (1.0, "mapped at 1 eV relative energy")):
        U = P.Ueff_yukawa(600, ER)
        Pr(f"   PdO (ii-lo) {lab}: U_eff {U:.0f} eV -> C3 events/day "
           f"{eff_pairs(c3, 'PdO', 'cp') * P.rate_pair(U, rates['PdO']['rho0']) * P.DAY:.2g}")

    # ---------------- figures
    fig, ax = plt.subplots(figsize=(8.2, 6.2))
    cmap = LinearSegmentedColormap.from_list("seq", ["#dce9f8", "#2a78d6", "#0b2d5c"])
    im = ax.imshow(M_U, cmap=cmap, vmin=150, vmax=320, aspect="auto")
    for i in range(M_U.shape[0]):
        for j in range(M_U.shape[1]):
            if np.isnan(M_U[i, j]):
                continue
            col = "white" if M_U[i, j] > 240 else "#0b0b0b"
            ax.text(j, i, f"{M_U[i, j]:.0f}{'*' if M_null[i, j] else ''}\n10$^{{{M_n[i, j]:.0f}}}$",
                    ha="center", va="center", fontsize=6.5, color=col)
    ax.set_xticks(range(len(C)), list(C))
    ax.set_yticks(range(len(CLASSES)), CLASSES, fontsize=8)
    cb = fig.colorbar(im, ax=ax, shrink=0.8)
    cb.set_label("U$_{eff}$ the class must have to be seen (eV)")
    ax.set_title("M1-C: minimal localised anomaly detectable at 5$\\sigma$ in 30 d\n"
                 "cell: required U$_{eff}$ (eV) and enhancement over Thomas-Fermi rate; "
                 "* = already excluded by prior nulls; blank = class absent", fontsize=8.5)
    fig.tight_layout()
    fig.savefig(os.path.join(P.OUT, "m1_matrix.png"), dpi=150)

    fig, ax = plt.subplots(figsize=(6.8, 4.4))
    Us = np.linspace(130, 330, 300)
    for j, (cn, c) in enumerate(C.items()):
        # best class not already excluded by prior nulls (lowest required U_eff)
        cand = []
        for k in c["inv"]:
            lm = lam_min(c, k)
            lim = S.null_limit_per_pair(S.SITE[k], c.get("null_host", "Pd"))
            if np.isfinite(lm) and (lim is None or lm <= lim):
                cand.append((P.Ue_required(lm, rates[k]["rho0"]), k))
        if not cand:
            cand = [(P.Ue_required(lam_min(c, k), rates[k]["rho0"]), k) for k in c["inv"]]
            tag = " (all excluded)"
        else:
            tag = ""
        Ur, kbest = min(cand)
        lm = lam_min(c, kbest)
        per_day = P.rate_pair(Us, rates[kbest]["rho0"]) / lm * (lam_min(c, kbest) * 0 + 1)
        thr = [ch["smin"] / 30 for n, ch in channel_thresholds(c).items()]
        thr_best = None
        for n, ch in channel_thresholds(c).items():
            Ne = eff_pairs(c, kbest, "cp" if n == "cp" else "x") * ch["scale"]
            if Ne > 0 and np.isclose(ch["smin"] / (T * Ne), lm):
                thr_best = ch["smin"] / 30
        ev = per_day * thr_best
        ax.semilogy(Us, ev, color=NT.COL[j], lw=2, label=f"{cn}: {kbest}{tag}")
        ax.plot([Ur], [thr_best], "o", color=NT.COL[j], ms=6, mec="white", zorder=4)
        lim = S.null_limit_per_pair(S.SITE[kbest], c.get("null_host", "Pd"))
        if lim is not None:
            Ul = P.Ue_required(lim, rates[kbest]["rho0"])
            if Us[0] < Ul < Us[-1]:
                ax.plot([Ul], [thr_best * lim / lm], "x", color=NT.COL[j], ms=7, mew=2, zorder=4)
    ax.set_ylim(1e-3, 1e6)
    ax.set_xlabel("U$_{eff}$ at the site class (eV)  [anomaly expressed as effective screening]")
    ax.set_ylabel("detected events per day (best channel)")
    ax.set_title("M1-C: yield ladder, most sensitive not-yet-excluded class per configuration\n"
                 "dot = 5$\\sigma$ in 30 d; x = ceiling from prior null experiments", fontsize=9)
    ax.legend(fontsize=7, frameon=False, loc="upper left")
    ax.grid(color="#e6e6e6", lw=0.6)
    fig.tight_layout()
    fig.savefig(os.path.join(P.OUT, "m1_ladder.png"), dpi=140)

    with open(os.path.join(P.OUT, "m1_ladder.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    return C, rates, M_U, M_n, M_null


if __name__ == "__main__":
    main()
