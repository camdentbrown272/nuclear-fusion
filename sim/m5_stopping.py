#!/usr/bin/env python3
"""M5 stopping-power / range engine for light ions (p, d, t, 3He, alpha).

Electronic stopping: Ziegler (SRIM-type) elemental parameterisation as implemented in
catima/pycatima (https://pypi.org/project/pycatima/ , https://github.com/hrosiak/catima),
combined with Bragg's additivity rule for compounds. Isotopes (d, t, 3He) use velocity
scaling (stopping per atom depends only on T/u and projectile Z).
Nuclear stopping: ZBL universal nuclear stopping (Ziegler, Biersack, Littmark 1985).
Independent check: Bethe-Bloch with Barkas-Berger shell correction and Bloch term,
implemented here from first principles (valid >~ 2 MeV/u).
Validation reference: NIST PSTAR/ASTAR electronic stopping (transcribed in Geant4
G4PSTARStopping.cc / G4ASTARStopping.cc, file sim/m5_nist_star.json).

Run directly to regenerate the validation figure/table and range tables:
    python3 sim/m5_stopping.py
"""
import json
import os
from functools import lru_cache

import numpy as np

try:
    import pycatima as _cat
except ImportError as e:  # pragma: no cover
    raise SystemExit("pip install pycatima  (SRIM-type stopping coefficients)") from e

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FIGS = os.path.join(ROOT, "docs", "models", "figs")
NA = 6.02214076e23
ME_C2 = 0.51099895  # MeV
AMU = 931.49410242  # MeV

# ---------------------------------------------------------------------------------------
# Projectiles: (Z, mass in u)  -- nuclear masses (CODATA 2018)
PROJ = {
    "p": (1, 1.007276),
    "d": (1, 2.013553),
    "t": (1, 3.015501),
    "h": (2, 3.014932),  # 3He
    "a": (2, 4.001506),  # 4He
}

# Elements: Z -> (symbol, A [g/mol], I [eV] ICRU-37/49 solid/condensed values)
# I-values: ICRU Report 37 (1984) / 49 (1993) as tabulated by NIST ESTAR/PSTAR
# https://physics.nist.gov/PhysRefData/Star/Text/method.html
ELEM = {
    1: ("H", 1.008, 19.2), 6: ("C", 12.011, 78.0), 7: ("N", 14.007, 82.0),
    8: ("O", 15.999, 95.0), 13: ("Al", 26.982, 166.0), 14: ("Si", 28.086, 173.0),
    18: ("Ar", 39.948, 188.0), 20: ("Ca", 40.078, 191.0), 26: ("Fe", 55.845, 286.0),
    28: ("Ni", 58.693, 311.0), 29: ("Cu", 63.546, 322.0), 46: ("Pd", 106.42, 470.0),
    47: ("Ag", 107.868, 470.0), 79: ("Au", 196.967, 790.0), 78: ("Pt", 195.08, 790.0),
}
M_D = 2.014  # deuterium atomic mass (g/mol)

# Materials: name -> (list of (Z, atoms per formula unit, atomic mass), density g/cm3, source)
MAT = {
    "Pd": ([(46, 1, 106.42)], 12.02, "CRC Handbook; https://www.webelements.com/palladium/physics.html"),
    "PdD0.9": ([(46, 1, 106.42), (1, 0.9, M_D)], 10.90,
               "4 f.u./(N_A a^3), a = 4.04 A for beta-PdD_x at x~0.9 (a = 4.02-4.025 A at beta-phase "
               "onset; e.g. https://inis.iaea.org/records/9hbv3-2bz85); +/-2 %"),
    "PdO": ([(46, 1, 106.42), (8, 1, 15.999)], 8.3, "CRC Handbook (PdO 8.3 g/cm3)"),
    "Au": ([(79, 1, 196.967)], 19.32, "CRC Handbook"),
    "Ni": ([(28, 1, 58.693)], 8.908, "CRC Handbook"),
    "Cu": ([(29, 1, 63.546)], 8.96, "CRC Handbook"),
    "CaO": ([(20, 1, 40.078), (8, 1, 15.999)], 3.34, "CRC Handbook"),
    "D2O": ([(1, 2, M_D), (8, 1, 15.999)], 1.104, "CRC Handbook (25 C)"),
    "H2O": ([(1, 2, 1.008), (8, 1, 15.999)], 0.997, "CRC Handbook (25 C)"),
    "Mylar": ([(6, 10, 12.011), (1, 8, 1.008), (8, 4, 15.999)], 1.397,
              "NIST PSTAR material 222 (polyethylene terephthalate)"),
    "Al": ([(13, 1, 26.982)], 2.699, "NIST"),
    "D2gas": ([(1, 2, M_D)], 1.654e-4, "ideal gas, 1 bar, 293 K"),
    "Si": ([(14, 1, 28.086)], 2.329, "NIST"),
    "CR39": ([(6, 12, 12.011), (1, 18, 1.008), (8, 7, 15.999)], 1.31,
             "TASTRAK PADC C12H18O7, https://www.tasl.co.uk/tastrak-padc.php"),
    "Ag": ([(47, 1, 107.868)], 10.50, "CRC"),
    "Pt": ([(78, 1, 195.08)], 21.45, "CRC"),
    "air": ([(7, 1.562, 14.007), (8, 0.42, 15.999), (18, 0.0093, 39.948)], 1.205e-3, "dry air 1 atm 20 C"),
    "HDPE": ([(6, 1, 12.011), (1, 2, 1.008)], 0.95, "typical HDPE"),
}


def mat_mass(name):
    comp, rho, _ = MAT[name]
    return sum(n * a for _, n, a in comp)


# ---------------------------------------------------------------------------------------
def _srim_e(zp, zt, t_per_u):
    """Electronic stopping cross-section, eV/(1e15 atoms/cm2), SRIM-type parameterisation."""
    return _cat.srim_dedx_e(int(zp), int(zt), float(t_per_u), False)


def zbl_nuclear(zp, mp, zt, mt, e_keV):
    """ZBL universal nuclear stopping, eV/(1e15 atoms/cm2). e_keV = lab kinetic energy."""
    zz = zp ** 0.23 + zt ** 0.23
    eps = 32.53 * mt * e_keV / (zp * zt * (mp + mt) * zz)
    if eps <= 30:
        sn = np.log(1 + 1.1383 * eps) / (2 * (eps + 0.01321 * eps ** 0.21226 + 0.19593 * eps ** 0.5))
    else:
        sn = np.log(eps) / (2 * eps)
    return 8.462 * zp * zt * mp * sn / ((mp + mt) * zz)


def stopping(part, mat, T, nuclear=True):
    """Total mass stopping power, MeV cm^2/g, for projectile `part` of kinetic energy T (MeV)."""
    zp, mp = PROJ[part]
    comp, rho, _ = MAT[mat]
    M = mat_mass(mat)
    s = 0.0  # eV/(1e15 formula units / cm2)
    for zt, n, at in comp:
        se = _srim_e(zp, zt, T / mp)
        sn = zbl_nuclear(zp, mp, zt, at, T * 1e3) if nuclear else 0.0
        s += n * (se + sn)
    return s * 1e-15 * 1e-6 * NA / M  # MeV cm2/g


def electronic(part, mat, T):
    return stopping(part, mat, T, nuclear=False)


# ---------------------------------------------------------------------------------------
# Bethe-Bloch (independent first-principles check)
K_BB = 0.307075  # MeV cm2/mol (4 pi N_A r_e^2 m_e c^2)


def shell_barkas_berger(eta, I_eV):
    """Barkas-Berger (1964) shell correction C (total, divide by Z). eta = beta*gamma >= ~0.13."""
    I = I_eV
    return ((0.422377 * eta ** -2 + 0.0304043 * eta ** -4 - 0.00038106 * eta ** -6) * 1e-6 * I ** 2
            + (3.858019 * eta ** -2 - 0.1667989 * eta ** -4 + 0.00157955 * eta ** -6) * 1e-9 * I ** 3)


def bloch(y, n=200):
    k = np.arange(1, n + 1)
    return -y ** 2 * np.sum(1.0 / (k * (k ** 2 + y ** 2)))


def bethe_bloch_element(part, Z, A, I_eV, T):
    """Electronic mass stopping power (MeV cm2/g) from Bethe-Bloch + shell + Bloch."""
    zp, mp = PROJ[part]
    M = mp * AMU
    g = 1 + T / M
    b2 = 1 - 1 / g ** 2
    eta = np.sqrt(b2) * g
    L0 = np.log(2 * ME_C2 * 1e6 * b2 * g ** 2 / I_eV) - b2
    C = shell_barkas_berger(eta, I_eV) / Z
    L2 = bloch(zp / 137.036 / np.sqrt(b2))
    return K_BB * zp ** 2 * Z / A / b2 * (L0 - C + L2)


# ---------------------------------------------------------------------------------------
EGRID = np.logspace(-3, np.log10(40.0), 900)  # MeV


@lru_cache(maxsize=None)
def table(part, mat):
    """Return (E grid MeV, S_lin MeV/um, CSDA range um) for projectile in material."""
    rho = MAT[mat][1]
    S = np.array([stopping(part, mat, e) for e in EGRID])  # MeV cm2/g
    S_lin = S * rho * 1e-4  # MeV/um
    # range: integrate dE/S from ~0 (assume S ~ sqrt(E) below 1 keV -> R(1keV)=2 E/S)
    R0 = 2 * EGRID[0] / S_lin[0]
    inv = 1.0 / S_lin
    R = R0 + np.concatenate([[0], np.cumsum(0.5 * (inv[1:] + inv[:-1]) * np.diff(EGRID))])
    return EGRID, S_lin, R


def rng(part, mat, E):
    Eg, S, R = table(part, mat)
    return np.interp(np.log(np.maximum(E, 1e-3)), np.log(Eg), R)


def e_after(part, mat, E, L_um):
    """Residual energy (MeV) after straight path L_um; 0 if stopped. Vectorised."""
    Eg, S, R = table(part, mat)
    E = np.asarray(E, float)
    Rres = np.interp(np.log(np.maximum(E, 1e-3)), np.log(Eg), R) - np.asarray(L_um, float)
    out = np.interp(Rres, R, Eg)
    return np.where(Rres > R[0], out, 0.0)


def dedx_lin(part, mat, E):
    Eg, S, R = table(part, mat)
    return np.interp(np.log(np.maximum(E, 1e-3)), np.log(Eg), S)


def bohr_sigma2(part, mat, L_um):
    """Bohr energy-straggling variance (MeV^2) for path L_um (upper bound at low velocity)."""
    zp, _ = PROJ[part]
    comp, rho, _ = MAT[mat]
    M = mat_mass(mat)
    ne = rho / M * NA * sum(n * z for z, n, _ in comp)  # electrons/cm3
    # 4 pi z^2 e^4 n_e L ; e^2 = 1.44e-13 MeV cm
    return 4 * np.pi * zp ** 2 * (1.44e-13) ** 2 * ne * (L_um * 1e-4)


# ---------------------------------------------------------------------------------------
def validate(out):
    """Compare engine vs NIST PSTAR/ASTAR and Bethe-Bloch; write figure + text."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    d = json.load(open(os.path.join(HERE, "m5_nist_star.json")))
    T = np.array(d["T_MeV"])
    nist_map = {"Al": "G4_Al", "Si": "G4_Si", "Cu": "G4_Cu", "Ag": "G4_Ag", "Au": "G4_Au",
                "Pt": "G4_Pt", "Mylar": "G4_MYLAR", "H2O": "G4_WATER", "air": "G4_AIR"}
    print("== Validation: engine electronic stopping / NIST (PSTAR p <=2 MeV, ASTAR alpha <=10 MeV)", file=out)
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
    worst = {}
    for part, key, Emax, Es in [("p", "PSTAR", 2.0, [0.1, 0.3, 0.5, 1.0, 2.0]),
                                ("a", "ASTAR", 10.0, [0.5, 1.0, 2.0, 5.0, 8.0, 10.0])]:
        n = 60 if key == "PSTAR" else 78
        print(f"  {part}: E(MeV) " + " ".join(f"{e:>6}" for e in Es), file=out)
        ax = axs[0 if part == "p" else 1]
        for m, g in nist_map.items():
            ref = np.array(d[key][g])
            Tn = T[:n]
            mine = np.array([electronic(part, m, e) for e in Tn])
            ratio = mine / ref
            sel = Tn >= (0.1 if part == "p" else 0.4)
            worst[(part, m)] = np.max(np.abs(ratio[sel] - 1))
            ax.semilogx(Tn, ratio, label=m)
            vals = [np.interp(e, Tn, ratio) for e in Es]
            print(f"   {m:6s}        " + " ".join(f"{v:6.3f}" for v in vals), file=out)
        ax.axhspan(0.95, 1.05, color="0.9")
        ax.set_ylim(0.8, 1.2)
        ax.set_xlabel(f"{'proton' if part == 'p' else 'alpha'} energy (MeV)")
        ax.set_ylabel("engine / NIST (electronic)")
        ax.set_title(f"{key} validation ({'p' if part == 'p' else 'α'})")
        ax.legend(fontsize=7, ncol=3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "m5_stopping_validation.png"), dpi=130)
    plt.close(fig)
    print("  max |ratio-1| (p: 0.1-2 MeV, alpha: 0.4-10 MeV):", file=out)
    for k, v in worst.items():
        print(f"    {k[0]} in {k[1]:6s}: {100 * v:5.1f} %", file=out)

    # Pd has no NIST table: Ag (Z=47, same I=470 eV) scaled by Z/A is the natural proxy
    print("\n  Pd check: engine(Pd) / [NIST(Ag) * (Z/A)_Pd/(Z/A)_Ag]", file=out)
    zfac = (46 / 106.42) / (47 / 107.868)
    for part, key, n, Es in [("p", "PSTAR", 60, [0.3, 1.0, 2.0]), ("a", "ASTAR", 78, [1.0, 3.0, 5.0, 8.0])]:
        ref = np.array(d[key]["G4_Ag"]) * zfac
        for e in Es:
            print(f"    {part} {e:4.1f} MeV: {electronic(part, 'Pd', e) / np.interp(e, T[:n], ref):.3f}", file=out)

    print("\n  Bethe-Bloch(+shell+Bloch) / engine at >= 5 MeV/u (Barkas-Berger shell formula needs"
          " beta*gamma >~ 0.1; below that Bethe is not usable for Pd/Au):", file=out)
    for part, Es in [("p", [5, 10, 15, 30]), ("a", [20, 40, 80])]:
        for m, Z, A, I in [("Si", 14, 28.086, 173), ("Al", 13, 26.982, 166), ("Cu", 29, 63.546, 322),
                           ("Pd", 46, 106.42, 470), ("Au", 79, 196.967, 790)]:
            r = [bethe_bloch_element(part, Z, A, I, e) / electronic(part, m, e) for e in Es]
            print(f"    {part} {m:3s} E={Es}: " + " ".join(f"{x:.3f}" for x in r), file=out)

    print("\n  Literature cross-check (Mosier-Boss et al. EPJ AP 46 30901 (2009), LET-based): "
          "6 um Mylar stops p<0.45, t<0.55, 3He<1.40, alpha<1.45 MeV", file=out)
    for part in "pth" + "a":
        # find energy whose range equals 6 um in Mylar
        Eg, S, R = table(part, "Mylar")
        print(f"    engine: {part}: E with R=6 um in Mylar = {np.interp(6.0, R, Eg):.2f} MeV", file=out)
    print("  Common benchmarks: 5.486 MeV alpha in Si R =", f"{rng('a', 'Si', 5.486):.1f} um (SRIM ~28 um);",
          f"3.02 MeV p in Si R = {rng('p', 'Si', 3.02):.1f} um (PSTAR CSDA 3 MeV: 0.0216 g/cm2 = 92.7 um);",
          f"5.486 MeV alpha in air R = {rng('a', 'air', 5.486) / 1e4:.2f} cm (ASTAR/textbook 4.1 cm)", file=out)


PRODUCTS = [("p", 3.02, "D(d,p)T proton"), ("t", 1.01, "D(d,p)T triton"), ("h", 0.82, "D(d,n)3He helion"),
            ("a", 2.0, "alpha 2 MeV"), ("a", 4.0, "alpha 4 MeV"), ("a", 5.3, "210Po alpha 5.30"),
            ("a", 8.78, "212Po alpha 8.78"), ("a", 12.0, "alpha 12 MeV (Lipson 11-16)"),
            ("a", 16.0, "alpha 16 MeV"), ("t", 4.75, "t 4.75 (3D claim)"), ("h", 4.75, "3He 4.75 (3D claim)"),
            ("p", 14.7, "D(3He,p) secondary p")]
MATS = ["Pd", "PdD0.9", "PdO", "Au", "Ni", "Cu", "CaO", "D2O", "H2O", "Mylar", "Al", "Si", "CR39"]


def range_tables(out):
    print("\n== CSDA ranges (um); D2 gas (1 bar, 293 K) in mm", file=out)
    print("   product                       " + " ".join(f"{m:>7s}" for m in MATS) + "  D2gas(mm)", file=out)
    for part, E, lab in PRODUCTS:
        r = [rng(part, m, E) for m in MATS]
        print(f"   {lab:28s}{E:5.2f} " + " ".join(f"{x:7.2f}" for x in r)
              + f"  {rng(part, 'D2gas', E) / 1e3:8.1f}", file=out)
    print("\n== Stopping powers at birth energy (keV/um)", file=out)
    for part, E, lab in PRODUCTS[:3]:
        print(f"   {lab:28s} " + " ".join(f"{m}:{1e3 * dedx_lin(part, m, E):.1f}" for m in MATS), file=out)


if __name__ == "__main__":
    os.makedirs(FIGS, exist_ok=True)
    import sys

    class Tee:
        def __init__(self, f):
            self.f = f

        def write(self, s):
            self.f.write(s)
            sys.stdout.write(s)

        def flush(self):
            pass

    with open(os.path.join(FIGS, "m5_stopping.txt"), "w") as f:
        t = Tee(f)
        validate(t)
        range_tables(t)
