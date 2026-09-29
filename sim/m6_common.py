"""M6 common utilities: paths, physical constants, optical constants.

Optical data are the refractiveindex.info database files (CC0), fetched from
https://github.com/polyanskiy/refractiveindex.info-database and stored in
sim/m6_data/ so that every script runs offline and deterministically.

Sources (all tabulated n,k):
  Pd_Johnson   P.B. Johnson & R.W. Christy, PRB 9, 5056 (1974) (0.19-1.94 um)
  Pd_Palm      K.J. Palm et al., ACS Photonics 5, 4677 (2018) (0.25-1.68 um, 200 nm film)
  Pd_Rakic-LD  A.D. Rakic et al., Appl. Opt. 37, 5271 (1998) Lorentz-Drude fit (0.25-12 um)
  Au_Johnson   Johnson & Christy, PRB 6, 4370 (1972)
  Au_Olmon-ev  R.L. Olmon et al., PRB 86, 235147 (2012) evaporated Au (0.3-25 um)
  Ni_Johnson, Cu_Johnson, Ag_Johnson  Johnson & Christy (1972/1974)
  H2O_Hale     G.M. Hale & M.R. Querry, Appl. Opt. 12, 555 (1973)
  D2O_Kedenburg S. Kedenburg et al., Opt. Mater. Express 2, 1588 (2012) (0.5-1.6 um)

PdH_x / PdD_x: no tabulated PdD dielectric function was retrievable in this
session (publisher sites blocked). We use eps_PdD = F_HYD * eps_Pd with
F_HYD = 0.75 (range 0.6-0.9), a [BK] reading of von Rottkay, Rubin & Duine,
J. Appl. Phys. 85, 408 (1999) and Palm et al. (2018), both of which report a
20-40 % reduction of |eps| on hydriding across the visible/NIR. Electronic
structure is isotope-independent to <1e-3, so PdH_x and PdD_x at equal x share
eps; the scaling is applied equally to both.
"""
import os
import re
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(HERE, "m6_data")
FIGS = os.path.join(ROOT, "docs", "models", "figs")
os.makedirs(FIGS, exist_ok=True)

# --- constants (CODATA 2018) ---
c0 = 2.99792458e8
h = 6.62607015e-34
hbar = h / (2 * np.pi)
e = 1.602176634e-19
eps0 = 8.8541878128e-12
kB = 1.380649e-23
amu = 1.66053906660e-27
me = 9.1093837015e-31
NA = 6.02214076e23

F_HYD = 0.75          # eps(PdD)/eps(Pd) scaling, [BK], range 0.6-0.9


def eV_to_um(E):
    return 1.23984198 / np.asarray(E, float)


def um_to_eV(l):
    return 1.23984198 / np.asarray(l, float)


def _read_yml(fname):
    txt = open(os.path.join(DATA, fname)).read()
    blocks = re.split(r"\n\s*- type:", txt)
    tab_nk, tab_k, formula = None, None, None
    for b in blocks[1:]:
        kind = b.split("\n")[0].strip()
        if kind.startswith("tabulated nk"):
            rows = [list(map(float, ln.split())) for ln in b.split("data: |")[1].strip().split("\n")
                    if ln.strip() and re.match(r"^\s*[0-9.eE+-]+\s", ln)]
            tab_nk = np.array(rows)
        elif kind.startswith("tabulated k"):
            rows = [list(map(float, ln.split())) for ln in b.split("data: |")[1].strip().split("\n")
                    if ln.strip() and re.match(r"^\s*[0-9.eE+-]+\s", ln)]
            tab_k = np.array(rows)
        elif kind.startswith("formula 2"):
            coef = re.search(r"coefficients:\s*([^\n]+)", b).group(1).split()
            formula = np.array(list(map(float, coef)))
    return tab_nk, tab_k, formula


_cache = {}


def nk(name, lam_um):
    """Complex refractive index n + i k of a tabulated material at lam_um (um)."""
    lam = np.asarray(lam_um, float)
    if name not in _cache:
        _cache[name] = _read_yml(name + ".yml")
    tab_nk, tab_k, formula = _cache[name]
    if tab_nk is not None:
        l, n, k = tab_nk.T
        lo, hi = l.min(), l.max()
        if np.any(lam < lo * 0.999) or np.any(lam > hi * 1.001):
            raise ValueError(f"{name}: {lam.min()}-{lam.max()} um outside {lo}-{hi}")
        return np.interp(lam, l, n) + 1j * np.interp(lam, l, k)
    # formula 2 (Sellmeier): n^2 - 1 = C1 + sum Bi l^2/(l^2 - Ci)
    C = formula
    l2 = lam ** 2
    n2 = 1 + C[0] + C[1] * l2 / (l2 - C[2]) + C[3] * l2 / (l2 - C[4])
    kk = np.interp(lam, tab_k[:, 0], tab_k[:, 1]) if tab_k is not None else 0 * lam
    return np.sqrt(n2) + 1j * kk


def eps(material, lam_um):
    """Relative permittivity (e^{-i w t} convention: Im eps > 0 = loss)."""
    lam = np.asarray(lam_um, float)
    m = material
    if m == "vac":
        return np.ones_like(lam) + 0j
    if m == "D2O":
        # Kedenburg formula valid 0.5-1.6 um; outside, use Hale-Querry H2O shifted by
        # the D2O/H2O index offset (-0.005) and H2O k (overestimates NIR loss; conservative).
        out = np.empty(lam.shape, complex)
        inr = (lam >= 0.5) & (lam <= 1.6)
        out[inr] = nk("D2O_Kedenburg", lam[inr]) ** 2
        n_h = nk("H2O_Hale", lam[~inr])
        out[~inr] = (n_h.real - 0.005 + 1j * n_h.imag) ** 2
        return out
    if m == "H2O":
        return nk("H2O_Hale", lam) ** 2
    if m in ("Pd", "Pd_JC"):
        return _metal_joined("Pd_Johnson", "Pd_Rakic-LD", lam)
    if m == "Pd_Palm":
        return nk("Pd_Palm", lam) ** 2
    if m in ("PdD", "PdH"):
        return F_HYD * eps("Pd", lam)
    if m == "Au":
        return _metal_joined("Au_Johnson", "Au_Olmon-ev", lam)
    if m == "Ni":
        return nk("Ni_Johnson", lam) ** 2
    if m == "Cu":
        return nk("Cu_Johnson", lam) ** 2
    if m == "Ag":
        return nk("Ag_Johnson", lam) ** 2
    raise KeyError(m)


def _metal_joined(primary, ir, lam):
    """Johnson & Christy up to 1.9 um, then an IR dataset (for mid-IR SPP estimates)."""
    lam = np.asarray(lam, float)
    out = np.empty(lam.shape, complex)
    a = lam <= 1.90
    if np.any(a):
        out[a] = nk(primary, lam[a]) ** 2
    if np.any(~a):
        out[~a] = nk(ir, lam[~a]) ** 2
    return out


def savefig(fig, name):
    p = os.path.join(FIGS, name)
    fig.savefig(p, dpi=130, bbox_inches="tight")
    return p


class Tee:
    """Collect text output and write it to docs/models/figs/<name>.txt."""

    def __init__(self, name):
        self.name = name
        self.lines = []

    def __call__(self, *a):
        s = " ".join(str(x) for x in a)
        print(s)
        self.lines.append(s)

    def save(self):
        with open(os.path.join(FIGS, self.name), "w") as f:
            f.write("\n".join(self.lines) + "\n")
