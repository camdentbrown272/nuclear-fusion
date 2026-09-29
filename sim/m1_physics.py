"""
M1 shared physics library: cross-sections, screening models, stopping powers,
thick-target yields, detection thresholds.

Reuses the M0 constants and the Yukawa-WKB routine from sim/m0_rate_budget.py
(imported, not copied).  Every literature value carries a tag:
  [BK]   standard value from background knowledge (web access was blocked in this
         session); the reference given is the standard source to check.
  [R3]/[R5]/[R7]... value taken from the project research digests.
  [calc] computed here.
  [guess] an engineering estimate with no direct source; flagged in the write-up.
"""
import os
import sys
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(__file__))
import m0_rate_budget as m0  # noqa: E402  (reuse M0 constants and functions)

OUT = m0.OUT
E_G = m0.E_G                     # eV, D+D Gamow energy (985.8 keV)
A_BP = m0.A_cm3_s                # cm^3/s, bound-pair constant (M0)
e2 = m0.e2_eV_A                  # eV*A
hbarc = m0.hbarc_eV_A            # eV*A
mu_c2 = m0.mu_c2                 # eV, reduced mass D+D
m_d_c2 = m0.m_d_c2               # eV
a0_A = 0.529177                  # Bohr radius, A
Ha_eV = 27.211386                # Hartree, eV
kB = 8.617333e-5                 # eV/K
DAY = 86400.0
T30 = 30 * DAY                   # s, the ADR-001 30-day run

# ------------------------------------------------------------------ Bosch-Hale
# Bosch & Hale, Nucl. Fusion 32, 611 (1992), Table IV [BK]; S in keV*mb, E_cm in keV.
# Valid 0.5-5000 keV (n) / 0.5-5000 keV (p); used below 0.5 keV as S(E)~S(0)+... (the
# polynomial is smooth to E->0; the Gamow factor carries all the energy dependence).
_BH = {
    "p": dict(BG=31.3970, A=[5.5576e4, 2.1054e2, -3.2638e-2, 1.4987e-6, 1.8181e-10],
              B=[0.0, 0.0, 0.0, 0.0]),
    "n": dict(BG=31.3970, A=[5.3701e4, 3.3027e2, -1.2706e-1, 2.9327e-5, -2.5151e-9],
              B=[0.0, 0.0, 0.0, 0.0]),
}


def S_factor(E_keV, branch):
    """Bosch-Hale S-factor, keV*b, E_cm in keV."""
    c = _BH[branch]
    A, B = c["A"], c["B"]
    E = np.asarray(E_keV, float)
    num = A[0] + E * (A[1] + E * (A[2] + E * (A[3] + E * A[4])))
    den = 1 + E * (B[0] + E * (B[1] + E * (B[2] + E * B[3])))
    return num / den * 1e-3          # keV*mb -> keV*b


def sigma_dd(E_cm_eV, Ue=0.0, branch="tot", K=1.0):
    """
    D+D cross-section (cm^2) at centre-of-mass energy E (eV) with constant-shift
    screening Ue (eV) (sigma_b(E) * exp(pi*eta*Ue/E) form written as the
    M0/R3 convention S/E * exp(-sqrt(E_G/(E+Ue)))) and an optional S-factor
    multiplier K (resonance scenario).
    """
    E = np.asarray(E_cm_eV, float)
    br = ["p", "n"] if branch == "tot" else [branch]
    S = sum(S_factor(E / 1e3, b) for b in br) * 1e3 * 1e-24      # eV*cm^2
    return K * S / E * np.exp(-np.sqrt(E_G / (E + Ue)))


# ------------------------------------------------------------------ screening
def n_au_to_cm3(n_au):
    return n_au / (a0_A * 1e-8) ** 3


def thomas_fermi_length(n_au):
    """TF screening length (A) for a free-electron gas of density n (e/bohr^3)."""
    kF = (3 * np.pi ** 2 * n_au) ** (1 / 3)
    kTF = np.sqrt(4 * kF / np.pi)                # bohr^-1
    return a0_A / kTF


def U_TF(n_au):
    """Constant-shift (r->0) screening energy of a TF-Yukawa potential, eV."""
    return e2 / thomas_fermi_length(n_au)


def _lindhard_F(x):
    x = np.atleast_1d(np.asarray(x, float))
    small = np.abs(x - 1) < 1e-9
    xs = np.where(small, 0.5, x)
    out = 0.5 + (1 - xs ** 2) / (4 * xs) * np.log(np.abs((1 + xs) / (1 - xs)))
    out = np.where(small, 0.5, out)
    return out if out.size > 1 else float(out[0])


class LindhardPotential:
    """
    Statically screened D-D potential with the full Lindhard (RPA) dielectric
    function of a free-electron gas of density n (e/bohr^3):
        V(r) = e^2/r - (e^2/r)(2/pi) int_0^inf dq sin(qr)/q * (1 - 1/eps(q))
    Tabulated on a log grid in r and interpolated (in r*V) for WKB use.
    """

    def __init__(self, n_au, r_min=1e-5, r_max=8.0, npts=260):
        self.n = n_au
        kF = (3 * np.pi ** 2 * n_au) ** (1 / 3) / a0_A           # A^-1
        kTF = 1 / thomas_fermi_length(n_au)                      # A^-1
        self.kF, self.kTF = kF, kTF
        qmax = 60 * max(kF, kTF)

        def one_minus_inv_eps(q):
            eps = 1 + (kTF / q) ** 2 * _lindhard_F(q / (2 * kF))
            return 1 - 1 / eps

        # constant shift U_L = e^2 (2/pi) int (1-1/eps) dq
        self.U0 = e2 * 2 / np.pi * quad(one_minus_inv_eps, 1e-9, qmax, limit=400,
                                         points=[2 * kF])[0]
        rs = np.geomspace(r_min, r_max, npts)
        shift = []
        for r in rs:
            f = lambda q: one_minus_inv_eps(q) / q if q > 0 else 0.0
            val = quad(f, 1e-9, qmax, weight="sin", wvar=r, limit=800)[0]
            shift.append(e2 / r * 2 / np.pi * val)
        self.r = rs
        self.rV = rs * (e2 / rs - np.array(shift))                # r*V(r), eV*A

    def V(self, r):
        return np.interp(np.log(r), np.log(self.r), self.rV) / r


def wkb_penetration(Vfun, E, Rn_A=1e-5, r_hi=8.0):
    """exp(-2 int kappa dr) from Rn to the outer turning point of V(r)=E."""
    r_tp = brentq(lambda r: Vfun(r) - E, Rn_A, r_hi)
    k = lambda r: np.sqrt(max(2 * mu_c2 * (Vfun(r) - E), 0.0)) / hbarc
    pts = np.geomspace(Rn_A, r_tp, 60)
    G = sum(quad(k, a, b, limit=200)[0] for a, b in zip(pts[:-1], pts[1:]))
    return np.exp(-2 * G), r_tp


def Ueff_from_P(P):
    """M0 definition: P = exp(-sqrt(E_G/Ue_eff))."""
    return E_G / np.log(P) ** 2


def Ueff_yukawa(Ue, E=0.04):
    """Map a Yukawa-shaped screening of (r->0) strength Ue to the thermal Ue_eff (M0)."""
    P, _ = m0.wkb_yukawa(Ue, E=E)
    return Ueff_from_P(P)


def Ue_required(rate_per_pair, rho0):
    return m0.Ue_required(rate_per_pair, rho0)


def rate_pair(Ue_eff, rho0):
    return m0.rate_bound_pair(Ue_eff, rho0)


def rho0_zero_point(hw_eV):
    """
    Relative-coordinate probability density at the pair's mean separation for two
    deuterons each in a harmonic well of quantum hw: sigma_rel = sqrt(2)*sigma_1,
    sigma_1 = sqrt(hbar/(2 m_D omega)).  Returns cm^-3.
    """
    s1 = hbarc / np.sqrt(2 * m_d_c2 * hw_eV)                   # A
    srel = np.sqrt(2) * s1
    return (2 * np.pi * srel ** 2) ** -1.5 * 1e24


# ------------------------------------------------------------------ stopping
def _Se_LS(E_keV, Z1, M1, Z2):
    """Lindhard-Scharff electronic stopping, eV cm^2/atom (velocity-proportional)."""
    Z = (Z1 ** (2 / 3) + Z2 ** (2 / 3)) ** 1.5
    return 1.91e-14 * Z1 ** (7 / 6) * Z2 / Z * np.sqrt(E_keV / (24.8 * M1))


def _Sn_ZBL(E_keV, Z1, M1, Z2, M2):
    """ZBL universal nuclear stopping, eV cm^2/atom (Ziegler-Biersack-Littmark 1985) [BK]."""
    zz = Z1 ** 0.23 + Z2 ** 0.23
    eps = 32.53 * M2 * E_keV / (Z1 * Z2 * (M1 + M2) * zz)
    if np.ndim(eps) == 0:
        eps = float(eps)
    sn = np.where(eps <= 30,
                  np.log(1 + 1.1383 * eps) / (2 * (eps + 0.01321 * eps ** 0.21226
                                                   + 0.19593 * eps ** 0.5)),
                  np.log(eps) / (2 * eps))
    return 8.462e-15 * Z1 * Z2 * M1 * sn / ((M1 + M2) * zz)


# host data: Z, M, electronic-stopping correction vs SRIM-like data [BK, R3 used 1.4 for Pd],
# host atoms per cm^3, D per host atom (x), density g/cm^3
HOSTS = {
    # name: (Z_host, M_host, corr, n_host cm^-3, x_D)
    "PdD0.7": (46, 106.4, 1.4, 6.8e22, 0.7),
    "PdD0.9": (46, 106.4, 1.4, 6.8e22, 0.9),
    "TiD2": (22, 47.9, 1.3, 5.7e22, 2.0),      # TiH2 fluorite a=4.45 A -> 4 Ti/cell [BK]
    "ZrD2": (40, 91.2, 1.3, 4.3e22, 2.0),      # [BK]
    "LiD": (3, 6.94, 1.2, 5.5e22, 1.0),        # rho 0.82 g/cm3, M 8.95 [BK]
    "D2O": (8, 16.0, 1.2, 3.33e22, 2.0),       # per O; rho 1.107 g/cm3 [BK]
}


I_EXC = {1: 19.0, 3: 40.0, 8: 95.0, 22: 233.0, 40: 393.0, 46: 470.0}   # mean excitation, eV [BK, ICRU 49]


def _Se_bethe(E_keV, Z1, M1, Z2):
    """Bethe electronic stopping (no shell/Barkas corrections), eV cm^2/atom."""
    E = np.asarray(E_keV, float)
    Mc2 = M1 * 931494.0                      # keV
    gam = 1 + E / Mc2
    b2 = 1 - 1 / gam ** 2
    arg = 2 * 511e3 * b2 * gam ** 2 / I_EXC[Z2]      # eV / eV
    return 5.099e-19 * Z1 ** 2 * Z2 / b2 * np.log(1 + arg)     # ln(1+x): positive at all E


def _Se(E_keV, Z1, M1, Z2, corr):
    """
    Electronic stopping: LS (x corr) joined to Bethe harmonically (Andersen-Ziegler style).
    Checks [calc vs BK]: 3.02 MeV p range in water 148 um (NIST PSTAR 146 um); in PdD 30 um,
    1.01 MeV t 6 um, 0.82 MeV 3He 1.8 um (R3: 33/7/2 um); R3 Table C reproduced within +20%.
    """
    lo = corr * _Se_LS(E_keV, Z1, M1, Z2)
    hi = _Se_bethe(E_keV, Z1, M1, Z2)
    return lo * hi / (lo + hi)


def stopping_ion(E_lab_keV, host, Z1=1, M1=2.014):
    """Stopping cross-section of ion (Z1, M1) per host formula unit, eV cm^2."""
    Z2, M2, corr, _, x = HOSTS[host]
    S = _Se(E_lab_keV, Z1, M1, Z2, corr) + _Sn_ZBL(E_lab_keV, Z1, M1, Z2, M2)
    # D atoms in the host (Bragg additivity); H stopping ~ LS x 1.2 at low E
    S += x * (_Se(E_lab_keV, Z1, M1, 1, 1.2) + _Sn_ZBL(E_lab_keV, Z1, M1, 1, 2.014))
    if host == "D2O":                         # the O atom is the "host", 2 D per O
        pass
    return S


def stopping_per_host(E_lab_keV, host):
    """Stopping cross-section of a deuteron per host formula unit, eV cm^2."""
    return stopping_ion(E_lab_keV, host)


def ion_range_um(E0_keV, host, Z1=1, M1=2.014, E_end_keV=1e-3):
    """CSDA range (um) of an ion from E0 down to E_end in host."""
    n = HOSTS[host][3]
    E = np.geomspace(E_end_keV, E0_keV, 400)
    return np.trapezoid(1 / (n * stopping_ion(E, host, Z1, M1)), E * 1e3) * 1e4


def thick_target_yield(E0_lab_eV, host, Ue=0.0, branch="tot", K=1.0, npts=400):
    """
    Fusions per incident deuteron of lab energy E0 stopping in `host`:
      Y = int_0^E0 x * sigma(E/2) / S_host(E) dE      (target D at rest)
    """
    x = HOSTS[host][4]
    if E0_lab_eV <= 1.0:
        return 0.0
    E = np.geomspace(1.0, E0_lab_eV, npts)
    f = x * sigma_dd(E / 2, Ue, branch, K) / stopping_per_host(E / 1e3, host)
    return np.trapezoid(f, E)


def range_um(E0_lab_eV, host, rho_host_cm3=None):
    """Projected path length (um) of a deuteron (CSDA) in host."""
    n = HOSTS[host][3]
    E = np.geomspace(1.0, E0_lab_eV, 300)
    return np.trapezoid(1 / (n * stopping_per_host(E / 1e3, host)), E) * 1e4


# ------------------------------------------------------------------ statistics
def s_min_5sigma(b, rel_sys=0.0, Z=5.0):
    """
    Minimal signal counts s for a Z-sigma discovery over expected background b
    (counts), known background (Asimov, Cowan et al. 2011) with an optional
    fractional background systematic: Z_eff^-2 = Z_asimov^-2 + (rel_sys*b/s)^2 / ... we use
    the conservative combination  s >= max(s_asimov, Z*rel_sys*b)  and require s >= 3.
    """
    def Zas(s):
        return np.sqrt(2 * ((s + b) * np.log(1 + s / b) - s))
    s_as = brentq(lambda s: Zas(s) - Z, 1e-9, 1e12) if b > 0 else 0.0
    return max(s_as, Z * rel_sys * b, 3.0)
