"""
M2 common: physical constants, parameters (with sources) and the Pd-D loading model.

Everything that maps a local cathode current density i -> overpotential eta ->
effective D fugacity f -> loading x = D/Pd lives here so that the current-
distribution script (m2_current.py) and the cell/operating-window script
(m2_cell.py) use exactly the same physics.

Units: SI unless stated. Current density in A/cm^2 (electrochemistry convention),
fugacity in atm (relative to the standard state 1 atm D2), lengths in cm inside
the field solver.

Sources: see PARAMS below and docs/models/M2-electrochemistry.md section 3.
Several primary papers could not be opened from this session (publisher/
repository sites blocked by the network proxy). Values marked "recalled" are
taken from the cited paper as remembered by the author of this model and are
treated as uncertain; they are bracketed by the sensitivity analysis.
"""
import os
import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(HERE, "..", "docs", "models", "figs")
os.makedirs(FIGS, exist_ok=True)

# ------------------------------------------------------------------ constants
F = 96485.332          # C/mol
R = 8.314462           # J/mol/K
T0 = 298.15            # K
Vm_gas = 0.024465      # m^3/mol ideal gas at 25 C, 1 atm
kB_eV = 8.617333e-5

# ------------------------------------------------------------------ materials
c_Pd = 12.02 / 106.42  # mol Pd / cm^3  (rho_Pd = 12.02 g/cm^3)  = 0.1129
rho_D2O = 1104.5       # kg/m^3 at 25 C
mu_D2O = 1.095e-3      # Pa s at 25 C (D2O ~ 1.23x H2O viscosity)
D_D_Pd = 2.5e-7        # cm^2/s, D in beta-PdD near 300 K (H: ~4e-7; D/H ~ 0.6-0.7), +/- x2

# Thermochemistry of electrolysis
E_rev_D2O = 243.44e3 / (2 * F)   # V, from Delta_f G(D2O,l) = -243.44 kJ/mol  -> 1.2615 V
E_tn_D2O = 294.60e3 / (2 * F)    # V, Delta_f H(D2O,l) = -294.60 kJ/mol       -> 1.5267 V
E_tn_H2O = 285.83e3 / (2 * F)    # V, Delta_f H(H2O,l) = -285.83 kJ/mol       -> 1.4812 V
E_rev_H2O = 237.13e3 / (2 * F)   # V                                          -> 1.2288 V


def kappa_LiOD(c_molL, T_C=25.0):
    """Conductivity of LiOD in D2O (S/cm). Empirical fit
    kappa[mS/cm] = 60.56c - 14.25c^2 + 2.514cT - 0.5459c^2 T  (c mol/L, T in C).
    Source: JCMNS article 'Conductivity and Molar Conductivity of LiOD Heavy
    Water Solution' https://jcmns.org/api/v1/articles/72600-conductivity-and-molar-conductivity-of-liod-heavy-water-solution.pdf
    (formula quoted in search abstract; 25 C form 123.41c - 27.90c^2)."""
    c = c_molL
    return 1e-3 * (60.56 * c - 14.25 * c**2 + 2.514 * c * T_C - 0.5459 * c**2 * T_C)


# ====================================================================== ISOTHERM
# beta-PdD_x at high loading. Mean-field lattice gas (Lacher/Fowler-Guggenheim
# form) with the half D2 chemical potential:
#   (1/2) R T ln f = dh + W (x - 0.7) - T ds + R T ln(x/(1-x))
# dh = partial molar enthalpy of D at x=0.7 relative to 1/2 D2 (J/mol D);
# W  = effective D-D repulsion that makes the enthalpy less negative at high x.
# Two anchors at 298 K fix (dh - T ds) and W:
#   A1: x = 0.67 at f = 1 atm  (PdH ~0.70 at 1 atm [wiki/Palladium_hydride];
#       PdD slightly lower because its plateau pressure is higher)
#   A2: x = 0.90 at f = F90 atm, central F90 = 1e4, band 3e3 - 3e4 atm
#       (the widely used Baranowski-based "x=0.9 needs ~10^4 atm" figure).
# Checks (not fitted): x ~ 0.80 at ~1e2 atm (Baranowski, Filipek & Raczynski:
# electrolytic charging ~ equivalent to 80-150 bar excess D2 pressure, giving the
# commonly observed x~0.8) and x -> ~1 at GPa pressures (f ~ 1e9 atm).
DH_ISO = -17.5e3       # J/mol D  (PdD plateau enthalpy ~ -35 kJ/mol D2), +/- 3 kJ


def _iso_consts(F90=1e4):
    x1, x2 = 0.67, 0.90
    L = lambda x: np.log(x / (1 - x))
    # 0.5 ln f = a + L(x) + b (x - 0.7)   at 298 K
    b = (0.5 * np.log(F90) - L(x2) + L(x1)) / ((x2 - 0.7) - (x1 - 0.7))
    a = -L(x1) - b * (x1 - 0.7)
    W = b * R * T0                 # J/mol
    ds = (DH_ISO / T0 - a * R)     # J/mol/K   (a = dh/RT0 - ds/R)
    return a, b, W, ds


def half_ln_f(x, T=T0, F90=1e4):
    a, b, W, ds = _iso_consts(F90)
    return (DH_ISO + W * (x - 0.7)) / (R * T) - ds / R + np.log(x / (1 - x))


def x_of_f(f, T=T0, F90=1e4):
    """Loading x = D/Pd in equilibrium with D fugacity f (atm). Vectorised."""
    f = np.atleast_1d(np.asarray(f, float))
    out = np.empty_like(f)
    for k, fk in enumerate(f):
        g = 0.5 * np.log(max(fk, 1e-30))
        out[k] = brentq(lambda x: half_ln_f(x, T, F90) - g, 1e-6, 1 - 1e-12)
    return out if out.size > 1 else out[0]


def f_of_x(x, T=T0, F90=1e4):
    return np.exp(2 * half_ln_f(np.asarray(x, float), T, F90))


# ====================================================================== KINETICS
# Volmer-Tafel-Heyrovsky on Pd in alkaline D2O (Langmuir adsorption, y = theta/(1-theta)).
#   Volmer    D2O + e + *  <-> D_ad + OD-        v_V = kV/(1+y) [e^{-b phi} - (y/K) e^{(1-b)phi}]
#   Heyrovsky D2O + e + D_ad -> D2 + OD- + *     v_H = kH/(1+y) [y e^{-bH phi} - K p e^{(1-bH)phi}]
#   Tafel     2 D_ad -> D2 + 2*                  v_T = kT/(1+y)^2 [y^2 - K^2 p]
# phi = F eta / RT (eta vs the reversible D2 electrode in the same solution,
# eta<0 cathodic). Detailed balance: equilibrium at phi = -0.5 ln p with y = K sqrt(p).
# Surface-subsurface equilibrium (fast absorption step): f = y^2 / K^2 (atm).
# Steady state with a net absorption (permeation) flux J_abs (mol/cm^2/s):
#       v_V = v_H + 2 v_T + J_abs ;  i = F (v_V + v_H)
# Mechanism references: Volmer-Tafel at low i, Volmer-Heyrovsky at high i, with a
# possible loading maximum depending on symmetry factors (Zhang et al.,
# https://www.lenr-canr.org/acrobat/ZhangWSthemaximum.pdf ; HER mechanism on Pd:
# https://www.researchgate.net/publication/245146846 ;
# D vs H kinetics in alkaline solution: https://www.sciencedirect.com/science/article/pii/0022072896045895).
#
# Parameterisation by observables (so that each knob maps to a measurable):
#   i0V : Volmer exchange current density (A/cm^2)  -> sets |eta| and the Tafel line
#   cT  : Tafel-regime loading coefficient, f ~ i / cT (A/cm^2/atm) -> surface
#         recombination activity (poisons such as thiourea lower cT)
#   fmax: saturation fugacity of the Volmer-Heyrovsky regime, (kV/(kH K))^2 -> the
#         surface-state ceiling on loading
#   K   : adsorption constant (theta at eta=0 ~ K); only sets how quickly theta
#         saturates, weak influence.
SURFACES = {
    # name        i0V     cT      fmax   comment
    "poor":     dict(i0V=3e-5, cT=1e-3, fmax=3e3, beta=0.5, betaH=0.5, K=3e-3),
    "typical":  dict(i0V=3e-5, cT=1e-4, fmax=3e4, beta=0.5, betaH=0.5, K=3e-3),
    "good":     dict(i0V=3e-5, cT=1e-5, fmax=3e5, beta=0.5, betaH=0.5, K=3e-3),
    "exceptional": dict(i0V=3e-5, cT=1e-6, fmax=3e6, beta=0.5, betaH=0.5, K=3e-3),
}
# Class meaning (validation in m2_loading.py):
#  poor        untreated / cracked / recombination-active surface: x ~ 0.73-0.83
#  typical     as-received, electrolyte-cleaned:                    x ~ 0.79-0.88
#  good        vacuum-annealed + etched (0.91-0.93 routinely reported): x 0.90 at ~0.1 A/cm2
#  exceptional recombination strongly poisoned (additives); upper bound, x ~ 0.95 at ~0.1 A/cm2


class Kinetics:
    def __init__(self, i0V=3e-5, cT=1e-4, fmax=3e4, beta=0.5, betaH=0.5, K=3e-3,
                 T=T0, p=1.0):
        self.T, self.p, self.K = T, p, K
        self.beta, self.betaH = beta, betaH
        self.kV = i0V / F * (1 + K)                 # i0V = F kV (1-theta0), theta0=K/(1+K)
        self.kT = cT / (2 * F * K**2)               # v_T ~ kT K^2 f  for y<<1
        self.kH = self.kV / (K * np.sqrt(fmax))     # y_sat = kV/kH  -> f_max
        self.f_T = F / (R * T)

    def rates(self, y, phi):
        K, p = self.K, self.p
        vV = self.kV / (1 + y) * (np.exp(-self.beta * phi) - (y / K) * np.exp((1 - self.beta) * phi))
        vH = self.kH / (1 + y) * (y * np.exp(-self.betaH * phi) - K * p * np.exp((1 - self.betaH) * phi))
        vT = self.kT / (1 + y) ** 2 * (y**2 - K**2 * p)
        return vV, vH, vT

    def y_of_eta(self, eta, J_abs=0.0):
        phi = eta * self.f_T
        g = lambda ly: (lambda r: r[0] - r[1] - 2 * r[2] - J_abs)(self.rates(np.exp(ly), phi))
        lo, hi = np.log(self.K * 1e-6), np.log(1e8)
        if g(lo) * g(hi) > 0:
            return np.nan
        return np.exp(brentq(g, lo, hi, xtol=1e-12))

    def state(self, eta, J_abs=0.0):
        """Return dict with i (A/cm^2, cathodic +), f (atm), theta, fractions."""
        y = self.y_of_eta(eta, J_abs)
        vV, vH, vT = self.rates(y, eta * self.f_T)
        i = F * (vV + vH)
        return dict(eta=eta, i=i, y=y, theta=y / (1 + y), f=(y / self.K) ** 2,
                    frac_H=vH / vV if vV else np.nan, frac_T=2 * vT / vV if vV else np.nan)

    def table(self, eta_grid=None, J_abs=0.0):
        if eta_grid is None:
            eta_grid = -np.linspace(0.002, 1.2, 600)
        rows = [self.state(e, J_abs) for e in eta_grid]
        eta = np.array([r["eta"] for r in rows])
        i = np.array([r["i"] for r in rows])
        f = np.array([r["f"] for r in rows])
        th = np.array([r["theta"] for r in rows])
        ok = np.isfinite(i) & (i > 0)
        return eta[ok], i[ok], f[ok], th[ok]


class Polarization:
    """Fast interpolants i(eta), eta(i), f(i), x(i) built from a Kinetics object.
    Used as the nonlinear cathode boundary condition by the field solver."""

    def __init__(self, kin: Kinetics, F90=1e4, T=T0):
        eta, i, f, th = kin.table()
        o = np.argsort(-eta)                       # eta decreasing magnitude order -> i increasing
        self.eta, self.i, self.f = eta[o], i[o], f[o]
        self.li = np.log(self.i)
        self.x = x_of_f(self.f, T, F90)
        self.kin = kin

    def i_of_eta(self, eta):
        """eta<0 -> cathodic i>0; returns (i, di/deta)."""
        e = np.clip(eta, self.eta[-1], self.eta[0])
        # eta array is decreasing (0 -> -1.2); interpolate on reversed arrays
        li = np.interp(-e, -self.eta, self.li)
        i = np.exp(li)
        de = 1e-4
        li2 = np.interp(-(e - de), -self.eta, self.li)
        didE = -(np.exp(li2) - i) / de            # d i / d eta  (negative)
        # linear continuation beyond table (keeps Newton sane)
        return i, didE

    def eta_of_i(self, i):
        return np.interp(np.log(i), self.li, self.eta)

    def f_of_i(self, i):
        return np.exp(np.interp(np.log(i), self.li, np.log(self.f)))

    def x_of_i(self, i):
        return np.interp(np.log(i), self.li, self.x)
