"""
M3 common physics: Pd–H/D thermodynamics, diffusion, exit-face kinetics, and a
1-D finite-volume solver (slab / cylinder / sphere) with a non-linear isotherm
that contains the alpha/beta miscibility gap.

Every parameter carries its source in the comment next to it.  Where the source
could not be re-read online during this session (the egress proxy blocked
osti.gov, arxiv.org, jcmns.org, lenr-canr.org, PMC), the value is flagged
[mem] = "value quoted from the cited primary source from memory; verify".
Values confirmed from search-result abstracts/snippets are flagged [web].

Units: SI unless stated.  Fugacity f in bar.  Loading x = D/Pd atom ratio.

Import only; the scripts m3_loading.py, m3_membrane.py, m3_mechanics.py,
m3_cycling.py produce all numbers and figures.
"""
import numpy as np
from scipy.linalg import solve_banded

# ------------------------------------------------------------------ constants
kB = 1.380649e-23          # J/K
kB_eV = 8.617333e-5        # eV/K
R = 8.314462               # J/mol/K
NA = 6.02214076e23
e = 1.602176634e-19
amu = 1.66053907e-27
BAR = 1e5                  # Pa
ATM = 1.01325              # bar

# Pd lattice
n_Pd = 6.80e28   # Pd atoms / m^3 ; rho = 12.02 g/cm3, M = 106.42 g/mol (CRC Handbook) -> 6.80e22 cm^-3
Omega_Pd = 8.85e-6   # m^3/mol molar volume of Pd = M/rho (same source)
N_s = 1.53e19    # surface sites / m^2 on Pd(111) (a = 3.89 A; 2/(sqrt3 a^2)) - geometry

# ------------------------------------------------------------------ isotopes
# T_c, x_c : critical point of the Pd-H(D) miscibility gap.
#   H: 566 K  [mem] Flanagan & Oates, Annu. Rev. Mater. Sci. 21 (1991) 269,
#      https://doi.org/10.1146/annurev.ms.21.080191.001413
#   D: 549 K, P_cr 3.6 MPa  [web] (search snippet, Pd-D critical point)
#      https://www.researchgate.net/figure/Pressure-composition-isotherms-for-the-Pd-D-system-Solid-lines-are-for-this-work_fig1_29457632
# alpha_max / beta_min at 298 K for H: 0.017 / 0.59  [mem] Flanagan & Oates 1991 (same URL); +-0.01
# Plateau (absorption/desorption geometric mean), van 't Hoff per mol X2:
#   H: dH = -39.0 kJ/mol H2, p_pl(298 K) = 0.015 bar  [mem] Flanagan & Oates 1991; Lasser & Klatt,
#      PRB 28 (1983) 748, https://doi.org/10.1103/PhysRevB.28.748 ; hysteresis ~x2 in p
#   D: dH = -35.5 kJ/mol D2, p_pl(298 K) = 0.045 bar (about 3x H: inverse isotope effect in solubility)
#      [mem] same sources; also Luo/Flanagan calorimetry https://doi.org/10.1016/0022-5088(91)90431-3
# Diffusion in the alpha (dilute) phase, Voelkl & Alefeld (1978) in "Hydrogen in Metals I",
#   Topics Appl. Phys. 28, 321, https://doi.org/10.1007/3540087052_52   [mem]
#   H: D0 = 2.90e-7 m2/s, Ea = 0.230 eV ;  D: D0 = 1.73e-7 m2/s, Ea = 0.206 eV (inverse isotope effect:
#   D diffuses faster than H at RT; lower Ea for D confirmed [web] Majorowski & Baranowski,
#   J. Phys. Chem. Solids 43 (1982) 1119, https://www.sciencedirect.com/science/article/abs/pii/0022369782901408)
ISO = {
    "H": dict(Tc=566.0, xc=0.27, a298=0.017, b298=0.590, dH=-39.0e3, p298=0.015,
              D0=2.90e-7, Ea=0.230, mass=2 * 1.00794 * amu),
    "D": dict(Tc=549.0, xc=0.27, a298=None, b298=None, dH=-35.5e3, p298=0.045,
              D0=1.73e-7, Ea=0.206, mass=2 * 2.01410 * amu),
}
GAP_EXP = 0.40   # exponent of the gap-width law w ~ (1 - T/Tc)^GAP_EXP (fit to alpha_max(90 C) ~ 0.035-0.045)

# beta-branch slope: 1/2 ln(f/f_pl) = kb (x - beta_min) - g ln((1-x)/(1-beta_min)),
# kb calibrated so that PdD at 298 K has f(x = 0.90) = 1e4 atm (the "10^3-10^4 atm" often quoted,
# e.g. M2 brief; McKubre et al. EPRI TR-104195) and giving f(0.95) ~ 1e5 atm, consistent with
# x = 0.95 at p ~ 0.92 GPa of D2 [web] (Baranowski high-pressure data as summarised in
# Hagelstein, J. Condensed Matter Nucl. Sci. 17 (2015) 35, https://jcmns.org/article/72362.pdf ;
# fugacity coefficient of D2 at 0.9 GPa is ~10-30).  Checks: x(PdD, 1 bar) = 0.68, x(PdH, 1 bar) = 0.69.
G_DIV = 0.10     # weak divergence as x -> 1 (octahedral sites exhausted)
F90_ATM = 1.0e4  # anchor


def _gap_w(T, Tc):
    return max(1.0 - T / Tc, 1e-6) ** GAP_EXP


def gap(T, iso="D"):
    """alpha_max, beta_min at temperature T (K). Law-of-corresponding-states scaling of the H gap."""
    p = ISO[iso]
    h = ISO["H"]
    w_ref = _gap_w(298.15, h["Tc"])
    w = _gap_w(T, p["Tc"])
    a = p["xc"] - (p["xc"] - h["a298"]) * w / w_ref
    b = p["xc"] + (h["b298"] - p["xc"]) * w / w_ref
    return a, b


def p_plateau(T, iso="D"):
    """Plateau fugacity (bar) from van 't Hoff."""
    p = ISO[iso]
    dS = p["dH"] / 298.15 - R * np.log(p["p298"])
    return np.exp(p["dH"] / (R * T) - dS / R)


def _kb298():
    a, b = gap(298.15, "D")
    rhs = 0.5 * np.log(F90_ATM * ATM / p_plateau(298.15, "D")) + G_DIV * np.log((1 - 0.9) / (1 - b))
    return rhs / (0.9 - b)


KB298 = _kb298()


def kbeta(T):
    # interaction treated as enthalpic: slope of mu/kT scales as 1/T
    return KB298 * 298.15 / T


def half_ln_f(x, T, iso="D"):
    """0.5*ln(f/1 bar) = mu/kT (per atom, rel. to 1/2 X2 at 1 bar). Maxwell-constructed (flat in gap)."""
    x = np.asarray(x, dtype=float)
    a, b = gap(T, iso)
    lp = 0.5 * np.log(p_plateau(T, iso))
    xa = np.clip(x, 1e-12, a)
    out = np.where(x < a, lp + np.log(xa / (1 - xa)) - np.log(a / (1 - a)), lp)
    xb = np.clip(x, b, 1 - 1e-9)
    beta = lp + kbeta(T) * (xb - b) - G_DIV * np.log((1 - xb) / (1 - b))
    return np.where(x > b, beta, out)


def fugacity(x, T, iso="D"):
    return np.exp(2 * half_ln_f(x, T, iso))


def x_of_f(f_bar, T, iso="D"):
    """Inverse isotherm on the stable branches (returns beta branch above plateau)."""
    from scipy.optimize import brentq
    a, b = gap(T, iso)
    if f_bar <= p_plateau(T, iso):
        return brentq(lambda x: 2 * half_ln_f(x, T, iso) - np.log(f_bar), 1e-12, a)
    return brentq(lambda x: 2 * half_ln_f(x, T, iso) - np.log(f_bar), b + 1e-12, 1 - 1e-9)


def D_alpha(T, iso="D"):
    p = ISO[iso]
    return p["D0"] * np.exp(-p["Ea"] / (kB_eV * T))


def D_chem(x, T, iso="D", blocking=True):
    """Chemical (Fick) diffusivity = D* x d(mu/kT)/dx, with D* = D_alpha (1-x) (site blocking,
    default) or D* = D_alpha (blocking=False, sensitivity case).
    alpha phase: exactly D_alpha. Gap: 0 (Maxwell construction -> Stefan problem)."""
    x = np.asarray(x, dtype=float)
    a, b = gap(T, iso)
    Da = D_alpha(T, iso)
    xb = np.clip(x, b, 1 - 1e-9)
    blk = (1 - xb) if blocking else 1.0
    beta = Da * blk * xb * (kbeta(T) + G_DIV / (1 - xb))
    return np.where(x < a, Da, np.where(x > b, beta, 0.0))


class Material:
    """Tabulated isotherm/Kirchhoff potential for fast evaluation.
    Phi(x) = n_Pd * int_0^x D_chem dx   [atoms m^-1 s^-1];  flux J = -dPhi/dz exactly."""

    def __init__(self, T=298.15, iso="D", smooth=0.002, ideal=None, blocking=True):
        self.T, self.iso = T, iso
        self.xg = np.linspace(0.0, 0.995, 39801)
        if ideal is None:
            d = D_chem(self.xg, T, iso, blocking)
        else:
            d = ideal(self.xg)
        if smooth > 0:
            dx = self.xg[1] - self.xg[0]
            k = int(3 * smooth / dx)
            kern = np.exp(-0.5 * (np.arange(-k, k + 1) * dx / smooth) ** 2)
            kern /= kern.sum()
            dp = np.pad(d, k, mode="edge")
            d = np.convolve(dp, kern, mode="valid")
        self.Dg = d
        self.Phig = n_Pd * np.concatenate([[0], np.cumsum(0.5 * (d[1:] + d[:-1]) * np.diff(self.xg))])
        self.alpha, self.beta = gap(T, iso)

    def Phi(self, x):
        return np.interp(x, self.xg, self.Phig)

    def dPhi(self, x):
        return n_Pd * np.interp(x, self.xg, self.Dg)

    def x_of_Phi(self, P):
        # Phi is non-decreasing; flat in the gap -> return the upper (beta) edge convention
        return np.interp(P, self.Phig, self.xg)


# ------------------------------------------------------------------ exit-face laws
def c_eff(x, T, iso="D"):
    """Activity form of the Pick concentration: c_eff = n alpha_max sqrt(f/f_pl) (== c in alpha phase)."""
    a, _ = gap(T, iso)
    return n_Pd * a * np.exp(half_ln_f(x, T, iso) - 0.5 * np.log(p_plateau(T, iso)))


def J_pick(x, T, kr, iso="D", Jsat=None):
    """Recombination (Pick) law in activity form J = kr c_eff^2 = kr (n a)^2 f/f_pl; optional
    surface-saturation cap J_sat (series combination)."""
    Jl = kr * c_eff(x, T, iso) ** 2
    if Jsat is None:
        return Jl
    return Jl * Jsat / (Jl + Jsat)


def J_barrier(x, T, G, iso="D"):
    """Overlayer permeation law: J = G sqrt(f[Pa]),  G = 2 NA Phi_mol / d (atoms m^-2 s^-1 Pa^-1/2)."""
    return G * np.sqrt(fugacity(x, T, iso) * BAR)


# Real-surface recombination coefficient of Pd (Pick convention):
# k_r = 1.5e-27 exp(-0.48 eV/kT) m^4/s  [web] "Asymmetric surface recombination of hydrogen on palladium
# exposed to plasma", https://www.osti.gov/etdeweb/biblio/20369213 (plasma-facing side value).
def kr_Pd_lit(T):
    return 1.5e-27 * np.exp(-0.48 / (kB_eV * T))


def kr_Pd_baskes(T, iso="D", s0=0.5):
    """Baskes upper bound for an ideally clean surface: k_r = 2 s0 / (K^2 sqrt(2 pi m kT)),
    K = Sieverts constant in atoms m^-3 Pa^-1/2.  Baskes, J. Nucl. Mater. 92 (1980) 318,
    https://www.sciencedirect.com/science/article/abs/pii/0022311580901178  [web title]"""
    a, _ = gap(T, iso)
    K = n_Pd * a / np.sqrt(p_plateau(T, iso) * BAR)
    m = ISO[iso]["mass"]
    return 2 * s0 / (K ** 2 * np.sqrt(2 * np.pi * m * kB * T))


def Jsat_desorb(T, Ed=0.75, nu=1e13):
    """Saturated-surface second-order desorption cap: J = 2 nu N_s exp(-Ed/kT) (theta -> 1).
    Ed(theta->1) ~ 0.75 eV is the upper bound implied by Iwamura's 2 sccm through a bare Pd exit
    face at 343 K (see m3_membrane.py); the low-coverage D2/Pd(111) value is 1.0 eV [web]
    (JPCC 126 (2022) 14500, https://pubs.acs.org/doi/10.1021/acs.jpcc.2c04567)."""
    return 2 * nu * N_s * np.exp(-Ed / (kB_eV * T))


# Overlayer permeabilities (mol X2 m^-1 s^-1 Pa^-1/2), H values; D ~ H/sqrt(2) (classical mass scaling)
#   Ni: 5.94e-5 exp(-51.5 kJ/RT) cm3(NTP) cm^-1 s^-1 Pa^-1/2 [web] = 2.65e-7 mol/(m s Pa^.5)
#       https://www.osti.gov/etdeweb/biblio/6006582
#   Cu: 2.8e-6 exp(-85 kJ/RT) [web] tritium-tracer near RT,
#       https://www.sciencedirect.com/science/article/abs/pii/S0925838813008761
#   Au: "lowest permeability of all metals, many orders below others" [web]; no RT number retrieved.
#       Treated as <= Cu (upper bound), Ishikawa & McLellan J. Phys. Chem. Solids 46 (1985) 445.
PERM = {
    "Ni": (2.65e-7, 51.5e3),
    "Cu": (2.8e-6, 85.0e3),
    "Au_upper": (2.8e-6, 85.0e3),
}


def G_overlayer(metal, d, T, iso="D"):
    P0, E = PERM[metal]
    P = P0 * np.exp(-E / (R * T))
    if iso == "D":
        P /= np.sqrt(2)
    return 2 * NA * P / d


# ------------------------------------------------------------------ 1-D FV solver
class FV1D:
    """Implicit (backward Euler + Newton) finite-volume solver for
         n dx/dt = r^-k d/dr ( r^k dPhi(x)/dr ),  k = 0 slab, 1 cylinder, 2 sphere.
    Domain r in [0, L].  Left (r=0): 'sym' (no flux), 'dirichlet' (value) or callable flux-in law.
    Right (r=L): 'dirichlet' value or flux-out law J(x) (+ derivative by finite difference).
    Grid geometrically refined toward boundaries that carry a flux/dirichlet condition."""

    def __init__(self, mat, L, N=200, geom=0, refine=(True, True), ratio=1.04):
        self.mat, self.L, self.k = mat, L, geom
        # build faces
        if refine[0] and refine[1]:
            h = ratio ** np.arange(N // 2)
            h = np.concatenate([h, h[::-1]])
        elif refine[1]:
            h = ratio ** np.arange(N)[::-1]
        elif refine[0]:
            h = ratio ** np.arange(N)
        else:
            h = np.ones(N)
        faces = np.concatenate([[0], np.cumsum(h)])
        faces *= L / faces[-1]
        self.faces = faces
        self.rc = 0.5 * (faces[1:] + faces[:-1])
        self.dr = np.diff(faces)
        k = geom
        self.V = (faces[1:] ** (k + 1) - faces[:-1] ** (k + 1)) / (k + 1)   # per unit (angle) measure
        self.A = faces ** k
        self.N = len(self.rc)

    def mean(self, x):
        return np.sum(x * self.V) / np.sum(self.V)

    def step(self, x_old, dt, left, right, tol=1e-10, maxit=30):
        """left/right: ('sym',) | ('dir', value) | ('flux', f) with f(x_boundary_cell)-> outward flux."""
        mat, N = self.mat, self.N
        x = x_old.copy()
        dcen = np.diff(self.rc)                      # centre spacings, N-1
        Ai = self.A[1:-1]
        for it in range(maxit):
            P = mat.Phi(x)
            dP = mat.dPhi(x)
            # interior fluxes F_{i+1/2} (outward +r direction) = -A (P_{i+1}-P_i)/dcen
            Fint = -Ai * (P[1:] - P[:-1]) / dcen
            res = n_Pd * self.V * (x - x_old) / dt
            res[:-1] += Fint
            res[1:] -= Fint
            main = n_Pd * self.V / dt
            up = np.zeros(N)   # coefficient of x_{i+1} in eq i
            lo = np.zeros(N)   # coefficient of x_{i-1} in eq i
            g = Ai / dcen
            main[:-1] += g * dP[:-1]
            up[:-1] = -g * dP[1:]
            main[1:] += g * dP[1:]
            lo[1:] = -g * dP[:-1]
            # right boundary
            AL = self.A[-1]
            if right[0] == "dir":
                gR = AL / (0.5 * self.dr[-1])
                FR = -gR * (mat.Phi(right[1]) - P[-1])
                res[-1] += FR
                main[-1] += gR * dP[-1]
            elif right[0] == "flux":
                f = right[1]
                J = f(x[-1])
                h = 1e-7
                dJ = (f(x[-1] + h) - f(x[-1] - h)) / (2 * h)
                res[-1] += AL * J
                main[-1] += AL * dJ
            # left boundary (r = 0 face)
            A0 = self.A[0]
            if left[0] == "dir":
                gL = A0 / (0.5 * self.dr[0])
                FL = -gL * (P[0] - mat.Phi(left[1]))     # flux in +r direction through face 0
                res[0] -= FL
                main[0] += gL * dP[0]
            elif left[0] == "flux":   # outward (-r) flux law
                f = left[1]
                J = f(x[0])
                h = 1e-7
                dJ = (f(x[0] + h) - f(x[0] - h)) / (2 * h)
                res[0] += A0 * J
                main[0] += A0 * dJ
            ab = np.zeros((3, N))
            ab[0, 1:] = up[:-1]
            ab[1] = main
            ab[2, :-1] = lo[1:]
            dx = solve_banded((1, 1), ab, -res)
            # damping: limit step to 0.05 in x
            m = np.max(np.abs(dx))
            if m > 0.05:
                dx *= 0.05 / m
            x = np.clip(x + dx, 0.0, 0.994)
            if m < tol:
                return x, True
        return x, m < 1e-7

    def run(self, x0, t_end, left, right, dt0=None, dtmax=None, stop=None, record=None, t0=0.0):
        """Adaptive stepping. left/right may be callables of t returning BC tuples.
        stop(t, x) -> True ends run. record(t, x) called each accepted step."""
        x = x0.copy()
        t = t0
        dt = dt0 or 1e-6 * t_end
        dtmax = dtmax or t_end / 50
        while t < t_end:
            dt = min(dt, t_end - t, dtmax)
            lb = left(t + dt) if callable(left) else left
            rb = right(t + dt) if callable(right) else right
            xn, ok = self.step(x, dt, lb, rb)
            if not ok:
                dt *= 0.3
                if dt < 1e-14 * t_end:
                    raise RuntimeError("time step collapse")
                continue
            chg = np.max(np.abs(xn - x))
            x, t = xn, t + dt
            if record:
                record(t, x)
            if stop and stop(t, x):
                break
            dt *= 1.5 if chg < 0.02 else (1.0 if chg < 0.05 else 0.6)
        return t, x
