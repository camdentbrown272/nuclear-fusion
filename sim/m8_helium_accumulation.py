"""
M8 — Helium-4 channel (hypothesis H2: heat + 4He, "no radiation") for the
detector-facing permeation membrane (DFM, ADR-002).

Question: if a D+D -> 4He reaction deposits 23.85 MeV per 4He in (or near) a
10-25 um Pd membrane, how much 4He reaches a static, NEG-pumped, all-metal UHV
exit chamber, how well can it be measured against D2 and backgrounds, and what
geometry follows?

Sections (each prints to docs/models/figs/m8_output.txt):
  1  production budget (He/s per W; claimed W/cm^2 mapped onto a 2-5 cm^2 membrane)
  2  He stopping/ranges in Pd (CATIMA electronic + ZBL nuclear), validated vs R5
  3  He-in-Pd trapping physics -> trap-limited escape length lambda
  4  release fractions (exit face / entry face / retained) vs birth depth and
     birth energy; Monte Carlo verification; SRI-M4 inversion
  5  He accumulation in the reaction skin (bubble / TEM / accelerated release)
  6  detection: atoms per mbar, D2 gas load and NEG sizing, m/z-4 interferences,
     required residual D2, ion-pumping losses
  7  background budget (atoms/day), recommended vs "naive" build
  8  minimum detectable 4He production (5 sigma, 1/7/30 d) -> equivalent power
  9  geometry numbers (volume, aliquot, NEG, permeator)

Literature values quoted from memory are tagged [BK] in the markdown doc.
Run:  python3 sim/m8_helium_accumulation.py
Deterministic (fixed RNG seed).  Needs numpy, scipy, matplotlib; uses pycatima
if installed, else an embedded table produced by pycatima 1.982 (CATIMA 1.7).
"""
import os
import numpy as np
from scipy.optimize import brentq

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "docs", "models", "figs")
os.makedirs(OUT, exist_ok=True)
LOG = []


def say(s=""):
    print(s)
    LOG.append(s)


def hdr(s):
    say("")
    say("=" * 78)
    say(s)
    say("=" * 78)


RNG = np.random.default_rng(20260929)

# ------------------------------------------------------------------ constants
MeV_J = 1.602176634e-13
e_C = 1.602176634e-19
kB = 1.380649e-23
NA = 6.02214076e23
T_ROOM = 295.0
kT_eV = 8.617333e-5 * 300.0          # 300 K for kinetics
ATOMS_PER_MBAR_L = 0.1 / (kB * T_ROOM)   # 1 mbar L = 0.1 J -> molecules at 295 K
LOSCHMIDT = 2.686780e19               # cm^-3 (273.15 K, 1 atm) -> atoms per cm3 STP
X_HE_AIR = 5.24e-6                    # vol. fraction of He in air
X_AR_AIR = 9.34e-3
CMHG_PER_ATM = 76.0
Q_HE = 23.847                         # MeV per d+d -> 4He
HE_PER_J = 1.0 / (Q_HE * MeV_J)       # 4He per joule
DAY = 86400.0
M_HE, M_D2, M_AR = 4.002603, 4.028204, 39.948
RHO_PD, M_PD = 12.02, 106.42
N_PD = RHO_PD / M_PD * NA             # Pd atoms / cm3

# nominal geometry (ADR-002 platform)
L_MEM_UM = 25.0          # membrane thickness (um); 10 um checked too
A_MEM = 3.0              # membrane area (cm^2), ~19.5 mm diameter (range 2-5)
V_E = 0.5                # exit (accumulation) chamber volume, L
V_A = 0.15               # analysis manifold volume, L
V_S = 0.05               # external sample cylinder, L
F_A = V_A / (V_E + V_A)  # aliquot fraction expanded into manifold
F_S = V_S / (V_E + V_A + V_S)

# plotting palette (dataviz reference palette, light mode)
C = dict(blue="#2a78d6", orange="#eb6834", aqua="#1baf7a", yellow="#eda100",
         magenta="#e87ba4", green="#008300", violet="#4a3aa7", red="#e34948")
SERIES = [C["blue"], C["orange"], C["aqua"], C["yellow"], C["magenta"],
          C["green"], C["violet"], C["red"]]
INK, INK2, MUTED, GRID, AXIS, SURF = ("#0b0b0b", "#52514e", "#898781",
                                      "#e1e0d9", "#c3c2b7", "#fcfcfb")
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "xtick.color": MUTED,
    "ytick.color": MUTED, "text.color": INK, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.6, "grid.linestyle": "-",
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2,
    "font.size": 9.5, "axes.titlesize": 10.5, "axes.titleweight": "bold",
    "legend.frameon": False, "legend.fontsize": 8.5,
})


def savefig(fig, name):
    p = os.path.join(OUT, name)
    fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig)
    say(f"  [figure] docs/models/figs/{name}")


def eng(x, p=2):
    return f"{x:.{p}e}"


# =============================================================== 1. production
hdr("1. PRODUCTION BUDGET (d+d -> 4He, 23.85 MeV/He)")
say(f"4He per joule          = {HE_PER_J:.4e} /J  (R1: 2.62e11)")
say(f"4He per W*day          = {HE_PER_J*DAY:.4e}  = {HE_PER_J*DAY/LOSCHMIDT*1e3:.3f} uL(STP)"
    f"  (R5: 2.26e16, 0.84 uL)")
say(f"1 mJ of 24 MeV/He heat = {HE_PER_J*1e-3:.3e} 4He atoms")

# (label, P_low W, P_high W, area cm^2, source)
CLAIMS = [
    ("SRI Pd wire, ~W from ~0.94 cm2 (typ. best)", 0.5, 2.3, 0.94,
     "R1 Tab.2 (Ø1 mm; '~MJ over ~5 d' = 2.3 W)"),
    ("China Lake (Miles) rods", 0.05, 0.5, 2.5, "R1 Tab.2 (1-4 cm2)"),
    ("Letts-Cravens laser-triggered foil", 0.1, 1.0, 0.5, "R1 Tab.2"),
    ("Takahashi (Osaka) plate", 1.0, 10.0, 12.5, "R1 Tab.2"),
    ("Tohoku/Iwamura Ni-Cu films (gas, R2)", 1.0, 6.0, 12.5, "R2 (H and D alike)"),
    ("F&P 1993 boil-off (episodic outlier)", 145.0, 145.0, 0.85, "R1 Tab.2"),
]
say("")
say(f"{'claim':46s} {'W/cm2 lo':>9s} {'W/cm2 hi':>9s} | {'P on 2 cm2 (W)':>15s} "
    f"{'P on 5 cm2 (W)':>15s} | {'4He/s (3 cm2, geo-mean)':>22s}")
claim_rows = []
for lab, plo, phi, area, src in CLAIMS:
    flo, fhi = plo / area, phi / area
    gm = np.sqrt(flo * fhi)
    claim_rows.append((lab, flo, fhi, gm))
    say(f"{lab:46s} {flo:9.3f} {fhi:9.3f} | {2*flo:6.3f}-{2*fhi:<8.3f} "
        f"{5*flo:6.3f}-{5*fhi:<8.3f} | {gm*A_MEM*HE_PER_J:22.2e}")
say("")
say("Planning levels on the membrane (f_release = 1), accumulation in V_E = "
    f"{V_E} L:")
say(f"{'P (W)':>8s} {'4He/s':>10s} {'4He/day':>10s} {'dP_He per day in V_E (mbar)':>28s}")
for P in [1e-9, 1e-6, 1e-3, 1e-2, 0.1, 1.0]:
    r = P * HE_PER_J
    say(f"{P:8.0e} {r:10.2e} {r*DAY:10.2e} {r*DAY/ATOMS_PER_MBAR_L/V_E:28.2e}")

# ================================================================ 2. stopping
hdr("2. He STOPPING AND RANGES IN Pd (CATIMA electronic + ZBL nuclear)")
# embedded fallback: pycatima 1.982 (CATIMA 1.7, SRIM-85 low-energy), 4He in Pd
_E_TAB = np.array([0.004, 0.005566, 0.007746, 0.01078, 0.015002, 0.020877,
                   0.029053, 0.04043, 0.056263, 0.078297, 0.10896, 0.151631,
                   0.211013, 0.29365, 0.40865, 0.568685, 0.791394, 1.10132,
                   1.532619, 2.132824, 2.968081, 4.130442, 5.748007, 7.999044,
                   11.131633, 15.491008, 21.557603, 30.0])            # MeV
_SE_TAB = np.array([57.93, 70.56, 86.16, 105.62, 129.75, 159.19, 194.56,
                    236.46, 285.77, 343.86, 408.39, 453.65, 487.15, 520.9,
                    559.15, 595.48, 618.26, 616.56, 585.7, 530.45, 462.26,
                    392.35, 327.19, 269.1, 218.68, 176.08, 140.96, 112.49])  # MeV cm2/g
try:
    import pycatima as catima
    _PD_MAT = catima.get_material(46)

    def _catima_se(E_MeV):
        p = catima.Projectile(4.001506, 2)
        p.T(E_MeV / 4.001506)
        return catima.dedx(p, _PD_MAT)

    def _catima_range_um(E_MeV):
        p = catima.Projectile(4.001506, 2)
        p.T(E_MeV / 4.001506)
        return catima.range(p, _PD_MAT) / RHO_PD * 1e4
    HAVE_CATIMA = True
except Exception:          # pragma: no cover
    HAVE_CATIMA = False
say(f"pycatima available: {HAVE_CATIMA} (embedded table used for integration either way)")

CONV_EV15 = 1e6 * M_PD / NA * 1e15       # MeV cm2/g -> eV/(1e15 atoms/cm2)
EV15_TO_EV_NM = N_PD * 1e-7 / 1e15       # eV/(1e15 at/cm2) -> eV/nm


def se_eV15(E_keV):
    """Electronic stopping of 4He in Pd, eV/(1e15 atoms/cm2)."""
    E = np.atleast_1d(np.asarray(E_keV, float)) / 1e3
    out = np.exp(np.interp(np.log(np.clip(E, _E_TAB[0], _E_TAB[-1])),
                           np.log(_E_TAB), np.log(_SE_TAB))) * CONV_EV15
    lo = E < _E_TAB[0]
    out[lo] = _SE_TAB[0] * CONV_EV15 * np.sqrt(E[lo] / _E_TAB[0])   # velocity-proportional
    return out


def sn_zbl_eV15(E_keV, Z1=2, M1=M_HE, Z2=46, M2=M_PD):
    """ZBL universal nuclear stopping, eV/(1e15 atoms/cm2) (Ziegler-Biersack-Littmark 1985)."""
    E = np.atleast_1d(np.asarray(E_keV, float))
    zz = Z1 ** 0.23 + Z2 ** 0.23
    eps = 32.53 * M2 * E / (Z1 * Z2 * (M1 + M2) * zz)
    sn = np.where(eps <= 30,
                  np.log1p(1.1383 * eps) / (2 * (eps + 0.01321 * eps ** 0.21226 + 0.19593 * np.sqrt(eps))),
                  np.log(eps) / (2 * eps))
    return 8.462 * Z1 * Z2 * M1 * sn / ((M1 + M2) * zz)


def lindhard_se_eV15(E_keV, Z1=2, M1=M_HE, Z2=46):
    """Lindhard-Scharff (first-principles) electronic stopping, for comparison only."""
    zfac = Z1 ** (7 / 6) * Z2 / (Z1 ** (2 / 3) + Z2 ** (2 / 3)) ** 1.5
    return 19.15 * zfac * np.sqrt(np.asarray(E_keV) / (25.0 * M1))


def ion_ranges(E0_keV, n=4000):
    """Path range (nm), electronic-only range (nm), nuclear energy fraction, projected range."""
    E = np.logspace(-3, np.log10(E0_keV), n)                    # keV
    se, sn = se_eV15(E), sn_zbl_eV15(E)
    st = (se + sn) * EV15_TO_EV_NM                              # eV/nm
    dE = np.diff(E) * 1e3
    mid = lambda y: 0.5 * (y[1:] + y[:-1])
    Rpath = np.sum(dE / mid(st))
    Rel = np.sum(dE / mid(se * EV15_TO_EV_NM))
    Enuc = np.sum(dE * mid(sn / (se + sn)))
    fn = Enuc / (E0_keV * 1e3)
    Rp = Rpath / (1 + M_PD / (3 * M_HE) * fn)                   # LSS-type path->projected
    return Rpath, Rel, fn, Rp, Enuc


say("")
say("Validation vs R5 (ATIMA/pycatima) alpha ranges in Pd (um):")
say(f"{'E (MeV)':>8s} {'R5':>7s} {'table-integ (el.)':>18s} {'catima direct':>14s}")
for E, r5 in [(3.7, 6.1), (5.49, 10.1), (7.69, 16.2)]:
    _, rel, _, _, _ = ion_ranges(E * 1e3)
    cd = _catima_range_um(E) if HAVE_CATIMA else float("nan")
    say(f"{E:8.2f} {r5:7.1f} {rel/1e3:18.2f} {cd:14.2f}")

say("")
say("Low-energy He in Pd (birth energies relevant to H2):")
say(f"{'E0':>10s} {'Se':>7s} {'Sn':>6s} {'Se_LS':>6s} {'R_path':>8s} {'f_nuc':>6s} "
    f"{'R_proj':>8s} {'NRT vac':>8s}")
RP = {}
for lab, E0 in [("5 keV", 5.0), ("20 keV", 20.0), ("76 keV", 76.3), ("23.8 MeV", 23.8e3)]:
    Rpath, Rel, fn, Rp, Enuc = ion_ranges(E0)
    nrt = 0.8 * Enuc / (2 * 40.0)
    RP[lab] = Rp
    say(f"{lab:>10s} {se_eV15(E0)[0]:7.1f} {sn_zbl_eV15(E0)[0]:6.2f} "
        f"{lindhard_se_eV15(E0):6.1f} {Rpath:7.0f}nm {fn:6.3f} {Rp:7.0f}nm {nrt:8.0f}")
say("  (stopping in eV/(1e15 at/cm2); Se_LS = Lindhard-Scharff; NRT with E_d = 40 eV)")
E_R_GAMMA = Q_HE ** 2 / (2 * 3727.379) * 1e3
say(f"4He recoil from 4He+gamma: E_R = E_g^2/(2Mc^2) = {E_R_GAMMA:.1f} keV")
say(f"Max 4He recoil from e+e- pair (p_pair <= 22.8 MeV/c): "
    f"{22.8**2/(2*3727.379)*1e3:.0f} keV")

# ================================================================ 3. trapping
hdr("3. He-IN-Pd TRAPPING PHYSICS -> escape length lambda = 1/sqrt(sum k^2)")
R_CAP = 0.4e-7   # cm, capture radius ~ lattice parameter (0.389 nm)


def k2_vac(cv):
    return 4 * np.pi * R_CAP * cv * N_PD


def k2_bub(nb, rb_nm):
    return 4 * np.pi * rb_nm * 1e-7 * nb


say(f"{'microstructure':52s} {'k2 (cm^-2)':>11s} {'lambda':>10s}")
MICRO = [
    ("annealed foil: rho_d=1e8 cm-2 only", 1e8),
    ("annealed foil: rho_d=1e8 cm-2 + C_v=1e-7", 1e8 + k2_vac(1e-7)),
    ("annealed foil: rho_d=1e9 cm-2", 1e9),
    ("alpha/beta-cycled PdDx: rho_d=1e10", 1e10),
    ("alpha/beta-cycled PdDx: rho_d=1e11", 1e11),
    ("vacancy traps C_v=1e-6", k2_vac(1e-6)),
    ("vacancy traps C_v=1e-5", k2_vac(1e-5)),
    ("superabundant vacancies C_v=1e-4", k2_vac(1e-4)),
    ("superabundant vacancies C_v=1e-3", k2_vac(1e-3)),
    ("He bubbles N_b=1e18 cm-3, r_b=0.5 nm", k2_bub(1e18, 0.5)),
    ("He bubbles N_b=1e19 cm-3, r_b=0.75 nm", k2_bub(1e19, 0.75)),
]
for lab, k2 in MICRO:
    lam = 1 / np.sqrt(k2)
    say(f"{lab:52s} {k2:11.2e} {lam*1e7:8.1f}nm")
say("=> lambda spans ~1 nm (bubble-filled / vacancy-rich) to ~1 um (pristine annealed);")
say("   nominal 30 nm (cycled PdDx), bracket 3-300 nm.")

say("")
say("Interstitial He diffusivity D = D0 exp(-E_m/kT), D0 = 1e-3 cm2/s, T = 300 K:")
for Em in [0.1, 0.2, 0.4]:
    D = 1e-3 * np.exp(-Em / kT_eV)
    say(f"  E_m = {Em:.1f} eV: D = {D:.1e} cm2/s; time to diffuse 30 nm = {(3e-6)**2/D:.1e} s")
say("  -> release (if any) is prompt on all experimental time scales; trapping, not")
say("     diffusion speed, sets the release fraction (lambda is independent of D).")
say("")
say("Detrapping time tau = 1/(1e13 s^-1 * exp(-E_diss/kT)) at 300 K:")
for Ed in [0.8, 1.0, 1.2, 1.5, 2.5]:
    tau = 1 / (1e13 * np.exp(-Ed / kT_eV))
    say(f"  E_diss = {Ed:.1f} eV: tau = {tau:.1e} s = {tau/DAY:.1e} d")
say("  -> vacancy/bubble traps (>=2 eV) are permanent; dislocation traps are permanent")
say("     only if E_diss >~ 1.2 eV (weak-trap slow release is a flagged uncertainty).")
say("")
G_1W = 1.0 * HE_PER_J / (A_MEM * 100e-7)     # He/cm3/s, 1 W into 100 nm skin
D_slow = 1e-3 * np.exp(-0.4 / kT_eV)
Ci = G_1W / (D_slow * 1.1e11)
say(f"Self-trapping check (1 W into a 100 nm skin, D={D_slow:.1e}, lambda=30 nm): "
    f"G={G_1W:.1e} /cm3/s, interstitial C_He={Ci:.1e} /cm3, "
    f"k2_self={8*np.pi*R_CAP*Ci:.1e} cm-2 << 1e11 -> negligible")

# ================================================================= 4. release
hdr("4. RELEASE FRACTIONS vs BIRTH DEPTH AND BIRTH ENERGY")
L_NM = L_MEM_UM * 1e3
LAMBDA_SELF = 2.0    # nm, trap length inside the stopping cascade of a keV recoil
ETA_REFL = 0.10      # fraction of ballistically escaping recoil He that ends as gas


def _sinh_ratio(a, b):
    """sinh(a)/sinh(b) for 0<=a<=b, overflow-safe."""
    a = np.asarray(a, float)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        r = np.exp(a - b) * (-np.expm1(-2 * a)) / (-np.expm1(-2 * b))
    return np.where(a <= 0, 0.0, r)


def p_thermal(u, L, lam):
    """Thermal (zero-recoil) He born at distance u from the EXIT face of a slab of
    thickness L with trap-limited escape length lam.  Returns (exit, entry)."""
    u = np.clip(np.asarray(u, float), 0, L)
    return _sinh_ratio((L - u) / lam, L / lam), _sinh_ratio(u / lam, L / lam)


def fate_thermal(u, L, lam):
    pe, pi = p_thermal(u, L, lam)
    return dict(exit_gas=pe, wall=np.zeros_like(pe), entry=pi, retained=1 - pe - pi)


def fate_recoil(u, L, lam, Rp, strag=0.35, n=1, rng=RNG):
    """Isotropic recoil with straight path s ~ N(Rp, strag*Rp); ballistic escape if
    the path crosses a face, else thermal diffusion from the stopping point with
    lambda_eff (own-cascade traps).  u: array of birth depths from exit face."""
    u = np.repeat(np.asarray(u, float), n)
    mu = rng.uniform(-1, 1, u.size)
    s = np.clip(rng.normal(Rp, strag * Rp, u.size), 0.05 * Rp, None)
    dz = s * mu                         # >0 toward exit face
    ball_exit = dz >= u
    ball_entry = -dz >= (L - u)
    stop = ~(ball_exit | ball_entry)
    lam_eff = 1 / np.sqrt(lam ** -2 + LAMBDA_SELF ** -2)
    pe, pi = p_thermal(np.clip(u - dz, 0, L), L, lam_eff)
    exit_gas = ETA_REFL * ball_exit + stop * pe
    wall = (1 - ETA_REFL) * ball_exit
    entry = ball_entry + stop * pi
    retained = stop * (1 - pe - pi)
    return dict(exit_gas=exit_gas, wall=wall, entry=entry, retained=retained)


# --- verification 1: thermal analytic vs lattice random walk with trapping
say("Verification: 1-D lattice random walk with trapping vs analytic sinh ratio")
lam_v, a_v, Lv = 30.0, 1.5, 300.0
p_trap = a_v ** 2 / (2 * lam_v ** 2)
say(f"  lambda={lam_v} nm, step={a_v} nm, p_trap/step={p_trap:.2e}, slab {Lv} nm")
say(f"  {'u0 (nm)':>8s} {'MC exit':>9s} {'analytic':>9s} {'MC entry':>9s} {'analytic':>9s}")
max_dev = 0
for u0 in [0.0 + a_v, 15, 30, 60, 150]:
    nw = 40000
    pos = np.full(nw, float(u0))
    alive = np.ones(nw, bool)
    res = np.zeros(nw, int)      # 1 exit, 2 entry, 3 trapped
    for _ in range(200000):
        idx = np.nonzero(alive)[0]
        if idx.size == 0:
            break
        tr = RNG.random(idx.size) < p_trap
        res[idx[tr]] = 3
        alive[idx[tr]] = False
        idx = idx[~tr]
        pos[idx] += np.where(RNG.random(idx.size) < 0.5, a_v, -a_v)
        ex = pos[idx] <= 0
        en = pos[idx] >= Lv
        res[idx[ex]] = 1
        res[idx[en]] = 2
        alive[idx[ex | en]] = False
    pe, pi = p_thermal(u0, Lv, lam_v)
    mce, mci = np.mean(res == 1), np.mean(res == 2)
    max_dev = max(max_dev, abs(mce - pe))
    say(f"  {u0:8.1f} {mce:9.4f} {float(pe):9.4f} {mci:9.4f} {float(pi):9.4f}")
say(f"  max |MC-analytic| (exit) = {max_dev:.4f} (stat. sigma ~0.0025) -> PASS"
    if max_dev < 0.012 else f"  max dev {max_dev:.4f} -> CHECK")

# --- verification 2: ballistic escape vs (1-u/R)/2
say("Verification: ballistic planar escape (fixed path R) vs (1-u/R)/2:")
for uR in [0.0, 0.25, 0.5, 0.9]:
    mu = RNG.uniform(-1, 1, 400000)
    say(f"  u/R={uR:4.2f}: MC {np.mean(mu >= uR):.4f}  analytic {(1-uR)/2:.4f}")

# --- depth scenarios
SCEN = [  # label, (u_lo, u_hi) in nm from EXIT face
    ("exit surface (0-1 nm)", (0.0, 1.0)),
    ("exit skin 0-10 nm", (0.0, 10.0)),
    ("exit skin 0-100 nm", (0.0, 100.0)),
    ("exit layer 0-1 um", (0.0, 1000.0)),
    ("bulk (uniform)", (0.0, L_NM)),
    ("entry skin 0-100 nm", (L_NM - 100.0, L_NM)),
    ("entry surface (0-1 nm)", (L_NM - 1.0, L_NM)),
]
BIRTH = [("thermal", None), ("20 keV", "20 keV"), ("76 keV", "76 keV")]
LAMS = [3.0, 30.0, 300.0, 1000.0]


def scenario_fates(lo, hi, L, lam, birth, n_mc=60000):
    if birth is None:
        u = np.linspace(lo, hi, 200001) if hi - lo > 1e-9 else np.array([lo])
        f = fate_thermal(u, L, lam)
    else:
        u = RNG.uniform(lo, hi, n_mc)
        f = fate_recoil(u, L, lam, RP[birth])
    return {k: float(np.mean(v)) for k, v in f.items()}


RES = {}
for L_um in [25.0, 10.0]:
    L = L_um * 1e3
    say("")
    say(f"Membrane L = {L_um:.0f} um.  Fractions: exit->gas | recoil implanted in walls | "
        "entry (electrolyte) | retained in Pd")
    for lam in LAMS:
        say(f" lambda = {lam:g} nm")
        for bl, bkey in BIRTH:
            if L_um == 10.0 and bkey is not None:
                continue
            for sl, (lo, hi) in SCEN:
                if sl.startswith("entry"):
                    lo, hi = (L - (L_NM - lo), L - (L_NM - hi))
                if sl.startswith("bulk"):
                    hi = L
                f = scenario_fates(lo, hi, L, lam, bkey)
                RES[(L_um, lam, bl, sl)] = f
                say(f"   {bl:7s} {sl:24s} exit {f['exit_gas']:8.2e} | wall {f['wall']:7.2e} | "
                    f"entry {f['entry']:8.2e} | ret {f['retained']:7.3f}")

# SRI M4 inversion
say("")
say("SRI M4 comparison (thick wire cathode, reaction zone at the electrolyte surface):")
f_obs = 0.60
x_u = brentq(lambda x: (1 - np.exp(-x)) / x - f_obs, 1e-6, 50)
say(f"  uniform source 0..delta: f = (lambda/delta)(1-exp(-delta/lambda)) = 0.60 -> delta = {x_u:.2f} lambda")
say(f"  exponential source (scale delta): f = lambda/(lambda+delta) = 0.60 -> delta = {1/f_obs-1:.2f} lambda")
for lam in [3, 30, 300]:
    say(f"  lambda = {lam:4d} nm -> SRI source confined to delta <= {x_u*lam:.0f} nm of the wetted surface")
say(f"  bulk production in a 1 mm wire would give f ~ lambda/(r/2) ~ {3e-6/0.025:.0e} (lambda=30 nm)")
say("  => if SRI's ~60 % is real, the source is within ~lambda of the WETTED surface; in the")
say("     DFM that is the ENTRY face, whose He goes to the electrolyte, not the exit chamber.")

# --- Figure 1: release vs depth
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
u = np.logspace(-1, np.log10(L_NM), 800)
for i, lam in enumerate([3.0, 30.0, 300.0, 1000.0]):
    pe, _ = p_thermal(u, L_NM, lam)
    ax[0].plot(u, pe, color=SERIES[i], label=f"thermal birth, λ = {lam:g} nm")
for j, (bl, key) in enumerate([("20 keV recoil (λ = 30 nm)", "20 keV"),
                               ("76 keV recoil, 4He+γ (λ = 30 nm)", "76 keV")]):
    ub = np.logspace(-1, np.log10(L_NM), 120)
    fr = fate_recoil(ub, L_NM, 30.0, RP[key], n=4000)
    eg = fr["exit_gas"].reshape(ub.size, -1).mean(1)
    ax[0].plot(ub, eg, color=SERIES[4 + j], ls="-", lw=1.6, label=bl + " → gas")
ax[0].set_xscale("log")
ax[0].set_xlabel("birth depth below the EXIT face (nm)")
ax[0].set_ylabel("fraction reaching exit chamber as free gas")
ax[0].set_title("(a) He reaching the static exit chamber, L = 25 µm")
ax[0].set_ylim(-0.02, 1.02)
ax[0].legend(loc="upper right")
z = np.linspace(0, L_NM, 1500)
pe, pi = p_thermal(z, L_NM, 30.0)
ax[1].plot(z / 1e3, pe, color=SERIES[0], label="to exit face (vacuum)")
ax[1].plot(z / 1e3, pi, color=SERIES[1], label="to entry face (electrolyte)")
ax[1].plot(z / 1e3, 1 - pe - pi, color=SERIES[2], label="retained in Pd (melt)")
ax[1].set_xlabel("birth depth below the exit face (µm)")
ax[1].set_ylabel("fraction")
ax[1].set_title("(b) Fate of thermal-birth He, λ = 30 nm")
ax[1].legend(loc="center right")
ax[1].set_ylim(-0.02, 1.02)
savefig(fig, "m8_release_vs_depth.png")

# --- Figure 2: scenario summary (thermal lambda band + 76 keV)
fig, ax = plt.subplots(figsize=(8.6, 4.0))
labs = [s for s, _ in SCEN]
yy = np.arange(len(labs))
FLOOR = 1e-7
th = np.array([RES[(25.0, 30.0, "thermal", s)]["exit_gas"] for s in labs])
lo_ = np.array([min(RES[(25.0, l, "thermal", s)]["exit_gas"] for l in [3, 300]) for s in labs])
hi_ = np.array([max(RES[(25.0, l, "thermal", s)]["exit_gas"] for l in [3, 300]) for s in labs])
r76 = np.array([RES[(25.0, 30.0, "76 keV", s)]["exit_gas"] for s in labs])
cl = lambda a: np.clip(a, FLOOR, None)
ax.barh(yy + 0.2, cl(th), height=0.36, color=SERIES[0], label="thermal birth, λ = 30 nm (whisker: λ = 3–300 nm)")
ax.errorbar(cl(th), yy + 0.2, xerr=[cl(th) - cl(lo_), cl(hi_) - cl(th)], fmt="none",
            ecolor=INK2, elinewidth=1, capsize=2)
ax.barh(yy - 0.2, cl(r76), height=0.36, color=SERIES[5], label="76 keV recoil (4He+γ), λ = 30 nm")
ax.set_xscale("log")
ax.set_xlim(FLOOR, 1.5)
ax.set_yticks(yy)
ax.set_yticklabels(labs)
ax.invert_yaxis()
ax.set_xlabel(f"fraction of produced 4He reaching the exit chamber as gas (floor {FLOOR:.0e})")
ax.set_title("Release to the static exit chamber by reaction-zone location (L = 25 µm)")
ax.legend(loc="lower right")
savefig(fig, "m8_release_scenarios.png")

# ======================================================== 5. skin accumulation
hdr("5. He ACCUMULATION IN THE REACTION SKIN (if retained)")
say("Time for retained He to reach a given He/Pd in a skin of thickness delta over "
    f"A = {A_MEM} cm2:")
say(f"{'P (W)':>7s} {'delta':>7s} | {'He/Pd=1e-3 (TEM-visible bubbles)':>33s} "
    f"{'He/Pd=0.3 (accel. release/blister)':>36s}")
for P in [1e-3, 1e-2, 0.1, 1.0]:
    for d_nm in [10, 100, 1000]:
        nat = N_PD * d_nm * 1e-7 * A_MEM
        t1 = 1e-3 * nat / (P * HE_PER_J)
        t2 = 0.3 * nat / (P * HE_PER_J)
        say(f"{P:7.0e} {d_nm:5d}nm | {t1/DAY:30.2e} d {t2/DAY:33.2e} d")
say("=> at >=0.1 W from a <=100 nm exit skin, He bubbles are TEM-visible within hours-days")
say("   and accelerated/burst release (He/M ~0.3) sets in within weeks: post-run TEM/TDS")
say("   of each skin sector is an independent, spatially resolved He test.")

# ============================================================= 6. detection
hdr("6. DETECTION: atoms <-> pressure, D2 gas load, NEG, interferences")
say(f"1 mbar*L at {T_ROOM:.0f} K = {ATOMS_PER_MBAR_L:.4e} atoms")
say(f"{'MDPP (mbar)':>12s} " + " ".join(f"{'V='+str(v)+' L':>12s}" for v in [0.15, 0.5, 1.0, 2.0]))
for mdpp in [1e-14, 1e-13, 1e-12, 1e-11]:
    say(f"{mdpp:12.0e} " + " ".join(f"{mdpp*v*ATOMS_PER_MBAR_L:12.2e}" for v in [0.15, 0.5, 1.0, 2.0]))
say("  (atoms at the minimum detectable partial pressure; SEM ~1e-14-1e-13, Faraday ~1e-11)")
say("Static noble-gas MS (magnetic sector, R~600-700 at m/z 4): floor ~1e5-1e6 atoms,")
say("  line blank ~1e6-1e7 atoms [BK]; used here: sigma_MS = 3e5, blank 3e6 +- 1e6.")

# --- D2 gas load
say("")
say(f"D2 gas load from permeation through A = {A_MEM} cm2 (current-density equivalent j_perm):")
S_H2_LIST = [100.0, 400.0, 2000.0]
Q_MAX_TORRL_G = 20.0                          # St 707 / St 172 H2 embrittlement limit [BK]
Q_MAX = Q_MAX_TORRL_G * 1.33322               # mbar L / g
say(f"  NEG capacity limit q_max = {Q_MAX_TORRL_G} Torr L/g = {Q_MAX:.1f} mbar L/g (H2; D2 taken equal per molecule)")
say(f"{'j (mA/cm2)':>10s} {'D2/s':>9s} {'Q (mbar L/s)':>12s} | "
    + " ".join(f"{'P_D2 S='+str(int(s)):>13s}" for s in S_H2_LIST)
    + f" | {'NEG g: 1 d':>10s} {'7 d':>8s} {'30 d':>8s}")
JLIST = [0.1, 1.0, 10.0, 100.0]
NEG_MASS = {}
for j in JLIST:
    nd2 = j * 1e-3 * A_MEM / e_C / 2
    Q = nd2 / ATOMS_PER_MBAR_L
    ps = [Q / (0.71 * s) for s in S_H2_LIST]
    ms = [Q * t * DAY / Q_MAX for t in [1, 7, 30]]
    NEG_MASS[j] = ms
    say(f"{j:10.1f} {nd2:9.2e} {Q:12.2e} | " + " ".join(f"{p:13.2e}" for p in ps)
        + f" | {ms[0]:10.2g} {ms[1]:8.3g} {ms[2]:8.3g}")
say("  (S = NEG H2 speed in L/s, D2 speed = 0.71 S; NEG mass at the embrittlement limit, no margin)")
for T in [295.0, 473.0]:
    lp = 4.8 + 2 * np.log10(Q_MAX_TORRL_G) - 6116.0 / T
    say(f"  NEG (St 707 Sieverts fit [BK]) equilibrium H2 pressure at q = 20 Torr L/g, "
        f"T = {T:.0f} K: {10**lp*1.333:.1e} mbar")
say("  -> operate the NEG at room temperature during static windows.")
vbar4 = np.sqrt(8 * kB * T_ROOM / (np.pi * M_D2 * 1.66054e-27)) / 4 * 100 * 1e-3   # L/s/cm2
S_PERM = 0.1 * 0.5 * vbar4 * 15.0
say(f"  Pd-Ag permeator: D2 impingement conductance {vbar4:.1f} L/s/cm2; with sticking 0.1,"
    f" downstream fraction 0.5, 15 cm2 hot (300-400 C) tube -> S ~ {S_PERM:.0f} L/s, no capacity limit, He-tight")
Q_turbo = 100e-3 * A_MEM / e_C / 2 / ATOMS_PER_MBAR_L
say(f"  dynamic mode at j = 100 mA/cm2 with 250 L/s turbo: P_D2 = {Q_turbo/250:.1e} mbar")

# --- Figure 4: NEG sizing
fig, ax = plt.subplots(1, 2, figsize=(11, 4.0))
jj = np.logspace(-1, 2, 100)
Qj = jj * 1e-3 * A_MEM / e_C / 2 / ATOMS_PER_MBAR_L
for i, t in enumerate([1, 7, 30]):
    ax[0].plot(jj, Qj * t * DAY / Q_MAX, color=SERIES[i], label=f"{t}-day static window")
ax[0].axhspan(1, 100, color=SERIES[0], alpha=0.07, lw=0)
ax[0].text(0.12, 30, "practical NEG cartridge\n(1–100 g)", color=INK2, fontsize=8.5)
ax[0].set_xscale("log"); ax[0].set_yscale("log")
ax[0].set_xlabel("permeation flux, current-density equivalent (mA/cm²)")
ax[0].set_ylabel("NEG mass at H₂-embrittlement limit (g)")
ax[0].set_title(f"(a) NEG mass needed, A = {A_MEM:g} cm², no margin")
ax[0].legend(loc="upper left")
for i, s in enumerate(S_H2_LIST):
    ax[1].plot(jj, Qj / (0.71 * s), color=SERIES[i], label=f"NEG {int(s)} L/s (H₂)")
ax[1].plot(jj, Qj / (0.71 * 400 + S_PERM), color=SERIES[3], label=f"NEG 400 + Pd–Ag permeator")
ax[1].set_xscale("log"); ax[1].set_yscale("log")
ax[1].set_xlabel("permeation flux (mA/cm²)")
ax[1].set_ylabel("steady D₂ pressure in static exit chamber (mbar)")
ax[1].set_title("(b) Residual D₂ during static accumulation")
ax[1].legend(loc="upper left")
savefig(fig, "m8_neg_sizing.png")

# --- interferences
say("")
say("m/z-4 species and resolving power needed vs 4He (4.002603 u):")
SPEC = [("D2+", 4.028204), ("HT+", 4.023874), ("H2D+ (ion-molecule)", 4.029752),
        ("3HeH+ (only with 3He spike)", 4.024190), ("12C3+", 4.000000 - 3 * 0.000549 / 1)]
for lab, m in SPEC:
    say(f"  {lab:30s} m = {m:.6f}  dm = {abs(m-M_HE):.5f}  R = {M_HE/abs(m-M_HE):7.0f}")
say("  12C3+ needs >=83.5 eV (sum of C ionisation energies 11.26+24.38+47.89) plus bond energy;")
say("  O4+ needs >=181 eV -> both energetically closed at the standard 70 eV electron energy.")
T_OVER_D = 10.0 / 60 / (np.log(2) / (12.32 * 3.156e7)) / (1.107 / 20.03 * 2 * NA)
say(f"  HT: T/D in D2O at 10 dpm/mL = {T_OVER_D:.1e} -> HT+/D2+ ~ {2*T_OVER_D:.0e}, negligible")
K_IM, TAU_SRC = 2e-9, 1e-6    # cm3/s, s  [BK]
say("  Ion-molecule XY+ + H2 -> XYH+ in the ion source: ratio ~ k n tau (k=2e-9 cm3/s, tau=1 us):")
for P in [1e-11, 1e-9, 1e-7, 1e-6, 1e-5]:
    n = P * 100 / (kB * T_ROOM) * 1e-6
    say(f"    P_H2/D2 = {P:.0e} mbar: H2D+/HD+ (or H3+/H2+) ~ {K_IM*n*TAU_SRC:.1e}")
say("  -> H2D+ matters only for the H2O control at >=1e-6 mbar H2 (natural HD 3.1e-4 of H2).")
R_S = 2.5   # sensitivity ratio D2 (m/z 4) / He (m/z 4) on a QMS [BK]
say(f"  D2/He sensitivity ratio at m/z 4 used: r_S = {R_S} (He RSF 0.14-0.18, D2 ~0.35-0.45) [BK]")


def peak_transmission(R, a_tail, dm=M_D2 - M_HE, m=4.0):
    x = dm * R / m
    return np.exp(-4 * np.log(2) * x ** 2) + a_tail / max(x, 1.0) ** 2


ANALYZERS = [  # label, R, tail coefficient, delta_sub (fractional knowledge of residual)
    ("unit-resolution QMS", 4.0, 1e-3, 1.0),
    ("HR-QMS, R = 500 (e.g. fusion-grade)", 500.0, 1e-3, 0.3),
    ("magnetic-sector static MS, R = 700", 700.0, 1e-5, 0.3),
]
say("")
say(f"{'analyzer':40s} {'T(dm=0.0256)':>12s} | required P_D2 (mbar) so that residual < MDPP/3:")
say(f"{'':40s} {'':>12s} | {'MDPP=1e-14':>11s} {'MDPP=1e-13':>11s}")
REQ = {}
for lab, R, at, dsub in ANALYZERS:
    T = peak_transmission(R, at)
    req = [m / (3 * dsub * R_S * T) for m in [1e-14, 1e-13]]
    REQ[lab] = req
    say(f"{lab:40s} {T:12.2e} | {req[0]:11.1e} {req[1]:11.1e}")

# residual D2 after getter cleanup of an isolated aliquot
Q_WALL = 1e-12       # mbar L/s/cm2, D2/HD outgassing of baked 316L after D2 exposure [BK-scaled]
A_WALL_A = 300.0     # cm2 manifold wall
S_CLEAN = 7.0        # L/s D2, small RT NEG (NP10/GP50 class)
P_RES = Q_WALL * A_WALL_A / S_CLEAN
say(f"Getter-cleaned manifold: P_D2 ~ q A / S = {Q_WALL:.0e}*{A_WALL_A:.0f}/{S_CLEAN} = {P_RES:.1e} mbar "
    f"(vacuum-fired 316L: {P_RES/10:.1e}); cleanup time constant V/S = {V_A/S_CLEAN:.3f} s")
P_EXIT_1MA = 1.0e-3 * A_MEM / e_C / 2 / ATOMS_PER_MBAR_L / (0.71 * 400)
say(f"Exit chamber during permeation (1 mA/cm2, 400 L/s NEG): P_D2 = {P_EXIT_1MA:.1e} mbar")
say(f"  -> unit-res QMS directly on the exit chamber sees a D2 'He-equivalent' of "
    f"{R_S*P_EXIT_1MA:.1e} mbar = {R_S*P_EXIT_1MA*V_E*ATOMS_PER_MBAR_L:.1e} atoms")

# --- Figure 3: D2 interference
fig, ax = plt.subplots(figsize=(8.2, 4.4))
pd2 = np.logspace(-13, -4, 200)
for i, (lab, R, at, dsub) in enumerate(ANALYZERS):
    T = peak_transmission(R, at)
    ax.plot(pd2, dsub * R_S * T * pd2, color=SERIES[i], label=f"{lab} (δ_sub = {dsub:g})")
for m, ls_lab in [(1e-13, "He MDPP 1e-13 mbar (SEM, typical)"), (1e-14, "He MDPP 1e-14 mbar (SEM, best)")]:
    ax.axhline(m, color=MUTED, lw=1)
    ax.text(2e-13, m * 1.35, ls_lab, color=INK2, fontsize=8)
ax.axvspan(P_RES / 10, P_RES * 3, color=SERIES[2], alpha=0.10, lw=0)
ax.text(P_RES / 9, 3e-6, "getter-cleaned\naliquot", color=INK2, fontsize=8)
ax.axvspan(P_EXIT_1MA / 10, P_EXIT_1MA * 10, color=SERIES[1], alpha=0.10, lw=0)
ax.text(P_EXIT_1MA / 9, 3e-6, "exit chamber during\npermeation (0.1–10 mA/cm²)", color=INK2, fontsize=8)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_ylim(1e-19, 1e-4)
ax.set_xlabel("D₂ partial pressure at the analyzer (mbar)")
ax.set_ylabel("unremoved D₂ signal at 4.0026 u, He-equivalent (mbar)")
ax.set_title("m/z-4 D₂ interference vs analyzer resolution")
ax.legend(loc="lower right")
savefig(fig, "m8_d2_interference.png")

# --- ion pumping of He by analyzers / gauges in a static volume
say("")
say("He loss by ion pumping in a static volume, 1-exp(-S t / V):")
say(f"{'device':34s} {'S_He (L/s)':>11s} {'V=0.5 L, 1 d':>13s} {'7 d':>8s} {'V=0.15 L, 30 min':>17s}")
for lab, S in [("QMS ion source, low estimate", 1e-6), ("QMS ion source, high estimate", 1e-4),
               ("Bayard-Alpert gauge (on)", 3e-3), ("cold-cathode gauge (on)", 3e-2)]:
    f1 = 1 - np.exp(-S * DAY / 0.5)
    f7 = 1 - np.exp(-S * 7 * DAY / 0.5)
    fa = 1 - np.exp(-S * 1800 / V_A)
    say(f"{lab:34s} {S:11.0e} {f1:13.2%} {f7:8.2%} {fa:17.2%}")
say("  -> no ionising device may run on the accumulation volume; the QMS lives on the")
say("     manifold and sees the gas only for ~30 min per aliquot (loss corrected by spikes).")

# ============================================================ 7. backgrounds
hdr("7. BACKGROUND BUDGET (atoms/day entering the accumulation volume)")


def air_leak(L_std):
    """He in-leak (atoms/s) through leaks of total standard He leak rate L_std (mbar L/s)."""
    return L_std * X_HE_AIR * ATOMS_PER_MBAR_L


def ar_leak(L_std, regime="molecular"):
    f = np.sqrt(M_HE / M_AR) if regime == "molecular" else 1.0
    return L_std * f * X_AR_AIR * ATOMS_PER_MBAR_L


def glass_perm(K, A_cm2, d_mm):
    """Steady He permeation from air through glass (atoms/s). K in cm3STP mm/(s cm2 cmHg)."""
    return K * A_cm2 / d_mm * X_HE_AIR * CMHG_PER_ATM * LOSCHMIDT


def slab_outgas(V_cm3, S_sol, D, l_cm, t_s):
    """He outgassing (atoms/s) from an air-saturated polymer slab (thickness l, both faces)."""
    N0 = V_cm3 * S_sol * X_HE_AIR * LOSCHMIDT
    n = np.arange(0, 400)
    return N0 * np.sum(8 * D / l_cm ** 2 * np.exp(-(2 * n + 1) ** 2 * np.pi ** 2 * D * t_s / l_cm ** 2))


def radiogenic_per_g_ppm():
    lam238 = np.log(2) / (4.468e9 * 3.156e7)
    lam235 = np.log(2) / (7.04e8 * 3.156e7)
    lam232 = np.log(2) / (1.405e10 * 3.156e7)
    nU = 1e-6 / 238.03 * NA
    nTh = 1e-6 / 232.04 * NA
    return (nU * (0.992745 * lam238 * 8 + 0.0072 * lam235 * 7) * DAY,
            nTh * lam232 * 6 * DAY)


ALPHA_U, ALPHA_TH = radiogenic_per_g_ppm()
say(f"Radiogenic 4He production: {ALPHA_U:.2e} /day per g per ppm U; {ALPHA_TH:.2e} /day per g per ppm Th")
say(f"Air leak: std-He leak 1e-12 mbar L/s -> He in-leak {air_leak(1e-12):.1f}/s = {air_leak(1e-12)*DAY:.2e}/day")
say(f"  Ar tracer: Ar/He in-leak = {ar_leak(1,'molecular')/air_leak(1):.0f} (molecular) to "
    f"{ar_leak(1,'viscous')/air_leak(1):.0f} (viscous)")
ar_min_atoms = 3 * 1e-14 * V_E * ATOMS_PER_MBAR_L
L_ar_min = ar_min_atoms / (ar_leak(1.0) * DAY)
say(f"  40Ar at 3x MDPP (1e-14 mbar) in V_E after 1 day -> detects leaks down to {L_ar_min:.1e} mbar L/s (std He)")
say("Glass/sapphire viewport permeation (CF40: A = 10 cm2, d = 3 mm), steady state vs air:")
GLASS = [("fused silica", 1e-10), ("borosilicate 7740 (Pyrex)", 1e-11),
         ("borosilicate 7056/Kodial (std CF viewport)", 3e-12), ("soda-lime", 1e-13),
         ("aluminosilicate 1720", 3e-14), ("sapphire (single-crystal Al2O3)", 1e-18)]
for lab, K in GLASS:
    q = glass_perm(K, 10.0, 3.0)
    say(f"  {lab:44s} K={K:7.0e} -> {q:9.2e} /s = {q*DAY:9.2e} /day = {q/HE_PER_J:8.1e} W-eq")
say(f"  check vs R5 (Pyrex 100 cm2, 2 mm): {glass_perm(1e-11,100,2):.1e} /s (R5: ~5e6)")
say("He solubility-diffusion through the Pd membrane from air-saturated electrolyte:")
lam_th = 6.62607e-34 / np.sqrt(2 * np.pi * M_HE * 1.66054e-27 * kB * 300)
n_gas = X_HE_AIR * 101325 / (kB * 300)
for Es in [2.0, 2.9]:
    c_site = n_gas * lam_th ** 3 * np.exp(-Es / kT_eV)
    C = c_site * 2 * N_PD
    J = 1e-5 * C / (L_MEM_UM * 1e-4) * A_MEM * DAY
    say(f"  E_sol = {Es} eV: C_He = {C:.1e} /cm3 -> {J:.1e} atoms/day through the lattice (nil)")
say(f"  Dissolved He in 100 mL air-saturated electrolyte (Bunsen 0.0087): "
    f"{100*0.0087*X_HE_AIR*LOSCHMIDT:.2e} atoms (R5: 1.2e14)")


def budget(design):
    """Return list of (item, nominal atoms/day, sigma atoms/day, note)."""
    b = []
    if design == "recommended":
        L = 1e-13
        b.append(("air leaks, CF seals (1e-13 std-He, Ar-traced)", air_leak(L) * DAY, 0.25 * air_leak(L) * DAY))
        cellx = 5e-8
        b.append(("cell->exit seal leak (1e-12, 0.05 ppm He, Kr-traced)", 1e-12 * cellx * ATOMS_PER_MBAR_L * DAY,
                  0.3 * 1e-12 * cellx * ATOMS_PER_MBAR_L * DAY))
        q = glass_perm(1e-18, 10.0, 3.0) * DAY
        b.append(("sapphire CF40 viewport", q, q))
        b.append(("ceramic (alumina-brazed) feedthroughs", 1e2, 1e2))
        q = slab_outgas(0.05, 0.02, 5e-8, 5e-3, 14 * DAY) * DAY + 1e2
        b.append(("polymers: Kapton wire only, 14 d after vent", q, q))
        q = 5000 * (0.5e-3 * ALPHA_U + 1e-3 * ALPHA_TH) * 1e-3
        b.append(("316L radiogenic (5 kg, 0.5/1 ppb U/Th, rel.<=1e-3)", q, q))
        q = 50 * (ALPHA_U + ALPHA_TH) * 0.01
        b.append(("NEG radiogenic (50 g, 1 ppm U+Th, RT, rel. 1 %)", q, q))
        q = 1e9 * RHO_PD * L_MEM_UM * 1e-4 * A_MEM * 0.10 / 14
        b.append(("Pd stock 4He (1e9/g annealed, 10 % over 14 d)", q, q))
        b.append(("He standards (3-valve pipette, pumped interspace)", 1e0, 1e0))
        b.append(("He permeation through Pd lattice", 0.0, 0.0))
    else:
        L = 1e-12
        b.append(("air leaks, CF seals (1e-12 std-He, no tracer)", air_leak(L) * DAY, air_leak(L) * DAY))
        cellx = X_HE_AIR
        b.append(("cell->exit seal leak (1e-12, air-level He)", 1e-12 * cellx * ATOMS_PER_MBAR_L * DAY,
                  1e-12 * cellx * ATOMS_PER_MBAR_L * DAY))
        q = glass_perm(3e-12, 10.0, 3.0) * DAY
        b.append(("borosilicate (7056) CF40 viewport", q, 0.1 * q))
        q = 4 * glass_perm(3e-12, 0.2, 3.0) * DAY
        b.append(("4 glass-sealed SHV/BNC feedthroughs", q, 0.3 * q))
        q = slab_outgas(0.2, 0.02, 1e-8, 0.2, 14 * DAY) * DAY
        b.append(("epoxy detector mounts (0.2 cm3, 2 mm), 14 d", q, q))
        q = 5000 * (0.5e-3 * ALPHA_U + 1e-3 * ALPHA_TH) * 1e-3
        b.append(("316L radiogenic", q, q))
        q = 50 * (ALPHA_U + ALPHA_TH) * 1.0
        b.append(("NEG radiogenic, NEG heated (rel. 100 %)", q, q))
        q = 1e10 * RHO_PD * L_MEM_UM * 1e-4 * A_MEM * 0.10 / 14
        b.append(("Pd stock 4He (1e10/g as received)", q, q))
        pres = 1e9 / (0.2e-3 * ATOMS_PER_MBAR_L)
        q = pres * 1e-10 * ATOMS_PER_MBAR_L * DAY
        b.append(("He standard, 2-valve pipette, unpumped", q, q))
        b.append(("He permeation through Pd lattice", 0.0, 0.0))
    return b


BUD = {}
for design in ["recommended", "naive"]:
    b = budget(design)
    tot = sum(x[1] for x in b)
    sig = np.sqrt(sum(x[2] ** 2 for x in b))
    BUD[design] = (b, tot, sig)
    say("")
    say(f"--- {design} build ---")
    say(f"{'item':56s} {'atoms/day':>11s} {'sigma':>10s} {'W-eq (f=1)':>11s}")
    for it, v, s in b:
        say(f"{it:56s} {v:11.2e} {s:10.2e} {v/DAY/HE_PER_J:11.1e}")
    say(f"{'TOTAL':56s} {tot:11.2e} {sig:10.2e} {tot/DAY/HE_PER_J:11.1e}")

# --- Figure 5: budget
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6), sharex=True)
for k, design in enumerate(["recommended", "naive"]):
    b, tot, sig = BUD[design]
    items = [x for x in b if x[1] > 0]
    vals = np.array([max(x[1], 1e-1) for x in items])
    ax[k].barh(np.arange(len(items)), vals, color=SERIES[0] if k == 0 else SERIES[1], height=0.62)
    ax[k].set_yticks(np.arange(len(items)))
    ax[k].set_yticklabels([x[0] for x in items], fontsize=7.8)
    ax[k].invert_yaxis()
    ax[k].set_xscale("log")
    ax[k].set_xlim(1e-1, 1e11)
    ax[k].set_xlabel("4He atoms/day into the accumulation volume")
    ax[k].set_title(f"({'ab'[k]}) {design} build — total {tot:.1e}/day")
    for i, v in enumerate(vals):
        ax[k].text(v * 1.5, i, f"{v:.0e}", va="center", fontsize=7.5, color=INK2)
fig.tight_layout()
savefig(fig, "m8_background_budget.png")

# ====================================================== 8. minimum detectable
hdr("8. MINIMUM DETECTABLE 4He PRODUCTION (5 sigma) AND EQUIVALENT POWER")
MDPP_TYP = 1e-13
n_man = V_A * ATOMS_PER_MBAR_L / F_A          # total-equivalent atoms per mbar in manifold
S_EXIT = 0.71 * 400
MODES = {
    "U1 unit-res QMS on loaded exit chamber": dict(
        inst=MDPP_TYP * V_E * ATOMS_PER_MBAR_L,
        intf=1.0 * R_S * P_EXIT_1MA * V_E * ATOMS_PER_MBAR_L, bud="recommended"),
    "U2 unit-res QMS, getter-cleaned aliquot": dict(
        inst=MDPP_TYP * n_man, intf=1.0 * R_S * P_RES * n_man, bud="recommended"),
    "H1 HR-QMS (R=500), cleaned aliquot": dict(
        inst=MDPP_TYP * n_man,
        intf=0.3 * R_S * peak_transmission(500, 1e-3) * P_RES * n_man, bud="recommended"),
    "S1 external sector MS + 3He ID": dict(
        inst=np.hypot(3e5, 1e6) / F_S,
        intf=0.3 * R_S * peak_transmission(700, 1e-5) * P_RES * V_S * ATOMS_PER_MBAR_L / F_S,
        bud="recommended"),
    "N1 HR-QMS, naive build": dict(
        inst=MDPP_TYP * n_man,
        intf=0.3 * R_S * peak_transmission(500, 1e-3) * P_RES * n_man, bud="naive"),
}
say(f"Accumulation volume V_E = {V_E} L, manifold V_A = {V_A} L (aliquot fraction {F_A:.3f}),"
    f" sample cylinder {V_S} L (fraction {F_S:.3f}); QMS He MDPP = {MDPP_TYP:.0e} mbar")


def md_atoms(mode, t_s):
    m = MODES[mode]
    sigB = BUD[m["bud"]][2] / DAY * t_s
    return 5 * np.sqrt(m["inst"] ** 2 + m["intf"] ** 2 + sigB ** 2)


say(f"{'mode':42s} {'sig_inst':>9s} {'sig_D2':>9s} | {'MD power (W), f_release = 1':>30s}")
say(f"{'':42s} {'(atoms)':>9s} {'(atoms)':>9s} | {'1 d':>9s} {'7 d':>9s} {'30 d':>9s}")
MD_TABLE = {}
for mode in MODES:
    ps = [md_atoms(mode, t * DAY) / (t * DAY) / HE_PER_J for t in [1, 7, 30]]
    MD_TABLE[mode] = ps
    say(f"{mode:42s} {MODES[mode]['inst']:9.1e} {MODES[mode]['intf']:9.1e} | "
        + " ".join(f"{p:9.1e}" for p in ps))
say("  (MD rate [He/s] = MD power x 2.62e11; background sigma grows linearly with t, so")
say("   background-limited modes reach a time-independent floor = 5 sigma_B/day / 86400 s)")
for mode in MODES:
    r7 = MD_TABLE[mode][1] * HE_PER_J
    say(f"  {mode:42s}: MD rate (7 d) = {r7:.1e} He/s")

# melt channel
m_mem = RHO_PD * L_MEM_UM * 1e-4 * A_MEM
stock_sig = 0.5 * 1e9 * m_mem
melt_sig = np.hypot(stock_sig, 2e7)
say("")
say(f"Post-run melt (membrane {m_mem*1e3:.0f} mg, stock 1e9 He/g after anneal, 50 % coupon scatter,"
    f" furnace blank 1e8 +- 2e7): sigma = {melt_sig:.1e} atoms")
for t in [7, 30]:
    say(f"  run of {t:2d} d: MD retained-He production = {5*melt_sig/(t*DAY):.0f} He/s"
        f" = {5*melt_sig/(t*DAY)/HE_PER_J:.1e} W (divided by retained fraction)")

# scenario-resolved MD, recommended chain H1 at 7 d, lambda=30 (band 3-300)
say("")
say("MD power by reaction-zone scenario (H1 chain, 7 d; thermal birth; exit-chamber channel):")
say(f"{'scenario':26s} {'f_exit (l=30)':>13s} {'MD (W) l=30':>12s} {'MD band l=3..300 (W)':>24s} {'retained':>9s}")
base7 = MD_TABLE["H1 HR-QMS (R=500), cleaned aliquot"][1]
SCEN_MD = {}
for sl, _ in SCEN:
    f30 = RES[(25.0, 30.0, "thermal", sl)]["exit_gas"]
    fs = [RES[(25.0, l, "thermal", sl)]["exit_gas"] for l in [3.0, 300.0]]
    ret = RES[(25.0, 30.0, "thermal", sl)]["retained"]
    md30 = base7 / f30 if f30 > 1e-300 else np.inf
    band = sorted(base7 / max(f, 1e-300) for f in fs)
    SCEN_MD[sl] = (f30, md30, band, ret)
    say(f"{sl:26s} {f30:13.2e} {md30:12.1e} {band[0]:11.1e}-{band[1]:<11.1e} {ret:9.3f}")
CAL_R5, CAL_R6 = 5e-3, (20e-3, 60e-3)
say("")
say(f"Calorimetric floor: {CAL_R5*1e3:.0f} mW (R5 target) / {CAL_R6[0]*1e3:.0f}-{CAL_R6[1]*1e3:.0f} mW (R6 Seebeck)")
for mode in MODES:
    say(f"  {mode:42s}: calorimetry(5 mW)/MD(7 d) = {CAL_R5/MD_TABLE[mode][1]:.1e}")
say("")
say("ADR-002 arithmetic check: ADR says '~1e15 atoms per mJ'. Correct: 1 mJ = "
    f"{HE_PER_J*1e-3:.2e} He; the calorimetric floor 5-60 mW x 1 day = "
    f"{5e-3*DAY*HE_PER_J:.1e}-{60e-3*DAY*HE_PER_J:.1e} He.")

# --- Figure 6: MD vs time
fig, ax = plt.subplots(figsize=(8.6, 4.8))
tt = np.logspace(np.log10(0.1), np.log10(30), 120)
for i, mode in enumerate(MODES):
    p = [md_atoms(mode, t * DAY) / (t * DAY) / HE_PER_J for t in tt]
    ax.plot(tt, p, color=SERIES[i], label=mode)
pmelt = 5 * melt_sig / (tt * DAY) / HE_PER_J
ax.plot(tt, pmelt, color=SERIES[5], lw=1.6, label="post-run melt (retained He)")
ax.axhspan(CAL_R6[0], CAL_R6[1], color=SERIES[7], alpha=0.10, lw=0)
ax.axhline(CAL_R5, color=SERIES[7], lw=1)
ax.text(0.11, 7e-3, "calorimetry: 5 mW (R5), 20–60 mW band (R6)", color=INK2, fontsize=8)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("static accumulation time (days)")
ax.set_ylabel("5σ minimum detectable power at 23.85 MeV/He (W)")
ax.set_title("⁴He-equivalent detection floor, f_release = 1 (divide by f for real skins)")
ax.set_ylim(1e-11, 1e-1)
ax.legend(loc="upper right", fontsize=7.8)
savefig(fig, "m8_mdp_vs_time.png")

# ====================================================== 9. geometry numbers
hdr("9. GEOMETRY NUMBERS")
say("Instrument floor vs exit-chamber volume (H1 chain; background term unchanged):")
for ve in [0.3, 0.5, 1.0, 2.0]:
    fa = V_A / (ve + V_A)
    inst = MDPP_TYP * V_A * ATOMS_PER_MBAR_L / fa
    sigB = BUD["recommended"][2]
    p1 = 5 * np.hypot(inst, sigB) / DAY / HE_PER_J
    say(f"  V_E = {ve:3.1f} L: aliquot fraction {fa:.3f}, sigma_inst = {inst:.1e} atoms, MD(1 d) = {p1:.1e} W")
say("  -> volume matters little once background-limited; keep V_E <= 0.5 L for the in-situ floor")
say("     and so that a 1 cm3 split aliquot (dilution ~1/500) keeps W-level signals in QMS range.")
for P in [1e-3, 1.0]:
    n7 = P * HE_PER_J * 7 * DAY
    say(f"  {P:g} W for 7 d -> {n7:.1e} atoms -> {n7/(V_E*ATOMS_PER_MBAR_L):.1e} mbar in V_E;"
        f" 1 cm3 split -> {n7*1e-3/(V_E+1e-3)/(V_A*ATOMS_PER_MBAR_L):.1e} mbar in manifold")
say("")
say("NEG mass for a static window with 2x margin (A = 3 cm2):")
for j in JLIST:
    say(f"  j = {j:6.1f} mA/cm2: 1 d {2*NEG_MASS[j][0]:8.3g} g | 7 d {2*NEG_MASS[j][1]:8.3g} g | "
        f"30 d {2*NEG_MASS[j][2]:8.3g} g")
say("Recommended: ~50-100 g St 707/St 172 (CapaciTorr D400-D2000 class) supports 7-day windows")
say("  at <=1-2 mA/cm2; higher flux needs dynamic pumping, shorter windows or a Pd-Ag permeator.")

with open(os.path.join(OUT, "m8_output.txt"), "w") as fh:
    fh.write("\n".join(LOG) + "\n")
print("\nwrote docs/models/figs/m8_output.txt")
