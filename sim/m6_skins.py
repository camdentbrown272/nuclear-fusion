"""M6 Q5-Q6: active skins for the detector-facing membrane (C3) and the C1 cathode.

 A. Charged-particle ranges / escape depths in Pd and PdD (3.02 MeV p, 1.01 MeV t, 0.82 MeV 3He).
 B. Energy loss of p and t through each candidate skin (a)-(f); interface densities of the Iwamura
    Pd/CaO and Ni/Cu stacks; share of D sites inside the escape volume that are interface sites.
 C. Exit-face loading set by the skin: steady-state permeation with a detailed-balance desorption
    law J = 2 s Z1 f(x) (s = effective D2 sticking/recombination probability, f = D2 fugacity of PdD_x);
    required s to keep the exit face in the beta phase; PdO lifetime under D flux.

Stopping powers: pycatima 1.98 (ATIMA with SRIM-type low-energy electronic stopping), elemental
mass stopping combined by Bragg additivity (pycatima's compound input gave inconsistent results
for H-containing compounds in our tests, so we do not use it).  https://github.com/hrosiak/catima
Run: python3 sim/m6_skins.py -> docs/models/figs/m6_skins_*.png, m6_skins.txt
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pycatima as cat
from m6_common import savefig, Tee, amu, kB, NA

out = Tee("m6_skins.txt")

# ---------------- materials: (elements [(A, Z, atoms per formula unit)], density g/cm3) ----------------
# densities: CRC Handbook [BK]; PdD_x beta density from lattice expansion (a = 4.03 A, 4 Pd per cell)
def rho_pdd(x):
    a = 3.89e-8 + (4.03e-8 - 3.89e-8) * min(x, 1) / 0.7 if x < 0.7 else 4.03e-8 + 0.03e-8 * (x - 0.7) / 0.3
    return 4 * (106.42 + x * 2.014) / NA / a ** 3


MATS = {
    "Pd": ([(106.42, 46, 1)], 12.02),
    "PdD0.8": ([(106.42, 46, 1), (2.014, 1, 0.8)], None),
    "PdO": ([(106.42, 46, 1), (15.999, 8, 1)], 8.30),
    "CaO": ([(40.078, 20, 1), (15.999, 8, 1)], 3.34),
    "Au": ([(196.97, 79, 1)], 19.32),
    "Ni": ([(58.693, 28, 1)], 8.908),
    "Cu": ([(63.546, 29, 1)], 8.96),
    "C": ([(12.011, 6, 1)], 2.0),
}
MATS["PdD0.8"] = (MATS["PdD0.8"][0], rho_pdd(0.8))

PROJ = {"p 3.02 MeV": (1.00728, 1, 3.02), "t 1.01 MeV": (3.01550, 1, 1.01), "3He 0.82 MeV": (3.01493, 2, 0.82)}
_sp_cache = {}


def mass_stopping(proj, mat, E_MeV):
    """Electronic+nuclear mass stopping (MeV cm2/g) by Bragg additivity of pycatima elements."""
    A, Z, _ = proj
    els, rho = MATS[mat]
    mtot = sum(a * n for a, z, n in els)
    s = 0.0
    for a, z, n in els:
        key = (A, Z, z, round(E_MeV, 6))
        if key not in _sp_cache:
            p = cat.Projectile(A, Z, Z, E_MeV / A)
            m = cat.Material([[a, z, 1]], density=1.0)
            _sp_cache[key] = cat.dedx(p, m)
        s += (a * n / mtot) * _sp_cache[key]
    return s


def lin_stopping(proj, mat, E):
    return mass_stopping(proj, mat, E) * MATS[mat][1]     # MeV/cm


def transport(proj, stack, E0, cos_t=1.0, dz_nm=0.5):
    """Energy after crossing a stack [(mat, thickness_nm), ...] at angle acos(cos_t)."""
    E = E0
    for mat, t in stack:
        n = max(1, int(np.ceil(t / dz_nm)))
        dl = t * 1e-7 / cos_t / n
        for _ in range(n):
            E -= lin_stopping(proj, mat, E) * dl
            if E <= 0.005:
                return 0.0
    return E


def csda_range(proj, mat, E0, Emin=0.01, n=400):
    Es = np.geomspace(Emin, E0, n)
    S = np.array([lin_stopping(proj, mat, E) for E in Es])
    return np.trapezoid(1 / S, Es) * 1e4    # um


def depth_to_threshold(proj, mat, E0, Ethr, n=400):
    Es = np.geomspace(Ethr, E0, n)
    S = np.array([lin_stopping(proj, mat, E) for E in Es])
    return np.trapezoid(1 / S, Es) * 1e4    # um of path during which E stays above Ethr


# =====================================================================================
out("A. RANGES AND ESCAPE DEPTHS (CSDA, straight-line; M5 branch was not available, computed here)")
cat_check = cat.range(cat.Projectile(1.00728, 1, 1, 3.02 / 1.00728), cat.Material([[106.42, 46, 1]], density=12.02)) / 12.02 * 1e4
out(f"  cross-check pycatima direct range, p 3.02 MeV in Pd: {cat_check:.2f} um")
ranges = {}
for pn, pr in PROJ.items():
    for mat in ("Pd", "PdD0.8", "Au", "PdO", "Ni"):
        R = csda_range(pr, mat, pr[2])
        ranges[(pn, mat)] = R
    out(f"  {pn:13s}: R(Pd) = {ranges[(pn,'Pd')]:.2f} um, R(PdD0.8) = {ranges[(pn,'PdD0.8')]:.2f} um, "
        f"R(Au) = {ranges[(pn,'Au')]:.2f} um, R(PdO) = {ranges[(pn,'PdO')]:.2f} um, R(Ni) = {ranges[(pn,'Ni')]:.2f} um")
out("  Escape into a 2 pi detector hemisphere from depth z (isotropic, straight tracks) with E > E_thr:")
out("  P(z) = 0.5 (1 - z/R_thr); escape-equivalent thickness = integral P dz = R_thr/4.")
thr = {"p 3.02 MeV": [0.5, 1.0, 2.0], "t 1.01 MeV": [0.2, 0.5], "3He 0.82 MeV": [0.2]}
esc = {}
for pn, pr in PROJ.items():
    for Et in thr[pn]:
        Rt = depth_to_threshold(pr, "PdD0.8", pr[2], Et)
        esc[(pn, Et)] = Rt
        out(f"  {pn:13s} in PdD0.8, E_thr = {Et} MeV: R_thr = {Rt:.2f} um -> escape-equivalent thickness {Rt/4:.2f} um")

# =====================================================================================
out("\nB. SKINS: ENERGY LOSS OF 3.02 MeV p AND 1.01 MeV t (normal exit and 60 deg), INTERFACES")
skins = {
    "(a) bare Pd, annealed+etched (1.5 nm C/O contamination)": [("C", 1.5)],
    "(b) PdO 10 nm (thermal)": [("PdO", 10)],
    "(b') PdO 30 nm (thermal)": [("PdO", 30)],
    "(c) Au/Pd/PdO: exit face = PdO 10 nm (Au 20 nm on entry face)": [("PdO", 10)],
    "(d) Pd 40 / [CaO 2 / Pd 18]x5 (140 nm)": [("Pd", 40)] + [("CaO", 2), ("Pd", 18)] * 5,
    "(e) 6x[Cu 2 / Ni 14] (96 nm), Ni on top": [("Ni", 14), ("Cu", 2)] * 6,
    "(f) Au 20 nm": [("Au", 20)],
    "(f') Au 50 nm": [("Au", 50)],
}
res_skin = {}
out(f"  {'skin':62s} {'dE_p(0)':>8s} {'dE_p(60)':>8s} {'dE_t(0)':>8s} {'dE_t(60)':>8s} {'dE_3He(0)':>9s}   (keV)")
for sk, st in skins.items():
    vals = []
    for pn in ("p 3.02 MeV", "t 1.01 MeV", "3He 0.82 MeV"):
        pr = PROJ[pn]
        for ct in ((1.0, 0.5) if pn != "3He 0.82 MeV" else (1.0,)):
            vals.append((pr[2] - transport(pr, st, pr[2], ct)) * 1e3)
    res_skin[sk] = vals
    out(f"  {sk:62s} {vals[0]:8.2f} {vals[1]:8.2f} {vals[2]:8.2f} {vals[3]:8.2f} {vals[4]:9.2f}")
out("  Si-detector FWHM for p/t is ~15-25 keV (R5/M5): every skin <= 140 nm shifts the p peak by less than one FWHM;")
out("  t and 3He shifts from the 140 nm Pd/CaO stack are ~1-2 FWHM and must be included in the line-shape fit.")

# interface densities
nPd_surf = 1.53e15  # Pd(111) atoms/cm2 ; (100): 1.32e15   [BK, from a = 3.89 A]
out("\n  Interface densities (one monolayer of interface sites per interface, Pd(111)-like 1.5e15 cm^-2):")
for nm, nif in (("(d) Pd/CaO x5", 10), ("(e) Cu/Ni x6", 11)):
    out(f"   {nm}: {nif} buried interfaces -> {nif*nPd_surf:.1e} interface sites/cm^2")
nD = 0.8 * 4 / (4.03e-8) ** 3
Resc = esc[("p 3.02 MeV", 1.0)] * 1e-4
out(f"   D atoms within the p escape-equivalent thickness (R_thr/4 at E_thr=1 MeV, PdD0.8): {nD*Resc/4:.1e} /cm^2")
out(f"   => 10 interfaces hold {10*nPd_surf/(nD*Resc/4):.1e} of the D sites that can send a detectable proton;")
out("      an interface stack wins only if an interface site is >1e4-1e5 times more active than a bulk site.")

# =====================================================================================
out("\nC. EXIT-FACE LOADING SET BY THE SKIN (detailed balance, 300 K)")
# D2 fugacity isotherm of PdD_x at ~300 K (desorption branch), [BK] anchor points:
#   alpha: Sieverts, x_alpha,max ~0.015 at plateau; plateau ~0.04 bar (PdD, ~3x PdH);
#   beta: x = 0.65 @ ~1 bar; 0.70 @ ~10 bar; 0.80 @ ~3e2; 0.90 @ ~1e4; 0.95 @ ~1e5 bar (high-pressure
#   PdH/PdD data of Baranowski et al. and electrochemical-fugacity mapping; each anchor uncertain x10).
xs_iso = np.array([0.015, 0.60, 0.65, 0.70, 0.80, 0.90, 0.95, 1.00])
lnf_iso = np.log(np.array([0.04, 0.045, 1.0, 10.0, 3e2, 1e4, 1e5, 1e6]))


def fug(x):
    x = np.asarray(x, float)
    Ks = 0.015 / np.sqrt(0.04)
    return np.where(x < 0.015, (x / Ks) ** 2, np.exp(np.interp(x, xs_iso, lnf_iso)))


def x_of_f(f):
    f = np.asarray(f, float)
    Ks = 0.015 / np.sqrt(0.04)
    return np.where(f < 0.04, Ks * np.sqrt(f), np.interp(np.log(np.maximum(f, 1e-300)), lnf_iso, xs_iso))


mD2 = 4.028 * amu
T = 300.0
Z1 = 1e5 / np.sqrt(2 * np.pi * mD2 * kB * T) * 1e-4      # D2 cm^-2 s^-1 per bar
out(f"  D2 impingement rate at 1 bar, 300 K: Z1 = {Z1:.2e} cm^-2 s^-1; desorption flux J = 2 s Z1 f(x) D cm^-2 s^-1")
S_SKIN = {  # effective sticking (= recombination) probability, [BK] order of magnitude, with range
    "(a) clean Pd": (0.1, 1e-2, 0.5),
    "(a) Pd, contaminated (C/S/CO)": (1e-3, 1e-4, 1e-2),
    "(e) Ni top (Ni/Cu)": (1e-2, 1e-3, 0.1),
    "(d) Pd top (Pd/CaO)": (0.1, 1e-2, 0.5),
    "(f) Au 20 nm on Pd": (1e-7, 1e-10, 1e-5),
    "(b) PdO intact": (1e-6, 1e-9, 1e-4),
}
out("  Effective s by skin [BK]: clean Pd 0.01-0.5 (Conrad/Ertl/Latta, Surf. Sci. 41, 435 (1974)); Ni 1e-3-0.1;")
out("  Au(111) dissociation barrier ~1 eV -> s(bulk Au) < 1e-15, thin (<=20 nm) Au on Pd is leaky (pinholes,")
out("  intermixing): 1e-10-1e-5; intact PdO: no dissociative channel until reduced (1e-9-1e-4).")
Js = [1e14, 1e15, 1e16, 1e17]
out("  x_out (central s) with [range over the s range]; a = alpha phase, 2ph = two-phase plateau, b = beta phase")
out(f"  {'skin':32s} " + " ".join(f"{'x_out @J='+format(J,'.0e'):>14s}" for J in Js))
for sk, (s, lo, hi) in S_SKIN.items():
    row = []
    for J in Js:
        xo = x_of_f(J / (2 * s * Z1))
        xlo = x_of_f(J / (2 * hi * Z1))
        xhi = x_of_f(J / (2 * lo * Z1))
        ph = lambda v: "a" if v < 0.015 else ("2ph" if v < 0.6 else "b")
        row.append(f"{float(xo):.2g}{ph(float(xo))}[{float(xlo):.1g}-{float(xhi):.2g}]")
    out(f"  {sk:32s} " + " ".join(f"{r:>14s}" for r in row))
out("  s required to hold the exit face at x_out:")
for xo in (0.6, 0.8, 0.9, 0.95):
    out(f"   x_out = {xo}: f = {float(fug(xo)):.1e} bar; s <= " + ", ".join(f"{J/(2*Z1*float(fug(xo))):.1e} (J={J:.0e})" for J in Js))
Dd = 5e-7   # cm2/s, R1 sec 3.2
nPd = 12.02 / 106.42 * NA
for h_um in (25, 50, 100):
    Jmax = Dd * nPd * 0.9 / (h_um * 1e-4)
    out(f"  Diffusion-limited flux for x_in = 0.9, x_out = 0, h = {h_um} um: {Jmax:.1e} D cm^-2 s^-1 = {Jmax*1.602e-19:.2f} A/cm^2 equivalent")
out("  => with any bare-Pd (s >~ 1e-4) exit face the vacuum face sits at x_out < 0.02 (alpha phase) for every")
out("     achievable flux, and the membrane drains unless the entry current exceeds ~1 A/cm^2 of *absorbed* D.")
out("     The exit skin is therefore the loading valve: beta phase at the exit (x >= 0.6) needs s <~ 1e-7 at")
out("     J = 1e16, and x >= 0.8 needs s <~ 2e-11 -- achievable at best by a thick, pinhole-free Au or oxide cap.")
out("  Loading drop across the escape-equivalent depth (7 um) at J = 1e16: dx = J*7e-4/(D n_Pd) = "
    f"{1e16*7e-4/(Dd*nPd):.1e} -> the escape volume sits at x_out.")

# PdO lifetime under outgoing D
rho_pdo = 8.30
nO = rho_pdo / (106.42 + 15.999) * NA    # O per cm3
out("\n  PdO reduction by outgoing D (PdO + 2D -> Pd + D2O), lifetime tau = 2 n_O t / (eps_red J):")
for t_nm in (10, 30):
    NO = nO * t_nm * 1e-7
    row = []
    for J in (1e15, 1e16, 1e17):
        for er in (1e-4, 1e-2, 1.0):
            row.append(f"J={J:.0e},eps={er:g}: {2*NO/(er*J)/3600:.2g} h")
    out(f"   t = {t_nm} nm ({NO:.1e} O/cm^2): " + "; ".join(row))
out("  PdO is reduced by H/D at room temperature [BK]; eps_red is unknown for D arriving from the bulk. Measure")
out("  it: RGA mass 20 (D2O) from the PdO sector's vacuum volume and before/after XPS. Plan on hours-days.")

# =====================================================================================
out("\nD. ALTERNATIVE LOADING VALVE: D2 GAS BACKFILL ON THE DETECTOR SIDE")
out("  Equilibrium exit loading with D2 at pressure p (flux-free limit): x_out = x(f = p). Energy loss of products")
out("  over a 2 cm gas path (pycatima, D element, rho = p m_D2 / kT):")
for pbar in (0.01, 0.1, 1.0):
    rho_g = pbar * 1e5 / (kB * T) * mD2 * 1e-3      # g/cm3
    MATS["D2gas"] = ([(2.014, 1, 1)], rho_g)
    row = []
    for pn in ("p 3.02 MeV", "t 1.01 MeV", "3He 0.82 MeV"):
        pr = PROJ[pn]
        Eo = transport(pr, [("D2gas", 2e7)], pr[2], 1.0, dz_nm=2e5)
        row.append(f"{pn.split()[0]}: {(pr[2]-Eo)*1e3:.1f} keV")
    out(f"   p(D2) = {pbar:5.2f} bar: x_out(eq) = {float(x_of_f(pbar)):.2f}; losses over 2 cm: " + ", ".join(row))
out("  Si detectors at <=100 V bias sit below the Paschen minimum of H2/D2 (~270 V), so 10-1000 mbar is usable.")
out("  => a D2 backfill of ~0.1-1 bar holds the exit face at the plateau/beta boundary (x ~ 0.6-0.65) regardless")
out("     of the skin chemistry, at the cost of ~0.1-1 MeV of the triton/3He energy over 2 cm. Hand-off to M3/M5.")

# =====================================================================================
# figures
fig, ax = plt.subplots(1, 2, figsize=(12, 4.3))
tt = np.linspace(0, 200, 41)
for mat, ls in (("Pd", "-"), ("PdO", "--"), ("Au", ":"), ("CaO", "-."), ("Ni", (0, (3, 1, 1, 1)))):
    for pn, col in (("p 3.02 MeV", "C0"), ("t 1.01 MeV", "C3")):
        pr = PROJ[pn]
        ax[0].plot(tt, [(pr[2] - transport(pr, [(mat, t)], pr[2])) * 1e3 for t in tt], ls=ls, color=col,
                   label=f"{pn.split()[0]} in {mat}")
ax[0].axhline(20, color="k", lw=0.8)
ax[0].text(5, 21, "typical Si FWHM ~20 keV", fontsize=7)
ax[0].set(xlabel="skin thickness (nm), normal exit", ylabel="energy loss (keV)", title="Energy loss of 3.02 MeV p and 1.01 MeV t")
ax[0].legend(fontsize=6, ncol=2)
xx = np.linspace(0.001, 0.99, 400)
for sk, (s, lo, hi) in S_SKIN.items():
    Jx = 2 * s * Z1 * fug(xx)
    ax[1].loglog(Jx, xx, label=sk)
ax[1].axvspan(1e15, 1e17, color="grey", alpha=0.15, label="plausible C3 permeation flux")
ax[1].set(xlabel="permeation flux J (D cm^-2 s^-1)", ylabel="exit-face loading x_out", xlim=(1e10, 1e22), ylim=(1e-3, 1.0),
          title="Exit-face loading set by the skin (detailed balance)")
ax[1].legend(fontsize=6)
for a in ax:
    a.grid(alpha=0.3, which="both")
savefig(fig, "m6_skins.png")
out.save()
