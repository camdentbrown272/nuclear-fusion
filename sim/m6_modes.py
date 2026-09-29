"""M6 Q4: phonon and mechanical mode engineering.

 1. PdD / PdH optical phonons: rocksalt Born-von Karman (central NN Pd-Pd, Pd-D, D-D force constants)
    fitted to [BK] neutron-scattering anchors, dispersion along G-X-W-L-G, zero-group-velocity points,
    one- and two-phonon DOS; comparison with the Letts beat lines (8.3, 15.3, 20.4 THz).
 2. Symmetry and momentum selection rules for driving these modes with a two-laser beat.
 3. Energy cost of a coherent optical-phonon amplitude in the laser-accessible skin.
 4. Nanoparticle breathing (Lamb l=0) modes vs diameter.
 5. Membrane (C3) and wire (C1) flexural / longitudinal / thickness / radial modes, fluid loading, Q;
    pressure-load check of an unsupported membrane.
 6. Ultrasonic / megasonic drive: strain vs intensity, Hagelstein-type MHz proposal.
Run: python3 sim/m6_modes.py -> docs/models/figs/m6_phonons.png, m6_mech.png, m6_modes.txt
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq, least_squares
from scipy.special import j0, j1
from m6_common import savefig, Tee, amu, hbar, kB, e, h

out = Tee("m6_modes.txt")
THz = 1e12

# =====================================================================================
# 1. Rocksalt BvK model
# =====================================================================================
a0 = 4.03e-10            # PdD beta-phase lattice constant (m) [BK; R4 table: 4.03 A]
mPd, mD, mH = 106.42 * amu, 2.01410 * amu, 1.00794 * amu
# neighbour vectors (units of a0)
fcc_nn = np.array([[s1 * 0.5, s2 * 0.5, 0] for s1 in (1, -1) for s2 in (1, -1)] +
                  [[s1 * 0.5, 0, s2 * 0.5] for s1 in (1, -1) for s2 in (1, -1)] +
                  [[0, s1 * 0.5, s2 * 0.5] for s1 in (1, -1) for s2 in (1, -1)])
oct_nn = np.array([[0.5, 0, 0], [-0.5, 0, 0], [0, 0.5, 0], [0, -0.5, 0], [0, 0, 0.5], [0, 0, -0.5]])


def dynmat(q, phi, mh):
    """q in units of 2pi/a0 (cartesian). phi = (Pd-Pd, Pd-X, X-X) central force constants (N/m)."""
    f11, f12, f22 = phi
    D = np.zeros((6, 6), complex)

    def blk(vecs, f):
        S = np.zeros((3, 3), complex)
        S0 = np.zeros((3, 3))
        for r in vecs:
            rh = r / np.linalg.norm(r)
            P = np.outer(rh, rh) * f
            S += P * np.exp(2j * np.pi * np.dot(q, r))
            S0 += P
        return S, S0
    S11, S11_0 = blk(fcc_nn, f11)
    S22, S22_0 = blk(fcc_nn, f22)
    S12, S12_0 = blk(oct_nn, f12)
    D[:3, :3] = (S11_0 - S11 + S12_0) / mPd
    D[3:, 3:] = (S22_0 - S22 + S12_0) / mh
    D[:3, 3:] = -S12 / np.sqrt(mPd * mh)
    D[3:, :3] = D[:3, 3:].conj().T
    w2 = np.linalg.eigvalsh(D)
    return np.sqrt(np.clip(w2, 0, None)) / (2 * np.pi)   # Hz


G = np.array([0, 0, 0.]); X = np.array([1, 0, 0.]); W = np.array([1, 0.5, 0]); L = np.array([0.5, 0.5, 0.5])
K = np.array([0.75, 0.75, 0])
path = [G, X, W, L, G, K]
labels = ["G", "X", "W", "L", "G", "K"]

# [BK] anchors (with ranges carried into the text):
#   Pd LA(X) ~ 6.7 THz (Miiller & Brockhouse, Can. J. Phys. 49, 704 (1971))
#   PdD_0.63 optical: bottom ~ 8.0 THz near G, top ~ 11.5 THz at the zone boundary
#   (Rowe, Rush, de Graaf & Ferguson, PRL 29, 1250 (1972); INS DOS peak ~37 meV = 8.9 THz,
#    Ross et al., J. Phys.: Condens. Matter 10, 3219 (1998)).  Uncertainty: +/-0.7 THz on each.
anchors = dict(PdLAX=6.7e12, optG=8.0e12, optTop=11.5e12)


def fit(anch):
    def res(p):
        phi = np.exp(p)
        fX = dynmat(X, (phi[0], phi[1], phi[2]), mD)
        fG = dynmat(G, (phi[0], phi[1], phi[2]), mD)
        fL = dynmat(L, (phi[0], phi[1], phi[2]), mD)
        return [(fX[2] - anch["PdLAX"]) / 1e12, (fG[3] - anch["optG"]) / 1e12,
                (max(fX[5], fL[5]) - anch["optTop"]) / 1e12]
    r = least_squares(res, np.log([30., 20., 5.]))
    return np.exp(r.x), r.fun


phi, resid = fit(anchors)
out("1. ROCKSALT BvK FIT (central forces, N/m): Pd-Pd = %.1f, Pd-D = %.1f, D-D = %.2f; residuals (THz) %s"
    % (phi[0], phi[1], phi[2], np.round(resid, 3)))


def path_freqs(phi, mh, npts=60):
    qs, xs, x = [], [], 0.0
    ticks = [0.0]
    for a, b in zip(path[:-1], path[1:]):
        seg = np.linalg.norm(b - a)
        for t in np.linspace(0, 1, npts, endpoint=False):
            qs.append(a + t * (b - a))
            xs.append(x + t * seg)
        x += seg
        ticks.append(x)
    qs.append(path[-1]); xs.append(x)
    return np.array(xs), np.array([dynmat(q, phi, mh) for q in qs]), ticks


def dos(phi, mh, n=24, bins=None):
    g = (np.arange(n) + 0.5) / n
    fs = []
    b1, b2, b3 = np.array([-1, 1, 1.]), np.array([1, -1, 1.]), np.array([1, 1, -1.])   # fcc reciprocal (2pi/a)
    for i in g:
        for j in g:
            for k in g:
                fs.append(dynmat(i * b1 + j * b2 + k * b3, phi, mh))
    fs = np.array(fs)
    return fs


# isotope: harmonic mass scaling vs measured anharmonic ratio for the optical peak
# PdH optical peak ~56 meV (13.5 THz) vs PdD ~37 meV (8.9 THz): ratio ~1.5 > sqrt2 [BK: Ross 1998]
xsD, fD, ticks = path_freqs(phi, mD)
xsH, fH, _ = path_freqs(phi, mH)
anh = 1.50 / np.sqrt(2)              # extra H/D ratio from anharmonicity, [BK]
fsD = dos(phi, mD)
opt = fsD[:, 3:].ravel()
out(f"  PdD optical band (model): {opt.min()/THz:.2f} - {opt.max()/THz:.2f} THz "
    f"({opt.min()*h/e*1e3:.1f}-{opt.max()*h/e*1e3:.1f} meV)")
for nm, q in (("G", G), ("X", X), ("L", L), ("W", W)):
    fq = dynmat(q, phi, mD)
    out(f"  PdD at {nm}: acoustic {np.round(fq[:3]/THz,2)} THz; optical {np.round(fq[3:]/THz,2)} THz")
optH = dos(phi, mH)[:, 3:].ravel() * anh
out(f"  PdH optical band (model, harmonic x {anh:.3f} anharmonic factor): {optH.min()/THz:.2f} - {optH.max()/THz:.2f} THz")
twoD = (opt[:, None] + opt[None, ::97]).ravel()
out(f"  PdD two-optical-phonon (overtone/combination) band: {twoD.min()/THz:.1f} - {twoD.max()/THz:.1f} THz")
letts = [8.3, 15.3, 20.4]
for fl in letts:
    inD = opt.min() / THz <= fl <= opt.max() / THz
    inH = optH.min() / THz <= fl <= optH.max() / THz
    in2 = twoD.min() / THz <= fl <= twoD.max() / THz
    out(f"  Letts {fl:4.1f} THz ({fl*1e12*h/e*1e3:.1f} meV): in PdD 1-phonon band: {inD}; in PdH band: {inH}; in PdD 2-phonon band: {in2}")
# sensitivity: +/- 0.7 THz on the anchors
for dG, dT in ((-0.7, -0.7), (0.7, 0.7), (0, 1.5)):
    an = dict(anchors); an["optG"] += dG * THz; an["optTop"] += dT * THz
    p2, _ = fit(an)
    o2 = dos(p2, mD, n=12)[:, 3:].ravel()
    out(f"  sensitivity (optG{dG:+.1f}, top{dT:+.1f} THz): PdD optical {o2.min()/THz:.2f}-{o2.max()/THz:.2f} THz")
out("  Loading dependence [BK]: PdH_x optical peak falls from ~68-69 meV (alpha, x~0.02) to ~56 meV (beta, x~0.6-0.7)")
out("  and is ~flat to x~1 (Ross 1998, high-pressure PdH_0.99). Same trend in PdD scaled by ~1/1.5. So raising x")
out("  from 0.63 to 0.95 does not move PdD optical modes to 15.3 THz: that would need ~+35 % hardening.")

# =====================================================================================
# 2. Selection rules for two-laser beat excitation
# =====================================================================================
out("\n2. SELECTION RULES FOR A TWO-LASER BEAT")
out("  Rocksalt (Fm-3m): zone-centre optical mode is T1u -> IR-active, first-order Raman-INACTIVE (as in NaCl).")
out("  A difference-frequency (stimulated-Raman / impulsive) drive therefore cannot couple to the bulk G-mode")
out("  at first order; coupling exists only where inversion symmetry is broken (surface layer, defects).")
dk = 2 * np.pi * 1.33 * (1 / 0.670e-6 - 1 / 0.6827e-6)
kL = 2 * np.pi / a0 * np.sqrt(3) / 2
out(f"  Beat wavevector (co-propagating 670/682.7 nm in D2O): {dk:.2e} m^-1; counter-propagating max ~{2*2*np.pi*1.33/0.676e-6:.2e} m^-1;")
out(f"  L-point wavevector |q_L| = {kL:.2e} m^-1 -> photons supply <1e-3 of the momentum. L-point (zone-boundary)")
out("  modes can be driven only through atomic-scale disorder; Hagelstein's L-point assignment is kinematically")
out("  disfavoured for any optical drive.")

# =====================================================================================
# 3. Cost of a coherent optical-phonon amplitude
# =====================================================================================
out("\n3. COST OF A COHERENT OPTICAL-PHONON AMPLITUDE IN THE LASER-ACCESSIBLE SKIN")
nD = 0.9 * 4 / a0 ** 3                   # D per m^3 at x = 0.9
f = 8.3e12
w = 2 * np.pi * f
nth = 1 / np.expm1(h * f / (kB * 300))
u_th = np.sqrt(hbar / (2 * mD * w) * (2 * nth + 1))
out(f"  thermal occupation at 8.3 THz, 300 K: {nth:.2f}; rms D displacement per mode direction: {u_th*1e10:.3f} A")
for area_mm2, depth_nm in ((1.0, 20.0), (50.0, 20.0)):
    V = area_mm2 * 1e-6 * depth_nm * 1e-9
    ND = nD * V
    for u in (1e-12, 1e-11):
        for Q in (10, 50):
            U = 0.5 * ND * mD * w ** 2 * u ** 2
            P = U * w / Q
            out(f"  {area_mm2:4.0f} mm^2 x {depth_nm:.0f} nm skin ({ND:.1e} D): coherent amplitude {u*1e10:.2f} A, Q={Q}: stored {U:.1e} J, "
                f"drive power {P:.1e} W")
Pl = 0.03
tau = 50 / w
for eta in (1e-3, 1.0):
    Uc = Pl * 0.3 * eta * tau
    nco = Uc / (h * f)
    ND = nD * 1e-6 * 20e-9
    out(f"  30 mW laser, A=0.3, beat-to-mode efficiency {eta:g}, Q=50: coherent phonons {nco:.1e} vs thermal {3*ND*nth:.1e} in the skin")
out("  => a tens-of-mW beat cannot raise the coherent optical-phonon population above ~1e-9 of thermal.")

# =====================================================================================
# 4. Nanoparticle breathing modes (Lamb l = 0)
# =====================================================================================
out("\n4. NANOPARTICLE BREATHING (LAMB l=0, free sphere) MODES")
mats = {  # E (GPa), nu, rho (kg/m3) -- CRC Handbook / [BK]; PdD: moduli -10 %, rho from lattice expansion
    "Pd": (121e9, 0.39, 12020.0),
    "PdD(x~0.7)": (109e9, 0.39, 11070.0),
    "Ni": (200e9, 0.31, 8908.0),
    "Au": (79e9, 0.44, 19300.0),
}


def speeds(Em, nu, rho):
    cL = np.sqrt(Em * (1 - nu) / (rho * (1 + nu) * (1 - 2 * nu)))
    cT = np.sqrt(Em / (2 * rho * (1 + nu)))
    return cL, cT


def lamb0(beta):
    f = lambda x: x / np.tan(x) - 1 + x ** 2 / (4 * beta ** 2)
    return brentq(f, 0.5, 3.1)


d_nm = np.array([1, 2, 5, 10, 20, 50])
bm = {}
for m, (Em, nu, rho) in mats.items():
    cL, cT = speeds(Em, nu, rho)
    xi = lamb0(cT / cL)
    fb = xi * cL / (np.pi * d_nm * 1e-9)
    bm[m] = (xi, cL)
    out(f"  {m:11s}: c_L = {cL:.0f} m/s, c_T = {cT:.0f} m/s, xi = {xi:.3f}; f(d) = " +
        ", ".join(f"{d:g} nm: {x/1e9:.0f} GHz" for d, x in zip(d_nm, fb)))
    for ft in (8.3e12, 15.3e12):
        out(f"      diameter needed for {ft/1e12:.1f} THz: {xi*cL/(np.pi*ft)*1e9:.2f} nm (below continuum validity; ~atomic cluster)")
out("  => breathing modes of 2-50 nm particles lie at 0.07-1.8 THz: they cannot be tuned to the optical-phonon band.")

# =====================================================================================
# 5. Membrane (C3) and wire (C1) modes
# =====================================================================================
out("\n5. MEMBRANE (C3) AND WIRE (C1) MODES")
Em, nu, rho = mats["PdD(x~0.7)"]
cL, cT = speeds(Em, nu, rho)
cbar = np.sqrt(Em / rho)
lam2 = {"(0,1)": 10.2158, "(1,1)": 21.260, "(2,1)": 34.877, "(0,2)": 39.771}
rho_f = 1107.0     # D2O density
out("  Clamped circular plate, in vacuo and with D2O on one side (added-mass factor Gamma=0.6689, Kwak 1991):")
for a_mm in (1.0, 5.0, 10.0):
    for h_um in (25, 50, 100):
        a, hh = a_mm * 1e-3, h_um * 1e-6
        Dflex = Em * hh ** 3 / (12 * (1 - nu ** 2))
        f01 = lam2["(0,1)"] / (2 * np.pi * a ** 2) * np.sqrt(Dflex / (rho * hh))
        beta = 0.6689 * rho_f * a / (rho * hh)
        out(f"   a={a_mm:4.1f} mm, h={h_um:3d} um: f01 = {f01/1e3:9.2f} kHz (vacuum), {f01/np.sqrt(1+beta)/1e3:9.2f} kHz (one side D2O); "
            f"thickness mode c_L/2h = {cL/(2*hh)/1e6:.1f} MHz")
# pressure load (1 bar electrolyte vs vacuum) on an unsupported membrane
out("  Pressure check, dp = 1 bar across an unsupported clamped membrane (large-deflection membrane formulae,")
out("  Timoshenko; w0 = 0.662 a (p a/E h)^(1/3), sigma = 0.423 (E p^2 a^2/h^2)^(1/3)):")
p = 1e5
for a_mm in (0.5, 1.0, 2.0, 5.0, 10.0):
    for h_um in (25, 50, 100):
        a, hh = a_mm * 1e-3, h_um * 1e-6
        w0 = 0.662 * a * (p * a / (Em * hh)) ** (1 / 3)
        sig = 0.423 * (Em * p ** 2 * a ** 2 / hh ** 2) ** (1 / 3)
        out(f"   a={a_mm:4.1f} mm, h={h_um:3d} um: deflection {w0*1e6:7.1f} um, membrane stress {sig/1e6:6.1f} MPa")
sig_allow = 30e6
for h_um in (25, 50, 100):
    hh = h_um * 1e-6
    amax = hh / p * np.sqrt((sig_allow / 0.423) ** 3 / Em)
    out(f"   max unsupported span (diameter) for sigma <= 30 MPa (annealed-Pd yield ~50 MPa [BK]), h={h_um} um: {2*amax*1e3:.2f} mm")
out("  Wire (C1) PdD, free-free / clamped: longitudinal f_n = n c_bar/(2L); flexural f1 = 4.730^2 (d/4) c_bar/(2 pi L^2);")
out("  radial breathing from x J0(x)/J1(x) = 2 (c_T/c_L)^2, x = w r/c_L:")
xr = brentq(lambda x: x * j0(x) / j1(x) - 2 * (cT / cL) ** 2, 1.0, 2.4)
for d_mm in (0.25, 0.5, 1.0):
    for L_mm in (30.0, 100.0):
        dd, LL = d_mm * 1e-3, L_mm * 1e-3
        out(f"   d={d_mm} mm, L={L_mm:.0f} mm: longitudinal {cbar/(2*LL)/1e3:.1f} kHz; flexural {4.730**2*(dd/4)*cbar/(2*np.pi*LL**2):.0f} Hz; "
            f"radial {xr*cL/(2*np.pi*dd/2)/1e6:.2f} MHz")
out("  Q [BK]: bulk PdH_x internal friction Q^-1 ~ 1e-3-1e-2 (hydrogen relaxations, hydride-twin and dislocation")
out("  damping) -> Q_int ~ 1e2-1e3; with electrolyte on one face, acoustic radiation into the liquid dominates:")
Zl, Zs = 1107 * 1400, rho * cL
out(f"  thickness mode radiation-limited Q ~ (pi/2) Z_s/Z_l = {np.pi/2*Zs/Zl:.0f} (Z_s={Zs:.2e}, Z_l={Zl:.2e} Rayl);"
    " flexural modes in liquid Q ~ 10-50.")

# =====================================================================================
# 6. Ultrasonic / megasonic drive
# =====================================================================================
out("\n6. ULTRASONIC / MEGASONIC DRIVE")
for strain in (1e-6, 1e-5, 1e-4):
    I = 0.5 * rho * cL ** 3 * strain ** 2
    out(f"  travelling longitudinal wave, strain {strain:.0e}: intensity {I/1e4:.2e} W/cm^2, stress {Em*strain/1e6:.2f} MPa")
f_hag = 2.2e6
out(f"  Hagelstein/Metzler MHz proposal (~1-3 MHz, e.g. {f_hag/1e6} MHz): quanta per 23.85 MeV = {23.85e6*e/(h*f_hag):.1e};"
    f" thermal occupation {kB*300/(h*f_hag):.1e} per mode (classical regime).")
out("  Standing-wave energy density in a 50 um PdD membrane resonator (thickness mode, 45 MHz) at Q=20, 1 W/cm^2 input:")
hh = 50e-6
Uden = 1e4 * 20 / (2 * np.pi * cL / (2 * hh)) / hh
out(f"   {Uden:.2e} J/m^3 -> strain amplitude {np.sqrt(2*Uden/Em):.1e}; mean energy per D atom {Uden/nD/e:.1e} eV")
out("  Detector compatibility: Si surface-barrier/PIPS detectors and charge preamps are microphonic; MHz drive")
out("  couples as pickup. Use (i) drive on the electrolyte-side cell body, (ii) detectors on separate vacuum")
out("  flanges with elastomer isolation, (iii) gated acquisition (drive on/off blocks) to measure pickup.")

# =====================================================================================
# figures
# =====================================================================================
fig, ax = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw={"width_ratios": [2, 1]})
for b in range(6):
    ax[0].plot(xsD, fD[:, b] / THz, "b-", lw=1.2, label="PdD model" if b == 0 else None)
    ax[0].plot(xsH, fH[:, b] * (anh if b >= 3 else 1) / THz, "g--", lw=0.8, label="PdH (x anharm.)" if b == 0 else None)
for fl in letts:
    ax[0].axhline(fl, color="r", ls=":", lw=1)
    ax[0].text(ticks[-1] * 0.98, fl + 0.2, f"Letts {fl}", color="r", ha="right", fontsize=7)
ax[0].set_xticks(ticks, labels)
ax[0].set(ylabel="frequency (THz)", title="Rocksalt BvK model fitted to [BK] INS anchors", ylim=(0, 24))
ax[0].legend(fontsize=7, loc="upper left")
hist_b = np.linspace(0, 26, 131)
ax[1].hist(fsD.ravel() / THz, hist_b, orientation="horizontal", alpha=0.5, label="PdD 1-phonon", density=True)
ax[1].hist(optH / THz, hist_b, orientation="horizontal", alpha=0.4, label="PdH optical", density=True)
ax[1].hist(twoD / THz, hist_b, orientation="horizontal", alpha=0.3, label="PdD 2-optical", density=True)
for fl in letts:
    ax[1].axhline(fl, color="r", ls=":", lw=1)
ax[1].set(xlabel="DOS (norm.)", ylim=(0, 24))
ax[1].legend(fontsize=7)
savefig(fig, "m6_phonons.png")

fig, ax = plt.subplots(1, 1, figsize=(6, 4))
dd = np.logspace(0, 2, 100)
for m, (xi, cLm) in bm.items():
    ax.loglog(dd, xi * cLm / (np.pi * dd * 1e-9) / THz, label=m)
for fl in letts:
    ax.axhline(fl, color="r", ls=":", lw=1)
ax.set(xlabel="particle diameter (nm)", ylabel="breathing-mode frequency (THz)", title="Lamb l=0 breathing modes vs Letts lines")
ax.legend(fontsize=7)
ax.grid(alpha=0.3, which="both")
savefig(fig, "m6_breathing.png")
out.save()
