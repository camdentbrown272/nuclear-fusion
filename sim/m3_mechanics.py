"""
M3 / Q3: mechanics of the detector-facing membrane and of supported films.

  A. Hydrogen-induced (Vegard) strain and elastic limits.
  B. Pressure load across a vacuum-backed membrane: axisymmetric Foeppl-von Karman (FvK) clamped plate,
     solved with scipy solve_bvp; verification against small-deflection and Hencky membrane limits.
     Support-grid design: stress vs free-span radius, open fraction.
  C. Hydrogen misfit with an edge clamp: buckling and required radial compliance.
  D. Supported films: Stoney stress/curvature, buckle-delamination critical thickness.
  E. Cracking of a deloaded (alpha) skin on a beta substrate: channel-crack critical depth.

Run:  python3 sim/m3_mechanics.py
Writes docs/models/figs/m3_mech_plate.png, m3_mech_grid.png, m3_mech_films.png, m3_mechanics.txt
"""
import os
import numpy as np
from scipy.integrate import solve_bvp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from m3_common import gap, Omega_Pd

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "docs", "models", "figs")
LOG = []


def log(s=""):
    print(s)
    LOG.append(s)


# ------------------------------------------------------------------ mechanical parameters
E_Pd = 121e9      # Pa  [web] Wikipedia "Palladium" (Young's modulus 121 GPa) https://en.wikipedia.org/wiki/Palladium
NU_Pd = 0.39      # [web] same source (Poisson ratio 0.39)
# yield strength: annealed ~35-50 MPa, UTS ~170 MPa; hard (cold worked) yield ~200-300 MPa  [mem]
#   ASM Metals Handbook vol. 2 (Palladium); Vickers hardness of commercial Pd 400-600 MPa [web] (Wikipedia)
#   -> Tabor sigma_y ~ H/3 = 130-200 MPa for as-supplied metal.  Hydrogen alpha/beta cycling hardens Pd to
#   roughly cold-worked levels (dislocation densities 1e14-1e15 m^-2) [mem] Flanagan & Oates 1991.
SY = {"annealed": 40e6, "as-supplied/H-cycled": 150e6, "cold-worked": 250e6}
VBAR_D = 1.73e-6  # m^3/mol partial molar volume of H(D) in Pd [mem] Peisl, in Hydrogen in Metals I (1978)
#                   https://doi.org/10.1007/3540087052_50 (dV/Omega = 0.19);  brief quotes ~1.7 cm3/mol
ETA = VBAR_D / (3 * Omega_Pd)   # linear Vegard strain per unit x


def plate_bvp(P, eps_star=0.0, nu=NU_Pd, rho0=1e-4, guess=None):
    """Nondimensional FvK clamped circular plate, immovable edge.
    P = p a^4/(E h^4);  eps_star: uniform in-plane eigenstrain (hydrogen expansion) -> prestress.
    y = [U, U', Phi, Phi'];  U = u a/h^2, Phi = (dw/dr) a/h, rho = r/a."""
    s = (1 + nu) * eps_star  # multiplied by a^2/h^2 outside
    def f(r, y):
        U, Up, Ph, Php = y
        nr = Up + 0.5 * Ph ** 2 + nu * U / r - s_eff
        Upp = -Ph * Php - Up / r + U / r ** 2 - (1 - nu) * Ph ** 2 / (2 * r)
        Phpp = -Php / r + Ph / r ** 2 + 12 * nr * Ph - 6 * (1 - nu ** 2) * P * r
        return np.vstack([Up, Upp, Php, Phpp])

    def bc(ya, yb):
        # regular centre: U ~ rho U', Phi ~ rho Phi' (linear), clamped immovable edge
        return np.array([ya[0] - rho0 * ya[1], ya[2] - rho0 * ya[3], yb[0], yb[2]])
    global s_eff
    s_eff = plate_bvp.s_scaled
    r = np.linspace(rho0, 1, 400)
    if guess is None:
        Ph = 12 * (1 - nu ** 2) * P / 16 * r * (1 - r ** 2) * 0.5
        y0 = np.vstack([np.zeros_like(r), np.zeros_like(r), Ph, np.gradient(Ph, r)])
    else:
        y0 = guess.sol(r)
    sol = solve_bvp(f, bc, r, y0, tol=1e-6, max_nodes=200000)
    return sol


plate_bvp.s_scaled = 0.0


def plate_results(sol, nu=NU_Pd):
    r = np.linspace(sol.x[0], 1, 2000)
    U, Up, Ph, Php = sol.sol(r)
    W0 = np.trapezoid(Ph, r)
    nr = Up + 0.5 * Ph ** 2 + nu * U / r - plate_bvp.s_scaled
    nt = U / r + nu * (Up + 0.5 * Ph ** 2) - plate_bvp.s_scaled
    br = 0.5 * (Php + nu * Ph / r)
    bt = 0.5 * (Ph / r + nu * Php)
    # stresses in units of E/(1-nu^2) (h/a)^2
    vm = []
    for sgn in (1, -1):
        sr = nr + sgn * br
        st = nt + sgn * bt
        vm.append(np.sqrt(sr ** 2 - sr * st + st ** 2))
    vm = np.maximum(*vm)
    return dict(W0=W0, vm_max=vm.max(), vm_centre=vm[0], vm_edge=vm[-1], r=r, vm=vm, nr=nr)


def solve_pressure(p, a, h, nu=NU_Pd, E=E_Pd):
    """Returns max von Mises stress (Pa), centre deflection (m) for pressure p on clamped disc radius a."""
    P = p * a ** 4 / (E * h ** 4)
    plate_bvp.s_scaled = 0.0
    sol = None
    for Pk in np.geomspace(min(P, 1.0), P, max(2, int(np.log10(max(P, 1)) * 6) + 2)):
        sol = plate_bvp(Pk, nu=nu, guess=sol)
    res = plate_results(sol, nu)
    scale = E / (1 - nu ** 2) * (h / a) ** 2
    return res["vm_max"] * scale, res["W0"] * h, res, sol


def main():
    a_, b_ = gap(298.15, "D")
    log("=== A. Vegard strain and elastic limits (PdD) ===")
    log(f"eta = Vbar/(3 Omega) = {ETA:.4f} per unit x; alpha->beta jump dx = {b_ - a_:.3f} -> "
        f"linear misfit {ETA * (b_ - a_) * 100:.2f} % (measured lattice-parameter jump 3.89->4.02 A = 3.3-3.5 %)")
    log(f"0 -> 0.95 loading: linear strain {ETA * 0.95 * 100:.2f} % (linear Vegard), volume "
        f"{((1 + ETA * 0.95) ** 3 - 1) * 100:.1f} %; cross-check: a = 3.89 A (Pd) -> 4.09 A at x ~ 0.98 [web, "
        f"jcmns 17 (2015) 35] gives {(4.09 / 3.89 - 1) * 100:.1f} % linear, i.e. eta_eff = "
        f"{(4.09 / 3.89 - 1) / 0.977:.3f}; the linear law over-predicts high-x strain by ~20 % (conservative)")
    for k, sy in SY.items():
        dx_el = sy * (1 - NU_Pd) / (E_Pd * ETA)
        log(f"  {k:22s} sigma_y = {sy / 1e6:5.0f} MPa: biaxial elastic strain limit {sy * (1 - NU_Pd) / E_Pd * 1e4:.1f}e-4 "
            f"-> max constrained composition difference dx = {dx_el:.4f}")

    # ------------------------------------------------------------------ B. verification
    log("")
    log("=== B. FvK clamped plate: verification ===")
    # small deflection: W0 = 12(1-nu^2) P / 64
    plate_bvp.s_scaled = 0.0
    for nu in (0.3, NU_Pd):
        sol = plate_bvp(1e-3, nu=nu)
        r = plate_results(sol, nu)
        log(f"  nu={nu}: small-P W0/P = {r['W0'] / 1e-3:.5f} vs Kirchhoff 12(1-nu^2)/64 = {12 * (1 - nu ** 2) / 64:.5f}")
    # Hencky membrane limit (nu = 0.3): w0/a = 0.662 (pa/Eh)^(1/3); centre stress 0.423 (E p^2 a^2/h^2)^(1/3)
    a, h, E = 10e-3, 10e-6, 1.0
    for P in (1e5, 1e6, 1e7):
        sol = None
        for Pk in np.geomspace(1, P, 40):
            sol = plate_bvp(Pk, nu=0.3, guess=sol)
        r = plate_results(sol, 0.3)
        q = P * (1 / 1) ** 0  # p a/(E h) = P (h/a)^3 in nondim
        w_hencky = 0.662 * (P) ** (1 / 3)          # w0/h since (pa/Eh)^(1/3) a / h = P^(1/3)
        s_c = r["vm_centre"] / (1 - 0.3 ** 2)       # centre membrane stress in units E (h/a)^2
        s_h = 0.423 * P ** (2 / 3)                  # (E p^2 a^2/h^2)^(1/3) / (E h^2/a^2) = P^(2/3)
        log(f"  P={P:.0e}: W0 = {r['W0']:.3f} vs Hencky {w_hencky:.3f} ({(r['W0'] / w_hencky - 1) * 100:+.1f} %); "
            f"centre stress {s_c:.3f} vs Hencky {s_h:.3f} ({(s_c / s_h - 1) * 100:+.1f} %); max/centre = "
            f"{r['vm_max'] / r['vm_centre']:.2f} (edge bending)")
    log("  Clamped-plate buckling under uniform in-plane misfit: 12(1+nu) eps* a^2/h^2 = j_{1,1}^2 = 14.68 "
        "(Timoshenko & Gere, Theory of Elastic Stability, 1961)")

    # ------------------------------------------------------------------ B. pressure on membranes
    log("")
    log("=== B. 1 atm across a vacuum-backed membrane (clamped, immovable edge) ===")
    p = 1.01325e5
    hs = [10e-6, 25e-6, 50e-6, 100e-6]
    bs = np.geomspace(0.1e-3, 12.5e-3, 22)
    RES = {}
    for h in hs:
        for b in bs:
            s, w, _, _ = solve_pressure(p, b, h)
            RES[(h, b)] = (s, w)
    log("  unsupported discs (a = radius):")
    for a in (5e-3, 10e-3, 12.5e-3):
        for h in hs:
            s, w, _, _ = solve_pressure(p, a, h)
            log(f"    a={a * 1e3:4.1f} mm h={h * 1e6:4.0f} um: max von Mises {s / 1e6:7.1f} MPa, centre deflection "
                f"{w * 1e3:6.3f} mm")
    log("  largest free-span radius b keeping max von Mises below a stress allowable (1 atm):")
    for h in hs:
        sv = np.array([RES[(h, b)][0] for b in bs])
        out = []
        for allow in (20e6, 40e6, 100e6):
            ok = bs[sv <= allow]
            out.append(f"{allow / 1e6:.0f} MPa: b <= {ok.max() * 1e3:.2f} mm" if len(ok) else f"{allow / 1e6:.0f}: none")
        log(f"    h = {h * 1e6:4.0f} um: " + "; ".join(out))
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for h in hs:
        ax[0].loglog(bs * 1e3, [RES[(h, b)][0] / 1e6 for b in bs], "o-", ms=3, label=f"h = {h * 1e6:.0f} µm")
        ax[1].loglog(bs * 1e3, [RES[(h, b)][1] * 1e6 for b in bs], "o-", ms=3, label=f"h = {h * 1e6:.0f} µm")
    for k, sy in SY.items():
        ax[0].axhline(sy / 1e6, ls=":", color="grey")
        ax[0].text(0.11, sy / 1e6 * 1.1, f"σy {k}", fontsize=7)
    ax[0].set_xlabel("free-span radius b (mm)")
    ax[0].set_ylabel("max von Mises stress at 1 atm (MPa)")
    ax[0].legend(fontsize=7)
    ax[0].set_title("Clamped Pd disc/hole, 1 atm, FvK", fontsize=9)
    ax[1].set_xlabel("free-span radius b (mm)")
    ax[1].set_ylabel("centre deflection (µm)")
    ax[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "m3_mech_plate.png"), dpi=130)
    plt.close(fig)

    # grid: hexagonal array of round holes, diameter d = 2b, pitch s; open fraction = pi/(2 sqrt3) (d/s)^2
    log("")
    log("=== B. Support grid (vacuum side): hexagonal round holes, diameter d, web width w, pitch d+w ===")
    log("  open fraction f = 0.9069 (d/(d+w))^2 ; hole-wall shadowing for a grid of thickness t: particles")
    log("  with tan(theta) > d/t cannot pass; transmission at angle theta for a round hole:")
    log("  T(q) = (2/pi)(acos q - q sqrt(1-q^2)), q = t tan(theta)/d.")
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for h in (25e-6, 50e-6):
        for w in (0.1e-3, 0.2e-3):
            d = 2 * bs
            f = 0.9069 * (d / (d + w)) ** 2
            sv = np.array([RES[(h, b)][0] for b in bs])
            ax[0].plot(f, sv / 1e6, "o-", ms=3, label=f"h={h * 1e6:.0f} µm, web {w * 1e3:.1f} mm")
    ax[0].set_yscale("log")
    ax[0].set_xlabel("grid open fraction")
    ax[0].set_ylabel("membrane max stress at 1 atm (MPa)")
    ax[0].axhline(40, ls=":", color="grey")
    ax[0].legend(fontsize=7)
    ax[0].set_title("Open fraction vs membrane stress", fontsize=9)
    th = np.linspace(0, 80, 200)
    for ar in (0.25, 0.5, 1.0):
        q = np.clip(ar * np.tan(np.radians(th)), 0, 1)
        ax[1].plot(th, 2 / np.pi * (np.arccos(q) - q * np.sqrt(1 - q ** 2)), label=f"t/d = {ar}")
    ax[1].set_xlabel("emission angle from normal (deg)")
    ax[1].set_ylabel("transmission through one hole")
    ax[1].legend(fontsize=7)
    ax[1].set_title("Hole-wall shadowing (hand-off to M5)", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "m3_mech_grid.png"), dpi=130)
    plt.close(fig)
    for h, b, w, t in ((25e-6, 0.5e-3, 0.15e-3, 0.5e-3), (50e-6, 0.75e-3, 0.2e-3, 0.5e-3), (50e-6, 1.0e-3, 0.2e-3, 0.5e-3)):
        s, wdef, _, _ = solve_pressure(p, b, h)
        d = 2 * b
        f = 0.9069 * (d / (d + w)) ** 2
        # angle-averaged transmission for isotropic (per unit solid angle) emission into the forward hemisphere
        thg = np.linspace(0, np.pi / 2, 4000)
        q = np.clip(t / d * np.tan(thg), 0, 1)
        Tq = 2 / np.pi * (np.arccos(q) - q * np.sqrt(1 - q ** 2))
        Tavg = np.trapezoid(Tq * np.sin(thg), thg)
        Tcone = np.trapezoid((Tq * np.sin(thg))[thg < np.radians(30)], thg[thg < np.radians(30)]) / \
            np.trapezoid(np.sin(thg)[thg < np.radians(30)], thg[thg < np.radians(30)])
        log(f"  h={h * 1e6:.0f} um, hole d={d * 1e3:.1f} mm, web {w * 1e3:.2f} mm, grid t={t * 1e3:.1f} mm: "
            f"stress {s / 1e6:.1f} MPa, dimple {wdef * 1e6:.1f} um, open fraction {f:.2f}, "
            f"hemisphere transmission (x open) {Tavg * f:.2f}, within 30 deg cone {Tcone * f:.2f}")
    # grid plate itself: perforated plate, clamped, uniform load, radius 12.5 mm; ligament efficiency
    for mat, Eg, sy in (("Mo", 329e9, 550e6), ("316L", 193e9, 200e6), ("Ti (NOT recommended: absorbs D)", 116e9, 275e6)):
        for tg in (0.5e-3, 1.0e-3):
            A = 12.5e-3
            # thick clamped plate edge stress 0.75 p a^2/t^2 (Roark), perforated: divide by ligament efficiency
            eff = 0.2 / 1.7   # (s - d)/s for d=1.5, s=1.7 mm
            sg = 0.75 * p * A ** 2 / tg ** 2 / eff
            log(f"  grid {mat:30s} t={tg * 1e3:.1f} mm, R=12.5 mm: edge bending stress ~{sg / 1e6:.0f} MPa "
                f"(sigma_y {sy / 1e6:.0f} MPa; Roark clamped plate / ligament efficiency {eff:.2f})")

    # ------------------------------------------------------------------ C. misfit with edge clamp
    log("")
    log("=== C. Hydrogen expansion with an edge clamp (disc radius a, thickness h) ===")
    for a in (5e-3, 10e-3, 12.5e-3):
        for h in (25e-6, 50e-6):
            eps_cr = 14.68 / (12 * (1 + NU_Pd)) * (h / a) ** 2
            log(f"  a={a * 1e3:4.1f} mm h={h * 1e6:3.0f} um: buckling misfit eps* = {eps_cr:.2e} -> dx = "
                f"{eps_cr / ETA:.2e};  full loading 0->0.95 misfit {ETA * 0.95:.3f} = "
                f"{ETA * 0.95 / eps_cr:.0f} x critical; radial growth at edge {ETA * 0.95 * a * 1e3:.3f} mm")
    log("  Friction-limited stress for a floating membrane pressed on the grid by 1 atm: sigma ~ mu p a / h")
    for a in (10e-3, 12.5e-3):
        for h in (25e-6, 50e-6):
            log(f"    mu=0.3 a={a * 1e3:.1f} mm h={h * 1e6:.0f} um: {0.3 * p * a / h / 1e6:.1f} MPa")

    # ------------------------------------------------------------------ D. supported films
    log("")
    log("=== D. Films on rigid substrates (fully constrained in-plane) ===")
    s_el = E_Pd * ETA * 0.9 / (1 - NU_Pd)
    log(f"  elastic biaxial stress for x = 0.9: {s_el / 1e9:.1f} GPa (compressive) -> always capped by film flow "
        f"stress; measured hydrogen-loaded Pd film stresses are -1 to -3 GPa [mem] (Pundt & Kirchheim, "
        f"Annu. Rev. Mater. Res. 36 (2006) 555, https://doi.org/10.1146/annurev.matsci.36.090804.094451)")
    log("  steady-state buckle-delamination: G = (1-nu^2) h sigma^2 / (2E) >= Gamma_i  ->  h_c = 2 E Gamma_i / ((1-nu^2) sigma^2)")
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    sig = np.geomspace(0.1e9, 3e9, 100)
    for G in (0.5, 2, 10):
        hc = 2 * E_Pd * G / ((1 - NU_Pd ** 2) * sig ** 2)
        ax[0].loglog(sig / 1e9, hc * 1e9, label=f"Γ_i = {G} J/m²")
        log("    Gamma_i=%4.1f J/m2: h_c = %s" % (G, ", ".join(
            f"{2 * E_Pd * G / ((1 - NU_Pd ** 2) * s ** 2) * 1e9:.0f} nm @ {s / 1e9:.1f} GPa" for s in (0.5e9, 1e9, 2e9))))
    ax[0].set_xlabel("film compressive stress (GPa)")
    ax[0].set_ylabel("critical thickness for delamination (nm)")
    ax[0].legend(fontsize=7)
    ax[0].set_title("Supported Pd film: buckle-delamination", fontsize=9)
    # Stoney
    Es, nus = 130e9, 0.28   # Si(100) biaxial approx [mem]
    log("  Stoney curvature of a Si substrate (300 um) under a loaded Pd film (1 GPa): kappa = 6 sigma h_f (1-nu_s)/(E_s h_s^2)")
    for hf in (100e-9, 1e-6):
        kap = 6 * 1e9 * hf * (1 - nus) / (Es * (300e-6) ** 2)
        log(f"    h_f = {hf * 1e9:.0f} nm: kappa = {kap:.3g} 1/m (R = {1 / kap:.2f} m) -> in-situ loading monitor")

    # ------------------------------------------------------------------ E. channel cracking of alpha skin
    log("")
    log("=== E. Deloaded alpha skin (tension) on beta substrate: channel cracking ===")
    log("  sigma_skin = min(E eta dx/(1-nu), sigma_flow); crack when Z sigma^2 d / E >= Gamma_c, Z = 1.976 "
        "(Beuth 1992, matched film/substrate)")
    Z = 1.976
    for KIc in (1.0, 3.0, 10.0, 30.0):
        Gc = (KIc * 1e6) ** 2 * (1 - NU_Pd ** 2) / E_Pd
        out = []
        for sf in (150e6, 300e6, 600e6):
            dc = Gc * E_Pd / (Z * sf ** 2)
            out.append(f"{dc * 1e6:.2g} um @ {sf / 1e6:.0f} MPa")
        log(f"  K_c = {KIc:4.1f} MPa m^1/2 (Gamma = {Gc:.3g} J/m2): critical skin depth " + ", ".join(out))
        dd = np.geomspace(0.1e-6, 1e-3, 100)
        ax[1].loglog(dd * 1e6, np.sqrt(Gc * E_Pd / (Z * dd)) / 1e6, label=f"K_c = {KIc} MPa√m")
    ax[1].axhline(E_Pd * ETA * (b_ - a_) / (1 - NU_Pd) / 1e6, color="k", ls="--", lw=0.8)
    ax[1].text(0.12, E_Pd * ETA * (b_ - a_) / (1 - NU_Pd) / 1e6 * 0.6, "elastic α/β misfit stress", fontsize=7)
    for k, sy in SY.items():
        ax[1].axhline(sy / 1e6, color="grey", ls=":")
    ax[1].set_xlabel("depth of deloaded α skin (µm)")
    ax[1].set_ylabel("skin stress needed to channel-crack (MPa)")
    ax[1].set_ylim(10, 2e4)
    ax[1].legend(fontsize=7)
    ax[1].set_title("Crack onset: skin stress above the curve cracks", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "m3_mech_films.png"), dpi=130)
    plt.close(fig)
    open(os.path.join(OUT, "m3_mechanics.txt"), "w").write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
