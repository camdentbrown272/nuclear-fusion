"""
M3 / Q1 + Q5: loading time vs characteristic dimension (film, foil, wire, nanoparticle) with the
alpha/beta miscibility gap (moving phase boundary), temperature 20-90 C, isotope H vs D.

Run:  python3 sim/m3_loading.py
Writes docs/models/figs/m3_isotherm.png, m3_loading_verification.png, m3_loading_profiles.png,
       m3_loading_time.png, m3_loading.txt
"""
import os
import numpy as np
from scipy.special import erf, jn_zeros
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from m3_common import (Material, FV1D, gap, p_plateau, fugacity, x_of_f, D_alpha, D_chem, half_ln_f,
                       n_Pd, e, ATM, kbeta)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "docs", "models", "figs")
os.makedirs(OUT, exist_ok=True)
LOG = []


def log(s=""):
    print(s)
    LOG.append(s)


GEOM = {0: "slab (film on inert substrate: l = thickness; free foil: l = half-thickness)",
        1: "cylinder (wire: l = radius)", 2: "sphere (particle: l = radius)"}


# ---------------------------------------------------------------- analytic references
def frac_const_D(tau, k):
    """Fractional uptake M/Minf for constant D, surface step, uniform zero initial (Crank 1975)."""
    if k == 0:
        n = np.arange(0, 400)
        lam = (2 * n + 1) * np.pi / 2
        return 1 - np.sum(2 / lam ** 2 * np.exp(-lam ** 2 * tau))
    if k == 1:
        z = jn_zeros(0, 400)
        return 1 - np.sum(4 / z ** 2 * np.exp(-z ** 2 * tau))
    n = np.arange(1, 400)
    return 1 - 6 / np.pi ** 2 * np.sum(np.exp(-n ** 2 * np.pi ** 2 * tau) / n ** 2)


def tau90_const(k):
    return brentq(lambda t: frac_const_D(t, k) - 0.9, 1e-4, 5)


def t90_numeric(mat, k, xs, ell=1e-6, N=200, frac=0.9, ret_hist=False):
    s = FV1D(mat, ell, N=N, geom=k, refine=(False, True), ratio=1.02)
    x0 = np.zeros(s.N)
    target = frac * xs
    hist = []

    prev = [0.0, 0.0]
    cross = [None]

    def stop(t, x):
        m = s.mean(x)
        if m >= target:
            # linear interpolation of the crossing time between the last two accepted steps
            cross[0] = prev[0] + (t - prev[0]) * (target - prev[1]) / max(m - prev[1], 1e-30)
            return True
        prev[0], prev[1] = t, m
        return False

    def rec(t, x):
        if ret_hist:
            hist.append((t, s.mean(x), x.copy()))

    Dref = D_alpha(mat.T, mat.iso)
    tend = 1e3 * ell ** 2 / Dref
    t, x = s.run(x0, tend, ("sym",), ("dir", xs), dt0=1e-8 * ell ** 2 / Dref,
                 dtmax=0.004 * ell ** 2 / Dref, stop=stop, record=rec)
    t = cross[0] if cross[0] is not None else t
    return (t, s, hist) if ret_hist else t


def main():
    T0 = 298.15
    # ------------------------------------------------------------ isotherm figure + Q5 table
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    xs = np.linspace(0.001, 0.985, 2000)
    cols = {293.15: "#1f77b4", 333.15: "#2ca02c", 363.15: "#d62728"}
    for T, c in cols.items():
        for iso, ls in (("D", "-"), ("H", "--")):
            ax[0].semilogy(xs, fugacity(xs, T, iso) / ATM, ls, color=c,
                           label=f"{iso}, {T - 273.15:.0f} C")
            if iso == "D":
                ax[1].semilogy(xs, D_chem(xs, T, iso), ls, color=c, label=f"D, {T - 273.15:.0f} C")
                ax[1].semilogy(xs, D_chem(xs, T, iso, blocking=False), ":", color=c, lw=1)
    ax[0].axhline(1, color="grey", lw=0.5)
    ax[0].set_xlabel("x = D/Pd (H/Pd)")
    ax[0].set_ylabel("equilibrium fugacity (atm)")
    ax[0].set_ylim(1e-5, 1e6)
    ax[0].set_title("Model isotherm (Maxwell plateau; beta branch anchored f(0.90)=1e4 atm, PdD 25 C)",
                    fontsize=8)
    ax[0].legend(fontsize=7)
    ax[1].set_xlabel("x")
    ax[1].set_ylabel("chemical diffusivity D_chem (m$^2$/s)")
    ax[1].set_ylim(1e-12, 1e-8)
    ax[1].set_title("D_chem = D*·x·d(mu/kT)/dx; solid: D*=D_a(1-x), dotted: D*=D_a", fontsize=8)
    ax[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "m3_isotherm.png"), dpi=130)
    plt.close(fig)

    log("=== Q5: temperature and isotope dependence (model isotherm + Voelkl-Alefeld diffusion) ===")
    log(f"beta-branch slope kb(298 K) = {kbeta(298.15):.2f} (per unit x, in mu/kT)")
    log("iso  T(C)  alpha_max beta_min  p_plateau(bar)  x(1 bar)  f(x=.90) f(.95) f(.97) (atm)   "
        "D_alpha   D_chem(.90)  (m2/s)")
    for iso in ("D", "H"):
        for T in (293.15, 313.15, 333.15, 363.15):
            a, b = gap(T, iso)
            log(f" {iso}  {T - 273.15:4.0f}   {a:.3f}    {b:.3f}      {p_plateau(T, iso):.4f}        "
                f"{x_of_f(1.0, T, iso):.3f}    {fugacity(0.9, T, iso) / ATM:8.3g} "
                f"{fugacity(0.95, T, iso) / ATM:8.3g} {fugacity(0.97, T, iso) / ATM:8.3g}   "
                f"{D_alpha(T, iso):.3g}   {float(D_chem(0.9, T, iso)):.3g}")
    log("Isotope summary at 25 C: D_alpha(D)/D_alpha(H) = %.2f; plateau p(D)/p(H) = %.1f; "
        "fugacity for x=0.90: D %.3g atm vs H %.3g atm (x %.1f)" % (
            D_alpha(298.15, "D") / D_alpha(298.15, "H"), p_plateau(298.15, "D") / p_plateau(298.15, "H"),
            fugacity(0.9, 298.15, "D") / ATM, fugacity(0.9, 298.15, "H") / ATM,
            fugacity(0.9, 298.15, "D") / fugacity(0.9, 298.15, "H")))
    # x reached at a given fugacity for H and D
    log("Loading reached at equal fugacity (25 C):  f(atm)  x_D   x_H")
    for f in (1, 1e2, 1e3, 1e4, 3e4, 1e5):
        log(f"   {f:8.0e}   {x_of_f(f * ATM, 298.15, 'D'):.3f}  {x_of_f(f * ATM, 298.15, 'H'):.3f}")

    # ------------------------------------------------------------ verification
    log("")
    log("=== Verification V1: constant-D limit vs analytic series (tau90 = D t90 / l^2) ===")
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    Dc = 1e-10
    ideal = Material(ideal=lambda x: np.full_like(x, Dc), smooth=0)
    ideal_T = ideal.T
    for k in (0, 1, 2):
        ta = tau90_const(k)
        row = []
        for N in (25, 50, 100, 200):
            s = FV1D(ideal, 1e-6, N=N, geom=k, refine=(False, True), ratio=1.02)
            t, x = s.run(np.zeros(s.N), 10 * 1e-12 / Dc, ("sym",), ("dir", 0.5),
                         dt0=1e-6 * 1e-12 / Dc, dtmax=2e-4 * 1e-12 / Dc,
                         stop=lambda t, x, s=s: s.mean(x) >= 0.45)
            row.append(t * Dc / 1e-12)
        log(f"  geom {k}: analytic {ta:.4f}; FV N=25/50/100/200: " +
            " ".join(f"{r:.4f}" for r in row) + f"  (err N=200: {abs(row[-1] - ta) / ta * 100:.2f} %)")
    # V2: one-phase Stefan (Neumann) problem
    log("=== Verification V2: one-phase Stefan problem (gap 0->b, constant D above b) vs Neumann ===")
    bS, xsS = 0.5, 0.8
    stef = Material(ideal=lambda x: np.where(x > bS, Dc, 0.0), smooth=0.0005)
    St = (xsS - bS) / bS
    lam = brentq(lambda l: l * np.sqrt(np.pi) * np.exp(l ** 2) * erf(l) - St, 1e-4, 5)
    Ls = 10e-6
    s = FV1D(stef, Ls, N=400, geom=0, refine=(True, False), ratio=1.0)
    # put Dirichlet at r=0 (left), semi-infinite toward right (no flux far away)
    times, upt = [], []

    def rec(t, x):
        times.append(t)
        upt.append(np.sum(x * s.V))
    s.run(np.zeros(s.N), (0.6 * Ls / (2 * lam)) ** 2 / Dc, ("dir", xsS), ("sym",),
          dt0=1e-9, dtmax=2e-3, record=rec)
    times = np.array(times)
    upt = np.array(upt)
    a = 2 * np.sqrt(Dc * times)
    sfr = lam * a
    ana = xsS * sfr - (xsS - bS) / erf(lam) * (sfr * erf(lam) + a / np.sqrt(np.pi) * (np.exp(-lam ** 2) - 1))
    sel = times > times[-1] * 0.1
    err = np.max(np.abs(upt[sel] - ana[sel]) / ana[sel])
    log(f"  St = {St:.2f}, lambda = {lam:.4f}; max rel. error of uptake (t > 0.1 t_end): {err * 100:.2f} %")
    ax[0].plot(times, upt * 1e6, label="FV (N=400)")
    ax[0].plot(times, ana * 1e6, "--", label="Neumann exact")
    ax[0].set_xlabel("t (s)")
    ax[0].set_ylabel("uptake (x·µm)")
    ax[0].set_title("V2 Stefan problem", fontsize=9)
    ax[0].legend(fontsize=8)

    # V3: mesh convergence on the real isotherm (slab, xs=0.9)
    log("=== Verification V3: mesh convergence, real PdD isotherm, slab, x_s = 0.90, 25 C ===")
    mat = Material(T0, "D")
    Ns = (25, 50, 100, 200, 400)
    tv = []
    for N in Ns:
        t = t90_numeric(mat, 0, 0.9, ell=1e-6, N=N)
        tv.append(t)
    Da = D_alpha(T0, "D")
    for N, t in zip(Ns, tv):
        log(f"  N={N:4d}: tau90 = {t * Da / 1e-12:.4f}")
    ax[1].semilogx(Ns, np.array(tv) * Da / 1e-12, "o-")
    ax[1].set_xlabel("cells")
    ax[1].set_ylabel("tau90 = D_a t90 / l^2")
    ax[1].set_title("V3 mesh convergence (slab, x_s=0.9)", fontsize=9)

    # ------------------------------------------------------------ profiles
    t, s, hist = t90_numeric(mat, 0, 0.9, ell=1e-6, N=200, ret_hist=True)
    for fr in (0.1, 0.3, 0.5, 0.7, 0.9):
        tt, mm, xx = min(hist, key=lambda h: abs(h[1] - fr * 0.9))
        ax[2].plot(s.rc / 1e-6, xx, label=f"mean {mm:.2f}")
    ax[2].axhspan(*gap(T0, "D"), color="grey", alpha=0.15, label="miscibility gap")
    ax[2].set_xlabel("depth from inert back (l units)  -> loaded face at 1")
    ax[2].set_ylabel("x")
    ax[2].set_title("Profiles: shrinking alpha core, sharp alpha/beta front", fontsize=9)
    ax[2].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "m3_loading_verification.png"), dpi=130)
    plt.close(fig)

    # ------------------------------------------------------------ tau90 table
    log("")
    log("=== Q1: dimensionless loading time tau90 = D_alpha(T) t90 / l^2 (surface held at x_s; x0=0) ===")
    log("   (diffusion-limited: surface in equilibrium with the charging fugacity; see current/surface limits below)")
    TAU = {}
    for iso in ("D", "H"):
        for T in (298.15, 333.15, 363.15):
            m = Material(T, iso)
            for xs_ in (0.70, 0.90, 0.95):
                for k in (0, 1, 2):
                    TAU[(iso, T, xs_, k)] = t90_numeric(m, k, xs_, ell=1e-6, N=150) * D_alpha(T, iso) / 1e-12
            log(f" {iso} {T - 273.15:3.0f} C: " + "  ".join(
                f"x_s={xs_:.2f}: slab {TAU[(iso, T, xs_, 0)]:.3f} cyl {TAU[(iso, T, xs_, 1)]:.3f} "
                f"sph {TAU[(iso, T, xs_, 2)]:.3f}" for xs_ in (0.70, 0.90, 0.95)))
    log("  constant-D reference: slab %.3f, cyl %.3f, sph %.3f" % tuple(tau90_const(k) for k in (0, 1, 2)))
    # sensitivity: no site blocking
    mnb = Material(T0, "D", blocking=False)
    tnb = [t90_numeric(mnb, k, 0.9, ell=1e-6, N=150) * Da / 1e-12 for k in (0, 1, 2)]
    log("  sensitivity D*=D_alpha (no blocking), D, 25 C, x_s=0.9: slab %.3f cyl %.3f sph %.3f" % tuple(tnb))

    # ------------------------------------------------------------ dimension table
    log("")
    log("=== Q1: t90 (diffusion-limited) for PdD, 25 C, x_s = 0.90 [and 0.95] ===")
    dims = [5e-9, 20e-9, 100e-9, 1e-6, 5e-6, 25e-6, 50e-6, 100e-6, 250e-6, 500e-6, 1e-3]
    lab = ["5 nm", "20 nm", "100 nm", "1 um", "5 um", "25 um", "50 um", "100 um", "250 um", "500 um", "1 mm"]

    def fmt(t):
        if t < 1e-3:
            return f"{t * 1e6:7.2f} us"
        if t < 1:
            return f"{t * 1e3:7.2f} ms"
        if t < 120:
            return f"{t:7.1f} s "
        if t < 7200:
            return f"{t / 60:7.1f} min"
        if t < 2 * 86400:
            return f"{t / 3600:7.1f} h "
        return f"{t / 86400:7.1f} d "

    # current-limited time: Q = n x (V/A) e / (eta i);  V/A = l/(k+1)
    log("  l      | film/foil(slab) t90   | wire t90       | sphere t90     | charge-limited t (x=0.9, "
        "slab/wire/sphere) at eta*i = 10 mA/cm2")
    for d, lb in zip(dims, lab):
        tt = [TAU[("D", T0, 0.90, k)] * d ** 2 / Da for k in (0, 1, 2)]
        tq = [n_Pd * 0.9 * d / (k + 1) * e / (10e-3 * 1e4) for k in (0, 1, 2)]
        log(f"  {lb:7s}| {fmt(tt[0]):>12s}          | {fmt(tt[1]):>12s}   | {fmt(tt[2]):>12s}   | "
            f"{fmt(tq[0])} / {fmt(tq[1])} / {fmt(tq[2])}")
    log("  x_s = 0.95 multiplies t90 by: slab %.2f, wire %.2f, sphere %.2f" % tuple(
        TAU[("D", T0, 0.95, k)] / TAU[("D", T0, 0.90, k)] for k in (0, 1, 2)))
    log("  Temperature: t90(60 C)/t90(25 C) ~ D_a(25)/D_a(60) x tau ratio = %.2f ; t90(90 C)/t90(25 C) = %.2f "
        "(slab, x_s=0.9; note x_s=0.9 needs %.1fx / %.1fx the fugacity of 25 C)" % (
            TAU[("D", 333.15, 0.9, 0)] / D_alpha(333.15) / (TAU[("D", T0, 0.9, 0)] / Da),
            TAU[("D", 363.15, 0.9, 0)] / D_alpha(363.15) / (TAU[("D", T0, 0.9, 0)] / Da),
            fugacity(0.9, 333.15) / fugacity(0.9, T0), fugacity(0.9, 363.15) / fugacity(0.9, T0)))
    log("  Isotope: t90(H)/t90(D) at 25 C, slab x_s=0.9: %.2f" % (
        TAU[("H", T0, 0.9, 0)] / D_alpha(T0, "H") / (TAU[("D", T0, 0.9, 0)] / Da)))

    # figure
    fig, ax = plt.subplots(figsize=(7.5, 5))
    dd = np.logspace(-9, -2.5, 200)
    names = {0: "film / foil (l = thickness / half-thickness)", 1: "wire (l = radius)", 2: "particle (l = radius)"}
    for k, c in zip((0, 1, 2), ("C0", "C1", "C2")):
        ax.loglog(dd * 1e6, TAU[("D", T0, 0.9, k)] * dd ** 2 / Da, color=c, label=names[k] + ", diffusion-limited")
        ax.loglog(dd * 1e6, n_Pd * 0.9 * dd / (k + 1) * e / (10e-3 * 1e4), ":", color=c, lw=1)
    ax.loglog(dd * 1e6, TAU[("D", 363.15, 0.9, 0)] * dd ** 2 / D_alpha(363.15), "--", color="C0", lw=1,
              label="film/foil at 90 C")
    ax.loglog(dd * 1e6, TAU[("H", T0, 0.9, 0)] * dd ** 2 / D_alpha(T0, "H"), "-.", color="C0", lw=1,
              label="film/foil, H (25 C)")
    for y, l in ((1, "1 s"), (3600, "1 h"), (86400, "1 d"), (86400 * 30, "30 d")):
        ax.axhline(y, color="grey", lw=0.4)
        ax.text(1.2e-3, y * 1.2, l, fontsize=7)
    ax.set_xlabel("characteristic dimension l (µm)")
    ax.set_ylabel("t90 to 0.9·x_s (s)")
    ax.set_title("PdD loading time to 90 % of x_s=0.90, 25 C (dotted: charge-limited at η·i = 10 mA/cm²)",
                 fontsize=8)
    ax.legend(fontsize=7)
    ax.set_ylim(1e-10, 1e8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "m3_loading_time.png"), dpi=130)
    plt.close(fig)

    open(os.path.join(OUT, "m3_loading.txt"), "w").write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
