"""
M3 / Q2: detector-facing membrane (DFM, config C3).  Back face held at x_in (electrochemical, from M2),
front (exit) face in vacuum with a recombination / barrier law.  Steady state (exact via the Kirchhoff
transform), design chart thickness x k_r, comparison of exit-face finishes, and transients
(first loading, current interruption, current modulation).

Run:  python3 sim/m3_membrane.py
Writes docs/models/figs/m3_dfm_chart.png, m3_dfm_profiles.png, m3_dfm_transients.png, m3_membrane.txt
"""
import os
import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from m3_common import (Material, FV1D, gap, fugacity, c_eff, J_pick, J_barrier, kr_Pd_lit, kr_Pd_baskes,
                       Jsat_desorb, G_overlayer, n_Pd, e, NA, BAR, ATM, x_of_f, D_alpha, kB_eV)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "docs", "models", "figs")
LOG = []
T0 = 298.15


def log(s=""):
    print(s)
    LOG.append(s)


def mAcm2(J):
    """D atoms m^-2 s^-1 -> equivalent current density mA/cm^2."""
    return J * e * 0.1


def steady(mat, x_in, L, Jlaw):
    """Solve (Phi(x_in) - Phi(x_e))/L = J(x_e).  Returns x_e, J."""
    Pin = mat.Phi(x_in)
    g = lambda xe: (Pin - mat.Phi(xe)) / L - Jlaw(xe)
    if g(x_in - 1e-9) > 0:     # even at x_e = x_in the law passes less than zero-gradient flux -> x_e ~ x_in
        return x_in, Jlaw(x_in)
    if g(1e-9) < 0:
        return 0.0, Pin / L
    xe = brentq(g, 1e-9, x_in - 1e-9, xtol=1e-7)
    return xe, Jlaw(xe)


def profile(mat, x_in, L, J, z):
    return mat.x_of_Phi(mat.Phi(x_in) - J * z)


def main():
    mat = Material(T0, "D")
    a, b = gap(T0, "D")
    Da = D_alpha(T0)
    log("=== DFM exit-face kinetics: reference values at 25 C (PdD) ===")
    ce9 = c_eff(0.9, T0)
    log(f"c_eff(x=0.90) = {ce9:.3g} m^-3 (activity-form Pick concentration; c_eff = c in alpha phase)")
    log(f"k_r(Pd, real surface, lit.) = {kr_Pd_lit(T0):.3g} m^4/s ; Baskes clean-surface bound (s0=0.5) = "
        f"{kr_Pd_baskes(T0):.3g} m^4/s")
    # Iwamura calibration of the saturated desorption cap
    T_I = 343.15
    J_I = 2.0 / 22414.0 / 60 * 2 * NA / 6.25e-4      # 2 sccm D2 through 25 mm x 25 mm [web] (area [mem])
    x_I = x_of_f(1.0 * ATM, T_I)
    mI = Material(T_I, "D")
    Jdiff_I = mI.Phi(x_I) / 100e-6
    Ed_max = -kB_eV * T_I * np.log(J_I / (2 * 1e13 * 1.53e19))
    log(f"Iwamura (D2 1 atm, 70 C, 0.1 mm Pd, 2 sccm over 6.25 cm2): J = {J_I:.3g} D/m2/s = {mAcm2(J_I):.1f} mA/cm2")
    log(f"  model: x_in(1 atm, 70 C) = {x_I:.3f}; diffusion-limited capacity (exit x=0) = {Jdiff_I:.3g} "
        f"D/m2/s -> observed/capacity = {J_I / Jdiff_I:.2f}")
    log(f"  => a bare Pd exit face in vacuum passed >= {J_I:.2g} D/m2/s at 343 K; with a saturated-surface cap "
        f"2 nu N_s exp(-Ed/kT), Ed <= {Ed_max:.3f} eV; extrapolated J_sat(25 C) >= {Jsat_desorb(T0, Ed_max):.3g} "
        f"D/m2/s ({mAcm2(Jsat_desorb(T0, Ed_max)):.2f} mA/cm2)")

    # ---------------------------------------------------- capacity table: max flux with exit >= 0.90
    log("")
    log("=== Max steady flux compatible with x >= 0.90 everywhere: J_max = [Phi(x_in) - Phi(0.90)] / L ===")
    Ls = np.array([5, 10, 25, 50, 100, 200]) * 1e-6
    log("  x_in  | " + " | ".join(f"L={L * 1e6:5.0f} um" for L in Ls) + "   (D/m2/s ; mA/cm2)")
    JMAX = {}
    for xin in (0.85, 0.90, 0.95, 0.97):
        row = []
        for L in Ls:
            Jm = max(mat.Phi(xin) - mat.Phi(0.9), 0) / L
            JMAX[(xin, L)] = Jm
            row.append(f"{Jm:8.2g};{mAcm2(Jm):6.2f}")
        log(f"  {xin:.2f}  | " + " | ".join(row))
    log("  -> x_in = 0.90 or below leaves no flux budget: the exit face can only be >= 0.90 if x_in > 0.90.")
    # also for x_exit >= 0.85 and >= beta_min+0.05
    for xt in (0.85, 0.80):
        Jm = (mat.Phi(0.95) - mat.Phi(xt)) / 50e-6
        log(f"  (x_in=0.95, L=50 um: J_max for exit >= {xt:.2f} is {Jm:.3g} D/m2/s = {mAcm2(Jm):.1f} mA/cm2)")

    # ---------------------------------------------------- design chart
    J_target = 1.0e21
    Lg = np.logspace(np.log10(5e-6), np.log10(200e-6), 70)
    kg = np.logspace(-42, -28, 90)
    fig, axs = plt.subplots(2, 2, figsize=(12, 9.5), sharex=True, sharey=True)
    log("")
    log(f"=== Design chart: exit loading & flux vs L and k_r (Pick, activity form, no cap), 25 C. "
        f"Target: x_exit >= 0.90 and J >= {J_target:.0e} D/m2/s ({mAcm2(J_target):.0f} mA/cm2) ===")
    for ax, xin in zip(axs.ravel(), (0.85, 0.90, 0.95, 0.97)):
        XE = np.zeros((len(kg), len(Lg)))
        JJ = np.zeros_like(XE)
        for i, kr in enumerate(kg):
            law = lambda x, kr=kr: J_pick(x, T0, kr)
            for j, L in enumerate(Lg):
                XE[i, j], JJ[i, j] = steady(mat, xin, L, law)
        cs = ax.contourf(Lg * 1e6, kg, XE, levels=[0, 0.3, 0.587, 0.7, 0.8, 0.85, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98],
                         cmap="viridis")
        c2 = ax.contour(Lg * 1e6, kg, mAcm2(JJ), levels=[0.1, 1, 3, 10, 30, 100, 300], colors="w",
                        linewidths=0.7)
        ax.clabel(c2, fmt="%g mA/cm²", fontsize=6)
        ok = (XE >= 0.9) & (JJ >= J_target)
        if ok.any():
            ax.contour(Lg * 1e6, kg, ok.astype(float), levels=[0.5], colors="r", linewidths=2)
        for nm, G, c in (("Ni 2 nm", G_overlayer("Ni", 2e-9, T0), "cyan"), ("Ni 20 nm", G_overlayer("Ni", 20e-9, T0), "cyan")):
            keq = J_barrier(0.9, T0, G) / ce9 ** 2     # k_r giving the same flux at x_exit = 0.90
            ax.axhline(keq, color=c, ls=":", lw=1)
            ax.text(60, keq * 1.8, f"≈ {nm} overlayer", color=c, fontsize=7)
        ax.axhline(kr_Pd_lit(T0), color="orange", ls="--", lw=1)
        ax.text(5.5, kr_Pd_lit(T0) * 1.8, "bare Pd, lit. k_r (no cap)", color="orange", fontsize=7)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(f"x_in = {xin:.2f}   (fill: exit loading; white: flux; red: recommended)", fontsize=8)
        fig.colorbar(cs, ax=ax, label="x at exit face")
        # required k_r for target flux with exit at 0.90..
        if ok.any():
            kmin = kg[np.where(ok.any(axis=1))[0].min()]
            kmax = kg[np.where(ok.any(axis=1))[0].max()]
            Lmax = Lg[np.where(ok.any(axis=0))[0].max()]
            log(f"  x_in={xin:.2f}: recommended region exists; k_r in [{kmin:.2g}, {kmax:.2g}] m^4/s, "
                f"L <= {Lmax * 1e6:.0f} um")
        else:
            log(f"  x_in={xin:.2f}: NO (L, k_r) satisfies both x_exit >= 0.90 and J >= target")
    for ax in axs[1]:
        ax.set_xlabel("membrane thickness L (µm)")
    for ax in axs[:, 0]:
        ax.set_ylabel("exit-face k_r (m⁴/s), Pick activity form")
    fig.suptitle("DFM steady state, PdD 25 °C: exit loading and permeation flux vs thickness and exit k_r",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "m3_dfm_chart.png"), dpi=130)
    plt.close(fig)
    # required k_r formula
    log("  Required k_r to deliver J with exit at exactly 0.90: k_r* = J / c_eff(0.90)^2 = "
        f"{J_target / ce9 ** 2:.3g} m^4/s for J = 1e21; scale linearly with J.")
    log(f"  k_r* / k_r(bare Pd, lit.) = {J_target / ce9 ** 2 / kr_Pd_lit(T0):.2g}; "
        f"k_r* / k_r(Baskes) = {J_target / ce9 ** 2 / kr_Pd_baskes(T0):.2g}")

    # ---------------------------------------------------- exit-face finishes compared
    log("")
    log("=== Exit-face finishes, x_in = 0.95, 25 C (steady state) ===")
    Jsat_lo = Jsat_desorb(T0, Ed_max)
    finishes = {
        "bare Pd, lit k_r, no cap": lambda x: J_pick(x, T0, kr_Pd_lit(T0)),
        "bare Pd, lit k_r + J_sat cap (Ed=%.2f eV)" % Ed_max: lambda x: J_pick(x, T0, kr_Pd_lit(T0), Jsat=Jsat_lo),
        "bare Pd, lit k_r + J_sat cap x100": lambda x: J_pick(x, T0, kr_Pd_lit(T0), Jsat=100 * Jsat_lo),
        "PdO (k_r x1e-3 of bare), no cap": lambda x: J_pick(x, T0, 1e-3 * kr_Pd_lit(T0)),
        "Au 20 nm (perm <= Cu)": lambda x: J_barrier(x, T0, G_overlayer("Au_upper", 20e-9, T0)),
        "Cu 20 nm": lambda x: J_barrier(x, T0, G_overlayer("Cu", 20e-9, T0)),
        "Ni 2 nm": lambda x: J_barrier(x, T0, G_overlayer("Ni", 2e-9, T0)),
        "Ni 20 nm": lambda x: J_barrier(x, T0, G_overlayer("Ni", 20e-9, T0)),
        "Ni/Cu 5x(10+10 nm)": lambda x: J_barrier(x, T0, 1 / (5 / G_overlayer("Ni", 10e-9, T0) +
                                                               5 / G_overlayer("Cu", 10e-9, T0))),
        "Pd/CaO (Iwamura) >= G_I, cap": lambda x: min(J_barrier(x, T0, J_I / np.sqrt(1.013e5)), Jsat_lo * 100),
    }
    log("  finish                                      L(um)  x_exit   J (D/m2/s)  J (mA/cm2)  min x >= 0.9?")
    rows = {}
    for name, law in finishes.items():
        for L in (25e-6, 50e-6, 100e-6):
            xe, J = steady(mat, 0.95, L, law)
            rows[(name, L)] = (xe, J)
            log(f"  {name:42s} {L * 1e6:5.0f}  {xe:6.3f}  {J:10.3g}  {mAcm2(J):9.3g}   {'yes' if xe >= 0.9 else 'no'}")
    # PdO lifetime
    nO = 8.3e3 / 0.1224 * NA            # O atoms / m^3 in PdO (rho 8.3 g/cm3, M 122.4 g/mol)
    for d in (1e-9, 5e-9):
        for J, fr in ((1e19, 0.1), (1e21, 0.1)):
            tl = nO * d / (2 * fr * J)
            log(f"  PdO {d * 1e9:.0f} nm reduced by permeating D (2 D per O, {fr * 100:.0f}% of flux reacting) "
                f"at J={J:.0e}: lifetime ~ {tl:.3g} s")
    # tuneable partial-coverage Au mask: J ~ phi * J_bare  (exit law), required open fraction
    log("  Patterned Au mask (open fraction phi) scales the bare-Pd exit law by ~phi (valid when the opening")
    log("  pitch << L so lateral diffusion equalises): phi required for J = 1e21 at exit 0.90:")
    for nm, law in (("lit k_r no cap", lambda x: J_pick(x, T0, kr_Pd_lit(T0))),
                    ("lit k_r + cap x100", lambda x: J_pick(x, T0, kr_Pd_lit(T0), Jsat=100 * Jsat_lo))):
        log(f"    bare law {nm}: J_bare(0.90) = {law(0.9):.3g} -> phi = {min(1, J_target / law(0.9)):.2g}")

    # ---------------------------------------------------- profiles figure
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    L = 50e-6
    z = np.linspace(0, L, 300)
    for name in ("bare Pd, lit k_r, no cap", "bare Pd, lit k_r + J_sat cap x100", "Ni 2 nm", "Au 20 nm (perm <= Cu)"):
        xe, J = rows[(name, L)]
        ax[0].plot(z * 1e6, profile(mat, 0.95, L, J, z), label=f"{name}: J={mAcm2(J):.3g} mA/cm²")
    krs = J_target / ce9 ** 2
    xe, J = steady(mat, 0.95, L, lambda x: J_pick(x, T0, krs))
    ax[0].plot(z * 1e6, profile(mat, 0.95, L, J, z), "k--", lw=2, label=f"design k_r*={krs:.1e}: J={mAcm2(J):.3g}")
    ax[0].axhline(0.9, color="grey", lw=0.5)
    ax[0].axhspan(a, b, color="grey", alpha=0.15)
    ax[0].set_xlabel("depth from back (electrolyte) face (µm)")
    ax[0].set_ylabel("x")
    ax[0].set_title("Steady profiles, L = 50 µm, x_in = 0.95, 25 °C", fontsize=9)
    ax[0].legend(fontsize=6)
    # flux-thickness capacity plot
    Lp = np.logspace(np.log10(5e-6), np.log10(200e-6), 100)
    for xin in (0.92, 0.95, 0.97):
        ax[1].loglog(Lp * 1e6, mAcm2((mat.Phi(xin) - mat.Phi(0.9)) / Lp), label=f"x_in = {xin}")
    ax[1].axhline(mAcm2(J_I), color="r", ls=":", label="Iwamura flux (70 °C, 2 sccm)")
    ax[1].axhline(mAcm2(J_target), color="k", ls=":", label="target 1e21 D/m²/s")
    ax[1].set_xlabel("L (µm)")
    ax[1].set_ylabel("max permeation (mA/cm² equiv) with x ≥ 0.90 throughout")
    ax[1].legend(fontsize=7)
    ax[1].set_title("Flux-thickness budget: J·L ≤ Φ(x_in) − Φ(0.90)", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "m3_dfm_profiles.png"), dpi=130)
    plt.close(fig)

    # ---------------------------------------------------- verification: transient -> steady
    log("")
    log("=== Verification V4: FV transient converges to the Kirchhoff steady state ===")
    krs = J_target / ce9 ** 2
    for L in (25e-6, 50e-6):
        xe_s, J_s = steady(mat, 0.95, L, lambda x: J_pick(x, T0, krs))
        s = FV1D(mat, L, N=160, geom=0, refine=(True, True), ratio=1.04)
        t, x = s.run(np.zeros(s.N), 400 * L ** 2 / Da, ("dir", 0.95), ("flux", lambda x: J_pick(x, T0, krs)),
                     dt0=1e-6, dtmax=2 * L ** 2 / Da)
        log(f"  L={L * 1e6:.0f} um: steady x_e = {xe_s:.4f}, FV x(last cell) = {x[-1]:.4f}; "
            f"J steady {J_s:.4g}, FV exit flux {J_pick(x[-1], T0, krs):.4g}")

    # ---------------------------------------------------- transients
    log("")
    log("=== Transients (design exit law k_r* for 1e21 at 0.90; x_in = 0.95; 25 C) ===")
    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    law = lambda x: J_pick(x, T0, krs)
    for L, c in zip((10e-6, 25e-6, 50e-6, 100e-6), ("C0", "C1", "C2", "C3")):
        s = FV1D(mat, L, N=120, geom=0, refine=(True, True), ratio=1.05)
        th, xh, mh = [], [], []

        def rec(t, x):
            th.append(t)
            xh.append(x[-1])
            mh.append(x.min())
        s.run(np.zeros(s.N), 20 * L ** 2 / Da, ("dir", 0.95), ("flux", law), dt0=1e-6, dtmax=0.05 * L ** 2 / Da,
              record=rec)
        th, xh, mh = map(np.array, (th, xh, mh))
        tb = th[np.argmax(xh > b)] if (xh > b).any() else np.nan
        t9 = th[np.argmax(mh >= 0.9)] if (mh >= 0.9).any() else np.nan
        log(f"  first loading L={L * 1e6:4.0f} um: exit face enters beta at {tb:.3g} s; whole membrane >= 0.90 "
            f"at {t9:.3g} s")
        ax[0].semilogx(th, xh, color=c, label=f"exit, L={L * 1e6:.0f} µm")
        ax[0].semilogx(th, mh, ":", color=c)
    ax[0].axhline(0.9, color="grey", lw=0.5)
    ax[0].set_xlabel("t (s)")
    ax[0].set_ylabel("x at exit face (solid), min x (dotted)")
    ax[0].set_title("First loading, back face stepped to 0.95", fontsize=9)
    ax[0].legend(fontsize=7)

    # current interruption: back face -> no flux (open circuit, optimistic) and exit keeps desorbing
    log("  Current interruption (back face no-flux; exit desorbs):")
    for lawname, lw in (("design k_r*", law),
                        ("bare Pd lit, cap x100", lambda x: J_pick(x, T0, kr_Pd_lit(T0), Jsat=100 * Jsat_lo)),
                        ("bare Pd lit, cap (lower bound)", lambda x: J_pick(x, T0, kr_Pd_lit(T0), Jsat=Jsat_lo))):
        for L in (25e-6, 50e-6):
            xe_s, J_s = steady(mat, 0.95, L, lw)
            s = FV1D(mat, L, N=120, geom=0, refine=(True, True), ratio=1.05)
            z = s.rc
            x0 = profile(mat, 0.95, L, J_s, z)
            th, xe_h = [], []

            def rec(t, x):
                th.append(t)
                xe_h.append(x[-1])
            tend = 20 * n_Pd * L * 0.4 / max(lw(0.7), 1e10)
            s.run(x0, min(tend, 3e7), ("sym",), ("flux", lw), dt0=1e-4, dtmax=min(tend, 3e7) / 400,
                  record=rec, stop=lambda t, x: x[-1] < 0.3)
            th, xe_h = np.array(th), np.array(xe_h)
            t90 = th[np.argmax(xe_h < 0.9)] if (xe_h < 0.9).any() else np.nan
            tbm = th[np.argmax(xe_h < b)] if (xe_h < b).any() else np.nan
            log(f"    {lawname:30s} L={L * 1e6:3.0f} um: exit < 0.90 after {t90:.3g} s; exit enters "
                f"two-phase (x<{b:.3f}) after {tbm:.3g} s")
            if L == 50e-6:
                ax[1].loglog(th, xe_h, label=lawname)
    ax[1].axhline(b, color="grey", lw=0.5)
    ax[1].set_xlabel("t after current loss (s)")
    ax[1].set_ylabel("exit-face x")
    ax[1].set_title("Current interruption, L = 50 µm", fontsize=9)
    ax[1].legend(fontsize=7)

    # modulation: square wave on x_in between 0.90 and 0.97
    log("  Back-face square-wave modulation 0.90 <-> 0.97, L = 50 um, design exit law:")
    L = 50e-6
    tauD = L ** 2 / float(np.interp(0.93, mat.xg, mat.Dg))
    log(f"    diffusion time L^2/D_chem(0.93) = {tauD:.3g} s")
    per, amp_x, amp_J = [], [], []
    for P in (0.3, 1, 3, 10, 30, 100, 300):
        s = FV1D(mat, L, N=120, geom=0, refine=(True, True), ratio=1.05)
        xe0, J0 = steady(mat, 0.935, L, law)
        x0 = profile(mat, 0.935, L, J0, s.rc)
        rec_t, rec_x = [], []

        def rec(t, x):
            rec_t.append(t)
            rec_x.append(x[-1])
        left = lambda t, P=P: ("dir", 0.97 if (t % P) < P / 2 else 0.90)
        s.run(x0, 6 * P + 5 * tauD, left, ("flux", law), dt0=P / 400, dtmax=P / 60, record=rec)
        rt, rx = np.array(rec_t), np.array(rec_x)
        sel = rt > rt[-1] - 2 * P
        per.append(P)
        amp_x.append(rx[sel].max() - rx[sel].min())
        amp_J.append(law(rx[sel].max()) - law(rx[sel].min()))
        log(f"    period {P:6.1f} s: exit x swing {amp_x[-1]:.4f}; exit flux swing {mAcm2(amp_J[-1]):.3g} mA/cm2")
    ax[2].loglog(per, amp_x, "o-", label="exit x swing (peak-peak)")
    ax[2].axvline(tauD, color="grey", ls=":")
    ax[2].text(tauD * 1.1, min(amp_x) * 1.2, "L²/D", fontsize=7)
    ax[2].set_xlabel("modulation period (s)")
    ax[2].set_ylabel("exit-face x swing")
    ax[2].set_title("Transfer of back-face modulation to exit face (L=50 µm)", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "m3_dfm_transients.png"), dpi=130)
    plt.close(fig)

    # masked rim: lateral loading length
    for t in (3600, 86400, 7 * 86400):
        l = np.sqrt(float(np.interp(0.9, mat.xg, mat.Dg)) * t)
        log(f"  Lateral diffusion under a masked rim: sqrt(D_chem(0.9) t) = {l * 1e3:.2f} mm after {t / 3600:.0f} h")

    open(os.path.join(OUT, "m3_membrane.txt"), "w").write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
