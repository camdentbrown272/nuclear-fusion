"""
M2 — loading model: i -> eta -> f -> x, validation, sensitivities, temperature,
recombination poisons, and the entry-side boundary condition table for M3
(including the effect of a permeation drain through a DFM membrane).

Run:  python3 sim/m2_loading.py
Writes docs/models/figs/m2_loading.txt, m2_isotherm.png, m2_loading_vs_i.png,
m2_sensitivity.png, m2_permeation.png, m2_m3_interface.txt
"""
import os
import sys
import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from m2_common import (FIGS, F, R, T0, SURFACES, Kinetics, Polarization, x_of_f, f_of_x,
                       half_ln_f, c_Pd, D_D_Pd, kB_eV)

OUT = []


def log(s=""):
    print(s)
    OUT.append(s)


I_GRID = np.array([0.001, 0.003, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0])

# validation bands (i range A/cm2, x range, label, grade)
VALID = [
    ((0.01, 0.3), (0.75, 0.82), "untreated Pd, 0.1 M LiOD [a]", "B"),
    ((0.05, 0.3), (0.91, 0.93), "vacuum-annealed + etched [a]", "B"),
    ((0.07, 0.09), (0.84, 0.86), "80 mA/cm2, 40 h -> 0.85 [b]", "C"),
    ((0.1, 0.6), (0.85, 0.95), "SRI/McKubre excess-heat cathodes [c]", "B"),
    ((0.005, 0.05), (0.78, 0.82), "Baranowski et al.: 80-150 bar equiv. [d]", "B"),
]


def isotherm_section():
    log("=" * 78)
    log("ISOTHERM beta-PdD (298 K): fugacity needed for a given loading")
    log("=" * 78)
    log("   x      f(atm) band [F90=3e3 | 1e4 | 3e4]     Delta mu_D = RT/2 ln f (eV)   eta_H = -(RT/2F) ln f (V)")
    for x in [0.70, 0.75, 0.80, 0.85, 0.88, 0.90, 0.92, 0.95, 0.97]:
        fs = [f_of_x(x, T0, F90) for F90 in [3e3, 1e4, 3e4]]
        dmu = 0.5 * kB_eV * T0 * np.log(fs[1])
        log(f"  {x:.2f}   {fs[0]:9.3g} | {fs[1]:9.3g} | {fs[2]:9.3g}       {dmu:+.3f}                      {-dmu:+.3f}")
    log("")
    log("Checks: x(f=1 atm)=%.3f (anchor), x(130 atm)=%.3f (Baranowski-Filipek 80-150 bar -> ~0.8), "
        "x(800 atm)=%.3f ('1 V ~ 800 atm', Berlinguette 2025), x(1e9 atm)=%.4f (GPa limit ~1)"
        % (x_of_f(1.0), x_of_f(130.0), x_of_f(800.0), x_of_f(1e9)))
    log("Temperature at fixed f=1e4 atm: " + ", ".join(f"{T-273.15:.0f} C: x={x_of_f(1e4, T):.3f}" for T in [278.15, 298.15, 318.15, 338.15, 358.15]))
    log("Temperature needed-f for x=0.90: " + ", ".join(f"{T-273.15:.0f} C: f={f_of_x(0.9, T):.3g} atm" for T in [278.15, 298.15, 318.15, 338.15]))

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    xs = np.linspace(0.62, 0.985, 300)
    for F90, ls in [(3e3, "--"), (1e4, "-"), (3e4, "--")]:
        ax[0].semilogx(f_of_x(xs, T0, F90), xs, ls, color="k", lw=1.5 if F90 == 1e4 else 0.8,
                       label=f"298 K, x=0.9 at {F90:.0e} atm" if F90 != 1e4 else "298 K central")
    for T, col in [(278.15, "b"), (318.15, "orange"), (338.15, "r")]:
        ax[0].semilogx(f_of_x(xs, T), xs, color=col, label=f"{T-273.15:.0f} C")
    ax[0].axhline(0.9, color="g", lw=.6); ax[0].axhline(0.95, color="g", lw=.6, ls=":")
    ax[0].plot([80, 150], [x_of_f(80), x_of_f(150)], "s", color="m", label="Baranowski-Filipek-Raczynski electrolytic equiv. (80-150 bar)")
    ax[0].plot([800], [x_of_f(800)], "D", color="c", label="'1 V ~ 800 atm' (Berlinguette 2025)")
    ax[0].set_xlabel("D fugacity f (atm)"); ax[0].set_ylabel("x = D/Pd"); ax[0].set_xlim(0.1, 1e9)
    ax[0].legend(fontsize=7); ax[0].grid(alpha=.3); ax[0].set_title("beta-PdD isotherm used by M2")
    ax2 = ax[1]
    ax2.plot(xs, 0.5 * kB_eV * T0 * np.log(f_of_x(xs)), "k")
    ax2.set_xlabel("x"); ax2.set_ylabel("Δμ_D vs ½D₂(1 atm) (eV)  = −η_H (V)")
    ax2.axvline(0.9, color="g", lw=.6); ax2.grid(alpha=.3)
    ax2.set_title("Only ~0.1-0.16 V of the cathode overpotential 'loads' the lattice")
    fig.tight_layout(); fig.savefig(os.path.join(FIGS, "m2_isotherm.png"), dpi=130); plt.close(fig)


def kinetics_section():
    log("")
    log("=" * 78)
    log("KINETIC LOADING MODEL: x vs current density (25 C, 0.1 M LiOD, no permeation drain)")
    log("=" * 78)
    pols = {k: Polarization(Kinetics(**v)) for k, v in SURFACES.items()}
    hdr = "  i(mA/cm2) " + " ".join(f"{k:>30s}" for k in pols)
    log(hdr)
    log("            " + " ".join(f"{'eta(V)  f(atm)   theta   x':>30s}" for _ in pols))
    for i in I_GRID:
        row = f"  {i*1e3:8.1f}  "
        for k, p in pols.items():
            eta = p.eta_of_i(i)
            st = p.kin.state(eta)
            row += f"   {eta:+.3f} {st['f']:9.3g} {st['theta']:.3f} {float(x_of_f(st['f'])):.3f}"
        log(row)
    # mechanism shares for typical
    log("")
    log("Mechanism shares ('typical'): fraction of adsorbed D leaving by Tafel vs Heyrovsky")
    k = Kinetics(**SURFACES["typical"])
    p = pols["typical"]
    for i in [0.001, 0.01, 0.1, 0.5, 1.0]:
        st = k.state(p.eta_of_i(i))
        log(f"   i={i*1e3:7.1f} mA/cm2: Tafel {100*st['frac_T']:.0f}%  Heyrovsky {100*st['frac_H']:.0f}%")
    # Tafel slope
    e1, e2 = p.eta_of_i(0.01), p.eta_of_i(0.1)
    e3, e4 = p.eta_of_i(0.1), p.eta_of_i(1.0)
    log(f"Tafel slope 10-100 mA/cm2: {1e3*(e1-e2):.0f} mV/dec; 0.1-1 A/cm2: {1e3*(e3-e4):.0f} mV/dec")
    # Nernst comparison
    log("")
    log("Nernstian fugacity f_N = exp(-2F eta/RT) that would follow if the whole overpotential loaded the lattice:")
    for i in [0.01, 0.1, 0.5]:
        eta = p.eta_of_i(i)
        log(f"   i={i*1e3:5.0f} mA/cm2: eta={eta:+.3f} V -> f_N={np.exp(-2*F*eta/(R*T0)):.2g} atm vs model f={p.f_of_i(i):.3g} atm"
            f"  (effective 'loading' overpotential eta_H={-(R*T0/(2*F))*np.log(p.f_of_i(i)):+.3f} V)")

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    ii = np.logspace(-3.3, 0.3, 200)
    cols = dict(poor="C3", typical="C0", good="C2", exceptional="C4")
    for kname, pol in pols.items():
        ax[0].semilogx(ii * 1e3, pol.x_of_i(ii), color=cols[kname], label=kname)
        ax[1].loglog(ii * 1e3, pol.f_of_i(ii), color=cols[kname], label=kname)
    ax[2].semilogx(ii * 1e3, pols["typical"].eta_of_i(ii), "k", label="total η (model)")
    ax[2].semilogx(ii * 1e3, -(R * T0 / (2 * F)) * np.log(pols["typical"].f_of_i(ii)), "C0--", label="η_H = −(RT/2F) ln f ('loading' part)")
    ax[2].semilogx(ii * 1e3, -(R * T0 / (2 * F)) * np.log(pols["good"].f_of_i(ii)), "C2--", label="η_H, good surface")
    for (ir, xr, lab, gr) in VALID:
        ax[0].add_patch(plt.Rectangle((ir[0] * 1e3, xr[0]), (ir[1] - ir[0]) * 1e3, xr[1] - xr[0], alpha=.18, color="gray"))
        ypos = xr[0] - 0.012 if "Baranowski" in lab else xr[1] + 0.003
        ax[0].text(ir[0] * 1e3, ypos, lab + f" ({gr})", fontsize=6)
    ax[0].axhline(0.9, color="k", lw=.6); ax[0].axhline(0.95, color="k", lw=.6, ls=":")
    ax[0].set_xlabel("cathode current density (mA/cm²)"); ax[0].set_ylabel("x = D/Pd (entry side)"); ax[0].set_ylim(0.65, 1.0)
    ax[0].legend(fontsize=8, loc="lower right"); ax[0].grid(alpha=.3); ax[0].set_title("Loading vs current density (surface classes)")
    ax[1].axhspan(3e3, 3e4, color="g", alpha=.12, label="f needed for x=0.9 (band)")
    ax[1].set_xlabel("i (mA/cm²)"); ax[1].set_ylabel("effective D fugacity (atm)"); ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
    ax[1].set_title("Effective fugacity")
    ax[2].set_xlabel("i (mA/cm²)"); ax[2].set_ylabel("V vs RHE(D)"); ax[2].legend(fontsize=8); ax[2].grid(alpha=.3)
    ax[2].set_title("Overpotential: most of it does not load")
    fig.tight_layout(); fig.savefig(os.path.join(FIGS, "m2_loading_vs_i.png"), dpi=130); plt.close(fig)
    return pols


def sensitivity_section():
    log("")
    log("=" * 78)
    log("SENSITIVITY of x(i) ('good' baseline) to model assumptions")
    log("=" * 78)
    base = dict(SURFACES["good"])
    variants = [("baseline", {}, 1e4),
                ("betaH=0.40", dict(betaH=0.40), 1e4), ("betaH=0.60", dict(betaH=0.60), 1e4),
                ("i0V x10", dict(i0V=3e-4), 1e4), ("i0V /10", dict(i0V=3e-6), 1e4),
                ("K=3e-4", dict(K=3e-4), 1e4), ("K=3e-2", dict(K=3e-2), 1e4),
                ("cT x3", dict(cT=3e-5), 1e4), ("cT /3", dict(cT=3.3e-6), 1e4),
                ("fmax x10", dict(fmax=3e6), 1e4), ("fmax /10", dict(fmax=3e4), 1e4),
                ("isotherm F90=3e3", {}, 3e3), ("isotherm F90=3e4", {}, 3e4)]
    res = {}
    log("   variant            " + " ".join(f"{i*1e3:>7.0f}" for i in [0.01, 0.05, 0.1, 0.3, 0.5, 1.0, 3.0]) + "  (mA/cm2)")
    for name, mod, F90 in variants:
        kw = dict(base); kw.update(mod)
        kin = Kinetics(**kw)
        eta, i, f, th = kin.table(eta_grid=-np.linspace(0.002, 1.6, 500))
        xs = x_of_f(f, T0, F90)
        row = [np.interp(np.log(q), np.log(i), xs) for q in [0.01, 0.05, 0.1, 0.3, 0.5, 1.0, 3.0]]
        res[name] = (i, xs)
        log(f"   {name:18s} " + " ".join(f"{v:7.3f}" for v in row))
    fig, ax = plt.subplots(figsize=(7, 4.6))
    for name, (i, xs) in res.items():
        ax.semilogx(i * 1e3, xs, lw=2 if name == "baseline" else 1, label=name)
    ax.set_xlim(1, 3000); ax.set_ylim(0.75, 1.0); ax.axhline(0.9, color="k", lw=.6)
    ax.set_xlabel("i (mA/cm²)"); ax.set_ylabel("x"); ax.legend(fontsize=7, ncol=2); ax.grid(alpha=.3)
    ax.set_title("Sensitivity of loading curve ('good' surface)")
    fig.tight_layout(); fig.savefig(os.path.join(FIGS, "m2_sensitivity.png"), dpi=130); plt.close(fig)


def poison_section():
    log("")
    log("=" * 78)
    log("RECOMBINATION POISONS / PROMOTERS (thiourea, S, As, CN-, ...): model mapping")
    log("  A poison that slows Tafel and Heyrovsky recombination by factor P multiplies the")
    log("  Volmer-Tafel fugacity by P and the Volmer-Heyrovsky ceiling by P^2.")
    log("=" * 78)
    for P in [1, 3, 10, 30, 100]:
        kw = dict(SURFACES["typical"]); kw["cT"] /= P; kw["fmax"] *= P**2
        pol = Polarization(Kinetics(**kw))
        log(f"   P={P:4d}: x(10)={float(pol.x_of_i(0.01)):.3f} x(100)={float(pol.x_of_i(0.1)):.3f} x(300)={float(pol.x_of_i(0.3)):.3f} x(1000 mA/cm2)={float(pol.x_of_i(1.0)):.3f}")


def temperature_section():
    log("")
    log("=" * 78)
    log("TEMPERATURE: x at fixed current (isotherm shift only; kinetic f(i) held at 25 C value ->")
    log("upper bound on x at higher T because recombination accelerates with T)")
    log("=" * 78)
    for kname in ["typical", "good"]:
        pol = Polarization(Kinetics(**SURFACES[kname]))
        for i in [0.1, 0.3]:
            f = pol.f_of_i(i)
            log(f"   {kname:8s} i={i*1e3:.0f} mA/cm2 (f={f:.3g} atm): " +
                ", ".join(f"{T-273.15:.0f}C x={x_of_f(f, T):.3f}" for T in [278.15, 288.15, 298.15, 313.15, 333.15, 353.15]))


# ------------------------------------------------------------ M3 interface
def f_in_with_drain(kin, i, ip):
    """Entry-side fugacity at cathodic current density i when a net absorption flux
    J_abs = ip/F (ip in A/cm^2 equivalent) enters the metal. Solve for eta such that
    total i matches, with the steady surface balance v_V = v_H + 2v_T + J_abs."""
    J = ip / F

    def g(eta):
        st = kin.state(eta, J_abs=J)
        return np.log(st["i"] / i) if np.isfinite(st["i"]) and st["i"] > 0 else -50

    try:
        eta = brentq(g, -1.8, -1e-4)
    except ValueError:
        return np.nan
    return kin.state(eta, J_abs=J)["f"]


def m3_interface():
    log("")
    log("=" * 78)
    log("M3 INTERFACE: entry-side x_in vs current density and permeation drain")
    log("  ip = F * J_perm (A/cm2 equivalent of the D flux leaving through the membrane)")
    log("=" * 78)
    lines = ["# M2 -> M3 entry-side boundary condition. 25 C, 0.1 M LiOD, isotherm F90=1e4 atm.",
             "# columns: surface_class  i(A/cm2)  ip/i  f_in(atm)  x_in",
             "# ip = F*J_perm is the permeation (absorption) flux expressed as current density.",
             "# Use: M3 computes J_perm(x_in) from diffusion + exit-face kinetics and iterates with this table.",
             "# Standalone function: sim/m2_loading.py: f_in_with_drain(Kinetics(**SURFACES[c]), i, ip)"]
    for kname in ["poor", "typical", "good", "exceptional"]:
        kin = Kinetics(**SURFACES[kname])
        log(f"  {kname}:  x_in for ip/i = 0, 0.01, 0.1, 0.3, 0.5, 0.8")
        for i in [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0]:
            xs = []
            for fr in [0.0, 0.01, 0.1, 0.3, 0.5, 0.8]:
                f = f_in_with_drain(kin, i, fr * i)
                x = float(x_of_f(f)) if np.isfinite(f) else np.nan
                xs.append(x)
                lines.append(f"{kname:12s} {i:7.3f} {fr:5.2f} {f:11.4g} {x:.4f}")
            log(f"     i={i*1e3:6.0f} mA/cm2: " + " ".join(f"{x:.3f}" for x in xs))
    with open(os.path.join(FIGS, "m2_m3_interface.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


def permeation_section():
    """Illustrative coupled DFM steady state: entry kinetics + linear diffusion + exit
    desorption J = k_des * f_out (vacuum side). k_des spans clean Pd -> barrier-coated."""
    log("")
    log("=" * 78)
    log("DFM STEADY STATE WITH PERMEATION (illustrative; M3 owns the full model)")
    log("  J = D c_Pd (x_in - x_out)/L = k_des f_out ;  entry: f_in(i, ip=F J)")
    log(f"  D = {D_D_Pd:.1e} cm2/s (x2 uncertain), c_Pd = {c_Pd:.4f} mol/cm3")
    # clean-Pd scale for k_des: detailed balance with gas at fugacity f: J = 2 s Z(f)
    m = 4.028 * 1.66054e-27
    Z = 101325 / np.sqrt(2 * np.pi * m * 1.380649e-23 * 298.15) / 6.02214e23 * 1e-4   # mol D2 /cm2/s/atm
    log(f"  D2 impingement at 1 atm: Z = {Z:.3g} mol D2/cm2/s; clean-Pd k_des ~ 2 s Z with s~1e-2..1"
        f" -> {2*1e-2*Z:.2g}..{2*Z:.2g} mol D/cm2/s/atm (upper bound; high-coverage s is lower)")
    log("=" * 78)
    kin = Kinetics(**SURFACES["good"])
    kdes = np.logspace(-15, -3, 49)
    res = {}
    for L_um in [25, 50, 100]:
        for i in [0.1, 0.3]:
            L = L_um * 1e-4
            xin_l, xout_l, J_l = [], [], []
            for kd in kdes:
                def resid(lJ):
                    J = np.exp(lJ)
                    f_in = f_in_with_drain(kin, i, F * J)
                    if not np.isfinite(f_in):
                        return 50.0
                    x_in = float(x_of_f(f_in))
                    x_out = x_in - J * L / (D_D_Pd * c_Pd)
                    if x_out <= 1e-3:
                        return 50.0
                    return np.log(kd * f_of_x(x_out)) - lJ
                x0 = float(x_of_f(f_in_with_drain(kin, i, 0.0)))
                lo, hi = np.log(1e-16), np.log(min(0.999 * i / F, 0.999 * x0 * D_D_Pd * c_Pd / L))
                grid = np.linspace(lo, hi, 60)
                vals = np.array([resid(q) for q in grid])
                sc = np.where((vals[:-1] > 0) & (vals[1:] <= 0))[0]
                try:
                    lJ = brentq(resid, grid[sc[0]], grid[sc[0] + 1])
                except (ValueError, IndexError):
                    xin_l.append(np.nan); J_l.append(np.nan); xout_l.append(np.nan)
                    continue
                J = np.exp(lJ)
                f_in = f_in_with_drain(kin, i, F * J)
                x_in = float(x_of_f(f_in))
                xin_l.append(x_in); J_l.append(J); xout_l.append(x_in - J * L / (D_D_Pd * c_Pd))
            res[(L_um, i)] = (np.array(xin_l), np.array(xout_l), np.array(J_l))
    for (L_um, i), (xi, xo, J) in res.items():
        log(f"  L={L_um:3d} um, i={i*1e3:.0f} mA/cm2:")
        for kd in [1e-13, 1e-11, 1e-9, 1e-7, 1e-5, 1e-3]:
            k = np.argmin(np.abs(np.log(kdes / kd)))
            log(f"     k_des={kdes[k]:.0e}: x_in={xi[k]:.3f} x_out={xo[k]:.3f} J={J[k]:.3g} mol/cm2/s"
                f" (={F*J[k]/i*100:.1f}% of i, {J[k]*6.022e23:.2g} D/cm2/s)")
    # required k_des for x_out >= 0.85 / 0.90
    for (L_um, i), (xi, xo, J) in res.items():
        s = []
        for thr in [0.85, 0.90]:
            ok = np.where(xo >= thr)[0]
            s.append(f"x_out>={thr}: k_des<={kdes[ok[-1]]:.1e}, J={J[ok[-1]]:.2g}" if ok.size else f"x_out>={thr}: not reached")
        log(f"  L={L_um} um i={i*1e3:.0f}: " + " ; ".join(s))
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    for (L_um, i), (xi, xo, J) in res.items():
        ls = "-" if i == 0.3 else "--"
        ax[0].semilogx(kdes, xo, ls, label=f"x_out, L={L_um} µm, i={i*1e3:.0f}")
        ax[1].loglog(kdes, J * 6.022e23, ls, label=f"L={L_um} µm, i={i*1e3:.0f}")
    xi, xo, J = res[(50, 0.3)]
    ax[0].semilogx(kdes, xi, "k:", label="x_in (L=50, i=300)")
    ax[0].axhline(0.9, color="g", lw=.6); ax[0].axvspan(2e-2 * Z, 2 * Z, color="r", alpha=.1)
    ax[0].text(3e-4, 0.3, "clean Pd\n(est.)", fontsize=7)
    ax[0].set_xlabel("exit-face desorption coefficient k_des (mol D cm⁻² s⁻¹ atm⁻¹)"); ax[0].set_ylabel("loading")
    ax[0].legend(fontsize=6.5); ax[0].grid(alpha=.3); ax[0].set_title("DFM ('good' entry surface): exit-face loading")
    ax[1].set_xlabel("k_des"); ax[1].set_ylabel("permeation flux (D cm⁻² s⁻¹)"); ax[1].legend(fontsize=7); ax[1].grid(alpha=.3)
    ax[1].set_title("Flux to the vacuum/detector side")
    fig.tight_layout(); fig.savefig(os.path.join(FIGS, "m2_permeation.png"), dpi=130); plt.close(fig)


if __name__ == "__main__":
    isotherm_section()
    kinetics_section()
    sensitivity_section()
    poison_section()
    temperature_section()
    m3_interface()
    permeation_section()
    with open(os.path.join(FIGS, "m2_loading.txt"), "w") as fh:
        fh.write("\n".join(OUT) + "\n")
