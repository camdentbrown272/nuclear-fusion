"""
M3 / Q4: loading protocols for the DFM.
  A  : high-loading steady state that avoids cracking (first loading, hold, keep-alive).
  B1 : controlled non-equilibrium inside the beta phase (flux / gradient pumping, no phase change).
  B2 : deliberate alpha/beta cycling (fresh crack faces, dislocations, deformation vacancies).

Back face: galvanostatic charging/stripping as flux laws (with a fugacity ceiling x_c set by the current;
M2 supplies the true x_in(i)), exit face: the design recombination law from m3_membrane.py.
Stresses: free (floating) plate, eigenstrain eta*x(z): the mean and linear parts are relieved (expansion +
bending), the remainder is the self-equilibrated stress.  Elastic-perfectly-plastic shakedown gives the
plastic strain range per cycle; that drives crack initiation (Coffin-Manson), short-crack growth (Tomkins)
and deformation-vacancy production.

Run:  python3 sim/m3_cycling.py
Writes docs/models/figs/m3_cycling.png, m3_cycling_damage.png, m3_cycling.txt
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from m3_common import Material, FV1D, gap, c_eff, J_pick, n_Pd, e, kB_eV, D_alpha
from m3_mechanics import E_Pd, NU_Pd, ETA, SY

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "docs", "models", "figs")
LOG = []
T0 = 298.15


def log(s=""):
    print(s)
    LOG.append(s)


def mA(J):
    return J * e * 0.1


def plate_mismatch(x, z):
    """Self-equilibrated eigenstrain mismatch of a free plate: eps*(z) minus its best linear fit."""
    es = ETA * x
    w = np.gradient(z)
    zc = z - np.sum(z * w) / np.sum(w)
    m0 = np.sum(es * w) / np.sum(w)
    m1 = np.sum(es * zc * w) / np.sum(zc * zc * w)
    return es - m0 - m1 * zc


def run_protocol(name, L, left_of_t, t_end, exit_law, x0, N=120, dtmax=None, sample=None):
    mat = Material(T0, "D")
    s = FV1D(mat, L, N=N, geom=0, refine=(True, True), ratio=1.05)
    s.chg_lo, s.chg_hi = 0.01, 0.03
    if np.isscalar(x0):
        x0 = np.full(s.N, x0)
    H = dict(t=[], mean=[], back=[], exit=[], mis=[], Jin=[])

    def rec(t, x):
        H["t"].append(t)
        H["mean"].append(s.mean(x))
        H["back"].append(x[0])
        H["exit"].append(x[-1])
        H["mis"].append(plate_mismatch(x, s.rc))
        lb = left_of_t(t)
        H["Jin"].append(-lb[1](x[0]) if lb[0] == "flux" else np.nan)
    s.run(x0, t_end, left_of_t, ("flux", exit_law), dt0=1e-3, dtmax=dtmax or t_end / 2000, record=rec)
    for k in H:
        H[k] = np.array(H[k])
    H["z"] = s.rc
    return H


def damage(H, t0, t1, sy):
    """Plastic strain range per cycle from the mismatch history within [t0, t1] (one cycle)."""
    sel = (H["t"] >= t0) & (H["t"] <= t1)
    mis = H["mis"][sel]
    rng = mis.max(axis=0) - mis.min(axis=0)
    ey = sy * (1 - NU_Pd) / E_Pd
    dep = np.maximum(0, rng - 2 * ey)
    sig_el = E_Pd / (1 - NU_Pd) * np.abs(mis).max()
    return dep.max(), dep.mean(), sig_el


def main():
    a, b = gap(T0, "D")
    krs = 1e21 / c_eff(0.9, T0) ** 2
    exit_law = lambda x: J_pick(x, T0, krs)
    L = 25e-6
    Dch = float(np.interp(0.9, Material(T0).xg, Material(T0).Dg))
    log(f"DFM L = {L * 1e6:.0f} um, exit law k_r* = {krs:.3g} m^4/s (1e21 D/m2/s at x_exit = 0.90); "
        f"tau_D = L^2/D_chem(0.9) = {L ** 2 / Dch:.2f} s")

    def charge(Jc, xc):     # inward galvanostatic flux, tapering to zero at the ceiling x_c (fugacity limit)
        return lambda x: -Jc * np.clip((xc - x) / 0.004, 0, 1)

    def strip(Ja, xmin=0.01):   # outward anodic stripping flux, stops at x_min (potential limit)
        return lambda x: Ja * np.clip((x - xmin) / 0.004, 0, 1)

    fig, ax = plt.subplots(3, 1, figsize=(10, 10))
    # ---------------------------------------------------------------- Protocol A
    log("")
    log("=== Protocol A: steady high loading (crack-avoiding) ===")
    ramp = [(0, 5e-3), (1800, 20e-3), (3600, 50e-3), (5400, 100e-3), (7200, 300e-3)]   # (t, A/cm2)
    xc_of_i = lambda i: 0.85 + 0.10 * np.clip(np.log10(i / 0.01) / np.log10(30), 0, 1)  # placeholder until M2
    eta_c = 0.5

    def leftA(t):
        i = [c for tt, c in ramp if t >= tt][-1]
        return ("flux", charge(eta_c * i * 1e4 / e, xc_of_i(i)))
    HA = run_protocol("A", L, leftA, 9000, exit_law, 0.0, dtmax=5.0)
    for tt, i in ramp:
        sel = HA["t"] <= tt + 1790
        log(f"  t = {tt:5.0f}-{tt + 1800:5.0f} s: i = {i * 1e3:5.0f} mA/cm2 (x_c placeholder {xc_of_i(i):.3f}): "
            f"end mean x = {HA['mean'][sel][-1]:.3f}, back {HA['back'][sel][-1]:.3f}, exit {HA['exit'][sel][-1]:.3f}")
    tb = HA["t"][np.argmax(HA["exit"] > b)]
    log(f"  whole membrane beta after {tb:.0f} s (single alpha->beta transit at 5 mA/cm2, eta = {eta_c})")
    for k, sy in SY.items():
        dmax, dmean, sel_ = damage(HA, 0, 9000, sy)
        log(f"  first loading, sigma_y {sy / 1e6:.0f} MPa: max elastic mismatch stress {sel_ / 1e6:.0f} MPa; plastic "
            f"strain (one transit) max {dmax:.4f}, thickness-mean {dmean:.4f}")
    # steady gradient stress at the end
    mis = HA["mis"][-1]
    log(f"  steady state (300 mA/cm2): self-equilibrated stress max {E_Pd / (1 - NU_Pd) * np.abs(mis).max() / 1e6:.1f} "
        f"MPa (free plate) -> elastic for all tempers")
    ax[0].plot(HA["t"] / 60, HA["back"], label="back face")
    ax[0].plot(HA["t"] / 60, HA["mean"], label="mean")
    ax[0].plot(HA["t"] / 60, HA["exit"], label="exit (detector) face")
    ax[0].axhspan(a, b, color="grey", alpha=0.15)
    ax[0].set_ylabel("x")
    ax[0].set_title("Protocol A: galvanostatic ramp 5→20→50→100→300 mA/cm², L = 25 µm", fontsize=9)
    ax[0].legend(fontsize=7)

    # ---------------------------------------------------------------- Protocol B1
    log("")
    log("=== Protocol B1: beta-phase flux/gradient pumping (no phase change) ===")
    x0 = HA["mean"][-1]
    res_B1 = []
    for P, xlo in ((20.0, 0.80), (60.0, 0.80), (240.0, 0.70)):
        Jc, Ja = eta_c * 0.3 * 1e4 / e, 0.1 * 1e4 / e     # 300 mA/cm2 cathodic, 100 mA/cm2 anodic
        leftB1 = lambda t, P=P, xlo=xlo: ("flux", charge(Jc, 0.95)) if (t % P) < P / 2 else ("flux", strip(Ja, xlo))
        H = run_protocol("B1", L, leftB1, 8 * P, exit_law, x0, dtmax=P / 200)
        sel = H["t"] > 6 * P
        out = [f"P={P:.0f} s, back 0.95<->{xlo}: exit x {H['exit'][sel].min():.3f}-{H['exit'][sel].max():.3f}, "
               f"mean {H['mean'][sel].min():.3f}-{H['mean'][sel].max():.3f}"]
        for k, sy in SY.items():
            dmax, dmean, sel_ = damage(H, 6 * P, 7 * P, sy)
            out.append(f"sy={sy / 1e6:.0f}: dep_max {dmax:.4f} mean {dmean:.4f}")
        log("  " + "; ".join(out))
        res_B1.append((P, xlo, H))
    P, xlo, H = res_B1[0]
    ax[1].plot(H["t"], H["back"], label="back")
    ax[1].plot(H["t"], H["mean"], label="mean")
    ax[1].plot(H["t"], H["exit"], label="exit")
    ax[1].set_ylabel("x")
    ax[1].set_title(f"Protocol B1: 300 mA/cm² (η=0.5) / −100 mA/cm² square wave, P = {P:.0f} s, back floor x={xlo}",
                    fontsize=9)
    ax[1].legend(fontsize=7)

    # ---------------------------------------------------------------- Protocol B2
    log("")
    log("=== Protocol B2: deliberate alpha/beta cycling ===")
    ic, ia = 0.1, 0.02       # A/cm2
    Jc, Ja = eta_c * ic * 1e4 / e, ia * 1e4 / e
    tc = n_Pd * L * (0.95 - 0.05) / Jc * 1.3
    ta = n_Pd * L * (0.95 - 0.05) / Ja * 1.1
    Pc = tc + ta
    log(f"  cathodic {ic * 1e3:.0f} mA/cm2 (eta {eta_c}) for {tc:.0f} s, anodic {ia * 1e3:.0f} mA/cm2 for {ta:.0f} s; "
        f"period {Pc / 60:.1f} min")
    leftB2 = lambda t: ("flux", charge(Jc, 0.95)) if (t % Pc) < tc else ("flux", strip(Ja, 0.01))
    HB = run_protocol("B2", L, leftB2, 3 * Pc, exit_law, x0, dtmax=Pc / 600)
    sel = HB["t"] > 2 * Pc
    log(f"  cycle 3: mean x {HB['mean'][sel].min():.3f}-{HB['mean'][sel].max():.3f}; exit {HB['exit'][sel].min():.3f}-"
        f"{HB['exit'][sel].max():.3f}; back {HB['back'][sel].min():.3f}-{HB['back'][sel].max():.3f}")
    DEP = {}
    for k, sy in SY.items():
        dmax, dmean, sel_ = damage(HB, 2 * Pc, 3 * Pc, sy)
        DEP[k] = (dmax, dmean)
        log(f"  sigma_y {sy / 1e6:.0f} MPa: plastic strain range per cycle max {dmax:.4f}, thickness-mean {dmean:.4f}")
    ax[2].plot(HB["t"] / 60, HB["back"], label="back")
    ax[2].plot(HB["t"] / 60, HB["mean"], label="mean")
    ax[2].plot(HB["t"] / 60, HB["exit"], label="exit")
    ax[2].axhspan(a, b, color="grey", alpha=0.15)
    ax[2].set_xlabel("t (min) [panel 2: s]")
    ax[2].set_ylabel("x")
    ax[2].set_title(f"Protocol B2: {ic * 1e3:.0f} mA/cm² cathodic / {ia * 1e3:.0f} mA/cm² anodic, period "
                    f"{Pc / 60:.0f} min", fontsize=9)
    ax[2].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "m3_cycling.png"), dpi=130)
    plt.close(fig)

    # ---------------------------------------------------------------- damage model
    log("")
    log("=== Damage model for B2 (per full alpha/beta cycle) ===")
    dep = DEP["as-supplied/H-cycled"][0]
    log(f"  plastic strain range used: {dep:.4f} (H-cycled temper; misfit {ETA * (b - a):.4f})")
    log("  Crack initiation, Coffin-Manson dep/2 = ef' (2 N_i)^c, c = -0.6 [generic fcc; Manson 1965]:")
    Ni = {}
    for ef in (0.1, 0.3, 1.0):
        Ni[ef] = 0.5 * (dep / 2 / ef) ** (1 / -0.6)
        log(f"    ef' = {ef}: N_i = {Ni[ef]:.0f} cycles")
    log("  Short-crack growth (Tomkins 1968, strain-controlled): da/dN = B dep a, B = 1-5; a0 = 1 um.")
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    Ncyc = np.arange(0, 301)
    for B, ls in ((1, "-"), (2, "--"), (5, ":")):
        Ng = np.log(L / 1e-6) / (B * dep)
        log(f"    B = {B}: cycles from initiation to through-thickness crack (L = 25 um): {Ng:.0f}; "
            f"L = 10 um: {np.log(10) / (B * dep):.0f}; L = 50 um: {np.log(50) / (B * dep):.0f}")
        for ef, c in ((0.1, "C3"), (0.3, "C1"), (1.0, "C2")):
            aN = np.where(Ncyc > Ni[ef], 1e-6 * np.exp(B * dep * (Ncyc - Ni[ef])), 0)
            aN = np.minimum(aN, L)
            if B == 2:
                ax[0].plot(Ncyc, aN * 1e6, ls, color=c, label=f"ef'={ef}, B={B}")
    ax[0].axhline(L * 1e6, color="k", lw=0.6)
    ax[0].set_xlabel("α/β cycles")
    ax[0].set_ylabel("crack depth (µm)")
    ax[0].set_title("B2: crack depth vs cycles (L = 25 µm)", fontsize=9)
    ax[0].legend(fontsize=7)
    log("  New crack-face area per unit membrane area per cycle = 4/S * da/dN (two faces, two crack families),")
    log("  S = crack spacing ~ grain size (intergranular) or ~5x depth (shear-lag saturation):")
    for S in (20e-6, 50e-6, 100e-6):
        for aa in (1e-6, 5e-6, 20e-6):
            dA = 4 / S * 2 * dep * aa
            log(f"    S = {S * 1e6:4.0f} um, depth {aa * 1e6:4.0f} um: {dA:.3g} cm2 per cm2 per cycle; cumulative "
                f"crack face {4 * aa / S:.2f} cm2/cm2")
    # vacancies
    log("  Deformation-induced vacancies: dc_v/deps = chi sigma Omega / E_f^v (Militzer, Sun & Jonas, Acta Metall.")
    log("  Mater. 42 (1994) 133 [mem]; chi = 0.1, sigma = flow stress, Omega = 1.47e-29 m3, E_f^v(Pd) = 1.5 eV [mem]).")
    Om = 1.47e-29
    lo_rate = 1e-5
    for k, sy in SY.items():
        hi_rate = 0.1 * sy * Om / (1.5 * e)
        log(f"    {k:22s}: {hi_rate:.2g} per unit strain (lower bound used: {lo_rate:.0e}); per B2 cycle "
            f"{lo_rate * dep:.1e} - {hi_rate * dep:.1e}; after 30 cycles {30 * lo_rate * dep:.1e} - "
            f"{30 * hi_rate * dep:.1e}")
    # vacancy mobility at RT
    for T in (298.15, 363.15):
        for Em in (1.0, 1.2):
            jumps = 1e13 * np.exp(-Em / (kB_eV * T)) * 86400
            ldiff = 2.75e-10 * np.sqrt(jumps)
            log(f"    vacancy (E_m = {Em} eV) at {T - 273.15:.0f} C: {jumps:.3g} jumps/day, diffusion length/day "
                f"{ldiff * 1e9:.2g} nm")
    # SAV thermodynamics (illustrative)
    for eb, r in ((0.23, 6),):
        Ef = 1.5 - r * eb
        log(f"    SAV thermodynamics: E_f(vac + {r} D) = 1.5 - {r} x {eb} = {Ef:.2f} eV -> equilibrium c_v(25 C) "
            f"~ exp(-Ef/kT) = {np.exp(-Ef / (kB_eV * T0)):.2g}, but vacancies cannot arrive (kinetics above)")
    # cumulative area plot
    for S, c in ((20e-6, "C0"), (50e-6, "C1")):
        for ef in (0.3,):
            aN = np.where(Ncyc > Ni[ef], np.minimum(1e-6 * np.exp(2 * dep * (Ncyc - Ni[ef])), L), 0)
            ax[1].plot(Ncyc, 4 * aN / S, color=c, label=f"S = {S * 1e6:.0f} µm, ef'={ef}, B=2")
    ax[1].set_xlabel("α/β cycles")
    ax[1].set_ylabel("cumulative crack-face area (cm² per cm² membrane)")
    ax[1].set_title("B2: fresh crack-face area", fontsize=9)
    ax[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "m3_cycling_damage.png"), dpi=130)
    plt.close(fig)
    open(os.path.join(OUT, "m3_cycling.txt"), "w").write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
