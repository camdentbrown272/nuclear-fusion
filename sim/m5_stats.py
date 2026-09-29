#!/usr/bin/env python3
"""M5 statistics: Currie limits, Poisson/on-off discovery thresholds, Feldman-Cousins,
Bayesian limits, current-modulation (on/off) power, look-elsewhere correction, run time to 5 sigma.

Library functions are imported by the other m5_* scripts. Running this file regenerates
figs/m5_stats.txt and figs/m5_stats_runtime.png using channel parameters from the other models.
"""
import os
import sys

import numpy as np
from scipy import optimize, stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

P5 = stats.norm.sf(5.0)  # one-sided 5 sigma = 2.87e-7


# ------------------------------------------------------------------------------ Currie
def currie(B, paired=True, alpha=0.05, beta=0.05):
    """Currie (1968) critical level and detection limit in counts. paired=True: blank measured
    for the same time (variance 2B); False: blank known exactly."""
    k = stats.norm.isf(alpha)
    f = np.sqrt(2.0) if paired else 1.0
    Lc = k * f * np.sqrt(B)
    Ld = k ** 2 + 2 * Lc  # = 2.71 + 4.65 sqrt(B) for paired, k=1.645
    return Lc, Ld


# ------------------------------------------------------------------------------ Poisson discovery
def n_crit(b, sigma=5.0):
    """Smallest n with P(N >= n | b) <= p(sigma), known background b."""
    p = stats.norm.sf(sigma)
    n = 1
    while stats.poisson.sf(n - 1, b) > p:
        n += 1
    return n


def z_onoff(n_on, n_off, tau):
    """Li & Ma (1983) eq. 17 significance; tau = t_off/t_on (alpha = 1/tau)."""
    a = 1.0 / tau
    n_on = np.asarray(n_on, float)
    n_off = np.asarray(n_off, float)
    tot = n_on + n_off
    with np.errstate(divide="ignore", invalid="ignore"):
        t1 = np.where(n_on > 0, n_on * np.log((1 + a) / a * n_on / tot), 0)
        t2 = np.where(n_off > 0, n_off * np.log((1 + a) * n_off / tot), 0)
    z = np.sqrt(np.maximum(2 * (t1 + t2), 0))
    return np.where(n_on > a * n_off, z, -z)


def discovery_signal(b, alpha_sigma=5.0, power=0.5, tau=None, b_unc_rel=None):
    """Signal counts s needed for a 5 sigma discovery with probability `power`.
    tau=None: background known (Poisson counting, exact); tau given: on/off with t_off = tau t_on
    (Asimov median of Li-Ma). b_unc_rel is accepted for API compatibility (ignored when tau given)."""
    if tau is None:
        b = max(b, 1e-6)
        n = n_crit(b, alpha_sigma)
        # find s: P(N >= n | b+s) = power
        f = lambda s: stats.poisson.sf(n - 1, b + s) - power
        return optimize.brentq(f, 1e-9, 1e4)
    f = lambda s: z_onoff(s + b, tau * b, tau) - alpha_sigma
    return optimize.brentq(f, 1e-6, 1e7)


# ------------------------------------------------------------------------------ Feldman-Cousins
def fc_interval(n_obs, b, cl=0.9, mu_max=50.0, dmu=0.005):
    """Feldman & Cousins (1998) confidence belt for Poisson signal with known background."""
    mus = np.arange(0, mu_max, dmu)
    nmax = int(mu_max + b + 10 * np.sqrt(mu_max + b) + 20)
    n = np.arange(nmax)
    mubest = np.maximum(n - b, 0)
    lo, hi = None, None
    acc = []
    for mu in mus:
        p = stats.poisson.pmf(n, mu + b)
        r = p / stats.poisson.pmf(n, mubest + b)
        order = np.argsort(-r)
        c = np.cumsum(p[order])
        k = np.searchsorted(c, cl) + 1
        inset = np.zeros(nmax, bool)
        inset[order[:k]] = True
        acc.append(inset[n_obs])
    acc = np.array(acc)
    if acc.any():
        lo, hi = mus[acc][0], mus[acc][-1]
    return lo, hi


def bayes_ul(n_obs, b, cl=0.9):
    """Flat-prior (s>=0) Bayesian upper limit, known background."""
    # closed form: posterior CDF = 1 - P(N<=n | b+s)/P(N<=n | b)
    g = lambda s: stats.poisson.cdf(n_obs, b + s) / stats.poisson.cdf(n_obs, b) - (1 - cl)
    return optimize.brentq(g, 0, 1e4)


# ------------------------------------------------------------------------------ modulation
def modulation_contrast(T_mod, tau):
    """Mean(on-half) - mean(off-half) of a first-order-lag response (time constant tau) to a
    50 %-duty square wave of period T_mod, in units of the full on-state signal."""
    h = T_mod / (2 * tau)
    return 1 - 2 * np.tanh(h / 2) / h if h > 1e-9 else 0.0


def modulation_time_to_5sigma(r_sig, r_bkg, T_mod, tau, sigma=5.0):
    """Total run time (s) for the on-vs-off difference (equal time) to reach `sigma` (Asimov/Li-Ma),
    for a signal present only when the drive is on (with lag tau)."""
    c = modulation_contrast(T_mod, tau)
    # mean on-state signal fraction and off-state fraction
    h = T_mod / (2 * tau)
    m_off = (1 - c) / 2 if h > 0 else 0.5
    m_on = m_off + c
    f = lambda T: z_onoff((r_bkg + m_on * r_sig) * T / 2, (r_bkg + m_off * r_sig) * T / 2, 1.0) - sigma
    lo = 1e-3
    if f(1e12) < 0:
        return np.inf
    return optimize.brentq(f, lo, 1e12)


def time_to_5sigma(r_sig, r_bkg, tau_bkg=1.0, sigma=5.0, known_bkg=False):
    """Run time (s) for signal+background counting vs a background-only control observed for
    tau_bkg times as long (equal-time H-control: tau_bkg=1). Asimov median."""
    if known_bkg:
        # Poisson known-b: Asimov continuous formula
        g = lambda T: np.sqrt(2 * ((r_sig + r_bkg) * T * np.log(1 + r_sig / max(r_bkg, 1e-30)) - r_sig * T)) - sigma
        return optimize.brentq(g, 1e-6, 1e14)
    g = lambda T: z_onoff((r_sig + r_bkg) * T, tau_bkg * r_bkg * T, tau_bkg) - sigma
    return optimize.brentq(g, 1e-6, 1e14)


# ------------------------------------------------------------------------------ look-elsewhere
def local_sigma_for_global(n_trials, global_sigma=5.0):
    pg = stats.norm.sf(global_sigma)
    pl = 1 - (1 - pg) ** (1.0 / n_trials)
    return stats.norm.isf(pl)


# ------------------------------------------------------------------------------ main
def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import m5_stopping as st

    FIGS = st.FIGS
    out = open(os.path.join(FIGS, "m5_stats.txt"), "w")

    def P(*a):
        print(*a)
        print(*a, file=out)

    P("== Currie limits (counts), alpha=beta=5 %")
    for B in [0, 0.1, 1, 10, 100, 1000]:
        lc, ld = currie(B, paired=True)
        lc2, ld2 = currie(B, paired=False)
        P(f"   B={B:7g}: paired blank L_C={lc:7.2f} L_D={ld:7.2f} | known blank L_C={lc2:7.2f} L_D={ld2:7.2f}")
    P("   (Currie's Gaussian L_D fails for B<~10: use the Poisson rows below)")

    P("\n== 5 sigma discovery (50 % / 90 % power): signal counts needed")
    for B in [0, 0.01, 0.1, 0.3, 1, 3, 10, 100, 1000]:
        s50 = discovery_signal(B)
        s90 = discovery_signal(B, power=0.9)
        s50oo = discovery_signal(max(B, 1e-3), tau=1.0) if B > 0 else np.nan
        P(f"   B={B:6g}: n_crit={n_crit(max(B, 1e-6)):3d}  s(50%)={s50:7.2f}  s(90%)={s90:7.2f}"
          f"  | on/off equal-time (Asimov) s={s50oo:7.2f}")

    P("\n== Feldman-Cousins 90 % CL intervals (signal counts) vs Bayesian flat-prior 90 % UL")
    for b in [0.0, 0.5, 3.0]:
        for n in [0, 1, 3, 5]:
            lo, hi = fc_interval(n, b)
            P(f"   b={b:3.1f} n={n}: FC [{lo:.2f}, {hi:.2f}]   Bayes UL {bayes_ul(n, b):.2f}")

    P("\n== Look-elsewhere: local significance required for global 5 sigma (Sidak)")
    for N in [1, 3, 10, 30, 100, 1000]:
        P(f"   N={N:5d} independent windows/channels -> local {local_sigma_for_global(N):.2f} sigma")
    P("   Pre-registered primary: N=1 (C3 proton PID window, D vs H, current on vs off). Secondary channels "
      "(t, 3He, n, alpha 5-20 MeV bins, 14.7 MeV p) count as N~10 -> 5.5 sigma local.")

    P("\n== Modulation contrast c = 1 - (4 tau/T) tanh(T/(4 tau)) for a first-order lag tau")
    for ratio in [0.5, 1, 2, 5, 10, 20, 50]:
        P(f"   T_mod/tau = {ratio:5.1f}: c = {modulation_contrast(ratio, 1.0):.3f}")
    # numeric check of the closed form
    tau = 1.0
    T = 6.0
    dt = 1e-3
    t = np.arange(0, 40 * T, dt)
    u = ((t % T) < T / 2).astype(float)
    y = np.zeros_like(t)
    for i in range(1, len(t)):
        y[i] = y[i - 1] + dt / tau * (u[i - 1] - y[i - 1])
    last = t > 30 * T
    c_num = y[last & (u == 1)].mean() - y[last & (u == 0)].mean()
    P(f"   numeric check T/tau=6: closed form {modulation_contrast(6, 1):.4f} vs ODE {c_num:.4f}")

    # --------------------------------------------------------------- run time to 5 sigma
    import m5_si_telescope as si
    import m5_neutron as nn
    import m5_cr39 as cr

    DAY = 86400.0
    si_res = si.main()
    n_res = nn.channel_summary()
    cr_res = cr.channel_summary()
    chans = {
        # name: (detected counts per fusion, background counts/s)
        "C3 Si telescope, PID window (U 0-5 um)": (0.5 * si_res["eff"]["U(0-5um)"][1], si_res["B_design"] / DAY),
        "C3 Si, single detector 2.6-3.1 (no PID)": (0.5 * si_res["eff"]["U(0-5um)"][0], si_res["B_single"] / DAY),
        "3He bank (C1/C3), 2.45 MeV n": (0.5 * n_res["eff"], n_res["bkg"]),
        "EJ-309 pair, PSD, 2.45 MeV n": (0.5 * n_res["ej_eff"], n_res["ej_bkg"]),
        "CR-39 behind 6 um Mylar (C1/C2)": (0.5 * cr_res["eff"], cr_res["bkg"]),
    }
    rates = np.logspace(-4, 0, 9)
    P("\n== Run time to 5 sigma (median, equal-duration H-control as background reference), per channel")
    P("   rate(fus/s) " + " ".join(f"{r:9.0e}" for r in rates))
    fig, ax = plt.subplots(figsize=(7.5, 5))
    rt = {}
    for name, (e, b) in chans.items():
        Ts = np.array([time_to_5sigma(r * e, b, tau_bkg=1.0) for r in rates])
        rt[name] = Ts
        P(f"   {name:42s} eff/fusion={e:.2e} bkg={b:.2e}/s")
        P("      days: " + " ".join(f"{T / DAY:9.2g}" for T in Ts))
        ax.loglog(rates, Ts / DAY, "o-", label=name)
    # modulation example on C3 Si
    e, b = chans["C3 Si telescope, PID window (U 0-5 um)"]
    Tm = np.array([modulation_time_to_5sigma(r * e, b, T_mod=3600, tau=300) for r in rates])
    ax.loglog(rates, Tm / DAY, "k--", label="C3 Si, current on/off (T=1 h, τ=5 min)")
    P("   C3 Si, on/off modulation T_mod=1 h, tau=5 min (signal only while drive on):")
    P("      days: " + " ".join(f"{T / DAY:9.2g}" for T in Tm))
    ax.axhline(30, color="0.5", ls=":")
    ax.text(1.2e-4, 33, "30 days (ADR-001)", fontsize=8)
    ax.set_xlabel("D–D fusion rate in the active layer (s⁻¹)")
    ax.set_ylabel("run time to 5σ (days, median)")
    ax.set_ylim(1e-3, 1e5)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "m5_stats_runtime.png"), dpi=130)
    plt.close(fig)
    # rate reachable in 30 days
    P("\n== Minimum rate reaching 5 sigma (median) in 1 / 14 / 30 days, equal-time control:")
    for name, (e, b) in chans.items():
        rr = []
        for days in [1, 14, 30]:
            T = days * DAY
            g = lambda lr: z_onoff((10 ** lr * e + b) * T, b * T, 1.0) - 5.0
            rr.append(10 ** optimize.brentq(g, -12, 14))
        P(f"   {name:42s} " + " ".join(f"{x:9.2e}" for x in rr))
    out.close()


if __name__ == "__main__":
    main()
