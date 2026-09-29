#!/usr/bin/env python3
"""Red-team 1C: statistics, inference and credibility checks on docs/design/iteration-1.md.

Sections (all written to docs/design/figs/rt1c_stats.txt):
  1. Trials factor implied by the design and the global significance of a local 5 sigma.
  2. T1 (Si PID window) reach in the real DFM geometry: per skin, per control topology,
     with live-time losses and the +-x3 background uncertainty; local-6 sigma thresholds.
  3. Fragility of a low-count "5 sigma" to the background model (true p if B is x3 / bursts).
  4. D-specific cosmic-neutron background (d(n,np) breakup, n-d recoil deuterons) that the
     H2O twin does not share. Order-of-magnitude only.
  5. "Absent in the matched control": power of the T-H null given its exposure.
  6. Sample count: binomial and lot-correlated (beta-binomial) P(>=k active); allocation of
     16 quadrant slots over skins; claim rule "replication in >=2 membranes".
  7. Bayesian decision analysis: posterior on "effect real" after all-null / single-membrane
     positive / full claim-rule positive.
  8. Sec. 10 ratio audit (reach / claimed magnitude).
Figures: docs/design/figs/rt1c_reach.png, docs/design/figs/rt1c_posterior.png
Fixed inputs only; no random numbers. Run: python3 sim/rt1c_stats.py
"""
import os
import sys

import numpy as np
from scipy import optimize, special, stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from m5_stats import n_crit, z_onoff  # noqa: E402

FIGS = os.path.join(HERE, "..", "docs", "design", "figs")
os.makedirs(FIGS, exist_ok=True)
DAY = 86400.0
P5 = stats.norm.sf(5.0)

_out = open(os.path.join(FIGS, "rt1c_stats.txt"), "w")


def P(*a):
    print(*a)
    print(*a, file=_out)


def sig(p):
    return stats.norm.isf(p)


# =============================================================================== 1. trials
def trials():
    P("== 1. Trials factor implied by iteration-1 (enumeration of analyses the text invites)")
    rows = [
        # name, count, comment
        ("T1 per skin (5) x phase (P2,P3a,P3b,P4,P5) x window (peak, PID) x contrast (vs T-H, vs off, lock-in)",
         5 * 5 * 2 * 3),
        ("T1 per membrane-quadrant (16) x phase (5) x window (2)  [claim rule needs per-membrane results]", 16 * 5 * 2),
        ("secondary lines: t 1.01, 3He 0.82 per quadrant (16) x phase (5)", 2 * 16 * 5),
        ("T2 511-511: rail position (6 cells) x phase (5) x window (511-511, 3-25 MeV)", 6 * 5 * 2),
        ("T3 He: channel (front, headspace) x A-cell (4) x static window (~10 per run)", 2 * 4 * 10),
        ("T3 melts: 16 quadrants + 4 headspace/entry-side", 20),
        ("T4 C3-G: candidate product masses (~10) x 2 labs", 20),
        ("secondary: neutron, heat, AE-correlation x 5 phases", 15),
        ("flux lock-in: period grid (~5) x lag grid (~4) x 16 quadrants", 5 * 4 * 16),
    ]
    naive = 0
    for name, n in rows:
        naive += n
        P(f"   {n:5d}  {name}")
    P(f"   naive total N = {naive}")
    P("   Correlations (shared data between windows/contrasts/phases) reduce this; an effective")
    P("   N_eff of 1/10 to 1/3 of the naive count is typical for overlapping counting analyses.")
    P("   N_eff   global sigma of a local 5.0 sigma   local sigma needed for global 5 sigma (Sidak)")
    for N in [1, 4, 20, 50, 100, 300, 1000]:
        pg = 1 - (1 - P5) ** N
        pl = 1 - (1 - P5) ** (1.0 / N)
        P(f"   {N:5d}   {sig(pg):5.2f}                               {sig(pl):5.2f}")
    return naive


# =============================================================================== 2. T1 reach in DFM
def reach_counts(b, tau=None, sigma=5.0, power=0.5):
    """Signal counts for `sigma` discovery at `power`; known background (tau=None) uses the exact
    Poisson n_crit; on/off with t_off = tau t_on uses Asimov Li-Ma."""
    if tau is None:
        n = n_crit(max(b, 1e-9), sigma)
        f = lambda s: stats.poisson.sf(n - 1, b + s) - power
        return optimize.brentq(f, 1e-9, 1e5)
    f = lambda s: z_onoff(s + b, tau * b, tau) - sigma
    return optimize.brentq(f, 1e-6, 1e7)


def t1_reach():
    P("\n== 2. T1 reach in the DFM geometry (median 5 sigma, fusions/s per quadrant of a given skin)")
    P("   Inputs: branch 0.5 to p; eps per emitted p per quadrant = 0.15 (iteration-1 risk 4) or 0.08")
    P("   (entry-face source, 2.05 MeV p, lower PID acceptance, from m5_si.txt U(0-20um) scaling);")
    P("   background per quadrant = cell B/4 with cell B = 0.01-0.03/day (M5), x3 uncertainty;")
    P("   live fraction 0.8 (EMI edge vetoes, keep-alive gaps, DAQ); P2 = 21 d nominal, 30 d = M5's number.")
    res = {}
    cases = []
    for eps in [0.15, 0.08]:
        for Bcell in [0.01, 0.03, 0.09]:
            for n_q in [1, 3, 4]:
                for days in [21, 30]:
                    cases.append((eps, Bcell, n_q, days))
    P("   eps   Bcell/d  n_q  days | known-B(10x)  | ctrl=T-H(1 quad/skin) | local 6.0 sigma known-B")
    for eps, Bcell, n_q, days in cases:
        T = days * 0.8 * DAY
        e_f = 0.5 * eps  # counts per fusion
        b = n_q * Bcell / 4 * days * 0.8  # background counts in A quadrants of that skin
        s_known = reach_counts(b, tau=10.0)
        # T-H has one quadrant per skin: t_off/t_on = 1/n_q
        s_th = reach_counts(b, tau=1.0 / n_q)
        s_6 = reach_counts(b, None, sigma=6.0)
        r = lambda s: s / (n_q * e_f * T)
        res[(eps, Bcell, n_q, days)] = (r(s_known), r(s_th), r(s_6))
        if days == 21 or (eps == 0.15 and Bcell == 0.01):
            P(f"   {eps:4.2f}  {Bcell:6.2f}  {n_q:3d}  {days:4d} | {r(s_known):8.2e} ({s_known:4.1f}c)"
              f" | {r(s_th):8.2e} ({s_th:5.1f}c)  | {r(s_6):8.2e} ({s_6:4.1f}c)")
    P("   M5 headline for comparison: 3.3e-5 (10x bkg) / 7.6e-5 (equal-time) fusions/s for the WHOLE")
    P("   Oe20 membrane in vacuum, eps=0.22, 30 d, 100 % live, surface source.")
    return res


# =============================================================================== 3. fragility
def fragility():
    P("\n== 3. Fragility of low-count 5 sigma: true significance if the background is mis-modelled")
    P("   b_assumed  n_crit(5s)  true sigma if b_true = x2 / x3 / x10 of assumed")
    for b in [0.05, 0.1, 0.3, 1.0, 3.0]:
        n = n_crit(b, 5.0)
        out = []
        for k in [2, 3, 10]:
            out.append(sig(stats.poisson.sf(n - 1, k * b)))
        P(f"   {b:8.2f}   {n:5d}       " + "  ".join(f"{x:5.2f}" for x in out))
    P("   Profiled over a +-x3 (log-uniform) background prior, the 5 sigma n_crit rises to:")
    for b in [0.1, 0.3, 1.0]:
        # worst case in the band = 3b
        P(f"     b_nom={b:4.2f}: n_crit(b)={n_crit(b)}, n_crit(3b)={n_crit(3 * b)}")
    P("   Burst artifacts: one unvetoed microdischarge/pile-up cluster that deposits k>=n_crit PID-")
    P("   window events is a '5 sigma' on its own. Required: burst test separate from rate test;")
    P("   rate test counts clusters (events within 10 s in one detector) as ONE event.")


# =============================================================================== 4. D-specific cosmic
def d_specific():
    P("\n== 4. D-specific cosmic-neutron backgrounds absent from the H2O twin (order of magnitude)")
    NA = 6.022e23
    nD_pd = 6.8e22 * 0.9  # D/cm3 in PdD0.9
    V_mem = np.pi * 1.0 ** 2 * 15e-4  # cm3
    ND_mem = nD_pd * V_mem
    nD_d2o = 1.107 / 20.03 * NA * 2  # D/cm3
    ND_film = nD_d2o * np.pi * 1.0 ** 2 * 100e-4  # 100 um electrolyte film reachable by >2 MeV d
    ND_gas = 0.5 * 2.69e19 * 2 * np.pi * 1.0 ** 2 * 0.4  # 4 mm of 0.5 bar D2
    ND = ND_mem + ND_film + ND_gas
    P(f"   D atoms within product range of the telescope: membrane {ND_mem:.1e}, electrolyte film"
      f" (100 um) {ND_film:.1e}, front gas {ND_gas:.1e}; total {ND:.1e}")
    # breakup: sigma ~0.1-0.2 b above ~5 MeV; flux >5 MeV inside house 1-4e-3 /cm2/s
    for lab, sig_b, phi, acc in [("low", 0.08e-24, 1e-3, 0.03), ("high", 0.2e-24, 4e-3, 0.08)]:
        r = ND * sig_b * phi * DAY
        P(f"   d(n,np) breakup [{lab}]: {r:.3f} breakups/day/cell; x acceptance*PID-window fraction"
          f" {acc} -> {r * acc:.4f} protons/day/cell in the T1 window")
    for lab, sig_e, phi, acc, leak in [("low", 1.0e-24, 5e-3, 0.05, 1e-3), ("high", 2.5e-24, 1.5e-2, 0.15, 1e-2)]:
        r = ND * sig_e * phi * DAY
        P(f"   n-d elastic recoil d [{lab}]: {r:.2f}/day/cell; reaching telescope with E>1.4 MeV x{acc}"
          f" = {r * acc:.3f}/day; p/d PID leakage {leak} -> {r * acc * leak:.5f}/day in T1 window")
    P("   H2O twin: n-p recoil protons 0.006/day (M5) -> the D-H difference is NOT signed-conservative:")
    P("   D-specific terms of 1e-3 to 3e-2/day are comparable to the whole modelled background.")
    P("   They are flux-independent: a within-cell D flux-on vs D flux-off (static-loaded) contrast")
    P("   cancels them; the D-vs-H contrast does not.")


# =============================================================================== 5. control absence
def control_absence():
    P("\n== 5. 'Absence in the matched control': how much does a T-H null say?")
    P("   A-cells observe s signal counts summed over n_q quadrants of one skin; if the effect is")
    P("   NOT D-specific (artifact common to D and H cells), T-H (1 quadrant) expects s/n_q.")
    P("   s   n_q  expected T-H  P(T-H = 0 | common artifact)  likelihood ratio for 'D-specific'")
    for s in [5, 10, 20]:
        for n_q in [3, 4]:
            mu = s / n_q
            p0 = np.exp(-mu)
            P(f"   {s:3d}  {n_q:3d}   {mu:6.2f}        {p0:6.3f}                        {1 / p0:7.1f}")
    P("   With all 4 T-H quadrants pooled against all 16 A quadrants (tau=1/4), s=10 gives LR ~12.")


# =============================================================================== 6. sample counts
def p_atleast(k, n, p):
    return stats.binom.sf(k - 1, n, p)


def bb_atleast(k, n, mean, rho):
    """Beta-binomial P(X>=k); rho = intra-class correlation (lot effect)."""
    if rho <= 0:
        return p_atleast(k, n, mean)
    if rho >= 1:
        return mean if k >= 1 else 1.0
    ab = 1 / rho - 1
    a, b = mean * ab, (1 - mean) * ab
    ks = np.arange(k, n + 1)
    pm = np.exp(special.gammaln(n + 1) - special.gammaln(ks + 1) - special.gammaln(n - ks + 1)
                + special.betaln(ks + a, n - ks + b) - special.betaln(a, b))
    return pm.sum()


def samples():
    P("\n== 6. Sample count. P(>=1) is what iteration-1 quotes; the claim rule needs >=2 membranes")
    P("   p     N   P(>=1)  P(>=2)   | lot-correlated rho=0.3: P(>=1) P(>=2) | rho=0.7: P(>=1) P(>=2)")
    for p in [0.05, 0.1, 0.2]:
        for N in [4, 6, 8, 12]:
            P(f"   {p:4.2f} {N:3d}   {p_atleast(1, N, p):5.3f}  {p_atleast(2, N, p):5.3f}   |"
              f"            {bb_atleast(1, N, p, 0.3):5.3f}  {bb_atleast(2, N, p, 0.3):5.3f}  |"
              f"         {bb_atleast(1, N, p, 0.7):5.3f}  {bb_atleast(2, N, p, 0.7):5.3f}")
    P("   Two lots x N/2 membranes (lot-level effect, within-lot p=0.2, lot 'good' w.p. 0.5 -> mean 0.1):")
    for N in [4, 8]:
        # each lot good w.p. 0.5; if good, per-membrane p=0.2
        pg, pin = 0.5, 0.2
        # P(>=2 active total) with 2 lots each of N/2
        m = N // 2
        dist_lot = [(1 - pg) + pg * stats.binom.pmf(0, m, pin)] + [pg * stats.binom.pmf(j, m, pin) for j in range(1, m + 1)]
        tot = np.convolve(dist_lot, dist_lot)
        one_lot = [(1 - pg) + pg * stats.binom.pmf(0, N, pin)] + [pg * stats.binom.pmf(j, N, pin) for j in range(1, N + 1)]
        P(f"     N={N}: two lots P(>=1)={1 - tot[0]:.3f} P(>=2)={1 - tot[0] - tot[1]:.3f} | "
          f"one lot P(>=1)={1 - one_lot[0]:.3f} P(>=2)={1 - one_lot[0] - one_lot[1]:.3f}")
    P("\n   Allocation of 16 quadrant slots (quadrant-level activity, p=0.2 per slot, independent):")
    allocs = {"5 skins (4,3,3,3,3) [iter-1]": [4, 3, 3, 3, 3], "4 skins (4,4,4,4)": [4, 4, 4, 4],
              "3 skins (6,5,5)": [6, 5, 5], "2 skins (8,8)": [8, 8]}
    for name, a in allocs.items():
        e1 = sum(p_atleast(1, n, 0.2) for n in a)
        e2 = sum(p_atleast(2, n, 0.2) for n in a)
        P(f"     {name:32s} sum_s P(>=1)={e1:4.2f}  sum_s P(>=2)={e2:4.2f}")
    P("   Membrane-level activity (entry face / bulk): skins do not matter; only N membranes does.")


# =============================================================================== 7. Bayes
def posterior(prior, lr):
    o = prior / (1 - prior) * lr
    return o / (1 + o)


def bayes():
    P("\n== 7. Bayesian decision analysis (H = effect real at >=1 % of a claimed magnitude in this regime)")
    P("   P(pos|H) = P(cond) x P(>=k active) x P(det); P(pos|not H) = artifact/fluke rate per campaign.")
    P("   P(cond)=0.5, P(det)=0.9, p=0.2, N=4 -> single-membrane P(pos|H)=0.5*0.59*0.9=0.27;")
    P("   claim-rule (>=2 membranes) P(pos|H)=0.5*0.18*0.9=0.08.")
    q1 = 0.5 * p_atleast(1, 4, 0.2) * 0.9
    q2 = 0.5 * p_atleast(2, 4, 0.2) * 0.9
    P("   prior   post|all-null   post|1-membrane +ve (art=0.05/0.01)   post|claim-rule +ve (art=1e-3/1e-4)")
    rows = []
    for pr in [1e-4, 1e-3, 1e-2, 3e-2, 0.1]:
        pn = posterior(pr, (1 - q1) / 1.0)
        p1a = posterior(pr, q1 / 0.05)
        p1b = posterior(pr, q1 / 0.01)
        p2a = posterior(pr, q2 / 1e-3)
        p2b = posterior(pr, q2 / 1e-4)
        rows.append((pr, pn, p1a, p1b, p2a, p2b))
        P(f"   {pr:6.0e}   {pn:8.2e}        {p1a:8.2e} / {p1b:8.2e}                   {p2a:8.2e} / {p2b:8.2e}")
    P(f"   Bayes factor of an all-null against H: {1 / (1 - q1):.2f} (N=4); N=8: {1 / (1 - 0.5 * p_atleast(1, 8, 0.2) * 0.9):.2f};"
      f" with P(cond)=0.9, N=8: {1 / (1 - 0.9 * p_atleast(1, 8, 0.2) * 0.9):.2f}")
    return rows, q1, q2


def fig_posterior(rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    r = np.array(rows)
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.loglog(r[:, 0], r[:, 1], "k-o", label="all null")
    ax.loglog(r[:, 0], r[:, 2], "C1--s", label="1-membrane positive, artifact rate 0.05")
    ax.loglog(r[:, 0], r[:, 3], "C1-s", label="1-membrane positive, artifact rate 0.01")
    ax.loglog(r[:, 0], r[:, 4], "C0--^", label="claim-rule positive, artifact rate 1e-3")
    ax.loglog(r[:, 0], r[:, 5], "C0-^", label="claim-rule positive, artifact rate 1e-4")
    ax.loglog(r[:, 0], r[:, 0], ":", color="0.5", label="posterior = prior")
    ax.set_xlabel("prior P(effect real in this regime)")
    ax.set_ylabel("posterior")
    ax.set_title("RT-1C: what the DFM-4 outcome is worth (N=4, p=0.2)", fontsize=9)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "rt1c_posterior.png"), dpi=130)
    plt.close(fig)


def fig_reach(res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    labels, vals = [], []
    labels.append("M5 headline\n(whole foil, vacuum,\n30 d, 10x bkg)")
    vals.append((3.3e-5, 3.3e-5))
    for n_q, lab in [(4, "skin L (4 quads)"), (3, "other skins (3 quads)"), (1, "one quadrant\n(per-membrane)")]:
        lo = res[(0.15, 0.01, n_q, 21)][0] * n_q
        hi = res[(0.08, 0.09, n_q, 21)][1] * n_q
        labels.append(lab + "\n21 d, total over quads")
        vals.append((lo, hi))
    x = np.arange(len(labels))
    for i, (lo, hi) in enumerate(vals):
        if hi > lo:
            ax.plot([i, i], [lo, hi], lw=8, color="C0", solid_capstyle="butt")
        else:
            ax.plot([i], [lo], "s", ms=9, color="0.4")
    ax.axhline(1e-4, color="C3", ls=":")
    ax.text(len(labels) - 0.6, 1.1e-4, "1 % of Lipson low end (1e-2/s)", fontsize=7, color="C3", ha="right")
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel("5σ median reach, fusions/s (sum over quadrants)")
    ax.set_title("T1 reach: M5 headline vs DFM-4 as specified\n(bar: best = eps 0.15, B 0.01/d, 10x bkg; "
                 "worst = eps 0.08, B 0.09/d, T-H as control)", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "rt1c_reach.png"), dpi=130)
    plt.close(fig)


# =============================================================================== 8. sec 10 audit
def audit():
    P("\n== 8. Sec. 10 audit: reach / claimed magnitude (design says <= 1e-2 everywhere)")
    rows = [
        ("Si vs Lipson (1e-2..1 /s charged particles)", "3e-5..8e-5 (M5, whole foil)", 3e-5 / 1, 8e-5 / 1e-2),
        ("Si vs Lipson, DFM per skin, sum over quads (sec. 2)", "6e-5..1e-3 per skin", 6e-5 / 1, 1e-3 / 1e-2),
        ("Si vs Tohoku-static 130/day", "2.2 events/day", 2.2 / 130, 2.2 / 130),
        ("Heat vs SRI 0.1..1 W (as quoted)", "4..20 mW", 4e-3 / 1, 20e-3 / 0.1),
        ("Heat vs SRI scaled to DFM volume (4.7e-3 vs 2.4e-2 cm3)", "4..20 mW", 4e-3 / (1 * 0.2), 20e-3 / (0.1 * 0.2)),
        ("He front vs SRI-scale (volume-scaled)", "0.16..0.39 nW (x10 if stock-He 10x)", 0.16e-9 / 0.2, 3.9e-9 / 0.02),
        ("C3-G products vs 1e12..1e14 cm-2", "1e10 cm-2 (not demonstrated)", 1e10 / 1e14, 1e10 / 1e12),
    ]
    for name, reach, lo, hi in rows:
        flag = "OK" if hi <= 1e-2 else ("BORDERLINE" if lo <= 1e-2 else "FAILS")
        P(f"   {name:58s} {reach:34s} ratio {lo:7.1e}..{hi:7.1e}  {flag}")


def main():
    trials()
    res = t1_reach()
    fragility()
    d_specific()
    control_absence()
    samples()
    rows, q1, q2 = bayes()
    audit()
    fig_reach(res)
    fig_posterior(rows)
    _out.close()


if __name__ == "__main__":
    main()
