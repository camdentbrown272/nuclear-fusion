#!/usr/bin/env python3
"""M5 CR-39 (PADC) stack design: filter thicknesses that separate p, t, 3He and alpha; efficiency and
background for in-situ (C1/C2) and vacuum (C3) use.

Registration model (assumption, documented in the write-up): a track is registered and etchable to an
optically resolvable pit when the ion arrives at the CR-39 surface with LET above ~LET_min = 12 keV/um
(protons <~ 4 MeV at a standard 6.25 N NaOH, 70 C, 6 h etch; longer etches extend to ~20 MeV, see
https://www.osti.gov/servlets/purl/1076447 and the Nature Sci. Rep. 7, 2152 (2017) calibration) and within
the critical-angle cone (theta <= 50 deg from normal for 0.5-3 MeV protons; ~70 deg for alphas).

    python3 sim/m5_cr39.py    # -> figs/m5_cr39.txt, figs/m5_cr39_filters.png
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m5_stopping as st  # noqa: E402
import m5_escape as es  # noqa: E402

FIGS = st.FIGS
RNG = np.random.default_rng(39)
LET_MIN = 0.012  # MeV/um in CR-39
THETA_C = {"p": 50.0, "d": 55.0, "t": 60.0, "h": 70.0, "a": 70.0}

SPECIES = [("p", 3.02, "p 3.02 (DD)"), ("t", 1.01, "t 1.01 (DD)"), ("h", 0.82, "3He 0.82 (DD)"),
           ("a", 5.49, "222Rn alpha 5.49"), ("a", 6.00, "218Po alpha 6.00"), ("a", 7.69, "214Po alpha 7.69"),
           ("a", 8.78, "212Po alpha 8.78"), ("a", 12.0, "alpha 12 (Lipson)"), ("a", 16.0, "alpha 16"),
           ("p", 14.7, "p 14.7 (D-3He)"), ("p", 1.7, "p 1.7 (Lipson)")]

# Recommended 4-channel stack (all behind a 6 um Mylar chemical barrier in C1/C2; bare in vacuum C3)
STACK = [
    ("A", []),                                   # all species (barrier only)
    ("B", [("Al", 12.0)]),                       # stops t(1.01), 3He, alpha < 3.5 MeV; passes p, most Rn alphas
    ("C", [("Al", 55.0)]),                       # stops every natural alpha (<= 8.78 MeV, with the 6 um Mylar); passes 3 MeV p at ~1.35 MeV
    ("D", [("Al", 100.0)]),                      # stops 3.02 MeV p: passes only p > 3.4 MeV and alpha > 14 MeV
]


def residual(part, E0, layers, mu=1.0):
    E = np.atleast_1d(np.asarray(E0, float))
    for m, t in layers:
        E = st.e_after(part, m, E, t / mu)
    return E


def registered(part, E):
    E = np.atleast_1d(np.asarray(E, float))
    let = st.dedx_lin(part, "CR39", np.maximum(E, 1e-3))
    return (E > 0.01) & (let >= LET_MIN)


def main(barrier=True, water_um=0.0):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = open(os.path.join(FIGS, "m5_cr39.txt"), "w")

    def P(*a):
        print(*a)
        print(*a, file=out)

    P("== Registration: LET at CR-39 surface (keV/um) for bare incidence; threshold", LET_MIN * 1e3)
    for part, E0, lab in SPECIES:
        P(f"   {lab:22s}: LET = {1e3 * st.dedx_lin(part, 'CR39', E0):7.1f} keV/um  registered={bool(registered(part, E0)[0])}")
    P("   (p above ~4-5 MeV are below threshold for a 6 h etch: 14.7 MeV p needs a 16-24 h etch or Si)")

    pre = [("Mylar", 6.0)] if barrier else []
    if water_um:
        pre = [("D2O", water_um)] + pre
    P(f"\n== Stack response (normal incidence / 45 deg): residual energy (MeV) at CR-39; pre-layers {pre}")
    P("   species                  " + "  ".join(f"{c:>13s}" for c, _ in STACK))
    table = {}
    for part, E0, lab in SPECIES:
        row = []
        for c, f in STACK:
            e0 = residual(part, E0, pre + f, 1.0)[0]
            e45 = residual(part, E0, pre + f, np.cos(np.radians(45)))[0]
            reg = registered(part, e0)[0]
            row.append((e0, e45, reg))
            table[(lab, c)] = (e0, e45, reg)
        P(f"   {lab:22s} " + "  ".join(f"{a:5.2f}/{b:4.2f}{'*' if r else ' '}" for a, b, r in row))
    P("   (* = registered at normal incidence)")

    P("\n== Species signature (which channels record it) -> identification logic")
    for part, E0, lab in SPECIES:
        sig = "".join(c if table[(lab, c)][2] else "-" for c, _ in STACK)
        P(f"   {lab:22s}: {sig}")
    P("   3.02 MeV p  = A B C - ; t/3He = A - - - (t and 3He separated by pit diameter/depth; 3He LET 3x t);")
    P("   natural alphas (Rn, U/Th, 210Po) = A B - - or A - - - : channel C is alpha-free by construction;")
    P("   energetic alphas >14 MeV or p >3.4 MeV (registered) = A B C D. Channel D is the 'blind' control for C.")

    # thickness tolerance: C must stop 8.78 MeV alpha with margin, pass 3.02 p with E>0.5 at 45 deg
    Ra = st.rng("a", "Al", 8.785)
    P(f"\n== Filter tolerance: R(8.785 MeV alpha, Al) = {Ra:.1f} um ; R(3.02 p, Al) = {st.rng('p', 'Al', 3.02):.1f} um")
    for t in [50, 55, 60, 65, 70, 75]:
        e0 = residual("p", 3.02, pre + [("Al", t)], 1.0)[0]
        e45 = residual("p", 3.02, pre + [("Al", t)], np.cos(np.radians(45)))[0]
        a0 = residual("a", 8.785, pre + [("Al", t)], 1.0)[0]
        P(f"   Al {t:3d} um: p residual {e0:.2f} (0 deg) / {e45:.2f} (45 deg) MeV; 8.78 alpha residual {a0:.2f} MeV")
    tm = st.rng("p", "Mylar", 3.02) - st.rng("p", "Mylar", residual("p", 3.02, pre + [("Al", 55)])[0])
    P(f"   all-Mylar equivalent of 6 um Mylar + 55 um Al (same p residual): {tm:.0f} um Mylar; it stops alphas up to "
      f"{np.interp(tm, st.table('a', 'Mylar')[2], st.table('a', 'Mylar')[0]):.2f} MeV")

    # figure: residual energy vs filter thickness
    fig, axs = plt.subplots(1, 2, figsize=(11, 4))
    for ax, m in zip(axs, ["Al", "Mylar"]):
        tt = np.linspace(0, 130 if m == "Al" else 170, 300)
        for part, E0, lab in SPECIES[:4] + SPECIES[6:9]:
            ax.plot(tt, [residual(part, E0, [(m, t)])[0] for t in tt], label=lab)
        ax.axhline(0, color="k", lw=0.5)
        ax.set_xlabel(f"{m} filter thickness (µm)")
        ax.set_ylabel("residual energy at CR-39 (MeV)")
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "m5_cr39_filters.png"), dpi=130)
    plt.close(fig)

    # efficiency in C1-style contact geometry and in C3 vacuum
    P("\n== Detection efficiency of channel C (55 um Al) per emitted 3.02 MeV p (4pi), source U(0-5 um) PdD")
    for film in [0.0, 10.0, 30.0, 100.0]:
        n = 200000
        z = RNG.uniform(0, 5.0, n)
        mu = RNG.uniform(0, 1, n)
        E = es.transport("p", 3.02, z, mu, [("PdO", 0.02)])
        layers = ([("D2O", film)] if film else []) + pre + [("Al", 55.0)]
        for m, t in layers:
            E = st.e_after("p", m, E, t / mu)
        ok = registered("p", E) & (mu >= np.cos(np.radians(THETA_C["p"])))
        P(f"   electrolyte film {film:5.1f} um: eff = {0.5 * ok.mean():.4f} (x CR-39 area fraction facing the cathode)")
    out.close()


def channel_summary():
    """In-situ C1/C2 geometry: CR-39 channel C strips covering 30 % of a cathode's projected surface,
    10 um electrolyte film + 6 um Mylar + 55 um Al; background 3 proton-like tracks/cm2 per 30-day
    exposure in the channel-C pit-size window over 4 cm2 (assumption: pre-etch mapped chips,
    N2-stored; literature fresh-CR-39 background ~30 tracks/cm2 in 28 d is mostly radon alpha, which
    channel C blocks)."""
    n = 200000
    z = RNG.uniform(0, 5.0, n)
    mu = RNG.uniform(0, 1, n)
    E = es.transport("p", 3.02, z, mu, [("PdO", 0.02)])
    for m, t in [("D2O", 10.0), ("Mylar", 6.0), ("Al", 55.0)]:
        E = st.e_after("p", m, E, t / mu)
    ok = registered("p", E) & (mu >= np.cos(np.radians(THETA_C["p"])))
    eff = 0.5 * ok.mean() * 0.30
    bkg = 3.0 * 4.0 / (30 * 86400.0)
    return dict(eff=eff, bkg=bkg)


if __name__ == "__main__":
    main()
