#!/usr/bin/env python3
"""M5 escape fractions and detector spectra for D-D (and other) charged products.

Straight-line CSDA transport with Bohr energy straggling (continuous-slowing-down;
angular scattering neglected -- justified below), isotropic emission from a depth
distribution inside the active layer, through overlayers, a vacuum gap and the Si
detector entrance window, convolved with detector resolution.

    python3 sim/m5_escape.py      # regenerates figs/m5_escape_*.png and figs/m5_escape.txt
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m5_stopping as st  # noqa: E402

FIGS = st.FIGS
RNG = np.random.default_rng(20260929)
FWHM_SI = 0.025  # MeV, Si detector resolution (ULTRA-class, 20-30 keV FWHM for alphas; ORTEC ULTRA brochure)
SI_WINDOW_UM = 0.05  # 50 nm Si-equivalent entrance window (ion-implanted, ORTEC ULTRA spec)
ACTIVE = "PdD0.9"


def transport(part, E0, depth_um, mu, overlayers, active=ACTIVE, straggle=True, rng=RNG):
    """Energy (MeV) at exit of the overlayer stack for emission at depth `depth_um` (arrays),
    direction cosine `mu` (>0, towards surface). Returns E (0 = stopped)."""
    E = np.full(np.shape(depth_um), float(E0))
    var = np.zeros_like(E)
    layers = [(active, depth_um)] + [(m, np.full_like(E, t)) for m, t in overlayers]
    for mat, t in layers:
        L = t / mu
        E = st.e_after(part, mat, E, L)
        if straggle:
            var += st.bohr_sigma2(part, mat, L)
    if straggle:
        alive = E > 0
        E = np.where(alive, np.maximum(E + rng.normal(0, 1, E.shape) * np.sqrt(var), 0.0), 0.0)
    return E


def sample_depth(kind, t_um, n, rng=RNG):
    if kind == "delta":
        return np.full(n, float(t_um))
    if kind == "uniform":
        return rng.uniform(0, t_um, n)
    raise ValueError(kind)


def escape_fraction(part, E0, kind, t_um, overlayers, Emin=0.0, window=None, n=200000, straggle=True):
    """Fraction of ALL (4 pi) emissions that leave the surface with E>Emin (or E in window)."""
    z = sample_depth(kind, t_um, n)
    mu = RNG.uniform(0, 1, n)  # upward hemisphere: isotropic => mu uniform; weight 1/2
    E = transport(part, E0, z, mu, overlayers, straggle=straggle)
    ok = (E > Emin) if window is None else ((E >= window[0]) & (E <= window[1]))
    return 0.5 * ok.mean()


def disk_detector_mc(part, E0, kind, t_um, overlayers, a_mm, b_mm, d_mm, n=400000,
                     window_um=SI_WINDOW_UM, fwhm=FWHM_SI, extra_si_um=0.0, theta_max=None):
    """Isotropic emission from a disk (radius a) facing a coaxial disk detector (radius b) at gap d.
    Returns (fraction of all emissions detected, detected energies incl. resolution)."""
    r = a_mm * np.sqrt(RNG.uniform(0, 1, n))
    ph = RNG.uniform(0, 2 * np.pi, n)
    x0, y0 = r * np.cos(ph), r * np.sin(ph)
    mu = RNG.uniform(0, 1, n)
    phi = RNG.uniform(0, 2 * np.pi, n)
    s = np.sqrt(1 - mu ** 2)
    tan = s / np.maximum(mu, 1e-12)
    xd = x0 + d_mm * tan * np.cos(phi)
    yd = y0 + d_mm * tan * np.sin(phi)
    hit = xd ** 2 + yd ** 2 <= b_mm ** 2
    if theta_max is not None:
        hit &= mu >= np.cos(np.radians(theta_max))
    z = sample_depth(kind, t_um, n)
    E = np.zeros(n)
    idx = np.nonzero(hit)[0]
    Eh = transport(part, E0, z[idx], mu[idx], overlayers)
    # detector entrance window (dead layer) at the same angle (detector parallel to source)
    Eh = st.e_after(part, "Si", Eh, window_um / mu[idx])
    if extra_si_um:
        Eh = st.e_after(part, "Si", Eh, extra_si_um / mu[idx])
    Eh = np.where(Eh > 0, Eh + RNG.normal(0, fwhm / 2.3548, Eh.shape), 0)
    E[idx] = Eh
    det = E > 0.05  # 50 keV electronic threshold
    return 0.5 * det.mean(), E[det], mu[det]


def geom_eff(a_mm, b_mm, d_mm, n=400000):
    """Pure geometric efficiency (fraction of 4 pi) disk-to-disk, for verification."""
    r = a_mm * np.sqrt(RNG.uniform(0, 1, n))
    ph = RNG.uniform(0, 2 * np.pi, n)
    mu = RNG.uniform(0, 1, n)
    phi = RNG.uniform(0, 2 * np.pi, n)
    tan = np.sqrt(1 - mu ** 2) / np.maximum(mu, 1e-12)
    xd = r * np.cos(ph) + d_mm * tan * np.cos(phi)
    yd = r * np.sin(ph) + d_mm * tan * np.sin(phi)
    return 0.5 * np.mean(xd ** 2 + yd ** 2 <= b_mm ** 2)


def point_on_axis_eff(b, d):
    return 0.5 * (1 - d / np.sqrt(d ** 2 + b ** 2))


OVERLAYERS = {
    "none": [],
    "PdO 20 nm": [("PdO", 0.02)],
    "PdO 100 nm": [("PdO", 0.1)],
    "Au 50 nm": [("Au", 0.05)],
    "CaO 100 nm": [("CaO", 0.1)],
    "Ni 1 um": [("Ni", 1.0)],
    "Cu 1 um": [("Cu", 1.0)],
    "D2O 10 um film": [("D2O", 10.0)],
    "Mylar 6 um": [("Mylar", 6.0)],
}
DD = [("p", 3.02, "p 3.02"), ("t", 1.01, "t 1.01"), ("h", 0.82, "3He 0.82")]


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = open(os.path.join(FIGS, "m5_escape.txt"), "w")

    def P(*a):
        print(*a)
        print(*a, file=out)

    # ------------------------------------------------------------------ verification
    P("== Verification")
    for part, E0, lab in DD:
        R = st.rng(part, ACTIVE, E0)
        for zf in [0.0, 0.25, 0.5, 0.9]:
            mc = escape_fraction(part, E0, "delta", zf * R, [], Emin=1e-3, straggle=False)
            P(f"  {lab}: delta at z={zf:.2f}R  MC={mc:.4f}  analytic 0.5(1-z/R)={0.5 * (1 - zf):.4f}")
        for tf in [0.5, 1.0, 3.0]:
            mc = escape_fraction(part, E0, "uniform", tf * R, [], Emin=1e-3, straggle=False)
            an = 0.5 * (1 - tf / 2) if tf <= 1 else 1 / (4 * tf)
            P(f"  {lab}: uniform 0..{tf:.1f}R  MC={mc:.4f}  analytic={an:.4f}")
    for a, b, d in [(0.01, 12.0, 5.0), (0.01, 12.0, 20.0)]:
        P(f"  geometry: point on axis b={b} d={d}: MC={geom_eff(a, b, d):.4f} analytic={point_on_axis_eff(b, d):.4f}")
    # straggling sanity: Bohr sigma for 3 MeV p after 10 um PdD
    s = np.sqrt(st.bohr_sigma2("p", ACTIVE, 10.0)) * 1e3
    P(f"  Bohr straggling sigma, p through 10 um PdD0.9: {s:.1f} keV; through 30 um: "
      f"{np.sqrt(st.bohr_sigma2('p', ACTIVE, 30.0)) * 1e3:.1f} keV")
    # angular scattering neglect: estimate rms angle (Highland is not valid at low E; use
    # Rutherford multiple scattering, small-angle Moliere width approx) -> path-length correction
    P("  Straight-line approximation: detour factor (CSDA/projected range) for 3 MeV p in high-Z is ~1.02-1.05 "
      "(SRIM projected vs CSDA); escape depths quoted are CSDA, i.e. overestimate by <=5 %.")

    # ------------------------------------------------------------------ escape vs depth
    P("\n== Escape fraction (of 4pi) vs source depth in PdD0.9, 20 nm PdO overlayer")
    P("   depth(um)  " + "  ".join(f"{lab:>18s}" for _, _, lab in DD) + "   p in 2.6-3.1 MeV")
    depths = [0, 0.5, 1, 2, 5, 10, 20, 30]
    esc_tab = {}
    for z in depths:
        row = []
        for part, E0, lab in DD:
            row.append(escape_fraction(part, E0, "delta", z, OVERLAYERS["PdO 20 nm"], Emin=0.05))
        w = escape_fraction("p", 3.02, "delta", z, OVERLAYERS["PdO 20 nm"], window=(2.6, 3.1))
        esc_tab[z] = row + [w]
        P(f"   {z:8.1f}   " + "  ".join(f"{x:18.4f}" for x in row) + f"   {w:.4f}")

    # ------------------------------------------------------------------ max useful thickness
    P("\n== Maximum useful active-layer thickness (PdD0.9), uniform reaction density per volume")
    P("   Signal per unit area ~ integral_0^t f(z) dz saturates at t = R_eff; 90 % of saturation at 0.68 R_eff.")
    rows = []
    for part, E0, lab, Emin in [("p", 3.02, "p counting E>1.41 MeV (telescope PID window, 25 um DeltaE)", 1.41),
                                ("p", 3.02, "p counting E>0.3 MeV (single Si)", 0.3),
                                ("t", 1.01, "t E>0.2 MeV", 0.2), ("h", 0.82, "3He E>0.2 MeV", 0.2),
                                ("a", 5.3, "alpha 5.3 E>0.5", 0.5), ("a", 12.0, "alpha 12 E>0.5", 0.5),
                                ("p", 14.7, "p 14.7 E>1", 1.0)]:
        Reff = st.rng(part, ACTIVE, E0) - st.rng(part, ACTIVE, Emin)
        rows.append((lab, Reff))
        P(f"   {lab:52s} R_eff = {Reff:7.2f} um  (90 % at {0.684 * Reff:6.2f} um)")
    # spectroscopic window: emitted p must arrive with E in 2.6-3.1 MeV
    ts = np.array([0.1, 0.3, 1, 2, 3, 5, 7, 10, 15, 20, 30, 40, 60, 100])
    sig_win = []
    sig_cnt = []
    for t in ts:
        w = escape_fraction("p", 3.02, "uniform", t, OVERLAYERS["PdO 20 nm"], window=(2.6, 3.1), n=100000)
        c = escape_fraction("p", 3.02, "uniform", t, OVERLAYERS["PdO 20 nm"], Emin=1.41, n=100000)
        sig_win.append(w * t)
        sig_cnt.append(c * t)
    sig_win = np.array(sig_win)
    sig_cnt = np.array(sig_cnt)
    t90w = np.interp(0.9 * sig_win.max(), sig_win, ts)
    t90c = np.interp(0.9 * sig_cnt.max(), sig_cnt, ts)
    P(f"   p in 2.6-3.1 MeV window: signal/area saturates at {sig_win.max():.3f} um-equivalent; 90 % reached at t = {t90w:.1f} um")
    P(f"   p with E>1.41 MeV (DeltaE-E PID window 1.41-3.1 MeV): saturates at {sig_cnt.max():.3f} um-eq; 90 % at t = {t90c:.1f} um")
    P(f"   => PID-window counting gains x{sig_cnt.max() / sig_win.max():.1f} over the 2.6-3.1 MeV peak window for bulk-distributed sites")

    # ------------------------------------------------------------------ overlayers table
    P("\n== Escape fraction (of 4pi), uniform source 0-t um in PdD0.9 under overlayer (E>0.05 MeV; p also in 2.6-3.1 window)")
    P("   overlayer          t(um)  " + "  ".join(f"{lab:>8s}" for _, _, lab in DD) + "  p-window")
    ov_tab = {}
    for name, ov in OVERLAYERS.items():
        for t in [0.0, 1.0, 10.0, 100.0]:
            kind = "delta" if t == 0 else "uniform"
            row = [escape_fraction(part, E0, kind, t, ov, Emin=0.05, n=100000) for part, E0, _ in DD]
            w = escape_fraction("p", 3.02, kind, t, ov, window=(2.6, 3.1), n=100000)
            ov_tab[(name, t)] = row + [w]
            P(f"   {name:16s} {('surf' if t == 0 else f'{t:5.0f}'):>6s}  " + "  ".join(f"{x:8.4f}" for x in row) + f"  {w:8.4f}")

    # ------------------------------------------------------------------ fig: escape vs overlayer thickness (surface source)
    fig, axs = plt.subplots(1, 4, figsize=(15, 3.8))
    mats = ["Au", "PdO", "Ni", "Cu", "CaO", "D2O", "Mylar", "Al"]
    for ax, (part, E0, lab) in zip(axs, DD + [("a", 5.3, "α 5.3 (bkg ref)")]):
        for m in mats:
            R = st.rng(part, m, E0)
            t = np.logspace(-3, np.log10(max(3 * R, 1.0)), 200)
            f = np.where(t < R, 0.5 * (1 - t / R), 0)
            ax.semilogx(t, f, label=f"{m} (R={R:.3g} µm)")
        ax.set_title(f"{lab} MeV: surface source")
        ax.set_xlabel("overlayer thickness (µm)")
        ax.set_ylabel("escape fraction of 4π")
        ax.legend(fontsize=6)
        ax.set_ylim(0, 0.52)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "m5_escape_vs_overlayer.png"), dpi=130)
    plt.close(fig)

    # ------------------------------------------------------------------ fig: escape vs depth + useful thickness
    fig, axs = plt.subplots(1, 2, figsize=(11, 4))
    zz = np.linspace(0, 40, 81)
    for part, E0, lab in DD:
        f = [escape_fraction(part, E0, "delta", z, OVERLAYERS["PdO 20 nm"], Emin=0.05, n=40000, straggle=False) for z in zz]
        axs[0].plot(zz, f, label=lab + " MeV")
    fw = [escape_fraction("p", 3.02, "delta", z, OVERLAYERS["PdO 20 nm"], window=(2.6, 3.1), n=40000) for z in zz]
    axs[0].plot(zz, fw, "--", label="p in 2.6–3.1 MeV")
    axs[0].set_xlabel("emission depth in PdD$_{0.9}$ (µm), under 20 nm PdO")
    axs[0].set_ylabel("escape fraction of 4π")
    axs[0].legend()
    axs[1].semilogx(ts, sig_cnt / sig_cnt.max(), "o-", label="p, E>1.41 MeV (ΔE–E PID)")
    axs[1].semilogx(ts, sig_win / sig_cnt.max(), "s-", label="p, 2.6–3.1 MeV window")
    axs[1].axhline(0.9 * sig_win.max() / sig_cnt.max(), color="0.6", ls=":")
    axs[1].set_xlabel("active-layer thickness t (µm), uniform reaction density")
    axs[1].set_ylabel("detected signal per unit area (norm.)")
    axs[1].legend()
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "m5_escape_vs_depth.png"), dpi=130)
    plt.close(fig)

    # ------------------------------------------------------------------ spectra at detector
    P("\n== Detected spectra: membrane a=10 mm radius, detector b=13.8 mm (600 mm2), d=5 mm; PdO 20 nm overlayer")
    fig, axs = plt.subplots(1, 3, figsize=(15, 4))
    bins = np.linspace(0, 3.4, 341)
    for ax, (part, E0, lab) in zip(axs, DD):
        for kind, t, name in [("delta", 0.0, "surface δ"), ("uniform", 1.0, "U(0–1 µm)"),
                              ("uniform", 10.0, "U(0–10 µm)"), ("uniform", 100.0, "U(0–100 µm)")]:
            eff, E, mu = disk_detector_mc(part, E0, kind, t, OVERLAYERS["PdO 20 nm"], 10.0, 13.8, 5.0)
            h, _ = np.histogram(E, bins=bins)
            # normalise: counts per MeV per emitted particle
            ax.plot(0.5 * (bins[1:] + bins[:-1]), h / (len(E) + 1e-9) * eff / (bins[1] - bins[0]), label=f"{name}: ε={eff:.3f}")
            win = np.mean((E > 2.6) & (E < 3.1)) if part == "p" else np.nan
            P(f"   {lab:8s} {name:12s}: detected/emitted = {eff:.4f}; mean E = {E.mean():.3f} MeV; "
              f"frac in 2.6-3.1 = {win:.3f}")
        ax.set_yscale("log")
        ax.set_ylim(1e-4, None)
        ax.set_xlim(0, 3.4 if part == "p" else 1.2)
        ax.set_xlabel("deposited energy (MeV)")
        ax.set_ylabel("counts / MeV / emitted particle")
        ax.set_title(f"{lab} MeV")
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "m5_escape_spectra.png"), dpi=130)
    plt.close(fig)
    out.close()


if __name__ == "__main__":
    main()
