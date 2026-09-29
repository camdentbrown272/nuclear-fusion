#!/usr/bin/env python3
"""M5 -- Si Delta-E/E telescope for the detector-facing membrane (C3).

1. Geometry scan: detected fraction vs detector area, distance and membrane size.
2. Delta-E/E particle identification (PID) Monte Carlo: p, d, t, 3He, alpha, muons.
3. Background budget in the proton windows (2.6-3.1 MeV peak and the PID window).
4. Vacuum requirements.
5. Minimum detectable D-D reaction rate vs run time.

    python3 sim/m5_si_telescope.py    # -> figs/m5_si_*.png, figs/m5_si.txt
"""
import os
import sys

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m5_stopping as st  # noqa: E402
import m5_escape as es  # noqa: E402
import m5_stats as ms  # noqa: E402

FIGS = st.FIGS
RNG = np.random.default_rng(7)
DAY = 86400.0

# --------------------------------------------------------------------------------------- baseline design
BASE = dict(
    a_mm=10.0,          # membrane active radius (3.14 cm2); set by M2/M3 -- scanned below
    d_mm=4.0,           # membrane front face -> Delta-E front face
    b_mm=13.8,          # Delta-E active radius (600 mm2)
    coll_mm=12.8,       # collimator aperture radius at the Delta-E plane (1 mm inside active edge)
    gap_mm=1.5,         # Delta-E back -> E front
    bE_mm=15.1,         # E active radius (~715 mm2 ; >= Delta-E + gap*tan(theta_max))
    tdE_um=25.0,        # Delta-E thickness
    tE_um=500.0,        # E thickness
    theta_max=50.0,     # software angular cut from strip positions (deg)
    fwhm_dE=0.050,      # MeV; 25 um x 150 mm2 quadrant ~0.6 nF -> ~40-50 keV FWHM (preamp noise slope)
    fwhm_E=0.025,       # MeV; ULTRA-class
    win_um=0.05,        # entrance windows (Si-equivalent), each detector face
    thr=0.10,           # MeV trigger threshold per detector
    overlayer=[("PdO", 0.02)],
)
MU_FLUX = 1.0 / 60.0  # muons cm^-2 s^-1 on horizontal area, sea level (PDG Cosmic-ray review,
#                       https://pdg.lbl.gov/2022/reviews/rpp2022-rev-cosmic-rays.pdf)


# --------------------------------------------------------------------------------------- telescope MC
def telescope_events(part, E0, depth_kind, depth_t, n, cfg=BASE, directions=None, E0_arr=None):
    """Isotropic emission from membrane; returns dict with dE, E deposits, hit mask, mu, weight norm.
    Detected fraction is relative to all 4pi emissions."""
    a, d, b, bc, g, bE = cfg["a_mm"], cfg["d_mm"], cfg["b_mm"], cfg["coll_mm"], cfg["gap_mm"], cfg["bE_mm"]
    r = a * np.sqrt(RNG.uniform(0, 1, n))
    ph = RNG.uniform(0, 2 * np.pi, n)
    x0, y0 = r * np.cos(ph), r * np.sin(ph)
    mu = RNG.uniform(0, 1, n)
    phi = RNG.uniform(0, 2 * np.pi, n)
    tan = np.sqrt(1 - mu ** 2) / np.maximum(mu, 1e-12)
    x1, y1 = x0 + d * tan * np.cos(phi), y0 + d * tan * np.sin(phi)
    x2, y2 = x0 + (d + g) * tan * np.cos(phi), y0 + (d + g) * tan * np.sin(phi)
    in_dE = x1 ** 2 + y1 ** 2 <= bc ** 2
    in_E = x2 ** 2 + y2 ** 2 <= bE ** 2
    ang_ok = mu >= np.cos(np.radians(cfg["theta_max"]))
    z = es.sample_depth(depth_kind, depth_t, n)
    E = np.full(n, float(E0)) if E0_arr is None else E0_arr
    Eex = np.zeros(n)
    idx = np.nonzero(in_dE)[0]
    # transport with per-event E0 (vectorised through e_after chain)
    Ei = E[idx]
    var = np.zeros_like(Ei)
    for mat, t in [(es.ACTIVE, z[idx])] + [(m, np.full_like(Ei, tt)) for m, tt in cfg["overlayer"]]:
        L = t / mu[idx]
        Ei = st.e_after(part, mat, Ei, L)
        var += st.bohr_sigma2(part, mat, L)
    Ei = np.where(Ei > 0, np.maximum(Ei + RNG.normal(0, 1, Ei.shape) * np.sqrt(var), 0), 0)
    Eex[idx] = Ei
    res = deposit(part, Eex, mu, in_dE, in_E, cfg)
    res.update(mu=mu, in_dE=in_dE, ang_ok=ang_ok, Eexit=Eex, n=n)
    return res


def deposit(part, Ein, mu, in_dE, in_E, cfg):
    t, tE, w = cfg["tdE_um"], cfg["tE_um"], cfg["win_um"]
    n = len(Ein)
    E1 = st.e_after(part, "Si", Ein, w / mu)                     # dead window of Delta-E
    E2 = st.e_after(part, "Si", E1, t / mu)                      # after Delta-E active layer
    dE = E1 - E2 + RNG.normal(0, 1, n) * np.sqrt(st.bohr_sigma2(part, "Si", t / mu)) * (E2 > 0)
    E3 = st.e_after(part, "Si", E2, 2 * w / mu)                  # back contact + E window
    E3 = np.where(in_E, E3, 0.0)
    E4 = st.e_after(part, "Si", E3, tE / mu)                     # punch-through of E
    Ed = E3 - E4
    dE = np.where(in_dE, dE + RNG.normal(0, cfg["fwhm_dE"] / 2.3548, n), 0)
    Ed = np.where(Ed > 0, Ed + RNG.normal(0, cfg["fwhm_E"] / 2.3548, n), 0)
    return dict(dE=dE, E=Ed, Etot=dE + Ed, punch=E4 > 0)


# --------------------------------------------------------------------------------------- PID
_pid_cache = {}


def dE_p_normal(Etot, cfg=BASE):
    """Expected Delta-E of a proton at normal incidence, as a function of total deposit."""
    key = (cfg["tdE_um"], cfg["win_um"])
    if key not in _pid_cache:
        E = np.linspace(0.3, 16, 4000)
        E1 = st.e_after("p", "Si", E, cfg["win_um"])
        E2 = st.e_after("p", "Si", E1, cfg["tdE_um"])
        dE = E1 - E2
        E3 = st.e_after("p", "Si", E2, 2 * cfg["win_um"])
        tot = dE + E3
        ok = E2 > 0
        _pid_cache[key] = (tot[ok], dE[ok])
    tot, dE = _pid_cache[key]
    return np.interp(Etot, tot, dE, left=np.nan, right=np.nan)


def pid_ratio(ev, cfg=BASE):
    return ev["dE"] / dE_p_normal(ev["Etot"], cfg)


def proton_cut(ev, cfg=BASE, lo=0.75, hi=None):
    """Coincidence (both > threshold) + angular cut + Delta-E within proton band."""
    if hi is None:
        hi = 1.0 / np.cos(np.radians(cfg["theta_max"])) + 0.2
    coinc = (ev["dE"] > cfg["thr"]) & (ev["E"] > cfg["thr"]) & ev["in_dE"] & ev["ang_ok"]
    with np.errstate(invalid="ignore"):
        r = pid_ratio(ev, cfg)
        return coinc & (r > lo) & (r < hi) & ~ev["punch"]


def p_window_low(cfg=BASE):
    """Lowest total proton energy reaching E above threshold at normal incidence."""
    E = np.linspace(0.2, 3.5, 3000)
    E1 = st.e_after("p", "Si", E, cfg["win_um"])
    E2 = st.e_after("p", "Si", E1, cfg["tdE_um"])
    E3 = st.e_after("p", "Si", E2, 2 * cfg["win_um"])
    return E[np.argmax(E3 > cfg["thr"] + 0.1)]


# --------------------------------------------------------------------------------------- muons
def muon_bkg(cfg=BASE, n=2_000_000):
    """Cosmic muons through horizontal Delta-E and E disks: Landau energy loss;
    probability that a muon mimics a proton in the PID window."""
    # zenith distribution for flux through horizontal area: I ~ cos^2 -> dN/dcos ~ cos^3
    c = RNG.uniform(0, 1, n) ** 0.25
    rE = cfg["bE_mm"] / 10
    area = np.pi * rE ** 2
    rate = MU_FLUX * area  # through E (approx.)
    out = {}
    dep = {}
    for key, t_um in [("dE", cfg["tdE_um"]), ("E", cfg["tE_um"])]:
        L = t_um * 1e-4 / c  # cm (ignore edge clipping; conservative)
        L = np.minimum(L, 2 * rE)
        xi = 0.1535 * (14 / 28.086) * 2.329 * L  # MeV (beta=1)
        bg = 30.0
        mpv = xi * (np.log(2 * 0.511e6 * bg ** 2 / 173.0) + np.log(xi * 1e6 / 173.0) + 0.2 - 1.0)
        lam = stats.landau.rvs(size=n, random_state=RNG)
        dep[key] = np.maximum(mpv + xi * (lam + 0.22278), 0)
    # half of muons crossing E also cross Delta-E (geometric overlap ~ area ratio)
    frac_both = (cfg["coll_mm"] / cfg["bE_mm"]) ** 2
    ev = dict(dE=dep["dE"] + RNG.normal(0, cfg["fwhm_dE"] / 2.3548, n),
              E=dep["E"] + RNG.normal(0, cfg["fwhm_E"] / 2.3548, n))
    ev["Etot"] = ev["dE"] + ev["E"]
    ev["in_dE"] = np.ones(n, bool)
    ev["ang_ok"] = np.ones(n, bool)  # muons carry no angle info from the membrane; worst case
    ev["punch"] = np.zeros(n, bool)
    wl = p_window_low(cfg)
    sel = proton_cut(ev, cfg) & (ev["Etot"] > wl) & (ev["Etot"] < 3.1)
    sel_pk = sel & (ev["Etot"] > 2.6)
    single_pk = (dep["E"] > 2.6) & (dep["E"] < 3.1)
    out["rate_s"] = rate
    out["p_pid"] = sel.mean() * frac_both
    out["p_peak"] = sel_pk.mean() * frac_both
    out["p_single_peak"] = single_pk.mean()
    out["mpv_E_keV"] = np.median(dep["E"]) * 1e3
    return out


# --------------------------------------------------------------------------------------- neutron recoils
def gammel_np(E):
    """n-p elastic cross-section (b), Gammel parameterisation (E in MeV), valid ~0-40 MeV."""
    return (3 * np.pi / (1.206 * E + (-1.86 + 0.09415 * E + 0.0001306 * E ** 2) ** 2)
            + np.pi / (1.206 * E + (0.4223 + 0.13 * E) ** 2))


def cosmic_fast_n(n):
    """Sample sea-level cosmic neutron energies 1-200 MeV from a piecewise model normalised to
    Gordon et al. IEEE TNS 51 (2004) 3427: flux(>10 MeV)=3.6e-3 cm-2 s-1 (NYC); 1-10 MeV taken
    as 2.5e-3 cm-2 s-1 (evaporation peak; +/-50 %). Returns (E, total flux >1 MeV)."""
    f1, f2 = 2.5e-3, 3.6e-3
    u = RNG.uniform(0, 1, n)
    lo = u < f1 / (f1 + f2)
    E = np.where(lo, 10 ** RNG.uniform(0, 1, n), 10 ** RNG.uniform(1, np.log10(200), n))  # ~1/E in each band
    return E, f1 + f2


def recoil_bkg(target, cfg=BASE, n=400_000):
    """Cosmic fast-neutron elastic recoils (p from PdH, d from PdD) emitted from the top 60 um
    of the membrane into the telescope. Returns rate (s^-1) in PID window & peak window."""
    En, ftot = cosmic_fast_n(n)
    if target == "p":
        sig = gammel_np(En)
        kmax = 1.0
        part = "p"
    else:
        sig = gammel_np(En)  # n-d elastic ~ n-p within ~30 % over 2-20 MeV (ENDF/B-VIII); assumption
        kmax = 8.0 / 9.0
        part = "d"
    # recoil energy uniform in [0, kmax En] (isotropic CM); direction isotropic (neutrons from all angles)
    Er = RNG.uniform(0, 1, n) * kmax * En
    depth = 60.0
    ev = telescope_events(part, 0.0, "uniform", depth, n, cfg, E0_arr=Er)
    nD = 0.9 * st.MAT["PdD0.9"][1] / st.mat_mass("PdD0.9") * st.NA  # per cm3
    # interactions per s in the (depth x area) volume: flux * sigma * N ; sample weights
    vol = np.pi * (cfg["a_mm"] / 10) ** 2 * depth * 1e-4
    w = ftot * np.mean(sig) * 1e-24 * nD * vol  # total recoils/s (isotropic emission => 4pi weights)
    sig_w = sig / np.mean(sig)
    # proton-hypothesis PID cut (d must be REJECTED; p recoils accepted)
    wl = p_window_low(cfg)
    sel = proton_cut(ev, cfg) & (ev["Etot"] > wl) & (ev["Etot"] < 3.1)
    pk = sel & (ev["Etot"] > 2.6)
    # isotropic emission: telescope_events samples only upward hemisphere -> factor 0.5
    return dict(total=w, pid=w * 0.5 * np.mean(sel * sig_w), peak=w * 0.5 * np.mean(pk * sig_w),
                single_peak=w * 0.5 * np.mean((((ev["E"] + ev["dE"]) > 2.6) & ((ev["E"] + ev["dE"]) < 3.1) & ev["in_dE"]) * sig_w))


# --------------------------------------------------------------------------------------- alpha backgrounds
U238 = [4.198, 4.775, 4.687, 4.784, 5.490, 6.002, 7.687, 5.304]   # U-238 chain (secular eq.)
TH232 = [4.012, 5.423, 5.685, 6.288, 6.778, 6.051 * 0.36 + 8.785 * 0.0, 8.785 * 0.64]  # approx.


def alpha_bkg(cfg=BASE, n=300_000):
    """Alpha backgrounds reaching the telescope: (i) bulk U/Th in the membrane (1 ppb each -> scale),
    (ii) 210Po on membrane surface, (iii) detector-intrinsic (ORTEC ULTRA-AS <=24/day, 3-8 MeV, 450 mm2).
    Returns dict of rates in single-detector peak window and after PID."""
    out = {}
    wl = p_window_low(cfg)
    # (i) bulk U/Th, uniform in top 30 um (alpha R <= 25 um in PdD) ; activity per cm3 for 1 ppb U, 1 ppb Th
    A_U = 12.4e-3 * 1e-9 * 1e6  # Bq per g for 1 ppb U-238 (12.4 Bq/mg)
    A_Th = 4.06e-3 * 1e-9 * 1e6
    rho = st.MAT["PdD0.9"][1]
    vol = np.pi * (cfg["a_mm"] / 10) ** 2 * 30e-4
    lines = [(e, A_U * rho * vol) for e in U238] + [(e, A_Th * rho * vol * (1 if e > 1 else 0)) for e in TH232]
    tot_single = tot_pid = 0.0
    for e, act in lines:
        if act == 0 or e < 1:
            continue
        ev = telescope_events("a", e, "uniform", 30.0, n // 8, cfg)
        et = ev["Etot"]
        tot_single += act * 0.5 * np.mean((et > 2.6) & (et < 3.1) & ev["in_dE"])
        tot_pid += act * 0.5 * np.mean(proton_cut(ev, cfg) & (et > wl) & (et < 3.1))
    out["bulkUTh_1ppb_single_peak"] = tot_single
    out["bulkUTh_1ppb_pid"] = tot_pid
    # (ii) surface 210Po, 1 alpha/cm2/day emitted
    ev = telescope_events("a", 5.304, "delta", 0.0, n, cfg)
    A = np.pi * (cfg["a_mm"] / 10) ** 2 / DAY
    out["Po210_surf_1perday_cm2_single_peak"] = A * 0.5 * np.mean((ev["Etot"] > 2.6) & (ev["Etot"] < 3.1))
    out["Po210_surf_1perday_cm2_pid"] = A * 0.5 * np.mean(proton_cut(ev, cfg) & (ev["Etot"] > wl) & (ev["Etot"] < 3.1))
    out["Po210_surf_1perday_cm2_detected"] = A * 0.5 * np.mean(ev["dE"] > cfg["thr"])
    return out


# --------------------------------------------------------------------------------------- main
def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = open(os.path.join(FIGS, "m5_si.txt"), "w")

    def P(*a):
        print(*a)
        print(*a, file=out)

    # ---------------------------------------------------------------- 1. geometry scan
    P("== 1. Geometry: detected fraction of 4pi (3.02 MeV p, source U(0-5 um) under 20 nm PdO), telescope accept.")
    areas = [150, 300, 450, 600, 900, 1200]
    ds = [2, 3, 4, 6, 8, 12, 20]
    fig, axs = plt.subplots(1, 2, figsize=(11, 4))
    table = {}
    for a in [5.0, 10.0, 15.0]:
        P(f"  membrane radius a = {a} mm ({np.pi * a * a / 100:.2f} cm2):  rows=Delta-E area (mm2), cols=d (mm) {ds}")
        for A in areas:
            b = np.sqrt(A / np.pi)
            row = []
            for d in ds:
                cfg = dict(BASE, a_mm=a, b_mm=b, coll_mm=b - 1.0, d_mm=d, bE_mm=b + 1.5 * np.tan(np.radians(60)))
                ev = telescope_events("p", 3.02, "uniform", 5.0, 60000, cfg)
                sel = proton_cut(ev, cfg) & (ev["Etot"] > 0.5)
                row.append(0.5 * sel.mean())
            table[(a, A)] = row
            P(f"    {A:5d}: " + " ".join(f"{x:.3f}" for x in row))
            if a == 10.0:
                axs[0].plot(ds, row, "o-", label=f"{A} mm²")
    axs[0].set_xlabel("membrane–ΔE distance d (mm)")
    axs[0].set_ylabel("p detected / p emitted (4π)")
    axs[0].set_title("a = 10 mm membrane, θ ≤ 50°, PID-accepted")
    axs[0].legend(fontsize=7)
    # angular cut effect
    cuts = [30, 40, 50, 60, 70, 90]
    effs = []
    for tm in cuts:
        cfg = dict(BASE, theta_max=tm)
        ev = telescope_events("p", 3.02, "uniform", 5.0, 80000, cfg)
        effs.append(0.5 * (proton_cut(ev, cfg, hi=1 / np.cos(np.radians(min(tm, 80))) + 0.2) & (ev["Etot"] > 0.5)).mean())
    axs[1].plot(cuts, effs, "o-")
    axs[1].set_xlabel("θ_max cut (deg)")
    axs[1].set_ylabel("p detected / emitted (baseline geometry)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "m5_si_geometry.png"), dpi=130)
    plt.close(fig)
    P("  theta_max scan (baseline): " + ", ".join(f"{c}deg:{e:.3f}" for c, e in zip(cuts, effs)))

    # ---------------------------------------------------------------- 2. Delta-E thickness -> window
    P("\n== 2. Delta-E thickness vs lowest identifiable proton energy and alpha stopping")
    for t in [10, 15, 20, 25, 40, 65]:
        cfg = dict(BASE, tdE_um=t)
        P(f"   Delta-E {t:3d} um: p window low edge = {p_window_low(cfg):.2f} MeV; "
          f"max alpha stopped in Delta-E = {np.interp(t, st.table('a', 'Si')[2], st.table('a', 'Si')[0]):.2f} MeV; "
          f"t stopped up to {np.interp(t, st.table('t', 'Si')[2], st.table('t', 'Si')[0]):.2f} MeV; "
          f"dE(3.02 p) = {1e3 * (3.02 - st.e_after('p', 'Si', 3.02, t)):.0f} keV")
    wl = p_window_low()
    P(f"   baseline ({BASE['tdE_um']} um): PID proton window = {wl:.2f}-3.10 MeV")

    # ---------------------------------------------------------------- 3. PID scatter
    fig, ax = plt.subplots(figsize=(7, 5.5))
    evs = {}
    for part, E0, kind, t, col, lab in [("p", 3.02, "uniform", 20.0, "C0", "p 3.02, U(0–20 µm)"),
                                        ("d", 3.0, "uniform", 20.0, "C2", "d (recoil-like) 3 MeV"),
                                        ("t", 1.01, "uniform", 3.0, "C1", "t 1.01"),
                                        ("h", 0.82, "uniform", 1.0, "C4", "³He 0.82"),
                                        ("a", 5.304, "delta", 0.0, "C3", "²¹⁰Po α surface"),
                                        ("a", 8.785, "uniform", 30.0, "C5", "α 8.78 bulk"),
                                        ("a", 12.0, "uniform", 20.0, "C6", "α 12 MeV"),
                                        ("t", 4.75, "uniform", 20.0, "C7", "t 4.75"),
                                        ("h", 4.75, "uniform", 5.0, "C8", "³He 4.75")]:
        ev = telescope_events(part, E0, kind, t, 40000, BASE)
        m = (ev["dE"] > BASE["thr"]) & ev["in_dE"] & ev["ang_ok"]
        evs[lab] = ev
        ax.scatter(ev["Etot"][m][:3000], ev["dE"][m][:3000], s=1, color=col, label=lab)
    mu = muon_bkg(n=200000)
    Et = np.linspace(wl, 3.1, 50)
    ax.plot(Et, 0.75 * dE_p_normal(Et), "k--", lw=0.8)
    ax.plot(Et, (1 / np.cos(np.radians(BASE["theta_max"])) + 0.2) * dE_p_normal(Et), "k--", lw=0.8, label="proton PID band")
    ax.set_xlabel("E_total = ΔE + E (MeV)")
    ax.set_ylabel("ΔE (MeV), 25 µm")
    ax.set_xlim(0, 13)
    ax.set_ylim(0, 6)
    ax.legend(fontsize=7, markerscale=5)
    ax.set_title("C3 telescope: 25 µm ΔE / 500 µm E, θ ≤ 50°")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "m5_si_pid.png"), dpi=130)
    plt.close(fig)

    P("\n== 3. PID performance (fractions of detected-in-Delta-E events passing the proton cut in the PID window)")
    for lab, ev in evs.items():
        det = (ev["dE"] > BASE["thr"]) & ev["in_dE"] & ev["ang_ok"]
        sel = proton_cut(ev) & (ev["Etot"] > wl) & (ev["Etot"] < 3.1)
        P(f"   {lab:22s}: passes p-cut {sel.sum() / max(det.sum(), 1):.2e}  ({sel.sum()} of {det.sum()})")

    # proton signal efficiency for source distributions
    P("\n== 4. Signal efficiency (detected in window / emitted p, 4pi), baseline telescope")
    eff = {}
    for kind, t, name in [("delta", 0.0, "surface"), ("uniform", 1.0, "U(0-1um)"), ("uniform", 5.0, "U(0-5um)"),
                          ("uniform", 10.0, "U(0-10um)"), ("uniform", 20.0, "U(0-20um)"), ("uniform", 33.3, "U(0-33um)=R_p"),
                          ("uniform", 100.0, "U(0-100um)")]:
        ev = telescope_events("p", 3.02, kind, t, 200000, BASE)
        sel = proton_cut(ev)
        e_pid = 0.5 * np.mean(sel & (ev["Etot"] > wl) & (ev["Etot"] < 3.1))
        e_pk = 0.5 * np.mean(sel & (ev["Etot"] > 2.6) & (ev["Etot"] < 3.1))
        eff[name] = (e_pk, e_pid)
        P(f"   {name:16s}: peak(2.6-3.1) {e_pk:.4f}   PID window({wl:.2f}-3.1) {e_pid:.4f}")
    for part, E0, t, name in [("t", 1.01, 0.0, "triton surface (Delta-E only, 0.8-1.05 MeV)"),
                              ("h", 0.82, 0.0, "3He surface (Delta-E only, 0.6-0.85 MeV)")]:
        ev = telescope_events(part, E0, "delta", t, 100000, BASE)
        lo, hi = (0.8, 1.05) if part == "t" else (0.6, 0.85)
        e = 0.5 * np.mean(ev["in_dE"] & ev["ang_ok"] & (ev["E"] < BASE["thr"]) & (ev["dE"] > lo) & (ev["dE"] < hi))
        eff[name] = (e, e)
        P(f"   {name:42s}: {e:.4f}")

    # ---------------------------------------------------------------- 5. background budget
    P("\n== 5. Background budget per telescope (counts/day), baseline, sea level, no neutron shield")
    mu = muon_bkg()
    P(f"   muon rate through E: {mu['rate_s']:.3f}/s; MPV in E = {mu['mpv_E_keV']:.0f} keV")
    rp = recoil_bkg("p")
    rd = recoil_bkg("d")
    al = alpha_bkg()
    # detector-intrinsic alphas (ORTEC ULTRA-AS warranty <=24 counts/day in 3-8 MeV for 450 mm2;
    # https://www.ortec-online.com/-/media/ametekortec/brochures/u/ultra_ultra-as-a4.pdf) scaled to 600 mm2
    det_alpha_day = 24.0 * 600 / 450
    # Si(n,p)/(n,a) inside Delta-E: flux(>5 MeV) ~ 4.5e-3 cm-2 s-1 (Gordon 2004 spectrum), sigma ~0.45 b
    nSi = st.MAT["Si"][1] / 28.086 * st.NA
    vdE = np.pi * (BASE["b_mm"] / 10) ** 2 * BASE["tdE_um"] * 1e-4
    si_np = 4.5e-3 * 0.45e-24 * nSi * vdE
    budget = [
        # name, single-detector peak window (no PID), PID window after cuts (per day)
        ("cosmic muons (no veto)", mu["rate_s"] * mu["p_single_peak"] * DAY, mu["rate_s"] * mu["p_pid"] * DAY),
        ("detector-intrinsic alphas (spec max)", det_alpha_day * 0.5 / 5.0, det_alpha_day * 1e-4),
        ("membrane bulk U+Th, 1 ppb each", al["bulkUTh_1ppb_single_peak"] * DAY, max(al["bulkUTh_1ppb_pid"] * DAY, 1e-4 * al["bulkUTh_1ppb_single_peak"] * DAY * 10)),
        ("membrane surface 210Po, 1 /cm2/day", al["Po210_surf_1perday_cm2_single_peak"] * DAY, max(al["Po210_surf_1perday_cm2_pid"], 1e-4 * al["Po210_surf_1perday_cm2_detected"]) * DAY),
        ("n-p recoils from residual H (H/D=1 %)", 0.01 * rp["single_peak"] * DAY, 0.01 * rp["pid"] * DAY),
        ("n-d recoils (misID as p)", rd["single_peak"] * DAY, rd["pid"] * DAY),
        ("Si(n,p),(n,a) in Delta-E (upper bound)", si_np * DAY * 0.1, si_np * DAY * 0.3),
    ]
    tot_s = tot_p = 0
    P(f"   {'source':42s} {'single Si, 2.6-3.1':>20s} {'telescope PID window':>22s}")
    for name, s1, s2 in budget:
        tot_s += s1
        tot_p += s2
        P(f"   {name:42s} {s1:20.3g} {s2:22.3g}")
    P(f"   {'TOTAL':42s} {tot_s:20.3g} {tot_p:22.3g}")
    P(f"   n-p recoil rate in a PdH control (H/Pd=0.9) would be {rp['pid'] * DAY:.3g}/day in the PID window "
      f"(this is a REAL-proton background specific to the H control; total recoils {rp['total'] * DAY:.2g}/day)")
    P(f"   with muon veto (99 %) + 10 cm HDPE around chamber (fast-n x0.5 assumed): see MDA table")
    B_design = tot_p
    B_single = tot_s

    # ---------------------------------------------------------------- 6. vacuum
    P("\n== 6. Vacuum requirements")
    for part, E0 in [("p", 3.02), ("t", 1.01), ("h", 0.82)]:
        Sair = st.stopping(part, "air", E0) * 1e3  # keV cm2/g
        # energy loss over 10 mm path at pressure P (mbar), air
        for Pm in [1000, 10, 1, 1e-2]:
            rho = 1.205e-3 * Pm / 1013.25
            P(f"   {part} {E0} MeV: loss over 10 mm at {Pm:g} mbar air = {Sair * rho * 1.0:.3g} keV")
    # D2 gas load from permeation
    for J in [1e15, 1e16, 1e17]:
        Q = J * np.pi * (BASE["a_mm"] / 10) ** 2 / 2 * 1.380649e-23 * 293 * 1e4  # mbar L/s  (D2 molecules)
        P(f"   D permeation flux {J:.0e} D/cm2/s -> gas load {Q:.2e} mbar L/s -> P = {Q / 50:.1e} mbar with 50 L/s (D2) pumping")

    # ---------------------------------------------------------------- 7. MDA
    P("\n== 7. Minimum detectable D-D fusion rate (fusions/s; BR(p+t)=0.5) -- 5 sigma discovery at 50 % power,"
      " background known from equal-duration H-control/off data (see m5_stats)")
    scen = {
        "single Si, 2.6-3.1 (no PID), 30 d": (eff["U(0-5um)"][0], B_single),
        "telescope peak 2.6-3.1, 30 d": (eff["U(0-5um)"][0], B_design * 0.4),
        "telescope PID window, 30 d": (eff["U(0-5um)"][1], B_design),
        "telescope PID, +veto+HDPE, 30 d": (eff["U(0-5um)"][1], B_design * 0.3),
    }
    mda_rows = []
    for name, (e, B) in scen.items():
        for days in [1, 14, 30]:
            b = B * days
            s = ms.discovery_signal(b, alpha_sigma=5.0, power=0.5, b_unc_rel=1 / np.sqrt(max(b, 1.0)) if b > 0 else 0)
            rate = s / (e * 0.5 * days * DAY)
            mda_rows.append((name, days, e, B, s, rate))
            P(f"   {name:36s} {days:3d} d: eff={e:.3f}, B={B:.3g}/day -> s_5sig={s:.1f} counts -> {rate:.2e} fusions/s")
    out.close()
    return dict(eff=eff, B_design=B_design, B_single=B_single, budget=budget, mda=mda_rows)


if __name__ == "__main__":
    main()
