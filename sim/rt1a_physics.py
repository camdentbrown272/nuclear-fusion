#!/usr/bin/env python3
"""Red-team 1A physics checks on docs/design/iteration-1.md (DFM-4).

Sections (output -> docs/design/figs/rt1a_physics.txt, figure rt1a_tomography.png):
  1. Depth tomography: 3.02 MeV p, t, 3He, alphas, 6Li(n,t) tritons through loaded PdD0.9 + 4 mm of
     0.5 bar D2, at normal incidence AND over the real acceptance (emission angle smears depth).
  2. Quadrant (x, J) regimes with the M3 isotherm / Kirchhoff potential: F quadrant, 15 um,
     0.5 bar D2 exit, entry face charge-limited.
  3. Electrolyte-borne recoil backgrounds through the foil (n-p in H2O twin / D2O with H, n-d in D2O),
     and n-d recoils in the 0.5 bar D2 front gas, using M5's telescope deposit and PID code.
  4. He-4 / He-3 background inventory for the cell-headspace and front channels (order-of-magnitude).
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import m5_stopping as st  # noqa: E402
import m5_si_telescope as tel  # noqa: E402
import m3_common as m3  # noqa: E402

ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "docs", "design", "figs")
os.makedirs(OUT, exist_ok=True)
RNG = np.random.default_rng(11)
LINES = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LINES.append(s)


# ---------------------------------------------------------------- geometry of the design (iteration-1 §3)
L_PD0 = 15.0                      # um, unloaded Pd
L_LOADED = L_PD0 * (4.04 / 3.89)  # um, beta-PdD0.9 lattice expansion (a 3.89 -> 4.04 A), same Pd areal density
GAS_UM_EQ = 4000.0 * 0.5          # 4 mm of 0.5 bar D2 == 2 mm of 1 bar D2gas (m5_stopping density is 1 bar)
CFG = dict(tel.BASE)
CFG["overlayer"] = []


def through(part, E0, z_um, mu, gas=True):
    """Energy after depth z_um of PdD0.9 and (optionally) the D2 gap at direction cosine mu."""
    E = st.e_after(part, "PdD0.9", E0, np.asarray(z_um) / mu)
    if gas:
        E = st.e_after(part, "D2gas", E, GAS_UM_EQ / mu)
    return E


# ---------------------------------------------------------------- 1. tomography
def tomography():
    P("== 1. Depth tomography (loaded thickness %.2f um PdD0.9; gas = 4 mm x 0.5 bar D2)" % L_LOADED)
    P("   iteration-1 table: exit 3.00 / mid 2.55 / entry 2.05 MeV (p); 214Po alpha 1.26 MeV")
    rows = [("p", 3.02), ("t", 1.01), ("h", 0.82), ("a", 7.687), ("a", 8.785), ("a", 6.0), ("t", 2.73), ("a", 2.05)]
    for part, E0 in rows:
        vals = []
        for z in [0.0, L_LOADED / 2, 13.0, 15.0, L_LOADED, 17.0 * 4.04 / 3.89]:
            vals.append("z=%5.2f:%6.3f" % (z, float(through(part, E0, z, 1.0))))
        P("  %s %.3f MeV  " % (part, E0) + "  ".join(vals))
    # angle: acceptance of the 4 mm / r_coll 12.8 mm telescope from a 10 mm radius membrane
    n = 200_000
    r = 10.0 * np.sqrt(RNG.uniform(0, 1, n))
    ph = RNG.uniform(0, 2 * np.pi, n)
    mu = RNG.uniform(0, 1, n)
    phi = RNG.uniform(0, 2 * np.pi, n)
    tan = np.sqrt(1 - mu ** 2) / np.maximum(mu, 1e-9)
    x1 = r * np.cos(ph) + 4.0 * tan * np.cos(phi)
    y1 = r * np.sin(ph) + 4.0 * tan * np.sin(phi)
    hit = x1 ** 2 + y1 ** 2 <= 12.8 ** 2
    muh = mu[hit]
    P("  accepted direction cosines: median %.2f, 10th pct %.2f, 5th pct %.2f (theta_90%% = %.0f deg)"
      % (np.median(muh), np.percentile(muh, 10), np.percentile(muh, 5), np.degrees(np.arccos(np.percentile(muh, 10)))))
    # energy at the detector for entry-face vs mid-foil vs exit-face protons over the acceptance
    res = {}
    for lab, z in [("exit", 0.0), ("z=4", 4.0), ("mid", L_LOADED / 2), ("entry", L_LOADED)]:
        E = through("p", 3.02, z, muh)
        res[lab] = E
        P("  p from %-5s: E_det median %.3f, 10-90%% [%.3f, %.3f] MeV, frac in peak(2.6-3.1) %.2f"
          % (lab, np.median(E), np.percentile(E, 10), np.percentile(E, 90), np.mean((E > 2.6) & (E < 3.1))))
    # overlap: probability an entry-face proton is reconstructed shallower than mid-foil by energy alone
    thr = np.median(res["mid"])
    P("  P(entry-face p has E above median mid-foil E) = %.3f; P(mid-foil p below median entry E) = %.3f"
      % (np.mean(res["entry"] > thr), np.mean(res["mid"] < np.median(res["entry"]))))
    # triton marker (ADR-004 6LiF on entry face)
    Et = through("t", 2.73, L_LOADED, muh)
    P("  6LiF tritons (2.73 MeV) from entry face through %.1f um: emerge %.0f %% of accepted, median %.3f MeV "
      "(<0.10 MeV threshold for most); at 12 um nominal they would give %.3f MeV"
      % (L_LOADED, 100 * np.mean(Et > 0.1), np.median(Et), float(through("t", 2.73, 12.0 * 4.04 / 3.89, 1.0))))
    # figure
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6.4, 3.8))
        bins = np.linspace(1.0, 3.1, 85)
        for lab, c in [("exit", "C0"), ("z=4", "C1"), ("mid", "C2"), ("entry", "C3")]:
            ax.hist(res[lab], bins, histtype="step", color=c, label=f"p born at {lab}")
        ax.axvspan(2.6, 3.1, color="0.9", zorder=0)
        ax.set_xlabel("proton energy at dE (MeV)")
        ax.set_ylabel("accepted events (arb.)")
        ax.set_title("RT-1A: depth tomography over the real acceptance")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(OUT, "rt1a_tomography.png"), dpi=130)
    except Exception as e:  # pragma: no cover
        P("  (figure skipped: %s)" % e)


# ---------------------------------------------------------------- 2. regimes
def regimes():
    P("\n== 2. Quadrant regimes (M3 isotherm/Kirchhoff potential, PdD, 295 K)")
    T = 295.0
    mat = m3.Material(T=T, iso="D")
    for p in [0.2, 0.5, 1.0]:
        P("  equilibrium x(%.1f bar D2, 295 K) = %.3f" % (p, m3.x_of_f(p, T, "D")))
    xe = m3.x_of_f(0.5, T, "D")
    L = L_LOADED * 1e-6
    for xin in [0.90, 0.95]:
        J = (mat.Phi(xin) - mat.Phi(xe)) / L
        P("  F quadrant, x_in=%.2f, x_exit=%.3f clamped: J = %.2e D m-2 s-1 = %.0f mA cm-2 absorbed"
          % (xin, xe, J, J * 1.602e-19 / 10.0 * 1e3 / 1e3))
    q = 1.602e-19
    for i_mA, eta in [(300, 0.5), (300, 0.2), (100, 0.5), (20, 0.5)]:
        J = i_mA * 1e-3 * 1e4 / q * eta          # D m-2 s-1 absorbed (upper bound: all absorbed D permeates)
        xin = mat.x_of_Phi(mat.Phi(xe) + J * L)
        P("  F quadrant charge-limited: i=%d mA/cm2, eta=%.1f -> J=%.1e, x_in over F = %.3f"
          % (i_mA, eta, J, xin))
    P("  => over F/X/ED quadrants the whole foil (incl. the ENTRY face) sits at x ~ 0.64-0.66, not >= 0.9")
    P("     lateral coupling length (first Dirichlet/Neumann eigenmode): L/pi = %.1f um, 2L/pi = %.1f um"
      % (L_LOADED / np.pi, 2 * L_LOADED / np.pi))


# ---------------------------------------------------------------- 3. electrolyte recoils
def electrolyte_recoils(n=400_000):
    P("\n== 3. Cosmic fast-neutron recoils born in the ELECTROLYTE, crossing the foil (M5 PID + windows)")
    En, ftot = tel.cosmic_fast_n(n)
    wl = tel.p_window_low(CFG)
    out = {}
    nD2O = 2 * 1.104 / 20.03 * st.NA
    for lab, part, kmax, liquid, nden in [
        ("n-p, H2O twin (T-H)", "p", 1.0, "H2O", 2 * 0.997 / 18.015 * st.NA),
        ("n-p, D2O with H/D=0.5%", "p", 1.0, "D2O", 0.005 * nD2O),
        ("n-d, D2O (misID d->p)", "d", 8.0 / 9.0, "D2O", nD2O),
        ("D(n,np)n breakup p, D2O", "bu", None, "D2O", nD2O),
    ]:
        sig = tel.gammel_np(En)
        if part == "bu":
            # deuteron breakup: threshold 3.34 MeV; sigma ~0.2 b at 10 MeV, ~0.15 b at 20-50 MeV (ENDF/B-VIII
            # shape, from memory, +-50 %); proton energy ~U[0, 2/3 (En - 2.22)], isotropic. Order of magnitude.
            sig = np.where(En < 3.34, 0.0, np.where(En < 10, 0.2 * (En - 3.34) / 6.66, np.where(En < 30, 0.18, 0.13)))
            Er = RNG.uniform(0, 1, n) * np.maximum(En - 2.22, 0) * 2 / 3
            part = "p"
        else:
            Er = RNG.uniform(0, 1, n) * kmax * En
        zmax_um = 2000.0      # liquid layer considered (> range of every recoil that can cross the foil)
        zl = RNG.uniform(0, zmax_um, n)
        r = 10.0 * np.sqrt(RNG.uniform(0, 1, n))
        ph = RNG.uniform(0, 2 * np.pi, n)
        mu = RNG.uniform(0, 1, n)          # downward hemisphere
        phi = RNG.uniform(0, 2 * np.pi, n)
        tan = np.sqrt(1 - mu ** 2) / np.maximum(mu, 1e-9)
        d = 4.0 + L_LOADED * 1e-3
        x1 = r * np.cos(ph) + d * tan * np.cos(phi)
        y1 = r * np.sin(ph) + d * tan * np.sin(phi)
        x2 = r * np.cos(ph) + (d + CFG["gap_mm"]) * tan * np.cos(phi)
        y2 = r * np.sin(ph) + (d + CFG["gap_mm"]) * tan * np.sin(phi)
        in_dE = x1 ** 2 + y1 ** 2 <= CFG["coll_mm"] ** 2
        in_E = x2 ** 2 + y2 ** 2 <= CFG["bE_mm"] ** 2
        E = st.e_after(part, liquid, Er, zl / mu)
        E = st.e_after(part, "PdD0.9", E, L_LOADED / mu)
        E = st.e_after(part, "D2gas", E, GAS_UM_EQ / mu)
        E = np.where(in_dE, E, 0.0)
        ev = tel.deposit(part, E, mu, in_dE, in_E, CFG)
        ev.update(in_dE=in_dE, ang_ok=np.ones(n, bool))
        sel = tel.proton_cut(ev, CFG) & (ev["Etot"] > wl) & (ev["Etot"] < 3.1)
        pk = sel & (ev["Etot"] > 2.6)
        vol = np.pi * 1.0 ** 2 * zmax_um * 1e-4                    # cm3 of liquid over the 10 mm radius
        w = ftot * np.mean(sig) * 1e-24 * nden * vol               # recoils/s (isotropic)
        sw = sig / np.mean(sig)
        pid = w * 0.5 * np.mean(sel * sw) * 86400
        peak = w * 0.5 * np.mean(pk * sw) * 86400
        out[lab] = pid
        P("  %-26s PID window %.3g /day, peak window %.3g /day (unshielded cosmic flux %.1e cm-2 s-1 >1 MeV)"
          % (lab, pid, peak, ftot))
    # front-gas n-d recoils
    n2 = n
    En, _ = tel.cosmic_fast_n(n2)
    sig = tel.gammel_np(En)
    Er = RNG.uniform(0, 1, n2) * 8 / 9 * En
    zg = RNG.uniform(0, 4.0, n2)          # mm above the dE
    r = 12.8 * np.sqrt(RNG.uniform(0, 1, n2))
    ph = RNG.uniform(0, 2 * np.pi, n2)
    mu = RNG.uniform(0, 1, n2)
    phi = RNG.uniform(0, 2 * np.pi, n2)
    tan = np.sqrt(1 - mu ** 2) / np.maximum(mu, 1e-9)
    x1 = r * np.cos(ph) + zg * tan * np.cos(phi)
    y1 = r * np.sin(ph) + zg * tan * np.sin(phi)
    x2 = r * np.cos(ph) + (zg + CFG["gap_mm"]) * tan * np.cos(phi)
    y2 = r * np.sin(ph) + (zg + CFG["gap_mm"]) * tan * np.sin(phi)
    in_dE = x1 ** 2 + y1 ** 2 <= CFG["coll_mm"] ** 2
    in_E = x2 ** 2 + y2 ** 2 <= CFG["bE_mm"] ** 2
    E = np.where(in_dE, st.e_after("d", "D2gas", Er, zg * 1000 * 0.5 / mu), 0.0)
    ev = tel.deposit("d", E, mu, in_dE, in_E, CFG)
    ev.update(in_dE=in_dE, ang_ok=np.ones(n2, bool))
    sel = tel.proton_cut(ev, CFG) & (ev["Etot"] > wl) & (ev["Etot"] < 3.1)
    nD = 2 * 0.5e5 / (1.380649e-23 * 295) * 1e-6
    vol = np.pi * 1.28 ** 2 * 0.4
    w = tel.cosmic_fast_n(10)[1] * np.mean(sig) * 1e-24 * nD * vol
    P("  %-26s PID window %.3g /day" % ("n-d, 0.5 bar D2 front gas", w * 0.5 * np.mean(sel * sig / np.mean(sig)) * 86400))
    P("  M5 total telescope PID background (vacuum front, no electrolyte): 0.028 /day")
    P("  NOTE: shield factor for fast n inside 20 cm BPE + 5 cm Cu is ~0.3-1 (M5 uses 0.6); muon-spallation")
    P("        neutrons from the 5 cm Cu inner shield are NOT in the model.")
    return out


# ---------------------------------------------------------------- 4. helium
def helium():
    P("\n== 4. He-4 / He-3 background inventory (order of magnitude)")
    atoms_per_ccstp = 2.687e19
    p_he_atm = 5.24e-6
    barrer = 1e-10  # cm3STP cm /(cm2 s cmHg)
    dp_cmhg = p_he_atm * 76.0
    floor_front = 43.0            # He/s (M8: 0.16 nW)
    floor_1d = 0.39e-9 / (23.85e6 * 1.602e-19)
    P("  floors: 0.16 nW = %.0f He/s; 0.39 nW = %.0f He/s" % (floor_front, floor_1d))
    # PCTFE tube as air boundary (bore 20, wall 5 mm, height 90 mm)
    A = np.pi * 2.5 * 9.0
    for Pb in [1.0, 10.0, 30.0]:
        Q = Pb * barrer * A * dp_cmhg / 0.5 * atoms_per_ccstp
        P("  PCTFE wall %.0f cm2 x 5 mm, He permeability %4.0f Barrer -> %.1e He/s = %.1e x the 0.39 nW floor"
          % (A, Pb, Q, Q / floor_1d))
    # FFKM seal, air-exposed: M4 scaling from Viton 30 mm -> 2e8/s
    Q = 2e8 * 26.0 / 30.0
    P("  FFKM rim seal (Ø26), if air-exposed, scaled from M4 Viton 30 mm (2e8/s): %.1e He/s = %.1e x floor" % (Q, Q / floor_1d))
    # stored He in polymers (air-saturated), solubility 0.02 cm3STP/cm3/atm (M8 epoxy value as proxy)
    Vpoly = A * 0.5 + 1.0
    N = 0.02 * p_he_atm * Vpoly * atoms_per_ccstp
    P("  He stored in %.0f cm3 air-equilibrated polymer: %.1e atoms; outgassing over tau ~ 5 d -> ~%.0e He/s"
      % (Vpoly, N, N / (5 * 86400)))
    # cell -> front crossover through the FFKM seal, cell He partial pressure from that inventory
    n_head = 0.5e5 * 15e-6 / (1.380649e-23 * 295)
    for Ncell in [1e11, 1e13, 1e14]:
        frac = Ncell / n_head
        Q = 2e8 * 26 / 30 * frac / p_he_atm * (0.5 / 1.0)
        P("  cell He inventory %.0e -> x_He=%.1e in headspace; FFKM crossover to front ~%.1e He/s (%.1e x 43/s)"
          % (Ncell, frac, Q, Q / floor_front))
    # 6Li(n,alpha)t in electrolyte
    nLi6 = 1.0 * 0.012 * 0.0759 * st.NA
    for phi_th in [1e-3, 1e-2]:
        r = nLi6 * 940e-24 * phi_th
        P("  6Li(n,a)t in 12 mL 1 M LiOD at phi_th=%.0e: %.1e He-4/s (and same T/s)" % (phi_th, r))
    # tritium in D2O -> He-3
    m_d2o = 13.2e-3  # kg
    for a_bqkg in [1.5e2, 1e5, 1e7]:
        P("  D2O tritium %.1e Bq/kg -> He-3 ingrowth %.1e /s = %.1e /day (M8 assumes T/D 1.4e-15 = 1.5e2 Bq/kg)"
          % (a_bqkg, a_bqkg * m_d2o, a_bqkg * m_d2o * 86400))


def efficiency_vs_depth(n=300_000):
    P("\n== 1b. Telescope PID-window efficiency vs source depth, WITH the 4 mm / 0.5 bar D2 gap (M5 code)")
    cfg = dict(tel.BASE)
    cfg["overlayer"] = [("D2gas", GAS_UM_EQ)]
    cfg0 = dict(tel.BASE)
    cfg0["overlayer"] = []
    wl = tel.p_window_low(cfg)
    for z in [0.0, 4.0, L_LOADED / 2, 13.0, 15.0, L_LOADED, 17.66]:
        effs = []
        for c in (cfg0, cfg):
            ev = tel.telescope_events("p", 3.02, "delta", z, n, c)
            sel = tel.proton_cut(ev, c) & (ev["Etot"] > wl) & (ev["Etot"] < 3.1)
            effs.append(0.5 * sel.mean())
        P("  delta source at z=%5.2f um: eff(PID, vacuum) %.3f, eff(PID, 0.5 bar D2) %.3f" % (z, *effs))
    # Rn progeny plated on the entry face -> Delta-E-only events (no PID) in the t / 3He windows
    P("  Rn/Tn progeny alphas plated on the ENTRY face, through %.1f um (+ -2 um tolerance), Delta-E-only:" % L_LOADED)
    for Ea, lab in [(7.687, "214Po"), (8.785, "212Po"), (6.002, "218Po"), (6.05, "212Bi"), (5.304, "210Po")]:
        for Lz in [13.0 * 4.04 / 3.89, L_LOADED]:
            ev = tel.telescope_events("a", Ea, "delta", Lz, n, cfg)
            dEo = (ev["dE"] > cfg["thr"]) & (ev["E"] < cfg["thr"]) & ev["in_dE"]
            tw = dEo & (ev["dE"] > 0.85) & (ev["dE"] < 1.10)
            hw = dEo & (ev["dE"] > 0.55) & (ev["dE"] < 0.85)
            P("    %-6s L=%5.2f um: P(dE-only) %.3f, P(t window 0.85-1.10) %.4f, P(3He window 0.55-0.85) %.4f per decay"
              % (lab, Lz, 0.5 * dEo.mean(), 0.5 * tw.mean(), 0.5 * hw.mean()))


def main():
    tomography()
    efficiency_vs_depth()
    regimes()
    electrolyte_recoils()
    helium()
    with open(os.path.join(OUT, "rt1a_physics.txt"), "w") as f:
        f.write("\n".join(LINES) + "\n")


if __name__ == "__main__":
    main()
