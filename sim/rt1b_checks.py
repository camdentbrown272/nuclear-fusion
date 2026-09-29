"""
Red-team 1B (engineering/safety) quantitative checks on iteration-1 rev B.

Output: docs/design/figs/rt1b_checks.txt, rt1b_vdp_shunt.png, rt1b_gas_balance.png
Sections
  A  closed-cell gas and D2O balance (loading, permeation, recombiner)
  B  front pressurisation when the Pd-Ag exhaust stops
  C  thermal: V_cell with bubbles, liner conduction, house heat
  D  mechanics: Mo grid, lifted membrane, corrugation, Si dE window, transformation plasticity
  E  FX anodic half-cycle charge vs foil inventory
  F  He leak budget vs the 43-100 He/s floor; make-up D2 He load
  G  Si telescope: dE capacitance/noise, Paschen margin in D2
  H  rim van der Pauw shunted by the bonded Pt annulus (2-D finite differences)
All inputs are the rev-B values (iteration-1.md §3.2, §9) unless tagged [est.].
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

sys.path.insert(0, os.path.dirname(__file__))
import m3_common as m3  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), "..", "docs", "design", "figs")
os.makedirs(OUT, exist_ok=True)
lines = []


def p(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    lines.append(s)


F = 96485.0
R = 8.314
T0 = 295.0
kB = 1.380649e-23
A_cath = np.pi * 1.0 ** 2          # cm2, active dia 20 mm
V_head = 15e-6                      # m3
V_front = 30e-6                     # m3
P_fill = 0.5e5                      # Pa
t_pd = 12e-4                        # cm
n_Pd = 0.113                        # mol / cm3

# ------------------------------------------------------------------ A
p("== A. Closed-cell gas balance ==")
n_head = P_fill * V_head / (R * T0)
p(f"headspace D2 inventory: {n_head*1e6:.0f} umol")
V_pd = A_cath * t_pd * 1.10          # +10 % for the rim out to the bond
nD_load = V_pd * n_Pd * 0.90         # mol D atoms
O2_excess = nD_load / 4
D2_used = 2 * O2_excess
p(f"D absorbed to x=0.9: {nD_load*1e6:.0f} umol D -> excess O2 {O2_excess*1e6:.0f} umol, "
  f"recombiner burns {D2_used*1e6:.0f} umol D2 = {100*D2_used/n_head:.0f} % of headspace D2")
p(f"  -> cell pressure after loading without make-up: {P_fill*(n_head-D2_used)/n_head/100:.0f} mbar "
  f"(front stays 500): dP = {-(P_fill*D2_used/n_head)/100:.0f} mbar vs interlock -40 / relief foil -20")
res = {}
for name, j in (("H-M 10", 10e-3), ("FX 50", 50e-3), ("FX 150", 150e-3)):
    I = j * A_cath
    o2 = I / (4 * F)                   # mol/s excess O2 (permeated D never returns as D2)
    dD2 = 2 * o2
    t_exh = n_head / dD2
    dpdt = 3 * o2 * R * T0 / V_head / 100   # mbar/s while D2 lasts (2D2+O2 -> liquid)
    d2o_g_day = I / (2 * F) * 86400 * 20.03
    makeup_Lday = dD2 * 86400 * 22.4
    res[name] = (I, t_exh, dpdt, d2o_g_day, makeup_Lday)
    p(f"{name} mA/cm2 permeating: I_perm={I:.3f} A, headspace D2 gone in {t_exh:.0f} s "
      f"({t_exh/60:.1f} min), dP/dt={-dpdt:.2f} mbar/s, then O2 accumulates at "
      f"{o2*1e6:.2f} umol/s; D2O loss {d2o_g_day:.2f} g/day ({d2o_g_day*42:.0f} g over P2 vs 13 g); "
      f"make-up D2 needed {makeup_Lday:.2f} L(STP)/day")
# recombiner heat
for j in (0.2, 0.3, 0.5):
    p(f"recombiner heat at {j*1e3:.0f} mA/cm2: {1.527*j*A_cath:.2f} W")
gas_rate = 0.3 * A_cath / (2 * F) * 1.5 * R * T0 / P_fill * 1e6 * 60  # mL/min at 0.5 bar
p(f"gas generation at 300 mA/cm2 and 0.5 bar: {gas_rate:.1f} mL/min (headspace {V_head*1e6:.0f} mL "
  f"turned over every {V_head*1e6/gas_rate*60:.0f} s); 1 % recombiner lag -> "
  f"{0.01*gas_rate/60/ (V_head*1e6) *500:.3f} mbar/s")

# ------------------------------------------------------------------ B
p("\n== B. Front pressurisation with the Pd-Ag exhaust stopped ==")
for name, j in (("H-M 10", 10e-3), ("keep-alive 20 (FX)", 20e-3), ("FX 50", 50e-3), ("FX 150", 150e-3)):
    I = j * A_cath
    nd2 = I / (2 * F)
    dpdt = nd2 * R * T0 / V_front / 100
    p(f"{name}: front +{dpdt:.2f} mbar/s -> dP from +30 to -20 (relief foil bursts) in {50/dpdt:.0f} s")
# Pd-Ag area needed
Perm = 1.0e-8 / np.sqrt(2)   # mol/(m s Pa^0.5), Pd-23Ag 350 C, D isotope /sqrt2 [est.]
L = 75e-6
flux = Perm / L * np.sqrt(P_fill)
p(f"Pd-Ag (75 um wall, 350 C) flux at 0.5 bar -> vacuum: {flux:.3f} mol/m2/s; area for FX 150: "
  f"{0.15*A_cath/(2*F)/flux*1e4:.2f} cm2")

# ------------------------------------------------------------------ C
p("\n== C. Thermal ==")
kappa = 0.15      # S/cm, 1.0 M LiOD/D2O at 22 C [est.; LiOH/H2O 0.18]
g = 0.40          # cm
h_liq = 12.0 / A_cath   # cm
A_liner = np.pi * 2.0 * h_liq
G_liner = 0.25 * A_liner * 1e-4 / 1e-3
p(f"electrolyte column {h_liq*10:.0f} mm in a 20 mm bore; wetted liner {A_liner:.0f} cm2; "
  f"G_liner = {G_liner:.2f} W/K (PTFE k=0.25, 1 mm)")
rows = []
for j in (0.2, 0.3, 0.5):
    ug = j / (2 * F) * R * T0 / P_fill * 1e4 * 1e3   # mm/s superficial D2 at 0.5 bar
    for d_b in (50e-6, 100e-6):
        vr = 2 / 9 * 1100 * 9.81 * (d_b / 2) ** 2 / 1.25e-3 * 1e3   # mm/s Stokes
        eps = min(ug / vr, 0.6)
        Rohm = g / kappa * (1 - eps) ** -1.5
        V = 1.26 + 0.35 + 0.40 + j * Rohm     # E_rev(D2O) + eta_c + eta_a [est.] + IR
        Q = V * j * A_cath                     # all heat stays in the closed cell
        dT = Q / G_liner
        rows.append((j, d_b, eps, V, Q, dT))
        p(f"j={j*1e3:.0f} mA/cm2 bubbles {d_b*1e6:.0f} um: u_g={ug:.2f} mm/s eps={eps:.2f} "
          f"V_cell={V:.2f} V Q={Q:.1f} W  dT(liner)={dT:.1f} K  tau={60/G_liner:.0f} s")
p("free convection off an unjacketed 316L body (0.02 m2, h=5): G = 0.10 W/K -> "
  "dT = 36-80 K at 3.6-8 W: body must be liquid-jacketed")
Q_house = 8 * 3.6 + 8 * 20 + 40
G_house = 0.4 * 1.8 / 0.20
p(f"house: {Q_house:.0f} W (8 cells P2, 8 Pd-Ag heaters 20 W [est.], electronics 40 W) through "
  f"20 cm HDPE G={G_house:.1f} W/K -> +{Q_house/G_house:.0f} K without internal water cooling")
dvp = 1.4   # mbar/K D2O vapour near 22 C [est.]
p(f"cell dP sensitivity: gas 500/295={500/295:.2f} mbar/K + vapour {dvp} = {500/295+dvp:.1f} mbar/K; "
  f"a 500 mA/cm2 step (+~10 K) moves dP by ~{10*(500/295+dvp):.0f} mbar (interlock band 40)")

# ------------------------------------------------------------------ D
p("\n== D. Mechanics ==")
eta = 0.10 / 1.10            # ligament efficiency of 1.0 mm hex holes, 0.10 webs
for dP, lab in ((0.03, "normal +30 mbar"), (0.05, "+50 mbar"), (0.5, "loss-of-front 0.5 bar"),
                (1.0, "loss-of-front 1.0 bar (P4)")):
    s5 = 0.75 * dP * 0.1 * 2.5 ** 2 / 0.10 ** 2 / eta
    s10 = 0.75 * dP * 0.1 * 5.0 ** 2 / 0.10 ** 2 / eta
    p(f"Mo grid 0.10 mm, {lab}: sigma = {s5:.0f} MPa if ribs make 5 mm cells, "
      f"{s10:.0f} MPa if the only support is the cross septum (10 mm cells)")
p("  Mo sheet yield 550-700 MPa stress-relieved, 350-450 recrystallised; DBTT ~ RT recrystallised [est.]")
# lifted membrane (Hencky) scaled from ADR-005 54 MPa at 50 mbar
for dP in (0.02, 0.03, 0.05, 0.5):
    p(f"lifted 12 um membrane at reverse {dP*1e3:.0f} mbar: {54*(dP/0.05)**(2/3):.0f} MPa "
      f"(annealed yield 40, UTS 170)")
p(f"P3 FX front step to 1.0 bar with cell at 0.5 bar (unbalanced): {54*(0.5/0.05)**(2/3):.0f} MPa -> rupture")
# Si dE diaphragm (Hencky), a = 12 mm, h = 25 um
for dP in (1e3, 1e4, 5e4):
    s = 0.423 * (130e9 * dP ** 2 * 0.012 ** 2 / 25e-6 ** 2) ** (1 / 3) / 1e6
    p(f"25 um Si dE (a=12 mm) under {dP/100:.0f} mbar across it: {s:.0f} MPa (thin-wafer fracture ~50-150)")
# corrugation
p(f"Pt corrugation guided-leg strain for 0.45 mm radial travel (t=0.10, leg 3 mm): "
  f"{3*0.10*0.225/3**2*100:.2f} % (annealed Pt yield ~0.03 %): plastic on first load, "
  f"ratchets on each full load/deload")
for dP in (0.5, 1.0):
    p(f"corrugation as 3 mm diaphragm strip at {dP} bar: {dP*0.1*(3/0.1)**2/2:.0f} MPa (annealed Pt yield 35-50)")
# transformation plasticity (Greenwood-Johnson) during alpha<->beta under load
for sig, sy, lab in ((6, 40, "P1 transit, +30 mbar"), (28, 40, "P6a/pump-down 0.5 bar, annealed"),
                     (28, 150, "P6a 0.5 bar, hardened"), (44, 150, "P6a after P4 at 1.0 bar")):
    p(f"transformation plasticity {lab}: eps_tp ~ 5/6*0.10*{sig}/{sy} = {5/6*0.10*sig/sy*100:.1f} %")
# CTE mismatch at the bond
p(f"bond cool-down 850->22 C: (11.8-8.8)e-6*828 = {(11.8-8.8)*828e-4:.2f} % tensile misfit in Pd (yields, flat)")
p(f"hydride misfit on bonded rim: ~3.3 % linear at x=0.65-0.9 -> {121e3*0.033/0.7:.0f} MPa elastic "
  "equivalent vs 40-150 yield: the bond line is plastically cycled once per full load/deload")

# ------------------------------------------------------------------ E
p("\n== E. FX anodic half-cycle ==")
inv = t_pd * n_Pd * 0.65 * F   # C/cm2
q_an = 0.100 * 120              # 240 s period, square wave: 120 s anodic
p(f"FX foil D inventory at x=0.65: {inv:.1f} C/cm2; anodic half-cycle extracts {q_an:.0f} C/cm2 "
  f"({q_an/inv:.1f}x inventory) unless the front resupplies >= 100 mA/cm2")
p(f"x at 0.2 bar front, 22 C: {m3.x_of_f(0.2,295):.3f}; plateau {m3.p_plateau(295)*1e3:.0f} mbar; "
  f"at 60 C plateau {m3.p_plateau(333)*1e3:.0f} mbar")
p(f"P3 cycles: 14 d / 240 s = {14*86400/240:.0f} (M3 perforation life if each is an alpha/beta transit: 12-180)")

# ------------------------------------------------------------------ F
p("\n== F. He budget ==")
N_mbarL = 2.455e19
xHe_air = 5.24e-6
for L in (1e-10, 1e-12, 1e-13):
    p(f"total std-He leak {L:.0e} mbar L/s from air: {L*xHe_air*N_mbarL:.1f} He/s (floor 43-100 He/s)")
p(f"max total leak for <=10 % of 43 He/s: {4.3/(xHe_air*N_mbarL):.1e} mbar L/s std-He "
  "(below the ~5e-12 limit of a He leak detector; Ar-tracer route only)")
for name in ("H-M 10", "FX 50", "FX 150"):
    Ld = res[name][4]
    molec = Ld * 1e3 * 2.687e19
    for xHe in (1e-9, 1e-15):
        p(f"make-up D2 for {name}: {Ld:.2f} L/day at He fraction {xHe:.0e} -> {molec*xHe/86400:.1e} He/s")
V_manifold_spike = 1e-6  # mbar He partial in manifold after a spike [est.]
p(f"closed all-metal valve seat, 1e-10 std-He conductance, manifold He 1e-6 mbar: "
  f"{1e-10*1e-6*N_mbarL:.1f} He/s -> guard-vacuum interspace required on every cell valve")

# ------------------------------------------------------------------ G
p("\n== G. Si telescope ==")
eps0 = 1.04e-12  # F/cm Si
C_dE = eps0 * 1.5 / 25e-4
p(f"dE quadrant capacitance (150 mm2, 25 um): {C_dE*1e12:.0f} pF; series ENC ~3-5 e/pF -> "
  f"{C_dE*1e12*3*3.62e-3:.1f}-{C_dE*1e12*5*3.62e-3:.1f} keV rms")
for T in (22, 35, 60):
    I = 5e-9 * 2 ** ((T - 22) / 7.5)
    enc = np.sqrt(I * 1e-6 / 1.6e-19)
    p(f"leakage {I*1e9:.0f} nA/quadrant at {T} C [5 nA at 22 C est.]: shot ENC {enc:.0f} e = {enc*3.62e-3:.1f} keV rms")
A_p, B_p, gam = 5.0, 130.0, 0.01   # H2 Paschen constants (Torr cm), used for D2
x_min = np.e * np.log(1 + 1 / gam) / A_p
V_min = B_p * x_min / (np.log(A_p * x_min) - np.log(np.log(1 + 1 / gam)))
p(f"Paschen minimum (H2 constants, applied to D2): pd = {x_min:.1f} Torr cm, V_min = {V_min:.0f} V "
  f"(literature H2 ~270-300 V); 150 V bias -> margin {V_min/150:.1f}x, eroded by edge fields and Kr (Penning)")

# ------------------------------------------------------------------ H
p("\n== H. Rim vdP with the Pt annulus bonded under the Pd ==")


def vdp(rs_pd, rs_pt, n=161, R_out=17.0, r_bond=10.25, r_c=10.5):
    h = 2 * R_out / (n - 1)
    xs = np.linspace(-R_out, R_out, n)
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    r = np.hypot(X, Y)
    inside = r <= R_out
    # sheet conductance: Pd disk everywhere r<=r_Pd_edge(=12), Pt where r>=r_bond
    g = np.zeros_like(r)
    g[r <= 12.0] += 1 / rs_pd
    g[(r >= r_bond) & inside] += 1 / rs_pt
    idx = -np.ones(r.shape, int)
    idx[inside] = np.arange(inside.sum())
    N = inside.sum()
    A = lil_matrix((N, N))
    for i in range(n):
        for j in range(n):
            k = idx[i, j]
            if k < 0:
                continue
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                ii, jj = i + di, j + dj
                if 0 <= ii < n and 0 <= jj < n and idx[ii, jj] >= 0:
                    gg = 2 / (1 / max(g[i, j], 1e-12) + 1 / max(g[ii, jj], 1e-12))
                    A[k, k] += gg
                    A[k, idx[ii, jj]] -= gg

    def node(ang):
        return idx[np.argmin(np.abs(xs - r_c * np.cos(ang))), np.argmin(np.abs(xs - r_c * np.sin(ang)))]
    c = [node(a) for a in (0, np.pi / 2, np.pi, 3 * np.pi / 2)]
    A[c[3], :] = 0
    A[c[3], c[3]] = 1          # ground
    A = A.tocsr()
    b = np.zeros(N)
    b[c[0]] = 1.0
    b[c[1]] = -1.0
    b[c[3]] = 0
    v = spsolve(A.tocsc(), b)
    return (v[c[3]] - v[c[2]])   # V(D)-V(C) per unit current from A to B


rs_pt = 10.6e-8 / 1.0e-4                 # ohm/sq, 0.10 mm Pt
res_h = []
for rho_rel in (1.0, 1.4, 1.8):
    rs_pd = 10.5e-8 * rho_rel / 12e-6
    r_with = vdp(rs_pd, rs_pt)
    r_no = vdp(rs_pd, 1e12)
    res_h.append((rho_rel, r_with, r_no))
    p(f"rho_Pd/rho0={rho_rel}: R_vdP with Pt annulus {r_with*1e3:.4f} mOhm, Pd-only {r_no*1e3:.3f} mOhm")
s_with = (res_h[-1][1] / res_h[0][1] - 1)
s_no = (res_h[-1][2] / res_h[0][2] - 1)
p(f"R change for rho x1.8 (x 0 -> ~0.9): with annulus {s_with*100:.1f} %, Pd-only {s_no*100:.1f} %; "
  f"Pt TCR 0.39 %/K -> 1 K drift = {0.39/ (s_with*100) *100:.0f} % of the full loading signal")

fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
js = sorted(set(r[0] for r in rows))
for d_b, ls in ((50e-6, "--"), (100e-6, "-")):
    ax[0].plot([j * 1e3 for j in js], [r[5] for r in rows if r[1] == d_b], ls, marker="o",
               label=f"bubbles {d_b*1e6:.0f} µm")
ax[0].axhline(1, color="grey", lw=0.8)
ax[0].set_xlabel("cathodic current density (mA cm$^{-2}$)")
ax[0].set_ylabel("electrolyte − jacket ΔT across 1 mm PTFE (K)")
ax[0].set_title("Lined-cell thermal rise (jacketed body)")
ax[0].legend(frameon=False)
tt = np.linspace(0, 3600, 400)
for name in ("H-M 10", "FX 50", "FX 150"):
    I, t_exh, dpdt, _, _ = res[name]
    ax[1].plot(tt / 60, np.maximum(500 - dpdt * tt, 0), label=name + " mA cm$^{-2}$")
ax[1].axhline(500 - 40 - 30, color="r", lw=0.8, ls=":")
ax[1].text(1, 440, "interlock (ΔP = −40 mbar)", color="r", fontsize=8)
ax[1].set_xlabel("time without D$_2$ make-up (min)")
ax[1].set_ylabel("cell headspace D$_2$ (mbar)")
ax[1].set_title("Permeating cells consume headspace D$_2$")
ax[1].legend(frameon=False)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "rt1b_checks.png"), dpi=130)
with open(os.path.join(OUT, "rt1b_checks.txt"), "w") as f:
    f.write("\n".join(lines) + "\n")
