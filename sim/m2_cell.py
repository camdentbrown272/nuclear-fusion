"""
M2 — cell-level budget for the recommended C1 (coax) and C3 (DFM) cells:
cell voltage, ohmic heat, electrolyte temperature-rise bound, bubble void and
coverage, recombiner gas flow and heat, charge/time needed to load, and the
operating window vs LiOD concentration and temperature.

Run:  python3 sim/m2_cell.py
Writes docs/models/figs/m2_cell.txt, m2_cell_voltage.png, m2_bubbles_recombiner.png
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from m2_common import (FIGS, F, R, T0, Vm_gas, SURFACES, Kinetics, Polarization, kappa_LiOD,
                       E_rev_D2O, E_tn_D2O, E_tn_H2O, E_rev_H2O, c_Pd, D_D_Pd, rho_D2O, mu_D2O)

OUT = []


def log(s=""):
    print(s)
    OUT.append(s)


POL = Polarization(Kinetics(**SURFACES["typical"]))
POLG = Polarization(Kinetics(**SURFACES["good"]))

# anode: O2 evolution on Pt in alkaline solution, Tafel form (approximate; +-0.15 V)
B_OER, I0_OER = 0.12, 1e-6        # V/dec, A/cm2  -> eta_a(100 mA/cm2) ~ 0.6 V
K_TH_D2O = 0.595                  # W/m/K thermal conductivity of D2O near 25 C
CP_D2O = 4.21                     # J/g/K
U_B = 1.0                         # cm/s bubble-swarm rise velocity (0.3-3)


def eta_anode(i_a):
    return B_OER * np.log10(np.maximum(i_a, 1e-9) / I0_OER)


def void(i, ub=U_B):
    """Void fraction in the D2 column above an upward-facing cathode at current density i."""
    jg = i * Vm_gas * 1e6 / (2 * F)          # cm/s superficial gas velocity
    return min(jg / (jg + ub), 0.6)


def bruggeman(eps):
    return (1 - eps) ** 1.5


def bubble_coverage(i):
    """Vogt & Balzer (2005) stagnant-electrolyte correlation, coverage ~ j^0.3.
    Prefactor 0.023 (j in A/m^2) as recalled; exponent confirmed by abstract."""
    return min(0.023 * (i * 1e4) ** 0.3, 0.9)


# ------------------------------------------------------------ recommended cells
def c1_cell(i, c=0.1, T_C=25.0, a=0.05, L=3.0, Ra=1.0, anode_area_ratio=1.0, bubbles=True):
    """Coax wire between PTFE end plates (1-D, uniform)."""
    kap = kappa_LiOD(c, T_C)
    A_c = 2 * np.pi * a * L
    I = i * A_c
    R_ohm = np.log(Ra / a) / (2 * np.pi * kap * L)
    if bubbles:   # plume near wire: effective 5 % increase at 100 mA/cm2 (from m2_current A2b), scale ~ eps
        R_ohm *= 1 + 0.5 * void(i) * 1.0
    i_a = I / (2 * np.pi * Ra * L * anode_area_ratio)
    V = E_rev_D2O + abs(POL.eta_of_i(i)) + eta_anode(i_a) + I * R_ohm
    return dict(I=I, A=A_c, R=R_ohm, V=V, iR=I * R_ohm, eta_c=POL.eta_of_i(i), eta_a=eta_anode(i_a),
                P=I * V, Q_open=I * (V - E_tn_D2O), Q_rec=I * E_tn_D2O, kap=kap)


def c3_cell(i, c=0.1, T_C=25.0, Rd=1.0, g=1.5, anode="ring", bubbles=True):
    """DFM tube cell (bore = wetted disk radius Rd), ring anode at gap g (ring
    adds ~10 % constriction resistance, from m2_current C6/C7) or mesh at g."""
    kap = kappa_LiOD(c, T_C)
    A_c = np.pi * Rd**2
    I = i * A_c
    eps = void(i) if bubbles else 0.0
    R_ohm = g / (kap * bruggeman(eps) * A_c) * (1.10 if anode == "ring" else 1.0)
    A_an = (2 * np.pi * (Rd - 0.1) * np.pi * 0.05) if anode == "ring" else A_c   # 0.5 mm Pt wire loop / mesh
    i_a = I / A_an
    V = E_rev_D2O + abs(POL.eta_of_i(i)) + eta_anode(i_a) + I * R_ohm
    # 1-D conduction upper bound on electrolyte / membrane temperature rise:
    # uniform Joule heating q''' = i^2/kappa_eff in the gap, cathode (vacuum-backed) adiabatic,
    # anode plane held at bulk temperature; plus cathode interfacial heat i*(|eta_c|) at z=0.
    i_SI = i * 1e4
    kap_SI = kap * bruggeman(eps) * 100
    g_SI = g * 1e-2
    q3 = i_SI**2 / kap_SI
    # interfacial heat at the cathode: i*(|eta_c| - (E_tn - E_rev)) -> the part not stored as enthalpy
    q_int = i_SI * max(abs(POL.eta_of_i(i)) - (E_tn_D2O - E_rev_D2O), 0.0)
    dT = q3 * g_SI**2 / (2 * K_TH_D2O) + q_int * g_SI / K_TH_D2O
    return dict(I=I, A=A_c, R=R_ohm, V=V, iR=I * R_ohm, eta_c=POL.eta_of_i(i), eta_a=eta_anode(i_a), eps=eps,
                P=I * V, Q_open=I * (V - E_tn_D2O), Q_rec=I * E_tn_D2O, dT_cond=dT, kap=kap, i_a=i_a)


def thermo():
    log("=" * 78)
    log("THERMOCHEMISTRY")
    log("=" * 78)
    log(f"  E_rev(D2O) = {E_rev_D2O:.4f} V   (Delta_f G D2O(l) = -243.44 kJ/mol)")
    log(f"  E_tn (D2O) = {E_tn_D2O:.4f} V   (Delta_f H D2O(l) = -294.60 kJ/mol)")
    log(f"  E_rev(H2O) = {E_rev_H2O:.4f} V ; E_tn(H2O) = {E_tn_H2O:.4f} V (light-water control)")
    log(f"  Recombination heat at the recombiner = I * E_tn: {E_tn_D2O:.3f} W/A (D2O), {E_tn_H2O:.3f} W/A (H2O)")
    log(f"  => D2O vs H2O control differ by {1e3*(E_tn_D2O-E_tn_H2O):.1f} mW/A in open-cell enthalpy bookkeeping")
    for I in [0.1, 0.3, 1.0, 2.0]:
        nD2 = I / (2 * F)
        log(f"  I={I:4.1f} A: D2 {nD2*Vm_gas*1e6*60:6.2f} mL/min + O2 {0.5*nD2*Vm_gas*1e6*60:6.2f} mL/min (25 C, 1 atm);"
            f" recombiner heat {I*E_tn_D2O:.2f} W; D2O consumed if vented {nD2*20.03*3600:.2f} g/h")
    # headspace hazard
    for Vh in [5, 20, 50, 200]:
        n = Vh * 1e-6 / Vm_gas
        E = (2 / 3) * n * 249e3            # 2D2 + O2 -> 2 D2O(g), ~249 kJ/mol D2O(g)
        log(f"  stoichiometric D2/O2 headspace {Vh:4d} mL -> stored chemical energy {E:6.0f} J")


def c1_section():
    log("")
    log("=" * 78)
    log("C1 COAX (a=0.5 mm, L=30 mm, Ra=10 mm, between PTFE end plates): voltage and heat")
    log("=" * 78)
    for c in [0.1, 0.5, 1.0]:
        log(f"  LiOD {c} M (kappa={kappa_LiOD(c)*1e3:.1f} mS/cm at 25 C)")
        for i in [0.01, 0.05, 0.1, 0.2, 0.3, 0.5]:
            d = c1_cell(i, c)
            log(f"    i={i*1e3:4.0f} mA/cm2 I={d['I']:.3f} A: eta_c={d['eta_c']:+.2f} eta_a={d['eta_a']:.2f} iR={d['iR']:5.2f}"
                f" V_cell={d['V']:5.2f} V  P_in={d['P']:.2f} W  heat(open)={d['Q_open']:.2f} W  recomb={d['Q_rec']:.2f} W")
    # loading time
    a = 0.05
    for x in [0.9]:
        Q = x * c_Pd * np.pi * a**2 * F           # C per cm of wire
        log(f"  Charge to load 1 cm of wire to x={x}: {Q:.1f} C/cm -> at 10 mA/cm2 and 100 % efficiency {Q/(0.01*2*np.pi*a)/3600:.1f} h"
            f"; diffusion time a^2/D = {a**2/D_D_Pd/3600:.1f} h (90 % equilibration ~0.5 a^2/D = {0.5*a**2/D_D_Pd/3600:.1f} h)")
    for a_ in [0.0125, 0.025, 0.05, 0.1]:
        log(f"    wire diameter {20*a_:.2f} mm: 0.5 a^2/D = {0.5*a_**2/D_D_Pd/3600:5.2f} h ; I at 300 mA/cm2 over 30 mm = {0.3*2*np.pi*a_*3:.3f} A")


def c3_section():
    log("")
    log("=" * 78)
    log("C3 DFM (Rd=10 mm tube, ring anode g=15 mm unless noted): voltage, heat, temperature-rise bound")
    log("=" * 78)
    rows = {}
    for c in [0.1, 0.5, 1.0]:
        for anode, g in [("ring", 1.5), ("mesh", 0.5), ("mesh", 0.3)]:
            log(f"  LiOD {c} M, {anode} anode g={10*g:.0f} mm")
            for i in [0.01, 0.05, 0.1, 0.2, 0.3, 0.5]:
                d = c3_cell(i, c, g=g, anode=anode)
                rows[(c, anode, g, i)] = d
                log(f"    i={i*1e3:4.0f} mA/cm2 I={d['I']:.3f} A: eps={d['eps']:.3f} iR={d['iR']:5.2f} V_cell={d['V']:5.2f} V"
                    f" P_in={d['P']:6.2f} W heat(open)={d['Q_open']:6.2f} W  dT_cond(bound)={d['dT_cond']:6.1f} K"
                    f"  i_anode={d['i_a']:.2f} A/cm2")
    Q = 0.9 * c_Pd * 0.005 * F
    log(f"  Charge to load a 50 um membrane to x=0.9: {Q:.1f} C/cm2 -> {Q/0.02/60:.0f} min at 20 mA/cm2 (100 % eff.);"
        f" diffusion time L^2/D = {0.005**2/D_D_Pd:.0f} s")
    return rows


def window_section():
    log("")
    log("=" * 78)
    log("OPERATING WINDOW vs temperature (0.1 M and 1.0 M LiOD)")
    log("=" * 78)
    for T_C in [5, 10, 20, 25, 40, 60]:
        log(f"  T={T_C:3d} C: kappa(0.1 M)={kappa_LiOD(0.1, T_C)*1e3:5.1f} mS/cm, kappa(1.0 M)={kappa_LiOD(1.0, T_C)*1e3:5.1f} mS/cm;"
            f" C3 ring g=15mm @200 mA/cm2 V_cell(0.1 M)={c3_cell(0.2, 0.1, T_C)['V']:.1f} V, V_cell(1 M)={c3_cell(0.2, 1.0, T_C)['V']:.1f} V")


def bubbles_section():
    log("")
    log("=" * 78)
    log("BUBBLES: void fraction in the D2 column (u_b sensitivity) and adherent-bubble coverage")
    log("=" * 78)
    for i in [0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0]:
        e = [void(i, ub) for ub in [0.3, 1.0, 3.0]]
        th = bubble_coverage(i)
        log(f"  i={i*1e3:5.0f} mA/cm2: eps(u_b=0.3/1/3 cm/s) = {e[0]:.3f}/{e[1]:.3f}/{e[2]:.3f}; kappa factor {bruggeman(e[0]):.2f}/{bruggeman(e[1]):.2f}/{bruggeman(e[2]):.2f};"
            f" coverage Theta={th:.2f} -> local i on free area x{1/(1-th):.2f}, Delta x (good) = {float(POLG.x_of_i(i/(1-th))-POLG.x_of_i(i)):+.4f}")
    # Stokes rise velocities
    for d in [20e-6, 50e-6, 100e-6, 300e-6]:
        u = (rho_D2O - 0.2) * 9.81 * d**2 / (18 * mu_D2O)
        log(f"  Stokes rise velocity of a {d*1e6:.0f} um D2 bubble: {u*100:.2f} cm/s")


def figures(rows):
    ii = np.array([0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5])
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))
    for c, col in [(0.1, "C0"), (0.5, "C1"), (1.0, "C2")]:
        ax[0].plot(ii * 1e3, [c1_cell(i, c)["V"] for i in ii], "-o", color=col, ms=3, label=f"C1 coax, {c} M")
        ax[0].plot(ii * 1e3, [c3_cell(i, c)["V"] for i in ii], "--s", color=col, ms=3, label=f"C3 DFM ring g=15, {c} M")
        ax[1].semilogy(ii * 1e3, [c3_cell(i, c)["dT_cond"] for i in ii], "--s", color=col, ms=3, label=f"C3 ring g=15, {c} M")
        ax[1].semilogy(ii * 1e3, [c3_cell(i, c, g=0.3, anode='mesh')["dT_cond"] for i in ii], ":", color=col, label=f"C3 mesh g=3, {c} M")
    ax[0].set_xlabel("cathode current density (mA/cm²)"); ax[0].set_ylabel("cell voltage (V)"); ax[0].legend(fontsize=7); ax[0].grid(alpha=.3)
    ax[0].set_title("Required cell voltage (25 C)")
    ax[1].axhline(10, color="k", lw=.6); ax[1].set_xlabel("i (mA/cm²)"); ax[1].set_ylabel("ΔT bound (K)")
    ax[1].set_title("DFM electrolyte ΔT, conduction-only upper bound"); ax[1].legend(fontsize=7); ax[1].grid(alpha=.3)
    for I_, lab in [(1, "1 A")]:
        pass
    II = np.linspace(0, 2, 50)
    ax[2].plot(II, II * E_tn_D2O, label="recombiner heat D2O (I·1.527 V)")
    ax[2].plot(II, II * E_tn_H2O, "--", label="recombiner heat H2O (I·1.481 V)")
    ax2 = ax[2].twinx()
    ax2.plot(II, II / (2 * F) * Vm_gas * 1e6 * 60 * 1.5, "k:", label="D2+O2 gas flow")
    ax2.set_ylabel("D2+O2 flow (mL/min)")
    ax[2].set_xlabel("cell current (A)"); ax[2].set_ylabel("heat (W)"); ax[2].legend(fontsize=7, loc="upper left"); ax[2].grid(alpha=.3)
    ax[2].set_title("Internal recombiner load")
    fig.tight_layout(); fig.savefig(os.path.join(FIGS, "m2_cell_voltage.png"), dpi=130); plt.close(fig)


if __name__ == "__main__":
    thermo()
    c1_section()
    rows = c3_section()
    window_section()
    bubbles_section()
    figures(rows)
    with open(os.path.join(FIGS, "m2_cell.txt"), "w") as fh:
        fh.write("\n".join(OUT) + "\n")
