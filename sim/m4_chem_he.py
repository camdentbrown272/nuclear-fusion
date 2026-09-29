"""
M4 — electrochemical energy bookkeeping (recombination, loading, permeation)
and the 4He co-measurement budget for a sealed cell.

Run: python3 sim/m4_chem_he.py
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from m4_params import (F, R_GAS, K_B, E_TN_D2O, E_TN_H2O, DH_ABS_D, DHF_D2O, DH_VAP_D2O,
                       Q_DD_HE4, HE_AIR_FRAC, K_HE_PYREX_25C, EA_HE_PYREX, K_HE_VITON,
                       KH_HE_WATER, N_A)

OUT = os.path.join(os.path.dirname(__file__), "..", "docs", "models", "figs")
LOSCHMIDT = 2.6868e19   # molecules per cm^3 at STP


def p_vap_water(Tc):
    """Buck (1996) saturation vapour pressure of H2O, Pa."""
    return 611.21 * np.exp((18.678 - Tc / 234.5) * (Tc / (257.14 + Tc)))


def p_vap_d2o(Tc):
    # D2O vapour pressure is ~12 % below H2O at 25 C and ~7 % below at 80 C
    # (Hill & MacMillan, Ind. Eng. Chem. Res. 1988); linear interpolation of the ratio.
    ratio = 0.88 + (0.93 - 0.88) * (Tc - 25) / 55
    return p_vap_water(Tc) * ratio


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    t = ["M4 electrochemical energy bookkeeping and 4He budget", "=" * 60]
    t.append(f"Thermoneutral potentials: E_tn(D2O) = {E_TN_D2O:.4f} V, E_tn(H2O) = {E_TN_H2O:.4f} V, "
             f"difference = {1e3*(E_TN_D2O-E_TN_H2O):.1f} mV")
    t.append("")
    t.append("1. Recombination heat released at the recombiner of a closed cell, P_rec = E_tn * I")
    t.append("   I (A)   D2O (W)   H2O (W)   error if H2O control analysed with D2O E_tn (mW)")
    for I in [0.1, 0.2, 0.5, 1.0, 2.0]:
        t.append(f"   {I:4.1f}    {E_TN_D2O*I:6.3f}    {E_TN_H2O*I:6.3f}    {1e3*(E_TN_D2O-E_TN_H2O)*I:6.1f}")

    t.append("")
    t.append("2. OPEN cell: P_heat = (V - E_tn) I ; unaccounted recombination fraction phi gives")
    t.append("   apparent excess = phi * E_tn * I ;  vapour carried by the vent gas removes heat")
    t.append("   T(C)  phi=0.01 @1A (mW)   vapour loss @1A (mW, D2O, 1 atm)")
    for Tc in [25, 40, 60, 80, 95]:
        pv = p_vap_d2o(Tc)
        ngas = 0.75 * 1.0 / F                      # mol/s of D2 + O2/2 at 1 A
        Pvap = ngas * pv / (101325 - pv) * DH_VAP_D2O
        t.append(f"   {Tc:4d}   {1e3*0.01*E_TN_D2O:8.1f}             {1e3*Pvap:8.2f}")

    t.append("")
    t.append("3. CLOSED cell: recombination inefficiency (1-eta) -> heat deficit and pressure rise")
    Vh = 33e-6       # m^3 free headspace of the reference cell (30 mm x 40 mm dia minus basket)
    T = 308.0
    for I in [0.2, 0.5, 1.0, 2.0]:
        ngas = 0.75 * I / F
        for eta_def in [1e-2, 1e-3, 1e-4]:
            dpdt = eta_def * ngas * R_GAS * T / Vh
            t.append(f"   I={I:3.1f} A, 1-eta={eta_def:.0e}: heat deficit {1e3*eta_def*E_TN_D2O*I:7.3f} mW, "
                     f"dp/dt = {dpdt*3600/1e3:8.3f} kPa/h")
    # pressure-based bound
    sig_p = 20.0     # Pa, rms resolution of a 0-3.5 bar abs. capacitance/quartz gauge after 1 h fit
    W = 3600.0
    t.append(f"   Pressure bound: sigma_p = {sig_p} Pa on a {W/3600:.0f} h slope fit -> sigma(dp/dt) ~ "
             f"{sig_p*np.sqrt(12)/W:.2e} Pa/s")
    for I in [0.2, 0.5, 1.0, 2.0]:
        ngas = 0.75 * I / F
        s_eta = (sig_p * np.sqrt(12) / W) * Vh / (R_GAS * T) / ngas
        t.append(f"     I={I:3.1f} A: sigma(1-eta) = {s_eta:.2e}  ->  heat bound {1e6*s_eta*E_TN_D2O*I:.1f} uW")
    t.append("   (temperature drift of the headspace must be corrected: dp/p = dT/T -> 1 mK gives 0.3 Pa)")

    t.append("")
    t.append("4. Loading enthalpy: D2O(l) -> D(in Pd) + 1/4 O2 stays unrecombined (no D2 partner)")
    e_load = (DHF_D2O / 2 + DH_ABS_D) / F        # V per electron that ends up in the lattice
    t.append(f"   effective 'thermoneutral' for the loading current fraction: {e_load:.4f} V "
             f"(= [dHf(D2O)/2 + dH_abs(D)]/F)")
    cath = {
        "C1 wire 1 mm x 30 mm": np.pi * (0.05) ** 2 * 3.0,             # cm^3
        "C1 wire 2 mm x 20 mm": np.pi * (0.10) ** 2 * 2.0,
        "C3 membrane 50 um x 20 mm dia": np.pi * 1.0 ** 2 * 50e-4,
        "C3 membrane 100 um x 25 mm dia": np.pi * 1.25 ** 2 * 100e-4,
    }
    rho_pd, M_pd = 12.02, 106.42
    for name, v in cath.items():
        nD = v * rho_pd / M_pd * 0.9
        E = nD * F * e_load
        nO2 = nD / 4
        dp = nO2 * R_GAS * T / Vh
        t.append(f"   {name:32s}: V={v*1e3:7.2f} mm3, D stored at x=0.9 = {nD*1e3:7.4f} mmol, "
                 f"chemical energy {E:7.1f} J, excess O2 dp = {dp/1e3:6.2f} kPa")
        for tau_h in [1, 10]:
            t.append(f"        released over {tau_h:2d} h on deloading -> mean {1e3*E/(tau_h*3600):7.2f} mW")
    t.append("   During loading at current I with loading fraction chi: deficit = chi * I * "
             f"{e_load:.3f} V, e.g. chi=0.3, I=0.5 A -> {1e3*0.3*0.5*e_load:.0f} mW")

    t.append("")
    t.append("5. C3 (DFM) permeation: D leaving through the membrane to vacuum carries ")
    t.append("   dHf(D2O)/2 per D = E_tn(D2O) per electron; O2 partner accumulates in the cell.")
    for J in [1e-9, 1e-8, 1e-7]:            # mol D / cm^2 / s
        A = np.pi * 1.0 ** 2
        Iperm = J * A * F
        dp = J * A / 4 * R_GAS * T / Vh * 3600
        t.append(f"   J={J:.0e} mol D/cm2/s over 3.14 cm2: I_perm={1e3*Iperm:7.3f} mA, heat deficit "
                 f"{1e3*Iperm*E_TN_D2O:8.3f} mW, O2 build-up {dp/1e3:7.3f} kPa/h")

    # --------------------------------------------------------------- 4He
    t.append("")
    t.append("6. 4He budget (H2 hypothesis: 23.85 MeV per 4He)")
    he_per_J = 1 / Q_DD_HE4
    t.append(f"   {he_per_J:.3e} 4He per J  = {he_per_J*1e-3*86400:.3e} atoms per mW-day "
             f"= {he_per_J*1e-3*86400/LOSCHMIDT:.3e} cm3 STP per mW-day")
    n_head = 1.0e5 * Vh / (K_B * T)
    t.append(f"   headspace {Vh*1e6:.0f} mL at 1 bar: {n_head:.2e} molecules; 1 mW-day -> "
             f"{1e6*he_per_J*1e-3*86400/n_head:.4f} ppm (air = {HE_AIR_FRAC*1e6:.2f} ppm)")
    p_he_air = HE_AIR_FRAC * 101325   # Pa
    rows = []
    # (a) Viton O-ring, 50 mm seal diameter, compressed section: flux ~ K * dp * pi*D
    for D in [0.030, 0.050, 0.070]:
        Q = K_HE_VITON * (p_he_air / 100) * np.pi * D * 1e6       # cm3 STP / s
        rows.append((f"Viton O-ring, seal dia {D*1e3:.0f} mm", Q * LOSCHMIDT))
    # (b) Pyrex wall, area 150 cm2, 2 mm thick, 25 and 60 C
    for Tc in [25, 60, 90]:
        K = K_HE_PYREX_25C * np.exp(-EA_HE_PYREX / R_GAS * (1 / (Tc + 273.15) - 1 / 298.15))
        Q = K * 150 * (p_he_air / 1333.22) / 2.0                    # cmHg partial pressure, 2 mm
        rows.append((f"Pyrex wall 150 cm2 x 2 mm at {Tc} C (steady state)", Q * LOSCHMIDT))
    # (c) metal seals (CF copper gasket / VCR): He leak spec < 1e-10 mbar L/s total at 1 bar He;
    #     at air partial pressure scale by 5.24e-6 -> negligible
    Qm = 1e-10 * 1e-3 * 1e5 / 1.01325e5 * HE_AIR_FRAC   # cm3STP/s
    rows.append(("all-metal seals (<=1e-10 mbar L/s He-leak-tight, air side)", Qm * 1e3 * LOSCHMIDT))
    t.append("   Air-He in-leak sources (atoms/s) and equivalent H2 'excess power'")
    for name, r in rows:
        t.append(f"     {name:58s}: {r:9.2e} /s  = {1e3*r*Q_DD_HE4:9.4f} mW-equivalent")
    # (d) dissolved He in air-saturated electrolyte
    m_el = 0.066   # kg (60 mL D2O)
    n_diss = KH_HE_WATER * (p_he_air / 1e5) * m_el * N_A
    t.append(f"   He dissolved in 60 mL air-saturated D2O: {n_diss:.2e} atoms = {n_diss*Q_DD_HE4:.0f} J-equivalent "
             f"(= {1e3*n_diss*Q_DD_HE4/86400:.1f} mW for a day) -> must degas (sparge with He-free D2) ")
    # (e) initial headspace air
    t.append(f"   1 % residual air in headspace: {0.01*HE_AIR_FRAC*n_head:.2e} atoms = "
             f"{0.01*HE_AIR_FRAC*n_head*Q_DD_HE4:.0f} J-equivalent -> evacuate/purge to <1e-4 air fraction")
    # headspace volume trade-off
    t.append("")
    t.append("   Headspace volume trade-off at 10 mW excess for 3 days (all He released to gas):")
    for V in [10e-6, 20e-6, 33e-6, 60e-6, 100e-6]:
        n = 1e5 * V / (K_B * T)
        c = he_per_J * 0.010 * 3 * 86400 / n
        # O2 excess from loading the 1 mm x 30 mm wire
        nO2 = np.pi * 0.05 ** 2 * 3.0 * rho_pd / M_pd * 0.9 / 4
        t.append(f"     V={V*1e6:5.0f} mL: 4He = {c*1e6:7.3f} ppm ; loading O2 dp = {nO2*R_GAS*T/V/1e3:6.1f} kPa")

    with open(os.path.join(OUT, "m4_chem_he.txt"), "w") as fh:
        fh.write("\n".join(t) + "\n")
    print("\n".join(t))

    # figure: equivalent power of He backgrounds vs signal
    fig, ax = plt.subplots(figsize=(7, 3.8))
    names = [r[0] for r in rows] + ["dissolved He (60 mL, spread over 1 day)"]
    vals = [1e3 * r[1] * Q_DD_HE4 for r in rows] + [1e3 * n_diss * Q_DD_HE4 / 86400]
    y = np.arange(len(names))
    ax.barh(y, vals, color="C1")
    ax.axvline(5, color="k", ls=":", label="5 mW (≈5σ heat channel)")
    ax.set_xscale("log"); ax.set_yticks(y); ax.set_yticklabels(names, fontsize=7)
    ax.set_xlabel("air-⁴He background expressed as equivalent D+D→⁴He power (mW)")
    ax.legend(fontsize=7); ax.set_title("⁴He backgrounds vs heat-equivalent signal")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "m4_he_backgrounds.png"), dpi=130)


if __name__ == "__main__":
    main()
