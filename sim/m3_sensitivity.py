"""
M3 sensitivities: how the key design numbers move with the least-certain inputs.
  - beta-branch anchor f(x=0.90) for PdD at 25 C (1e3 / 1e4 / 1e5 atm)
  - site blocking in the tracer diffusivity (D* = D_a(1-x) vs D_a)
  - temperature (20-60 C) and isotope (H control)
Outputs: flux budget J_max (x_in = 0.95 -> exit 0.90, L = 25 um), required exit k_r*, alpha/beta
loading-time constant tau90 (slab, x_s = 0.9).

Run:  python3 sim/m3_sensitivity.py   -> docs/models/figs/m3_sensitivity.txt
"""
import os
import numpy as np
import m3_common as mc
from m3_loading import t90_numeric

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "docs", "models", "figs")
LOG = []


def log(s=""):
    print(s)
    LOG.append(s)


def metrics(T=298.15, iso="D", blocking=True, L=25e-6, xin=0.95):
    m = mc.Material(T, iso, blocking=blocking)
    Jmax = (m.Phi(xin) - m.Phi(0.90)) / L
    krs = 1e21 / mc.c_eff(0.9, T, iso) ** 2
    tau = t90_numeric(m, 0, 0.9, ell=1e-6, N=100) * mc.D_alpha(T, iso) / 1e-12
    t90_25 = tau * (25e-6) ** 2 / mc.D_alpha(T, iso)
    return Jmax, krs, tau, t90_25, mc.x_of_f(mc.ATM, T, iso), mc.fugacity(0.9, T, iso) / mc.ATM


def main():
    log("case                                   J_max(D/m2/s)  [mA/cm2]  k_r*(m4/s)  tau90   t90(25um)  x(1atm)  f(0.9) atm")
    base = mc.KB298

    def row(name, *a, **k):
        J, kr, tau, t9, x1, f9 = metrics(*a, **k)
        log(f"{name:38s} {J:10.3g}  [{J * mc.e * 0.1:7.1f}]  {kr:9.3g}  {tau:6.3f}  {t9:7.2f} s   {x1:.3f}   {f9:9.3g}")
    row("baseline (D, 25 C, blocking)")
    for F in (1e3, 1e5):
        mc.F90_ATM = F
        mc.KB298 = mc._kb298()
        row(f"anchor f(0.90) = {F:.0e} atm")
    mc.F90_ATM = 1e4
    mc.KB298 = base
    row("no site blocking (D* = D_alpha)", blocking=False)
    row("T = 20 C", T=293.15)
    row("T = 40 C", T=313.15)
    row("T = 60 C", T=333.15)
    row("H control, 25 C", iso="H")
    # exit recombination temperature scaling
    for T in (293.15, 313.15, 333.15):
        log(f"k_r(Pd lit.) at {T - 273.15:.0f} C = {mc.kr_Pd_lit(T):.3g} m4/s ; J_sat lower bound (Ed 0.75 eV) = "
            f"{mc.Jsat_desorb(T):.3g} D/m2/s ; ratio k_r*/k_r(lit) = "
            f"{1e21 / mc.c_eff(0.9, T) ** 2 / mc.kr_Pd_lit(T):.3g}")
    open(os.path.join(OUT, "m3_sensitivity.txt"), "w").write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
