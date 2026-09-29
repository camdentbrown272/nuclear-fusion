"""
M4 — Shanahan calibration-constant-shift (CCS) scenarios and the error budget
for each closed-cell calorimeter design.  Uses the 2D model (m4_thermal2d).

Estimator assumed for every design (pre-registered):
  P_out = [ Signal - (k_rec - k_cath) * E_tn * I ] / k_cath
i.e. calibrated with heaters at the cathode position AND at the recombiner,
assuming all recombination happens in the recombiner (verified by pressure).
CCS scenarios then move heat between locations without changing the total.

Run: python3 sim/m4_budget.py
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import m4_thermal2d as th  # noqa: E402
from m4_params import E_TN_D2O  # noqa: E402

OUT = th.OUT
DESIGNS = ["ISO", "FLOW", "SEEB0", "SEEB1", "SEEB1*"]
KW = {"SEEB1*": dict(f_top=0.85, t_shell=10 * th.mm)}      # recommended variant
OPS = [(0.2, 3.2), (0.5, 4.8), (1.2, 8.0)]                  # (A, V) -> 0.64, 2.4, 9.6 W


def kconst(design, **extra):
    d = "SEEB1" if design == "SEEB1*" else design
    kw = dict(KW.get(design, {})); kw.update(extra)
    r, _, _ = th.cal_constants(d, **kw)
    return r


def split(I, V):
    Prec = E_TN_D2O * I
    Pnr = (V - E_TN_D2O) * I
    return dict(recomb=Prec, cathode=0.2 * Pnr, anode=0.2 * Pnr, joule=0.6 * Pnr)


def measured(k, dist, ref="heater_cath"):
    return sum(P * k[p] for p, P in dist.items()) / k[ref]


def main():
    txt = ["M4 CCS scenarios and error budget", "=" * 60]
    K = {d: kconst(d) for d in DESIGNS}
    K3 = {d: kconst(d, k_el=3.0) for d in DESIGNS}
    KH = {d: kconst(d, k_head=2.0) for d in DESIGNS}
    ccs = {}
    txt.append("CCS scenarios: apparent excess (mW) and (% of input) with the two-heater estimator")
    txt.append("  A: 100 % of recombination moves from recombiner to cathode (O2 reduction at cathode)")
    txt.append("  B: 100 % of recombination moves to the electrolyte surface (flooded catalyst)")
    txt.append("  C: electrolyte stirring k_el 10 -> 3 W/mK (bubble regime change, lower current)")
    txt.append("  D: headspace effective conductivity 0.5 -> 2 W/mK (condensation/pressure change)")
    for I, V in OPS:
        P = I * V
        txt.append(f" Operating point I={I} A, V={V} V, P_in={P:.2f} W, recombiner heat {E_TN_D2O*I:.3f} W")
        dist = split(I, V)
        for d in DESIGNS:
            k = K[d]
            est0 = measured(k, dist) - (k["recomb"] - k["heater_cath"]) / k["heater_cath"] * dist["recomb"]
            def apparent(kk, dd):
                e = measured(kk, dd) - (k["recomb"] - k["heater_cath"]) / k["heater_cath"] * dist["recomb"]
                # calibration held at the baseline constants (k), physics uses kk
                e = sum(Pp * kk[p] for p, Pp in dd.items()) / k["heater_cath"] - \
                    (k["recomb"] - k["heater_cath"]) / k["heater_cath"] * dist["recomb"]
                return e - P
            dA = dict(dist); dA["cathode"] += dA.pop("recomb"); dA["recomb"] = 0.0
            dB = dict(dist); dB["surface"] = dB.pop("recomb"); dB["recomb"] = 0.0
            a = apparent(k, dA); b = apparent(k, dB); c = apparent(K3[d], dist); e = apparent(KH[d], dist)
            ccs[(d, I)] = (a, b, c, e, P)
            txt.append(f"   {d:6s} A {1e3*a:+9.2f} ({100*a/P:+.3f} %)  B {1e3*b:+9.2f} ({100*b/P:+.3f} %)  "
                       f"C {1e3*c:+9.2f} ({100*c/P:+.3f} %)  D {1e3*e:+9.2f} ({100*e/P:+.3f} %)   "
                       f"[baseline closure {1e3*(est0-P):+.2e} mW]")
    txt.append("")
    txt.append("Shanahan's reanalysis needs a +-2.5 % calibration-constant shift (Thermochim. Acta 382 (2002) 95).")
    txt.append("Maximum CCS each design can generate from ANY redistribution among modelled locations (% of the")
    txt.append("redistributed heat) = max |k_i/k_j - 1| over sources:")
    for d in DESIGNS:
        v = np.array([K[d][p] for p in th.POS])
        txt.append(f"   {d:6s}: {100*(v.max()/v.min()-1):.3f} %")

    # ------------------------------------------------------------------ error budgets
    # random (1 h window, 10 mK rms bath, tau_b 600 s) from m4_noise.txt
    noise = {"ISO": 0.08, "FLOW": 0.86, "SEEB0": 0.66, "SEEB1": 0.66, "SEEB1*": 0.66, "TWIN": 0.22}
    txt.append("")
    txt.append("Error budget at P_in = 2.4 W (0.5 A, 4.8 V) and 9.6 W (1.2 A, 8 V).  Systematic items in % of")
    txt.append("input unless marked; mW items converted.  CCS items use a plausible bound phi<=0.2 of the")
    txt.append("recombination relocating (worst-case phi=1 is in the CCS table above).")
    budget = {}
    for d in ["ISO", "FLOW", "SEEB1", "SEEB1*"]:
        for I, V in OPS[1:]:
            P = I * V
            a, b, c, e, _ = ccs[(d, I)]
            items = [
                ("cal. heater power (4-wire DC, 6.5 digit)", 0.010),
                ("input power, DC or modulated with 100 kS/s simultaneous V,I (<=100 Hz)", 0.005),
                ("calibration fit / repeatability", {"ISO": 0.10, "FLOW": 0.05, "SEEB1": 0.02, "SEEB1*": 0.02}[d]),
                ("CCS: recombination relocation phi<=0.2", 100 * 0.2 * max(abs(a), abs(b)) / P),
                ("CCS: stirring change (k_el 10->3)", 100 * abs(c) / P),
                ("CCS: headspace transport change", 100 * abs(e) / P),
                ("sensitivity drift (bath +-0.02 K; ISO radiative T^3; FLOW flow-meter)",
                 {"ISO": 0.10, "FLOW": 0.15, "SEEB1": 0.004, "SEEB1*": 0.004}[d]),
                ("heat-loss / bypass change (leads, insulation, coolant loop)",
                 {"ISO": 0.20, "FLOW": 0.20, "SEEB1": 0.01, "SEEB1*": 0.01}[d]),
                ("unrecombined gas (pressure-bounded, 50 uW)", 100 * 0.05e-3 / P),
                ("allowance: unmodelled 3D/contact/ageing effects (assumption)",
                 {"ISO": 0.30, "FLOW": 0.10, "SEEB1": 0.05, "SEEB1*": 0.05}[d]),
            ]
            if d == "FLOW":
                items += [("mass-flow calibration", 0.10), ("coolant c_p, density", 0.02),
                          ("inlet/outlet PRT calibration (dT ~ 2-8 K)", 0.05)]
            sys_tot = np.sqrt(sum(x[1] ** 2 for x in items))
            budget[(d, I)] = (items, sys_tot, noise[d])
            txt.append(f" {d:6s} P_in={P:.1f} W")
            for name, v in items:
                txt.append(f"     {name:70s} {v:7.4f} %  ({10*v*P:7.3f} mW)")
            txt.append(f"     {'TOTAL systematic (quadrature)':70s} {sys_tot:7.4f} %  ({10*sys_tot*P:7.3f} mW)")
            txt.append(f"     {'random sigma_P, 1 h, 10 mK bath (TWIN: ' + str(noise['TWIN']) + ' mW)':70s} "
                       f"          ({noise[d]:7.3f} mW)")
    with open(os.path.join(OUT, "m4_budget.txt"), "w") as fh:
        fh.write("\n".join(txt) + "\n")
    print("\n".join(txt))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4))
    x = np.arange(len(DESIGNS))
    for n, (lab, idx) in enumerate([("A: recomb→cathode", 0), ("B: recomb→surface", 1),
                                    ("C: stirring", 2), ("D: headspace", 3)]):
        ax.bar(x + (n - 1.5) * 0.2, [abs(100 * ccs[(d, 0.5)][idx] / ccs[(d, 0.5)][4]) + 1e-5 for d in DESIGNS],
               0.2, label=lab)
    ax.axhline(2.5, color="r", ls="--", lw=1, label="Shanahan CCS 2.5 %")
    ax.axhline(0.1, color="k", ls=":", lw=1, label="0.1 %")
    ax.set_yscale("log"); ax.set_xticks(x); ax.set_xticklabels(DESIGNS)
    ax.set_ylabel("|apparent excess| (% of 2.4 W input)")
    ax.set_title("Worst-case CCS by design (0.5 A, 4.8 V)")
    ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "m4_ccs.png"), dpi=130)


if __name__ == "__main__":
    main()
