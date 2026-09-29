"""
M4 — realistic calorimetric accuracy for heated gas-phase reactors (C4/C5:
200-300 C, 10-100 W heater).

Lumped steady-state network:
  powder/sample node (p) --G_pw(gas)--> reactor wall (w)
  wall --capture path--> calorimeter (jacket / coolant / thermopile box at T_j)
  wall --G_loss--> ambient, NOT captured (ports, flange, supports, insulation)
Heat sources: heater (wound on the wall or embedded in the bed) and any
excess (in the bed).  Calibration = heater power with an inert bed.
Apparent excess arises when a perturbation changes the fraction of heat lost
through the uncaptured path.  Designs:

  G1  oil mass-flow coil brazed on the reactor, mineral-wool insulation outside
      (Kitamura/Takahashi-type; measured recovery 0.88 +- 0.03)
  G2  reactor radiating inside an evacuated, black, water-cooled jacket
      (flow or Seebeck read-out of the jacket); losses only via ports/supports
  G3  G2 with every port/lead thermally anchored to the jacket (loss path
      starts at T_j, not at T_w) and a twin inert-bed reactor for differencing

Run: python3 sim/m4_gasphase.py
"""
import os
import sys
import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(__file__))
from m4_params import SIGMA_SB  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), "..", "docs", "models", "figs")

T_AMB = 293.15
T_J = 303.15


def solve(design, P_heat, P_ex=0.0, eps=0.40, G_pw=0.5, G_loss=None, heater="wall",
          A=0.0339, T_oil=None, G_coil=None, G_ins=None):
    """Return (captured W, T_w, T_p)."""
    P = P_heat + P_ex
    if design == "G1":
        T_oil = 423.15 if T_oil is None else T_oil          # oil at ~150 C
        G_coil = 0.50 if G_coil is None else G_coil        # W/K wall->oil
        G_ins = 0.030 if G_ins is None else G_ins          # insulation + ports to ambient
        # wall balance: P = G_coil (T_w - T_oil) + G_ins (T_w - T_amb)
        T_w = (P + G_coil * T_oil + G_ins * T_AMB) / (G_coil + G_ins)
        cap = G_coil * (T_w - T_oil)
    else:
        G_loss = (0.010 if design == "G2" else 0.010) if G_loss is None else G_loss
        T_ref_loss = T_AMB if design == "G2" else T_J   # anchored ports start at jacket temperature

        def bal(Tw):
            rad = eps * SIGMA_SB * A * (Tw ** 4 - T_J ** 4)
            gas = 0.002 * (Tw - T_J)       # residual gas conduction in jacket vacuum (1e-2 Pa)
            if design == "G2":
                loss = G_loss * (Tw - T_ref_loss)
            else:  # anchored: a small residual fraction (5 %) of the port conductance still bypasses
                loss = 0.05 * G_loss * (Tw - T_AMB)
            return P - rad - gas - loss

        T_w = brentq(bal, T_J + 1e-6, 2000)
        loss = (G_loss * (T_w - T_ref_loss)) if design == "G2" else 0.05 * G_loss * (T_w - T_AMB)
        cap = P - loss
    # bed temperature: heat generated in the bed crosses G_pw
    P_bed = P_ex + (P_heat if heater == "bed" else 0.0)
    T_p = T_w + P_bed / G_pw
    return cap, T_w, T_p


def apparent_excess(design, P_heat, pert, **base):
    """Calibrate with base params, then apply perturbation, return apparent excess (W)."""
    cap0, Tw0, _ = solve(design, P_heat, **base)
    k = cap0 / P_heat          # calibration constant (captured fraction)
    kw = dict(base); kw.update(pert)
    cap1, Tw1, _ = solve(design, P_heat, **kw)
    return cap1 / k - P_heat, Tw0, Tw1


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    txt = ["M4 gas-phase reactor calorimetry (C4/C5)", "=" * 60]
    # choose eps such that 50 W gives T_w ~ 250 C in G2
    base = {"G1": dict(), "G2": dict(eps=0.40), "G3": dict(eps=0.40)}
    txt.append("Baseline operating points (captured fraction, wall temperature)")
    for d in ["G1", "G2", "G3"]:
        for P in [10, 50, 100]:
            cap, Tw, Tp = solve(d, P, **base[d])
            txt.append(f"  {d} P={P:3d} W: captured {cap/P:.4f}, T_wall = {Tw-273.15:6.1f} C")
    txt.append("")
    perts = [
        ("emissivity +10 % (surface oxidation/deposition)", dict(eps=0.44)),
        ("emissivity -25 % (reduction/cleaning)", dict(eps=0.30)),
        ("insulation/port conductance +10 %", "ins"),
        ("oil inlet temperature +2 K (G1)", dict(T_oil=425.15)),
        ("coil-to-wall conductance -5 % (G1, oil viscosity/fouling)", dict(G_coil=0.475)),
    ]
    rows = []
    txt.append("Apparent excess (W and % of heater power) after a single perturbation, calibration held fixed")
    for name, pert in perts:
        for d in ["G1", "G2", "G3"]:
            if pert == "ins":
                pp = dict(G_ins=0.033) if d == "G1" else dict(G_loss=0.011)
            else:
                pp = pert
                if d == "G1" and "eps" in pp:
                    continue
                if d != "G1" and ("T_oil" in pp or "G_coil" in pp):
                    continue
            for P in [10, 50, 100]:
                ex, Tw0, Tw1 = apparent_excess(d, P, pp, **base[d])
                rows.append((name, d, P, ex))
                txt.append(f"  {name:55s} {d} P={P:3d} W: {ex*1e3:+9.1f} mW ({100*ex/P:+.3f} %)  "
                           f"T_w {Tw0-273.15:.1f}->{Tw1-273.15:.1f} C")
    txt.append("")
    txt.append("Source location: heater on wall vs excess in bed changes nothing in the lumped network")
    txt.append("(both leave through the wall), but the bed-to-wall gas conductance sets the powder temperature:")
    for gas, Gpw in [("D2 (k=0.14 W/mK)", 0.5), ("H2 (k=0.18 W/mK)", 0.5 * 0.187 / 0.141), ("vacuum", 0.05)]:
        _, Tw, Tp = solve("G2", 50, P_ex=0.0, heater="bed", G_pw=Gpw, eps=0.40)
        txt.append(f"  G2 50 W heater in bed, {gas:18s}: T_bed - T_wall = {Tp-Tw:6.1f} K")
    txt.append("  -> the D2 run and the H2 control run have different bed temperatures at equal heater power;")
    txt.append("     a thermally-activated artefact (e.g. desorption/absorption enthalpy) is NOT cancelled by the control.")
    txt.append("")
    # chemistry: absorption enthalpy transient scale
    nPd = 0.1 * 10 / 106.42      # 10 g composite with 10 % Pd (mol)  (PNZ-type, illustrative)
    for x in [0.5, 0.8]:
        E = nPd * x * 17.3e3      # J, heat of D absorption
        txt.append(f"  absorption heat, 10 g PNZ-type (1 g Pd), x={x}: {E:.0f} J -> {E/3600*1e3:.0f} mW for 1 h")
    txt.append("  (hydrogen isotope absorption/desorption and H/D exchange heats must be integrated and subtracted)")

    # realistic accuracy summary: quadrature of the perturbations that cannot be excluded
    txt.append("")
    txt.append("Realistic systematic accuracy (quadrature of listed perturbations at 50 W):")
    for d in ["G1", "G2", "G3"]:
        e = [r[3] for r in rows if r[1] == d and r[2] == 50]
        tot = np.sqrt(np.sum(np.square(e)))
        txt.append(f"  {d}: +-{tot*1e3:.0f} mW = +-{100*tot/50:.2f} % of input")
    with open(os.path.join(OUT, "m4_gasphase.txt"), "w") as fh:
        fh.write("\n".join(txt) + "\n")
    print("\n".join(txt))

    fig, ax = plt.subplots(figsize=(7.5, 4))
    labels = sorted(set((r[0], r[1]) for r in rows), key=lambda z: (z[1], z[0]))
    y = np.arange(len(labels))
    vals = [abs([r[3] for r in rows if (r[0], r[1]) == lab and r[2] == 50][0]) for lab in labels]
    ax.barh(y, vals, color=["C0" if l[1] == "G1" else "C1" if l[1] == "G2" else "C2" for l in labels])
    ax.set_yticks(y); ax.set_yticklabels([f"{l[1]}: {l[0]}" for l in labels], fontsize=6)
    ax.set_xscale("log"); ax.set_xlabel("|apparent excess| at 50 W heater (W)")
    ax.set_title("Gas-phase reactor: calibration-shift artefacts")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "m4_gasphase.png"), dpi=130)


if __name__ == "__main__":
    main()
