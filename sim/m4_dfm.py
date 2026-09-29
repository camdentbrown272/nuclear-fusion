"""
M4 — can a heat channel be added to the detector-facing membrane cell (C3)?

Lumped steady network for heat released at the vacuum-side (front) face of a
Pd membrane that is the cathode of a closed electrolytic half-cell inside a
SEEB1-type Seebeck envelope; the vacuum chamber with the Si detectors is
outside the calorimeter and is joined to the membrane flange by a vacuum
nipple that must cross the thermopile wall.

Nodes: membrane m, membrane clamp flange f, calorimeter (cell/shell) c,
vacuum chamber/ambient v (= sink temperature, not captured).
  m-c : membrane back face to electrolyte, h_e * A_m
  m-f : radial conduction in the membrane to the clamped rim, 8 pi k t
  m-v : radiation from front face to the detectors/chamber
  f-c : flange to cell body + wetted flange area
  f-v : flange to chamber through the thermal break (the design lever)
Captured fraction for membrane heat vs for the calibration heater in the
electrolyte gives the position-dependent error.

Run: python3 sim/m4_dfm.py
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from m4_params import SIGMA_SB, E_TN_D2O, F  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), "..", "docs", "models", "figs")


def network(G):
    """Solve for T_m, T_f with c and v held at 0 (c's rise is absorbed in the calibration).
    Returns captured fraction for unit heat at m and at f."""
    # unknowns Tm, Tf
    A = np.array([[G["mc"] + G["mf"] + G["mv"], -G["mf"]],
                  [-G["mf"], G["mf"] + G["fc"] + G["fv"]]])
    out = {}
    for name, q in [("membrane", [1.0, 0.0]), ("flange", [0.0, 1.0])]:
        Tm, Tf = np.linalg.solve(A, q)
        lost = G["mv"] * Tm + G["fv"] * Tf
        out[name] = 1 - lost
        out[name + "_T"] = Tm
    return out


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    txt = ["M4 DFM (C3) heat channel", "=" * 60]
    r_m, t_m, k_pd = 10e-3, 50e-6, 71.8
    A_m = np.pi * r_m ** 2
    eps_pd, T = 0.10, 303.0
    base = dict(mc=1000 * A_m, mf=8 * np.pi * k_pd * t_m, mv=4 * eps_pd * SIGMA_SB * T ** 3 * A_m, fc=0.30)
    txt.append("Conductances (W/K): " + ", ".join(f"{k}={v:.4g}" for k, v in base.items()))
    txt.append("")
    txt.append("Captured fraction of heat released at the membrane front face vs thermal-break conductance f->v")
    txt.append("(electrolyte calibration heater: captured fraction = 1 by definition of the reference)")
    breaks = [("316L weld, short DN40 nipple", 0.30), ("thin-wall 0.25 mm edge-welded bellows 50 mm", 0.03),
              ("ceramic (Al2O3) break + bellows", 0.010), ("PEEK/ceramic break + long bellows", 0.003)]
    res = []
    for name, Gfv in breaks:
        for h_e in [300, 1000, 3000]:
            G = dict(base); G["fv"] = Gfv; G["mc"] = h_e * A_m
            o = network(G)
            res.append((name, Gfv, h_e, o["membrane"], o["flange"], o["membrane_T"]))
            txt.append(f"  {name:45s} G_fv={Gfv:.3f}  h_e={h_e:5d}: membrane-heat captured {o['membrane']:.5f} "
                       f"(loss {100*(1-o['membrane']):.3f} %), flange-heat captured {o['flange']:.4f}, "
                       f"T_m rise {o['membrane_T']:.2f} K/W")
    txt.append("")
    txt.append("The bypass is itself calibratable only if a heater sits AT the membrane (thin-film heater on the")
    txt.append("flange rim, or ohmic heating of the membrane by a 4-wire current); the residual uncertainty is")
    txt.append("the variation of h_e (stirring) which moves the membrane/flange split.")
    for name, Gfv in breaks:
        caps = [r[3] for r in res if r[1] == Gfv]
        txt.append(f"  {name:45s}: captured varies {min(caps):.5f}-{max(caps):.5f} over h_e 300-3000 -> "
                   f"+-{50*(max(caps)-min(caps)):.3f} % of membrane heat")
    txt.append("")
    txt.append("Detection limit for membrane heat (P_in from loading current density):")
    vol_cm3 = A_m * t_m * 1e6
    for j in [0.05, 0.2, 0.5]:      # A/cm^2 on the back face
        I = j * A_m * 1e4
        V = 3.0 + 5.0 * I          # V, rough: E0+overpotentials + Rs(5 ohm) I
        Pin = V * I
        sig_sys = 1e-3 * Pin       # 0.1 % of input systematic (SEEB1 budget)
        sig_rnd = 0.3e-3           # W, 1 h window, twin, 10 mK bath (m4_noise)
        sig_perm = 0.5e-3          # W, uncertainty on permeation enthalpy leak (J known to +-10 % at 1e-8 mol/cm2/s)
        sig = np.sqrt(sig_sys ** 2 + sig_rnd ** 2 + sig_perm ** 2)
        txt.append(f"  j={j:4.2f} A/cm2: I={I:.3f} A, P_in~{Pin:.2f} W, sigma_P~{sig*1e3:.2f} mW, 5 sigma = "
                   f"{5*sig*1e3:.1f} mW = {5*sig/vol_cm3:.2f} W/cm3 of Pd (membrane {vol_cm3*1e3:.1f} mm3)")
    txt.append("")
    txt.append("Permeation enthalpy leak (D leaves to vacuum as D2; E_tn per electron):")
    for J in [1e-9, 1e-8, 1e-7]:
        P = J * A_m * 1e4 * F * E_TN_D2O
        txt.append(f"  J={J:.0e} mol D/cm2/s: {P*1e3:.3f} mW (must be measured: vacuum-side D2 flow or cell O2 rise)")
    with open(os.path.join(OUT, "m4_dfm.txt"), "w") as fh:
        fh.write("\n".join(txt) + "\n")
    print("\n".join(txt))

    fig, ax = plt.subplots(figsize=(6, 3.8))
    for h_e, c in zip([300, 1000, 3000], ["C0", "C1", "C2"]):
        xs = [r[1] for r in res if r[2] == h_e]
        ys = [100 * (1 - r[3]) for r in res if r[2] == h_e]
        ax.plot(xs, ys, "o-", color=c, label=f"h_e = {h_e} W/m²K")
    ax.axhline(0.1, color="k", ls=":", label="0.1 %")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("flange→chamber thermal-break conductance (W/K)")
    ax.set_ylabel("membrane heat lost to chamber (%)")
    ax.legend(fontsize=7); ax.set_title("C3 heat channel: bypass through the vacuum side")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "m4_dfm_bypass.png"), dpi=130)


if __name__ == "__main__":
    main()
