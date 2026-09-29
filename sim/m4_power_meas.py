"""
M4 — electrical input-power measurement with modulated (square-wave) cell current.

Cell model: galvanostat (first-order rise tau_g) driving a Randles cell
    V = E0 + eta + I*Rs ,   C_dl d(eta)/dt = I - 2 i0 sinh(eta/b)
Measurement methods compared (bias vs the true time average <V I>):
  M1  separate integrating DMMs: P = <V><I>              (misses cov(V,I))
  M2  simultaneous sampling of V and I with first-order anti-alias filters
      (bandwidths f_V, f_I) and inter-channel skew dt: P = <V_f(t) I_f(t+dt)>
  M3  level-gated: P = sum_levels I_level * <V>_level   (V averaged separately
      on the high and low half-cycles, with each level's current measured)

Run: python3 sim/m4_power_meas.py
"""
import os
import numpy as np

OUT = os.path.join(os.path.dirname(__file__), "..", "docs", "models", "figs")

# Cell parameters (0.1 M LiOD, 1 mm x 30 mm Pd cathode, coaxial Pt anode).
# Rs: coaxial gap ln(14/0.5)/(2 pi kappa L), kappa(0.1 M LiOH) ~ 2.2 S/m -> ~8 ohm; use 5 ohm
# (LiOD, 30-40 C).  Tafel b = 120 mV/dec -> asinh scale 52 mV; i0 = 1 mA;
# C_dl = 1 mF (roughened Pd, ~1 mF/cm^2 of real area x roughness; order-of-magnitude).
CELL = dict(E0=1.55, Rs=5.0, b=0.052, i0=1e-3, Cdl=1e-3, tau_g=5e-6)


def grid_period(T, n=2500):
    """Geometric grid clustered on both sides of each edge (t=0 and t=T/2)."""
    a = np.geomspace(1e-8, T / 4, n)
    g = np.unique(np.concatenate([[0.0], a, T / 2 - a, [T / 2]]))
    return np.concatenate([g[:-1], T / 2 + g])


def simulate(I0, m, f, cell=CELL, shape="square", n_periods=3):
    T = 1.0 / f
    tg = grid_period(T)
    t = np.concatenate([k * T + tg[:-1] for k in range(n_periods)] + [[n_periods * T]])
    if shape == "square":
        tgt = np.where((t % T) < T / 2, I0 * (1 + m), I0 * (1 - m))
        I = np.empty_like(t); I[0] = I0 * (1 - m)
        for k in range(1, len(t)):
            a = np.exp(-(t[k] - t[k - 1]) / cell["tau_g"])
            I[k] = tgt[k] + (I[k - 1] - tgt[k]) * a
    else:
        I = I0 * (1 + m * np.sin(2 * np.pi * f * t))
    eta = np.empty_like(t)
    eta[0] = cell["b"] * np.arcsinh(I[0] / (2 * cell["i0"]))
    for k in range(1, len(t)):
        dt = t[k] - t[k - 1]
        e = eta[k - 1]
        for _ in range(30):     # Newton on implicit Euler
            If = 2 * cell["i0"] * np.sinh(e / cell["b"])
            g = cell["Cdl"] * (e - eta[k - 1]) / dt - (I[k] - If)
            dg = cell["Cdl"] / dt + 2 * cell["i0"] / cell["b"] * np.cosh(e / cell["b"])
            de = g / dg
            e -= de
            if abs(de) < 1e-13:
                break
        eta[k] = e
    V = cell["E0"] + eta + I * cell["Rs"]
    # keep last period only (initial transient removed)
    s = t >= (n_periods - 1) * T
    return t[s] - t[s][0], V[s], I[s], T


def tavg(t, y, T):
    return np.trapezoid(y, t) / T


def lowpass(t, y, fc):
    if fc is None:
        return y
    tau = 1 / (2 * np.pi * fc)
    z = np.empty_like(y); z[0] = y[0]
    for k in range(1, len(t)):
        a = np.exp(-(t[k] - t[k - 1]) / tau)
        z[k] = y[k] + (z[k - 1] - y[k]) * a
    return z


def methods(I0, m, f, fV=100e3, fI=100e3, skew=1e-6, shape="square"):
    # run long enough for filters to settle: filter over 2 periods
    t, V, I, T = simulate(I0, m, f, shape=shape, n_periods=3)
    Ptrue = tavg(t, V * I, T)
    P1 = tavg(t, V, T) * tavg(t, I, T)
    # M2: filters applied over two concatenated copies of the (periodic) waveform
    t2 = np.concatenate([t, t[1:] + T]); V2 = np.concatenate([V, V[1:]]); I2 = np.concatenate([I, I[1:]])
    Vf = lowpass(t2, V2, fV); If = lowpass(t2, I2, fI)
    s = t2 >= T
    tt = t2[s] - T
    Ifs = np.interp((tt + skew) % T, tt, If[s])
    P2 = tavg(tt, Vf[s] * Ifs, T)
    # M3: level-gated (exclude nothing; average V on each half with that half's mean I)
    P3 = 0.0
    for sel in (t <= T / 2, t >= T / 2):
        ts = t[sel]
        P3 += np.trapezoid(V[sel], ts) * (np.trapezoid(I[sel], ts) / (ts[-1] - ts[0])) / T
    return Ptrue, P1, P2, P3


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    txt = ["M4 electrical power measurement with modulated current", "=" * 60,
           f"cell: {CELL}", ""]
    # verification: DC -> all methods identical; sine analytic check for M1: cov = Rs*m^2 I0^2/2 (+eta term)
    Pt, P1, P2, P3 = methods(0.5, 1e-9, 1.0)
    txt.append(f"verification DC (m=0): Ptrue={Pt:.6f} W, M1 bias={P1/Pt-1:+.1e}, M2 bias={P2/Pt-1:+.1e}")
    Pt, P1, P2, P3 = methods(0.5, 0.5, 1.0, shape="sine")
    # analytic: cov(V,I) = Rs*var(I) + cov(eta,I); var(I)=m^2 I0^2/2
    cov_R = CELL["Rs"] * 0.5 ** 2 * 0.5 ** 2 / 2
    txt.append(f"verification sine m=0.5 1 Hz: Ptrue-P1 = {Pt-P1:.5f} W; Rs*var(I) = {cov_R:.5f} W "
               f"(difference = eta-I covariance, {Pt-P1-cov_R:.5f} W)")
    txt.append("")
    I0 = 0.5
    fs = [0.01, 0.1, 1, 10, 100, 1000]
    ms = [0.1, 0.5, 0.9]
    res = {}
    txt.append("Square-wave modulation, I0 = 0.5 A.  Bias = (P_method - P_true)/P_true")
    txt.append(" m     f(Hz)   Ptrue(W)   M1 <V><I>      M2 (100k/100k, 1us)   M2 (1k/1.1k, 10us)   M3 gated")
    for m in ms:
        for f in fs:
            Pt, P1, P2, P3 = methods(I0, m, f)
            _, _, P2b, _ = methods(I0, m, f, fV=1e3, fI=1.1e3, skew=10e-6)
            res[(m, f)] = (Pt, P1, P2, P2b, P3)
            txt.append(f"{m:4.1f} {f:8.2f}   {Pt:8.4f}   {P1/Pt-1:+.3e}     {P2/Pt-1:+.3e}          "
                       f"{P2b/Pt-1:+.3e}        {P3/Pt-1:+.3e}")
    txt.append("")
    # skew / bandwidth mismatch scan at m=0.5, 10 Hz
    txt.append("M2 sensitivity at m=0.5, f=10 Hz and 1 kHz: skew and filter mismatch")
    for f in [10, 1000]:
        for fV, fI, sk in [(100e3, 100e3, 0), (100e3, 100e3, 1e-6), (100e3, 100e3, 10e-6),
                           (100e3, 90e3, 0), (10e3, 9e3, 0), (1e3, 0.9e3, 0)]:
            Pt, _, P2, _ = methods(I0, 0.5, f, fV=fV, fI=fI, skew=sk)
            txt.append(f"  f={f:5d} Hz fV={fV:8.0f} fI={fI:8.0f} skew={sk*1e6:5.1f} us  bias={P2/Pt-1:+.2e} "
                       f"({(P2-Pt)*1e3:+.3f} mW)")
    with open(os.path.join(OUT, "m4_power_meas.txt"), "w") as fh:
        fh.write("\n".join(txt) + "\n")
    print("\n".join(txt))

    fig, ax = plt.subplots(figsize=(6.5, 4))
    for m, c in zip(ms, ["C0", "C1", "C2"]):
        ax.plot(fs, [abs(res[(m, f)][1] / res[(m, f)][0] - 1) for f in fs], "o-", color=c,
                label=f"M1 <V><I>, m={m}")
        ax.plot(fs, [abs(res[(m, f)][3] / res[(m, f)][0] - 1) + 1e-9 for f in fs], "s--", color=c,
                label=f"M2 slow DAQ (1 kHz, 10% mismatch, 10 µs), m={m}")
        ax.plot(fs, [abs(res[(m, f)][2] / res[(m, f)][0] - 1) + 1e-9 for f in fs], "^:", color=c,
                label=f"M2 100 kHz, 1 µs skew, m={m}")
    ax.axhline(1e-4, color="k", lw=0.8, ls="-.", label="0.01 % budget")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(1e-9, 1)
    ax.set_xlabel("modulation frequency (Hz)"); ax.set_ylabel("|relative power bias|")
    ax.set_title("Input-power bias under square-wave current modulation (I0=0.5 A)")
    ax.legend(fontsize=6)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "m4_power_meas.png"), dpi=130)


if __name__ == "__main__":
    main()
