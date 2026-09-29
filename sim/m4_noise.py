"""
M4 — noise floor and environmental-drift model for ISO, FLOW and SEEB1 (single and twin).

Two-node lumped reduction of the 2D model (node 1 = electrolyte + internals,
node 2 = cell wall / Cu shell / jacket metal), parameters extracted from
m4_thermal2d.build().  The sink (bath / heat-sink block / coolant inlet)
temperature fluctuates as an Ornstein-Uhlenbeck process with rms sigma_b and
correlation time tau_b.  Each design's estimator is simulated at 1 s sampling
for 72 h with constant true power; the rms error of window-averaged power is
reported vs window length.

Estimators
  ISO   : P = K (T1 - Tb) + C dT1/dt                  (F&P lumped equation)
  FLOW  : P = mdot c (T_out - T_in) / eta + C dT2/dt (correction uses jacket thermistor)
  SEEB1 : P = V/S + C dT2/dt   (Tian correction with a shell thermistor)  or V/S only
  TWIN  : difference of two SEEB1 calorimeters in one bath (sink correlation rho,
          2 % mismatch in C and K), each corrected with its own calibration.

Run: python3 sim/m4_noise.py
"""
import os
import sys
import numpy as np
from scipy import signal

sys.path.insert(0, os.path.dirname(__file__))
import m4_thermal2d as th  # noqa: E402

OUT = th.OUT
RNG = np.random.default_rng(20260929)
DT = 1.0
NSTEP = 72 * 3600


def reduce_2node(design):
    m, info = th.build(design)
    S = th.sources(m, th.G)
    T = m.solve(S["uniform_el"])
    if design == "SEEB1":     # node 1 = whole cell (+ pad), node 2 = Cu shell under the thermopile
        inner = np.isin(m.label, ["el", "head", "rec", "rod", "cath", "an", "liner", "wall", "pad"])
    else:
        inner = np.isin(m.label, ["el", "head", "rec", "rod", "cath", "an", "liner"])
    outer_labels = {"ISO": ["wall"], "FLOW": ["wall", "pad"], "SEEB1": ["shell"]}[design]
    outer = np.isin(m.label, outer_labels) & ~m.dirichlet
    C = m.rhoc * m.vol
    C1, C2 = C[inner].sum(), C[outer].sum()
    T1 = (C * T)[inner].sum() / C1
    T2 = (C * T)[outer].sum() / C2
    K = 1.0 / T2
    G12 = 1.0 / (T1 - T2)
    return dict(C1=C1, C2=C2, K=K, G12=G12)


def ou(n, sigma, tau, rng):
    a = np.exp(-DT / tau)
    e = rng.normal(0, sigma * np.sqrt(1 - a * a), n)
    x = np.empty(n); x[0] = rng.normal(0, sigma)
    for k in range(1, n):
        x[k] = a * x[k - 1] + e[k]
    return x


def simulate_nodes(p, P, Ts):
    """States T1, T2 relative to 0; inputs P (W) and sink temperature Ts (K)."""
    C1, C2, K, G = p["C1"], p["C2"], p["K"], p["G12"]
    A = np.array([[-G / C1, G / C1], [G / C2, -(G + K) / C2]])
    B = np.array([[1 / C1, 0], [0, K / C2]])
    Cm = np.eye(2); D = np.zeros((2, 2))
    sysd = signal.cont2discrete((A, B, Cm, D), DT, method="zoh")
    Ad, Bd = sysd[0], sysd[1]
    u = np.vstack([P, Ts])
    x = np.zeros((2, len(Ts)))
    x0 = np.linalg.solve(A, -B @ u[:, 0])
    x[:, 0] = x0
    for k in range(1, len(Ts)):
        x[:, k] = Ad @ x[:, k - 1] + Bd @ u[:, k - 1]
    return x


def window_rms(err, W):
    n = int(W / DT)
    m = len(err) // n
    e = err[: m * n].reshape(m, n).mean(axis=1)
    return np.sqrt(np.mean(e ** 2))


def deriv(y):
    return np.gradient(y, DT)


def run(design, p, sigma_b, tau_b=600.0, P0=5.0, rng=RNG, twin=False, rho=0.95, correct=True):
    n = NSTEP
    Ts = ou(n, sigma_b, tau_b, rng)
    P = np.full(n, P0)
    x = simulate_nodes(p, P, Ts)
    T1, T2 = x
    C = p["C1"] + p["C2"]
    if design == "ISO":
        sT = 1e-4                 # K rms per 1 s reading, each thermistor (0.1 mK class bridge)
        T1m = T1 + rng.normal(0, sT, n); Tbm = Ts + rng.normal(0, sT, n)
        Kiso = 1.0 / (1 / p["K"] + 1 / p["G12"])         # calibrated steady-state constant
        est = Kiso * (T1m - Tbm) + (C * deriv(T1m) if correct else 0)
    elif design == "FLOW":
        mdc = 4181 * 1.2e-3                               # W/K for 1.2 g/s (72 g/min) water
        tr = 30                                           # s transit inlet->outlet sensor
        sT = 1e-3                                         # K rms, matched PRT pair per reading
        s_flow = 1e-3                                     # 0.1 % rms flow fluctuation (per 60 s)
        flow = 1 + np.repeat(rng.normal(0, s_flow, n // 60 + 1), 60)[:n]
        Q = p["K"] * (T2 - Ts)                           # heat into coolant
        Tin = Ts
        Tout = np.roll(Tin, tr) + Q / (mdc * flow)
        Tout[:tr] = Tout[tr]
        Tin_d = np.roll(Tin, tr); Tin_d[:tr] = Tin[0]      # inlet reading delayed by transit time
        dTm = Tout - Tin_d + rng.normal(0, sT * np.sqrt(2), n)
        est = mdc * dTm + (C * deriv(T2 + rng.normal(0, 1e-4, n)) if correct else 0)
    elif design in ("SEEB1", "TWIN"):
        S = 0.0988                                       # V/W (m4_thermal2d)
        sV = 50e-9                                       # V rms per 1 s reading (nanovoltmeter)
        V = S * p["K"] * (T2 - Ts) + rng.normal(0, sV, n)   # thermopile ~ heat flux through it
        est = V / S + (C * deriv(T2 + rng.normal(0, 1e-4, n)) if correct else 0)
        if design == "TWIN":
            Ts2 = rho * Ts + np.sqrt(1 - rho ** 2) * ou(n, sigma_b, tau_b, rng)
            p2 = dict(p); p2["C1"] *= 1.02; p2["C2"] *= 1.02; p2["K"] *= 1.02
            x2 = simulate_nodes(p2, P, Ts2)
            V2 = S * p2["K"] * (x2[1] - Ts2) + rng.normal(0, sV, n)
            est2 = V2 / S + (1.02 * C * deriv(x2[1] + rng.normal(0, 1e-4, n))
                                                   if correct else 0)
            est = est - est2 + P0          # differential channel; true difference = 0
    else:
        raise ValueError
    err = est - P0
    burn = 6 * 3600
    return err[burn:]


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    txt = ["M4 noise floor / environmental drift (72 h simulations, 1 s sampling, P = 5 W)", "=" * 60]
    params = {d: reduce_2node(d) for d in ["ISO", "FLOW", "SEEB1"]}
    params["TWIN"] = params["SEEB1"]
    for d, p in params.items():
        txt.append(f"{d:6s}: C1={p['C1']:.0f} J/K  C2={p['C2']:.0f} J/K  G12={p['G12']:.3f} W/K  K={p['K']:.3f} W/K")
    Ws = [60, 300, 900, 3600, 4 * 3600]
    sigmas = [1e-3, 1e-2, 1e-1]
    res = {}
    txt.append("")
    txt.append("rms error of window-averaged power (mW); sink OU fluctuation sigma_b, tau_b = 600 s")
    txt.append("design  corr  sigma_b(mK)  " + "  ".join(f"W={w:>5d}s" for w in Ws))
    for d in ["ISO", "FLOW", "SEEB1", "TWIN"]:
        for corr in [False, True]:
            for sb in sigmas:
                err = run(d, params[d], sb, correct=corr, rng=np.random.default_rng(1))
                r = [1e3 * window_rms(err, W) for W in Ws]
                res[(d, corr, sb)] = r
                txt.append(f"{d:6s}  {str(corr):5s}  {sb*1e3:8.0f}     " + "  ".join(f"{x:9.3f}" for x in r))
    # requirement: sigma_b giving 1 mW at 1 h (errors scale ~linearly with sigma_b once above sensor floor)
    txt.append("")
    txt.append("Required sink stability (rms, tau_b = 600 s) for sigma_P <= 1 mW at 1 h averaging, with correction:")
    for d in ["ISO", "FLOW", "SEEB1", "TWIN"]:
        r1, r2 = res[(d, True, 1e-2)][3], res[(d, True, 1e-1)][3]
        slope = (r2 - r1) / (0.1 - 0.01)          # mW per K
        floor = max(r1 - slope * 0.01, 0)
        need = (1.0 - floor) / slope if floor < 1 else float("nan")
        txt.append(f"  {d:6s}: floor {floor:.3f} mW, slope {slope/1e3:.3f} mW/mK -> sigma_b <= {need*1e3:.1f} mK")
    # systematic drift of sensitivity
    txt.append("")
    txt.append("Sensitivity drift with bath set-point (systematic, per K of long-term bath drift):")
    txt.append("  SEEB1: Bi2Te3 module S(T) ~ +0.2 %/K  -> 10 W x 0.2 %/K = 20 mW/K -> bath set-point +-0.02 K gives +-0.4 mW")
    txt.append("  ISO  : radiative K ~ T^3 -> 1.0 %/K of K; plus electrolyte-level/emissivity drifts (not modelled)")
    txt.append("  FLOW : c_p(T) and density known to <0.01 %; flow-meter calibration drift 0.1-0.2 %/yr (Coriolis spec)")
    with open(os.path.join(OUT, "m4_noise.txt"), "w") as fh:
        fh.write("\n".join(txt) + "\n")
    print("\n".join(txt))

    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    for d, c in zip(["ISO", "FLOW", "SEEB1", "TWIN"], ["C0", "C1", "C2", "C3"]):
        ax.plot(Ws, res[(d, True, 1e-2)], "o-", color=c, label=f"{d}, σ_b=10 mK, Tian/FP-corrected")
        ax.plot(Ws, res[(d, False, 1e-2)], "s:", color=c, alpha=0.6, label=f"{d}, σ_b=10 mK, uncorrected")
    ax.axhline(1, color="k", lw=0.8, ls="--")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("averaging window (s)"); ax.set_ylabel("rms power error (mW)")
    ax.legend(fontsize=6); ax.set_title("Noise floor with 10 mK rms sink fluctuations (τ_b=600 s)")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "m4_noise.png"), dpi=130)
    return params, res


if __name__ == "__main__":
    main()
