"""
M0 — First-principles D-D reaction-rate budget for a cold (room-temperature) lattice.

Purpose: turn "make cold fusion measurable" into numbers.
  1. How many reactions per second does each detection channel need?
  2. What effective barrier reduction ("screening energy" Ue_eff) would a
     population of N deuteron pairs need to reach that rate?
  3. How steep is the dependence (-> which sites dominate the rate)?
  4. Why heat-only claims imply a non-standard reaction channel.

Everything here is textbook nuclear physics (Gamow tunnelling, astrophysical
S-factor, bound-pair rate formula lambda = A*|psi(0)|^2). Constants in CGS/eV.
Run:  python3 sim/m0_rate_budget.py   (writes docs/models/figs/m0_*.png)
"""
import os
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(__file__), "..", "docs", "models", "figs")
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------- constants
alpha = 1 / 137.035999
hbarc_eV_A = 1973.2698          # eV * Angstrom
e2_eV_A = 14.399645             # e^2/(4 pi eps0) in eV * Angstrom
c_cm_s = 2.99792458e10
m_d_c2 = 1875.612e6             # eV
mu_c2 = m_d_c2 / 2              # reduced mass, D+D  (eV)
J_per_eV = 1.602176634e-19

# Gamow energy E_G = 2 mu c^2 (pi alpha Z1 Z2)^2
E_G = 2 * mu_c2 * (np.pi * alpha) ** 2          # eV  (~0.986 MeV)

# Low-energy S-factors (keV b) for the two main D+D branches
S_n = 55e3        # D(d,n)3He   eV*b  (~53-55 keV b at E->0)
S_p = 57e3        # D(d,p)T     eV*b
S_tot_eV_cm2 = (S_n + S_p) * 1e-24               # eV cm^2
Q_n, Q_p = 3.269e6, 4.033e6                      # eV
E_mean_per_fusion = (S_n * Q_n + S_p * Q_p) / (S_n + S_p)   # eV, standard branching

# Bound-pair reaction constant: lambda = A * rho(0),  A = S c / (pi alpha mu c^2)
A_cm3_s = S_tot_eV_cm2 * c_cm_s / (np.pi * alpha * mu_c2)

n_D_PdD = 6.8e22   # deuterons / cm^3 in PdD_x at x ~ 1


def gamow_factor(E_eff):
    """Barrier penetration exp(-sqrt(E_G/E)) for an effective energy E_eff (eV)."""
    return np.exp(-np.sqrt(E_G / E_eff))


# -------------------------------------------------------- WKB cross-check
def wkb_yukawa(Ue, E=0.025, Rn_A=1e-5):
    """
    Exact WKB exponent for a Yukawa-screened Coulomb potential
    V(r) = (e^2/r) exp(-r/lam),  lam = e^2/Ue (so V ~ e^2/r - Ue at r << lam).
    Returns penetration factor exp(-2 * int kappa dr) from r=Rn to turning point.
    """
    lam = e2_eV_A / Ue
    V = lambda r: e2_eV_A / r * np.exp(-r / lam)
    r_tp = brentq(lambda r: V(r) - E, Rn_A, 1e4 * lam)
    k = lambda r: np.sqrt(max(2 * mu_c2 * (V(r) - E), 0.0)) / hbarc_eV_A
    # split integral for accuracy near the 1/sqrt(r) singularity
    pts = np.geomspace(Rn_A, r_tp, 60)
    G = sum(quad(k, a, b, limit=200)[0] for a, b in zip(pts[:-1], pts[1:]))
    return np.exp(-2 * G), r_tp


# ----------------------------------------------------- rate models
def rate_free_gas(Ue, kT=0.0257, n=n_D_PdD):
    """
    'Optimistic' model: deuterons behave as a free gas at density n that can
    approach freely to the screening radius; sigma_s(E) = S/E exp(-sqrt(E_G/(E+Ue))).
    Maxwell average in the E << Ue limit: <E^-1/2> = 2/sqrt(pi kT).
    Returns fusions / s per deuteron (pair-counting factor 1/2 included).
    """
    sqrt_2_over_mu = np.sqrt(2 / mu_c2) * c_cm_s        # cm/s per sqrt(eV)
    sigv = S_tot_eV_cm2 * sqrt_2_over_mu * 2 / np.sqrt(np.pi * kT) * gamow_factor(Ue)
    return 0.5 * n * sigv


def rate_bound_pair(Ue_eff, rho0_cm3=1e25):
    """
    Bound-pair model: lambda = A * rho0 * exp(-sqrt(E_G/Ue_eff)).
    Ue_eff lumps together electron screening AND lattice confinement.
    rho0 ~ |psi|^2 at the start of the barrier; 1e24-1e26 cm^-3 spans
    lattice-site to molecular confinement.  Returns fusions / s per pair.
    """
    return A_cm3_s * rho0_cm3 * gamow_factor(Ue_eff)


def Ue_required(target_rate_per_pair, rho0_cm3=1e25):
    x = np.log(A_cm3_s * rho0_cm3 / target_rate_per_pair)
    return E_G / x ** 2


def main():
    lines = []
    P = lambda s="": (print(s), lines.append(s))

    P(f"Gamow energy E_G (D+D)              = {E_G/1e3:.1f} keV")
    P(f"Bound-pair constant A               = {A_cm3_s:.2e} cm^3/s")
    P(f"Mean energy per fusion (std branch) = {E_mean_per_fusion/1e6:.2f} MeV")

    # ---- calibration: D2 molecule (Koonin & Nauenberg 1989: ~3e-64 /s)
    hw = 0.371          # D2 vibrational quantum, eV
    x0 = hbarc_eV_A / np.sqrt(mu_c2 * hw)
    rho_D2 = (np.pi * x0 ** 2) ** -1.5 * 1e24      # cm^-3
    Ue_D2 = Ue_required(3e-64, rho_D2)
    P(f"D2 calibration: x0={x0:.3f} A, rho0={rho_D2:.1e} cm^-3 -> Ue_eff(D2) = {Ue_D2:.1f} eV "
      "(reproduces 3e-64 /s)")

    # ---- WKB cross-check of exp(-sqrt(E_G/Ue)) for a Yukawa potential
    P("\nWKB (Yukawa) vs constant-shift Gamow factor:")
    P("  Ue(eV)  lambda_s(pm)  r_tp(pm)  WKB-Yukawa   exp(-sqrt(EG/Ue))")
    for Ue in [30, 100, 300, 1000]:
        pw, rtp = wkb_yukawa(Ue)
        P(f"  {Ue:6.0f}  {e2_eV_A/Ue*100:10.2f}  {rtp*100:8.1f}  {pw:10.2e}   {gamow_factor(Ue):10.2e}")

    # ---- mapping: accelerator-style Yukawa Ue -> effective Ue at thermal energy
    P("\nYukawa Ue (accelerator-style) -> Ue_eff seen by a thermal pair (same penetration):")
    for Ue in [100, 300, 500, 800, 1000, 2000]:
        pw, _ = wkb_yukawa(Ue)
        P(f"  Ue={Ue:5d} eV -> Ue_eff = {E_G/np.log(pw)**2:6.0f} eV")

    # ---- detection thresholds (reactions / s needed)
    P("\nReactions/s needed for a 5-sigma detection:")
    def rate_needed(eff, bkg_cps, T_s, frac):
        # S*eff*frac*T >= 5*sqrt(bkg*T)  (background-dominated Gaussian limit)
        return 5 * np.sqrt(bkg_cps * T_s) / (eff * frac * T_s)
    T = 14 * 86400
    r_n = rate_needed(0.10, 0.05, T, 0.5)    # 3He bank, 10% eff, 0.05 cps, n-branch 50%
    r_p = rate_needed(0.05 * 0.3, 2e-5, T, 0.5)  # Si in vacuum: 5% solid angle, 30% escape, 3 MeV p window
    heat_W = 0.010                           # 10 mW best-case calorimetric resolution
    r_h = heat_W / (E_mean_per_fusion * J_per_eV)
    P(f"  neutrons (3He bank, eff 10%, bkg 0.05 cps, 14 d) : {r_n:.2e} fusions/s")
    P(f"  3 MeV protons (Si, 5% x 30% escape, bkg 2e-5 cps): {r_p:.2e} fusions/s")
    P(f"  heat (10 mW, standard D+D branching)             : {r_h:.2e} fusions/s")
    P(f"  -> heat is {r_h/r_n:.1e}x less sensitive than neutrons for standard D+D")

    # ---- neutron dose implied if watt-level heat came from standard D+D
    fus_per_W = 1 / (E_mean_per_fusion * J_per_eV)
    n_per_W = fus_per_W * S_n / (S_n + S_p)
    flux_1m = n_per_W / (4 * np.pi * 100 ** 2)
    h10 = 400e-12   # Sv cm^2, ICRP-74 H*(10) coefficient ~2.5 MeV neutrons
    P(f"\n1 W of standard D+D -> {n_per_W:.2e} n/s; at 1 m: {flux_1m:.2e} n/cm^2/s "
      f"-> {flux_1m*h10*3600:.1f} Sv/h  (lethal; so W-level heat claims REQUIRE a non-standard channel)")

    # ---- required Ue_eff vs number of active pairs
    P("\nRequired Ue_eff (eV) for the neutron threshold, vs number of active pairs N and rho0:")
    Ns = [1e9, 1e12, 1e15, 1e18, 1e21]
    P("  N_pairs   rho0=1e24   rho0=1e25   rho0=1e26")
    for N in Ns:
        vals = [Ue_required(r_n / N, rho) for rho in (1e24, 1e25, 1e26)]
        P(f"  {N:7.0e}   " + "   ".join(f"{v:8.0f}" for v in vals))

    # ---- free-gas model and null-limit constraint
    P("\nFree-gas model rate per deuteron in PdD (upper-bound picture):")
    for Ue in [25, 50, 100, 150, 200, 300, 800]:
        P(f"  Ue={Ue:4d} eV -> {rate_free_gas(Ue):.2e} /s per D")
    lim = 1e-25    # order of 1989-90 null limits (fusions per pair per s) - verify in R5
    Ue_lim = brentq(lambda U: rate_free_gas(U) - lim, 10, 1000)
    P(f"  A null limit of {lim:.0e} /pair/s implies Ue_eff(thermal) < {Ue_lim:.0f} eV in bulk PdD "
      "(free-gas picture) -> accelerator Ue~300-800 eV does NOT apply at thermal energies")

    # ---- exponent budget: what engineering can buy vs what physics must supply
    G_D2 = np.sqrt(E_G / Ue_D2)
    G_need_hi = np.sqrt(E_G / Ue_required(r_n / 1e21, 1e25))
    G_need_lo = np.sqrt(E_G / Ue_required(r_n / 1e12, 1e25))
    P("\nExponent budget (penetration = exp(-X)):")
    P(f"  molecular D2 today        X = {G_D2:5.1f}")
    P(f"  needed, N=1e21 pairs      X = {G_need_hi:5.1f}   gap = {G_D2-G_need_hi:5.1f} e-folds")
    P(f"  needed, N=1e12 pairs      X = {G_need_lo:5.1f}   gap = {G_D2-G_need_lo:5.1f} e-folds")
    P(f"  engineering levers: 1e9x more sites buys {np.log(1e9):.1f} e-folds; 100x better detection buys {np.log(100):.1f}")
    P("  -> geometry/engineering cannot close the gap alone; the design must maximise sensitivity to a")
    P("     real anomalous enhancement (if one exists) at the sites where theories/claims locate it.")

    # ---- steepness
    P("\nLocal steepness d ln(rate)/d ln(Ue) = 0.5*sqrt(E_G/Ue):")
    for Ue in [50, 100, 200, 300]:
        s = 0.5 * np.sqrt(E_G / Ue)
        P(f"  Ue={Ue:3d} eV: exponent {s:5.1f}  (+10% Ue -> x{1.1**s:.0f} rate)")

    # ---- tail dominance: sites with Gaussian-distributed Ue
    P("\nTail dominance: Ue ~ Normal(mu, sigma) across sites; share of total rate from top 0.1% of sites")
    rng = np.random.default_rng(1)
    for mu, sig in [(100, 5), (100, 10), (100, 20), (200, 20)]:
        U = np.clip(rng.normal(mu, sig, 2_000_000), 1, None)
        w = gamow_factor(U)
        k = int(0.001 * U.size)
        top = np.sort(w)[-k:].sum() / w.sum()
        P(f"  mu={mu}, sigma={sig}: top 0.1% of sites give {100*top:5.1f}% of all reactions")

    # ---- figures
    U = np.logspace(np.log10(20), np.log10(2000), 300)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for rho in (1e24, 1e25, 1e26):
        ax.loglog(U, rate_bound_pair(U, rho), label=f"bound pair, rho0={rho:.0e} cm$^{{-3}}$")
    ax.loglog(U, rate_free_gas(U), "k--", label="free gas in PdD (per D)")
    for N, ls in [(1e21, ":"), (1e15, "-."), (1e12, "--")]:
        ax.axhline(r_n / N, color="grey", ls=ls, lw=1)
        ax.text(22, r_n / N * 3, f"neutron-detectable, N={N:.0e}", fontsize=7, color="grey")
    ax.axvline(Ue_D2, color="C3", lw=1)
    ax.text(Ue_D2 * 1.05, 1e-60, "D$_2$ molecule", color="C3", fontsize=8, rotation=90)
    # accelerator Ue = 300-800 eV mapped to Ue_eff for freely approaching pairs (Yukawa shape)
    lo = E_G / np.log(wkb_yukawa(300)[0]) ** 2
    hi = E_G / np.log(wkb_yukawa(800)[0]) ** 2
    ax.axvspan(lo, hi, color="C2", alpha=0.12)
    ax.text(lo * 1.03, 1e-72, "Ue_eff IF accelerator\nUe=300-800 eV (Pd)\nacted on freely\napproaching pairs",
            fontsize=7, color="C2")
    ax.set_ylim(1e-90, 1e5)
    ax.set_xlabel("effective screening energy Ue_eff (eV)")
    ax.set_ylabel("fusions / s per pair")
    ax.set_title("M0: D-D rate vs effective barrier reduction")
    ax.legend(fontsize=7, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "m0_rate_vs_Ue.png"), dpi=140)

    with open(os.path.join(OUT, "m0_output.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
