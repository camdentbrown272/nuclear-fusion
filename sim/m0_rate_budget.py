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
Run:  python3 sim/m0_rate_budget.py   (writes docs/models/figs/m0_*; set M0_OUT=<dir> to write elsewhere)

Revision 2 (after red-team-0): D2 shell density, WKB from r=0, Bosch-Hale S(0),
matched-control thresholds, 30-day runs, ln-rate gap bookkeeping, 4He channel.
"""
import os
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.environ.get("M0_OUT", os.path.join(os.path.dirname(__file__), "..", "docs", "models", "figs"))
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
S_n = 53.7e3      # D(d,n)3He   eV*b  (Bosch-Hale 1992 S(0))
S_p = 55.6e3      # D(d,p)T     eV*b  (Bosch-Hale 1992 S(0))
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
def wkb_yukawa(Ue, E=0.025):
    """
    Exact WKB exponent for a Yukawa-screened Coulomb potential
    V(r) = (e^2/r) exp(-r/lam),  lam = e^2/Ue (so V ~ e^2/r - Ue at r << lam).
    Returns penetration factor exp(-2 * int_0^r_tp kappa dr).  The integral starts at
    r = 0 to stay consistent with the point-Coulomb Gamow factor that defines S(E);
    the substitution r = r_tp * s^2 removes both endpoint singularities.
    """
    lam = e2_eV_A / Ue
    V = lambda r: e2_eV_A / r * np.exp(-r / lam)
    r_tp = brentq(lambda r: V(r) - E, 1e-9, 1e4 * lam)
    def integrand(s):
        r = r_tp * s * s
        if r <= 0:
            return 2 * r_tp * np.sqrt(2 * mu_c2 * e2_eV_A / r_tp) / hbarc_eV_A
        return np.sqrt(max(2 * mu_c2 * (V(r) - E), 0.0)) / hbarc_eV_A * 2 * r_tp * s
    G = quad(integrand, 0, 1, limit=400, points=[1e-4, 1e-2, 0.1, 0.5, 0.9])[0]
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
    R_e = 0.7414        # D2 bond length, Angstrom
    # relative wavefunction is a radial shell at R_e with Gaussian width x0 (not centred at r=0)
    rho_D2 = 1 / (np.sqrt(np.pi) * x0 * 4 * np.pi * R_e ** 2) * 1e24   # cm^-3
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
        # Equal-time matched control (charter): S*eff*frac*T >= 5*sqrt(2*bkg*T).
        # Systematic floor ignored here: for neutrons a 1-2 % background drift sets
        # S >= 0.05-0.1 fusions/s regardless of T (R5 sec. 10.1).
        return 5 * np.sqrt(2 * bkg_cps * T_s) / (eff * frac * T_s)
    T = 30 * 86400
    r_n = rate_needed(0.10, 0.05, T, 0.5)    # 3He bank, 10% eff, 0.05 cps, n-branch 50%
    r_p = rate_needed(0.05 * 0.3, 2e-5, T, 0.5)  # Si in vacuum: 5% solid angle, 30% escape, 3 MeV p window
    heat_W = 0.010                           # 10 mW best-case calorimetric resolution
    r_h = heat_W / (E_mean_per_fusion * J_per_eV)
    He_limit_atoms = 1e10                    # 4He extraction/static-MS limit (red-team-0 sec. 2.2)
    r_he = He_limit_atoms / T                # reactions/s if every 4He is collected
    r_h2 = heat_W / (23.85e6 * J_per_eV)     # 10 mW via D+D -> 4He (23.85 MeV each)
    P(f"  neutrons (3He bank, eff 10%, bkg 0.05 cps, 30 d) : {r_n:.2e} fusions/s")
    P(f"  3 MeV protons (Si, 5% x 30% escape, bkg 2e-5 cps): {r_p:.2e} fusions/s")
    P(f"  heat (10 mW, standard D+D branching)             : {r_h:.2e} fusions/s")
    P(f"  -> heat is {r_h/r_n:.1e}x less sensitive than neutrons for standard D+D")
    P(f"  H2 channel: 10 mW of D+D->4He                    : {r_h2:.2e} reactions/s")
    P(f"  H2 channel: 4He, 1e10 atoms collected in 30 d    : {r_he:.2e} reactions/s "
      f"(= {r_he*23.85e6*J_per_eV*1e9:.0f} nW; x{r_h2/r_he:.0e} better than the calorimeter)")

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
    Ue_acc_lim = brentq(lambda U: np.sqrt(E_G / Ue_lim) + np.log(wkb_yukawa(U)[0]), 150, 2000)
    P(f"  A null limit of {lim:.0e} /pair/s implies Ue_eff(thermal) < {Ue_lim:.0f} eV in bulk PdD "
      f"(free-gas picture, S/E convention)")
    P(f"  -> equivalent Yukawa (accelerator-style) Ue < {Ue_acc_lim:.0f} eV: Tohoku's Pd ~310 eV is allowed,")
    P("     Bochum-class 500-800 eV cannot act on thermal pairs in bulk PdD")

    # ---- exponent budget: what engineering can buy vs what physics must supply
    gap_hi = np.log((r_n / 1e21) / 3e-64)
    gap_lo = np.log((r_n / 1e12) / 3e-64)
    P("\nRate budget in e-folds (ln of required rate per pair / molecular-D2 rate 3e-64 /s):")
    P(f"  molecular D2 exponent X(D2) = {np.sqrt(E_G / Ue_D2):5.1f}")
    P(f"  gap if all 1e21 pairs of a 0.5 g cathode are active : {gap_hi:5.1f} e-folds")
    P(f"  gap if only 1e12 special sites are active           : {gap_lo:5.1f} e-folds")
    P(f"  with ALL 1e21 pairs active AND 100x better detection: {gap_hi-np.log(100):5.1f} e-folds remain")
    P("  -> engineering cannot close the gap; a positive result needs a real anomaly. The design")
    P("     must reproduce the conditions where anomalies are claimed and detect them credibly.")

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
