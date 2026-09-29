"""
M1 Part D - deuterium permeation flux as a variable at the site level.

1. Bulk: flux is a ~1e-9 perturbation of site occupancies (drift velocity << hop velocity).
2. Exit face (vacuum side of a C3 membrane): steady state of arrival flux J and second-order
   recombinative desorption  J = 2 k2 (theta N_s)^2,  k2 = nu2 exp(-E_d/kT).
   -> coverage theta(J, T), the recombination-limited maximum flux J_max(T), and
      encounter rates per cm^2 per s:
        (a) recombination events (close D-D approach through the D2 transition state) = J/2
        (b) nearest-neighbour pair formations by surface hopping
        (c) transient subsurface occupancy
3. Fusion probability per encounter vs U_eff, and the minimum per-encounter probability that
   C3's Si channel would detect at 5 sigma in 30 d, as a function of J.
4. Lock-in on flux modulation: what modulation buys.

Run: python3 sim/m1_flux.py -> docs/models/figs/m1_flux.txt, m1_flux.png
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import m1_physics as P
import m1_ladder as L

lines = []


def Pr(s=""):
    print(s)
    lines.append(s)


# Surface kinetics of D on Pd(111)-like exit face [BK]:
NS = 1.53e15            # sites/cm^2
NU2 = 1e-2              # cm^2/s, second-order desorption prefactor (Christmann 1988; typical 1e-3..1e-1)
ED = (0.80, 0.90)       # eV, recombinative desorption barrier range for D/Pd (TPD peak ~300-350 K)
E_HOP = 0.15            # eV, surface diffusion barrier (H/Pd(111) ~0.1-0.15 eV)
NU_HOP = 1e13           # s^-1
E_SS = 0.35             # eV, subsurface -> surface barrier (DFT, 0.3-0.5)
TAU_TS = 1e-14          # s, dwell time of the pair at the D2 transition state (~1 vibrational period)
D_PD = 2e-7             # cm^2/s, D diffusivity in Pd at 300 K (H: 3.8e-7; D ~ H/sqrt(2)... ) [BK Völkl & Alefeld]


def coverage(J, T, Ed):
    k2 = NU2 * np.exp(-Ed / (P.kB * T))
    th = np.sqrt(J / (2 * k2)) / NS
    return np.minimum(th, 1.0), 2 * k2 * NS ** 2         # theta, J_max


def encounters(J, T, Ed):
    th, Jmax = coverage(J, T, Ed)
    Jeff = np.minimum(J, Jmax)
    rec = Jeff / 2                                        # recombination events /cm^2/s
    hop = NU_HOP * np.exp(-E_HOP / (P.kB * T))
    nn = NS * th ** 2 * (1 - th) * 3 * hop                # new NN pairs formed by hops, /cm^2/s
    tau_ss = 1 / (NU_HOP * np.exp(-E_SS / (P.kB * T)))
    n_ss = Jeff * tau_ss                                  # transient subsurface D per cm^2
    return dict(theta=th, Jmax=Jmax, rec=rec, nn=nn, n_ss=n_ss, tau_ss=tau_ss)


def main():
    Pr("M1 Part D - deuterium flux as a variable")
    # 1. bulk
    for J in (1e14, 1e16):
        ratio = J * 2.85e-8 / (6.1e22 * D_PD)
        Pr(f"  Bulk: J={J:.0e} D/cm2/s -> drift/hop-flux ratio J*a/(n D) = {ratio:.1e}: flux does not change bulk"
           " site occupancies (O, T, vacancy, dislocation) measurably.")
    Pr("  So any flux effect must live at the exit face, entry face or buried interfaces, where J sets coverage"
       " and the rate of recombination encounters.")
    # 2. exit face
    Pr("\n  Exit-face kinetics (Pd, vacuum side). J_max = recombination-limited flux (theta=1):")
    for T in (300, 350, 400, 500):
        Pr("   T=%d K: J_max = %s D/cm2/s" % (T, " .. ".join(f"{coverage(1, T, e)[1]:.1e}" for e in ED[::-1])))
    j3 = [coverage(1, 300, e)[1] for e in ED[::-1]]
    Pr(f"  -> at room temperature a bare Pd exit face into vacuum saturates at ~{j3[0]:.0e}-{j3[1]:.0e} D/cm2/s; above that"
       " D piles up subsurface and the permeation flux is recombination-limited.")
    Pr("     To run flux as a variable up to 1e16 D/cm2/s the exit face must be warm (>=350-400 K) or carry"
       " a recombination catalyst - but a Pt/Pd-black catalyst layer changes the very site population being"
       " tested. PdO also blocks recombination (lower J_max) - design conflict flagged for M3/M6/M7.")
    Pr("  Energetics: on Pd, 2 D(ads) -> D2(g) is ENDOTHERMIC by ~0.9 eV (adsorption heat, Christmann 1988)"
       " and D(bulk) -> 1/2 D2(g) is endothermic by ~0.2 eV/D; the 0.5-1 eV 'released per molecule' in the"
       " brief applies only to gas-phase atom recombination. Nascent D2 leaves near-thermal: no hot pairs.")

    Js = np.logspace(10, 17, 200)
    Pr("\n  Encounter rates at the exit face (Ed=0.85 eV):")
    Pr("   T(K)   J(D/cm2/s)  theta    recomb/cm2/s   NN-pair form./cm2/s   subsurf transient D/cm2")
    for T in (300, 450):
        for J in (1e12, 1e14, 1e16):
            e = encounters(J, T, 0.85)
            Pr(f"   {T:4d}  {J:9.0e}   {e['theta']:6.3f}   {e['rec']:10.2e}     {e['nn']:12.2e}         {e['n_ss']:10.2e}")
    Pr("  Static comparison: the exit face permanently holds 3*1.5e15 = 4.5e15 NN surface pairs/cm2.")
    Pr(f"  Flux-driven transition-state pair-time per cm2 = (J/2)*tau_TS = {1e16/2*TAU_TS:.0f} pair-s per s at"
       " J=1e16: flux adds the equivalent of ~50 permanent molecular-distance pairs per cm2.")

    # 3. fusion probability per encounter
    Pr("\n  Fusion probability per recombination encounter p = lambda_mol(U_eff) * tau_TS (rho0 of D2):")
    for U, lab in ((34, "D2-molecule (standard)"), (127, "accel. low, Pd"), (209, "accel. high, Pd surface"),
                   (247, "accel. PdO")):
        Pr(f"   U_eff={U:4d} eV ({lab:24s}): p = {P.rate_pair(U, 1.5e26)*TAU_TS:.1e}")
    rates = {s["k"]: L.S.site_models(s) for s in L.S.SITES}
    b5 = L.NT.part_B5()
    C3 = L.make_configs(b5)["C3"]
    ch = L.channel_thresholds(C3)["cp"]
    eps0 = L.eps_charged(0.0, C3["f_omega"])
    A = 2.0
    Rmin = ch["smin"] / (P.T30 * eps0)
    Pr(f"  C3 Si channel: s_min = {ch['smin']:.0f} counts/30 d, eps(surface) = {eps0:.3f} -> R_min = {Rmin:.1e}"
       " fusions/s")
    for T in (300, 450):
        for J in (1e12, 1e14, 1e16):
            rec = encounters(J, T, 0.85)["rec"] * A
            pmin = Rmin / rec
            Pr(f"   T={T} K, J={J:.0e}: {rec:.1e} encounters/s -> p_min = {pmin:.1e} per encounter"
               f" (x{pmin/(P.rate_pair(34,1.5e26)*TAU_TS):.0e} over D2-standard, U_eff-equiv"
               f" {P.Ue_required(pmin/TAU_TS, 1.5e26):.0f} eV)")
    Pr("  -> a flux-triggered anomaly is detectable only if a recombination encounter fuses with p >~ 1e-18..-16;"
       " static D2-molecule physics gives 3e-78. That is the M7 lock-in target.")

    # 4. lock-in
    Pr("\n  Flux modulation (M7): with square-wave on/off at 50 % duty, signal S appears only in 'on' bins;"
       " background (radon, cosmic, detector) is common-mode. Statistics are those of on-minus-off counting:"
       " s_min rises by sqrt(2) (~1.4x) vs a known-background count, but slow background systematics"
       " (barometric +-14 %, radon plate-out) cancel to first order if the period << their timescale (hours):"
       " choose 1-10 min periods (membrane diffusion time L^2/D = "
       f"{(25e-4)**2/D_PD:.0f} s for 25 um sets the upper frequency ~1/(2*that)).")

    # figure: two panels sharing the flux axis (no dual y-axis)
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(9.0, 3.9), sharex=True)
    for i, T in enumerate((300, 350, 400)):
        e = encounters(Js, T, 0.85)
        ax.loglog(Js, e["rec"], color=L.NT.COL[i], lw=2, label=f"{T} K")
        ax2.semilogx(Js, e["theta"], color=L.NT.COL[i], lw=2, label=f"{T} K")
    ax.set_xlabel("applied permeation flux J (D cm$^{-2}$ s$^{-1}$)")
    ax.set_ylabel("recombination encounters cm$^{-2}$ s$^{-1}$")
    ax.set_title("exit-face D-D recombination encounters", fontsize=9)
    ax2.set_xlabel("applied permeation flux J (D cm$^{-2}$ s$^{-1}$)")
    ax2.set_ylabel("exit-face D coverage $\\theta$")
    ax2.set_title("coverage; $\\theta$=1 means recombination-limited", fontsize=9)
    for a_ in (ax, ax2):
        a_.grid(color="#e6e6e6", lw=0.6)
        a_.legend(fontsize=7, frameon=False, loc="upper left")
    fig.suptitle("M1-D: flux at the exit face of a Pd membrane into vacuum (E$_d$ = 0.85 eV)", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(P.OUT, "m1_flux.png"), dpi=140)
    with open(os.path.join(P.OUT, "m1_flux.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
