"""
M1 Part B - non-thermal ("cold-compatible") energy sources for a D-D pair.

B0  How much pair kinetic energy is worth: e-folds of barrier per eV.
B1  Fracto-emission / crack-tip charge separation (PdD, oxide-skinned PdD/TiD2, LiD).
B2  Desorption / phase-transition transients (Lipson Pd/PdO claims).
B3  Electrochemical double-layer fields.
B4  Optical / plasmonic near fields.
B5  Cosmic-ray floor: knock-on d-d in flight, D(n,2n), D(gamma,n), muon catalysis,
    muon capture - for 1 cm3 PdD + 50 mL D2O at sea level.
B6  Czerski near-threshold 0+ resonance: S-factor multiplier K(E) (Breit-Wigner).
B7  Thermal-spike ("plateau") mechanism: what a beam-free device would need.

Run: python3 sim/m1_nonthermal.py -> docs/models/figs/m1_nonthermal.txt, m1_*.png
"""
import os
import warnings
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import m1_physics as P

warnings.filterwarnings("ignore")
e_C = 1.602176634e-19
eps0 = 8.8541878e-12
SIGMA_CX = 1e-15          # cm^2, D+ in D2 (charge exchange + elastic, 10 eV-1 keV) [BK, Phelps 1990 compilation]
COL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7"]

lines = []


def Pr(s=""):
    print(s)
    lines.append(s)


# ------------------------------------------------------------------ B0
def efolds_per_eV(U, E=0.0):
    """d X / d E where X = sqrt(E_G/(E+U)): e-folds of barrier gained per eV of pair energy."""
    return 0.5 * np.sqrt(P.E_G) * (E + U) ** -1.5


def rate_multiplier(dE, U):
    """Penetration ratio for adding dE (eV) of relative energy on top of screening U (constant-shift)."""
    return np.exp(np.sqrt(P.E_G / U) - np.sqrt(P.E_G / (U + dE)))


def part_B0():
    Pr("B0. Value of pair kinetic energy: e-folds per eV, and rate multiplier for +dE")
    Pr("    U_eff (eV) | e-folds/eV at E=0 | x rate for +1 eV | +10 eV | +100 eV | +1 keV")
    for U in (11, 34, 127, 250):
        Pr(f"    {U:6.0f}     | {efolds_per_eV(U):8.2f}          | {rate_multiplier(1, U):9.2e} |"
           f" {rate_multiplier(10, U):9.2e} | {rate_multiplier(100, U):9.2e} | {rate_multiplier(1000, U):9.2e}")
    Pr("    -> 1-2 eV transients matter only in relative terms; the absolute rate stays set by U_eff.")


# ------------------------------------------------------------------ B1
def fracto_yield_per_cm2(E_field, w_max_m, host, eta_D=1.0, gap_gas_cm3=0.0, Ue=0.0, branch="n"):
    """
    Fusion yield per cm^2 of freshly separated crack face.
      surface charge   sigma = eps0*E_field                 (C/m^2)
      D+ per cm^2      N = eta_D * sigma / e * 1e-4
      max ion energy   V_max = E_field * w_max              (vacuum gap: uniform spectrum on [0, eV_max])
      with gas (n cm^-3): collisional, exponential spectrum with mean E*lambda, capped at V_max.
    Returns (N_ions, mean energy eV, fusions per cm^2).
    """
    N = eta_D * eps0 * E_field / e_C * 1e-4
    Vmax = E_field * w_max_m
    Es = np.linspace(1.0, max(Vmax, 1.01), 200)
    if gap_gas_cm3 > 0:
        lam_m = 1 / (gap_gas_cm3 * SIGMA_CX) * 1e-2
        eps_l = E_field * lam_m
        w = np.exp(-Es / eps_l)
        w /= np.trapezoid(w, Es)
    else:
        w = np.ones_like(Es) / (Es[-1] - Es[0])
    Y = np.array([P.thick_target_yield(E, host, Ue, branch) for E in Es])
    Ybar = np.trapezoid(w * Y, Es)
    Ebar = np.trapezoid(w * Es, Es)
    return N, Ebar, N * Ybar


def part_B1():
    Pr("\nB1. Fracto-emission / crack-tip charge separation")
    Pr("    Physics limits: in a metal (PdD, TiD2) crack-face charge relaxes in ~eps0/sigma_el ~ 1e-18 s;"
       " the only sustained gap voltage is the contact-potential difference (<~1 V) [BK].")
    Pr("    Oxide-skinned faces (PdO 2-10 nm, TiO2 2-5 nm): V <= E_bd*t_ox ~ 1e9 V/m x 1e-8 m = 10 V [BK, thin-oxide breakdown].")
    Pr("    Ionic insulator (LiD): surface fields 1e7-1e9 V/m (brief/literature), vacuum gap (no D2 gas),"
       " V grows with opening w until field emission/flashover; w_max = 1-10 um [guess].")
    Pr("    D2 in PdD cracks at x~0.9: fugacity 1e3-1e4 bar -> n ~ 2e22 cm^-3, mean free path "
       f"{1/(2e22*SIGMA_CX)*1e7:.1f} nm -> collisional; ions gain E*lambda only.")
    Pr("    Case | E-field V/m | w_max | D+/cm^2 | <E> eV | n per cm^2 crack (U=0) | (U=site U_eff, empirical)")
    # last column: screening = thermal-mapped U_eff of the target site (m1_sites, empirical high/low);
    # using the keV accelerator U_e at eV impact energies would double-count (const-shift form).
    cases = [("PdD metal, contact pot.", 1e8, 1e-8, "PdD0.9", 2e22, 127),
             ("PdO-skinned PdD", 1e9, 1e-8, "PdD0.9", 0.0, 247),
             ("TiO2-skinned TiD2", 1e9, 5e-9, "TiD2", 0.0, 124),
             ("LiD, weak", 1e7, 1e-5, "LiD", 0.0, 25),
             ("LiD, typical", 1e8, 1e-5, "LiD", 0.0, 25),
             ("LiD, extreme", 1e9, 1e-5, "LiD", 0.0, 25),
             ("LiD, extreme, 1 um", 1e9, 1e-6, "LiD", 0.0, 25)]
    res = {}
    for name, E, w, host, gas, Ue in cases:
        N, Eb, Y0 = fracto_yield_per_cm2(E, w, host, 1.0, gas, 0.0)
        _, _, YU = fracto_yield_per_cm2(E, w, host, 1.0, gas, Ue)
        res[name] = (N, Eb, Y0, YU)
        Pr(f"    {name:24s} | {E:8.0e} | {w*1e6:6.2f} um | {N:8.1e} | {Eb:8.1f} | {Y0:9.1e} | {YU:9.1e}")
    # comparison with claims
    Pr("    Claims [BK]: Menlove et al. (J. Fusion Energy 9, 495, 1990) - Ti/D2 chips, bursts of up to ~1e2 counts"
       " (eff ~30%) during warm-up from 77 K, i.e. ~1e2-1e3 n per burst; not reproduced (R5).")
    Pr("    Klyuev/Lipson/Derjaguin (Sov. Tech. Phys. Lett. 1986; Nature 341, 492, 1989) - LiD / D2O-ice fracture,"
       " reported ~1e-2-1 n/s during vibro-milling; single group.")
    _, _, Yx, _ = res["LiD, extreme"]
    _, _, Yt, _ = res["LiD, typical"]
    Pr(f"    -> a 300-n burst needs {300/Yx:.1e} cm^2 of fresh LiD face at 1e9 V/m & 10 um opening, or "
       f"{300/max(Yt,1e-300):.1e} cm^2 at 1e8 V/m. For TiD2 (metal+oxide): {300/max(res['TiO2-skinned TiD2'][2],1e-300):.1e} cm^2.")
    Pr("    -> Fracto-fusion is physically possible ONLY in insulating hydrides (LiD) or thick insulating"
       " inclusions; in metallic PdD/TiD2 the achievable deuteron energy (<=10 eV) gives zero yield."
       " Menlove-type Ti bursts cannot be conventional fracto-fusion; they match cosmic spallation"
       " multiplets (R5 #8).")
    # per cycle, per configuration-relevant numbers
    Pr("    Per alpha/beta cycle of PdD (crack area 10 cm^2 per cm^3 per cycle [guess]):"
       f" {10*res['PdO-skinned PdD'][3]:.1e} n (oxide-skinned, screened), {10*res['PdD metal, contact pot.'][3]:.1e} n (bare, screened).")
    # figure
    Efs = np.logspace(6.5, 9.3, 18)
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for i, (w, lab) in enumerate([(1e-5, "LiD, opening 10 $\\mu$m"), (1e-6, "LiD, opening 1 $\\mu$m"),
                                  (1e-7, "LiD, opening 100 nm")]):
        Y = [fracto_yield_per_cm2(E, w, "LiD")[2] for E in Efs]
        ax.loglog(Efs, np.maximum(Y, 1e-40), color=COL[i], lw=2, label=lab)
    Y = [fracto_yield_per_cm2(E, 1e-8, "PdD0.9", Ue=247)[2] for E in Efs]
    ax.loglog(Efs, np.maximum(Y, 1e-40), color=COL[3], lw=2, ls="--",
              label="PdO-skinned PdD (10 nm oxide), U$_{eff}$=247 eV")
    ax.axhline(1, color="grey", lw=0.8, ls=":")
    ax.set_ylim(1e-40, 1e4)
    ax.set_xlabel("crack-face field (V/m)")
    ax.set_ylabel("D-D neutrons per cm$^2$ of fresh crack")
    ax.set_title("M1-B1: fracto-fusion yield (all charge carried by D$^+$)", fontsize=10)
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    ax.grid(color="#e6e6e6", lw=0.6)
    fig.tight_layout()
    fig.savefig(os.path.join(P.OUT, "m1_fracto.png"), dpi=140)
    return res


# ------------------------------------------------------------------ B2
def part_B2():
    Pr("\nB2. Desorption / phase-transition transients (Lipson Pd/PdO:Dx)")
    srcs = [("PdO/Pd Schottky built-in potential", 0.8),
            ("D+D -> D2 recombination energy (brief: 0.5-1 eV)", 1.0),
            ("Heyrovsky step at 1 V overpotential", 1.0),
            ("beta->alpha transition strain energy per D", 0.05)]
    for n, dE in srcs:
        Pr(f"    {n:52s}: <= {dE:4.2f} eV -> x{rate_multiplier(dE, 127):.2g} at U_eff=127 eV,"
           f" x{rate_multiplier(dE, 11):.2g} at U_eff=11 eV (TF)")
    Pr("    On Pd the recombinative desorption is endothermic (adsorption heat ~0.9 eV/D2, Christmann 1988 [BK]),"
       " so no hot D is produced; the 1 eV entries are upper bounds.")
    # Lipson claim
    claim = 0.1            # 3 MeV protons s^-1 cm^-2 [BK, order of magnitude of Lipson et al. 2005 claims]
    pairs_int = 7.8e15     # PdO/Pd interface pairs per cm^2 (m1_sites)
    pairs_1um = 6 * 6.1e22 * 1e-4
    for lab, N in (("interface layer", pairs_int), ("top 1 um", pairs_1um)):
        lam = 2 * claim / N
        Pr(f"    Lipson-scale claim {claim} p/s/cm^2 from the {lab} ({N:.1e} pairs/cm^2) needs"
           f" lambda = {lam:.1e} /pair/s -> static U_eff >= {P.Ue_required(lam, 5.4e24):.0f} eV")
    fint = P.rate_pair(247.5, 5.4e24) * pairs_int
    Pr(f"    Static accelerator-PdO screening mapped to thermal (U_eff 247 eV, m1_sites) predicts"
       f" {fint:.1e} fusions/s/cm^2 = {fint/2:.1e} p/s/cm^2 at the interface: {claim/(fint/2):.0f}x below the"
       " claim, but inside C3's reach (see Part C). Transients add <x1.5.")


# ------------------------------------------------------------------ B3
def part_B3():
    Pr("\nB3. Electrochemical double layer")
    V, L = 1.0, 3e-10                      # V across Helmholtz layer, m
    Ef = V / L
    grad = Ef / L
    for d in (0.74e-10, 2.85e-10):
        dE = grad * d ** 2 / 2             # eV (field gradient energy across pair, q=e)
        Pr(f"    field {Ef:.1e} V/m, gradient {grad:.1e} V/m^2: pair (d={d*1e10:.2f} A) relative-energy shift"
           f" {dE:.3f} eV -> x{rate_multiplier(dE, 34):.3f} (U=34), x{rate_multiplier(dE, 11):.2f} (U=11)")
    Pr("    (the 2.85 A row assumes a pair straddling the 3 A layer along the normal; an adsorbed pair lies in the"
       " plane, where the gradient gives no relative energy.)")
    Pr("    A uniform field exerts no relative force on two deuterons (same q/m). Ion energy between collisions in"
       " the diffuse layer (1e8 V/m x 0.3 nm) = 0.03 eV. Negligible.")


# ------------------------------------------------------------------ B4
def part_B4():
    Pr("\nB4. Optical / plasmonic near fields (CW, charter-compliant low power)")
    c, me, md = 2.998e8, 9.109e-31, 3.344e-27
    om = 2 * np.pi * c / 800e-9
    for I, enh in ((1e7, 10), (1e9, 100), (1e9, 1000)):
        E0 = np.sqrt(2 * I / (c * eps0))
        El = E0 * enh
        Up_e = e_C ** 2 * El ** 2 / (4 * me * om ** 2) / e_C
        Up_d = Up_e * me / md
        grad_dE = El / 10e-9 * (2.85e-10) ** 2 / 2          # hot-spot scale 10 nm
        Pr(f"    I={I:.0e} W/m^2, |E/E0|={enh}: E_loc={El:.1e} V/m, U_p(e)={Up_e:.1e} eV,"
           f" U_p(D)={Up_d:.1e} eV, pair gradient shift {grad_dE:.1e} eV")
    Pr("    -> 1e10-1e11 V/m would be needed for eV-scale effects; that is the damage/field-evaporation regime."
       " Optical fields cannot deliver relative D-D energy in a cold device.")


# ------------------------------------------------------------------ B5
NBANDS = [  # (label, flux cm^-2 s^-1, representative E_n MeV, sigma_nd elastic b, sigma_n2n b, frac above 3.34 MeV)
    ("1-10 MeV", 3.5e-3, 3.0, 2.0, 0.10, 0.4),
    ("10-100 MeV", 2.5e-3, 30.0, 0.6, 0.15, 1.0),
    (">100 MeV", 1.1e-3, 300.0, 0.1, 0.05, 1.0),
]
# Sea-level flux split [guess] consistent with R5: total 1.3e-2, >10 MeV 3.6e-3 (Gordon et al. 2004).
# sigma(n,d) elastic, D(n,2n) from ENDF/B-VIII.0 shapes [BK, +-30%].


def knockon_fusions(nD, V_cm3, host):
    """d-d fusions/s from cosmic-neutron elastic knock-on deuterons slowing in `host`."""
    tot = 0.0
    for _, phi, En, s_el, _, _ in NBANDS:
        recoils = phi * nD * V_cm3 * s_el * 1e-24
        Ed = np.linspace(0.02, 8 / 9, 25) * En * 1e6          # recoil spectrum ~uniform (isotropic CM)
        Ed = np.minimum(Ed, 6e6)                               # Bosch-Hale validity (E_cm<5 MeV) -> cap
        Y = np.mean([P.thick_target_yield(E, host) for E in Ed])
        tot += recoils * Y
    return tot


def part_B5():
    Pr("\nB5. Cosmic-ray 'cold floor' at sea level: 1 cm3 PdD0.9 cathode + 50 mL D2O")
    vol = {"PdD": (6.1e22, 1.0, "PdD0.9", 12.0), "D2O": (6.66e22, 50.0, "D2O", 55.4)}
    out = {}
    for k, (nD, V, host, mass) in vol.items():
        f_ko = knockon_fusions(nD, V, host)
        n2n = sum(phi * nD * V * s2 * 1e-24 * fr for _, phi, _, _, s2, fr in NBANDS)
        phi_g = 0.2           # 2.614 MeV 208Tl line flux in a concrete building, gamma/cm2/s [BK, 0.05-0.5]
        s_gn = 1.0e-27        # D(gamma,n) at 2.614 MeV [BK, Chadwick-Goldhaber classic, +-30%]
        gn = phi_g * nD * V * s_gn
        stop = 1.5e-5 * mass  # stopped mu per s (sea-level stopping rate 1.5e-5 /g/s [calc from PDG spectrum, x2])
        mum = stop * 0.44     # mu- fraction (mu+/mu- = 1.27) [BK]
        if k == "D2O":
            p_d, p_fus = 0.01, 3e-5   # initial capture on D; ddmu formation vs transfer to O [guess/BK rates]
            capfrac, nmult = 0.18, 1.0
        else:
            p_d, p_fus = 1e-3, 1e-6
            capfrac, nmult = 0.97, 1.6
        mucf = mum * p_d * p_fus
        mucap_n = mum * capfrac * nmult
        out[k] = dict(knock=f_ko, n2n=n2n, gn=gn, mucf=mucf, mucap=mucap_n)
        Pr(f"    [{k}] D atoms {nD*V:.1e}")
        Pr(f"       knock-on d-d fusions in flight   : {f_ko:.1e} /s  (half give 2.45 MeV n, half 3.02 MeV p)")
        Pr(f"       D(n,2n) breakup neutrons         : {n2n:.1e} n/s   <- D-specific, not fusion")
        Pr(f"       D(gamma,n) from 208Tl 2.614 MeV  : {gn:.1e} n/s (0.2 MeV n)  <- D-specific, not fusion")
        Pr(f"       muon-catalysed fusion            : {mucf:.1e} /s  (stopped mu- {mum:.1e}/s)")
        Pr(f"       mu- capture neutrons             : {mucap_n:.1e} n/s  (same in H control)")
    tot_fus = sum(v["knock"] + v["mucf"] for v in out.values())
    tot_Dn = sum(v["n2n"] + v["gn"] for v in out.values())
    Pr(f"    TOTAL real d-d fusion floor: {tot_fus:.1e} /s  ->  {tot_fus*86400:.1e} per day (undetectable)")
    Pr(f"    TOTAL D-specific non-fusion neutrons: {tot_Dn:.1e} n/s = {tot_Dn*86400:.0f} per day"
       " -> this, not fusion, is what the D-vs-H control must cancel;")
    Pr("       it is ~10-50% of a 5-sigma/30-day neutron threshold (M0: 2e-2/s for 14 d). Keep D2O volume small,"
       " shield 2.6 MeV gammas (5 cm Pb/Cu inner shield cuts it ~x10), or use PSD to veto 0.2 MeV photoneutrons.")
    Pr("    Recoil deuterons from n-d scattering in the top 20 um of a 2 cm2 foil reach a facing Si detector at"
       f" ~{3.5e-3*6.1e22*2e-24*2*20e-4*0.2:.0e} cps, spread over 0-8 MeV: negligible vs 2e-5 cps in the p window.")
    return out


# ------------------------------------------------------------------ B6
def resonance_K(E_eV, E_R, C):
    """S-factor multiplier from a narrow 0+ resonance (energy-independent reduced widths)."""
    Gam = 0.44                                  # eV: Gamma_p=40 meV [R3, v] x (1 + Gamma_ee/Gamma_p=10)
    return 1 + C / ((E_eV - E_R) ** 2 + Gam ** 2 / 4)


def part_B6():
    Pr("\nB6. Czerski 0+ threshold resonance: S-factor multiplier K(E)")
    Pr("    Breit-Wigner with Gamma_d(E) = 2 theta^2 gamma_W^2 (2 pi R/a) e^{4 sqrt(2R/a)} e^{-2 pi eta(E)}:")
    Pr("    for E << E_R the Gamow factors cancel and K = 1 + C/((E-E_R)^2 + Gamma^2/4) is energy independent.")
    hb2 = 20.76e6 * 1e-26                       # hbar^2/2mu in eV cm^2 (20.76 MeV fm^2)
    R, a = 5.0, 28.8                            # fm; nuclear Bohr radius for d+d
    gW2 = 3 * 20.76e6 / R ** 2                  # Wigner limit, eV
    pref = np.pi * hb2 * (2 / 9) / 1e-24        # eV*b
    Gout = 0.44
    Cmax = pref * 2 * gW2 * (2 * np.pi * R / a) * np.exp(4 * np.sqrt(2 * R / a)) * Gout / (110e3)
    Pr(f"    Unitarity-style ceiling (theta^2=1, R=5 fm, Gamma_out=0.44 eV): C_max = {Cmax:.2e} eV^2")
    Pr("    Calibration: the resonance must supply the excess at E_cm=2.5 keV between U_e=340 eV and the"
       f" resonance-fit U_e=105 eV (Szczecin): x{np.exp(0.5*np.sqrt(P.E_G/2500)*(340-105)/2500):.2f}.")
    Kcal = np.exp(0.5 * np.sqrt(P.E_G / 2500) * (340 - 105) / 2500)
    rows = []
    Pr("    E_R (eV above threshold) | theta^2 needed | K(2.5 keV) | K(thermal 0.04 eV) | equiv. dU_eff at U=127 / 34 eV")
    for E_R in (5000, 3000, 2000, 1000, 300, 30, 0.0):
        C = (Kcal - 1) * ((2500 - E_R) ** 2 + Gout ** 2 / 4)
        th2 = C / Cmax
        Kth = resonance_K(0.04, E_R, C)
        dU = [P.E_G / (np.sqrt(P.E_G / U) - np.log(Kth)) ** 2 - U for U in (127, 34)]
        rows.append((E_R, th2, Kth))
        ok = "" if th2 <= 1 else "  (theta^2>1: unphysical)"
        Pr(f"    {E_R:8.0f}                 | {th2:9.2e}      | {Kcal:6.2f}     | {Kth:10.2e}          |"
           f" +{dU[0]:.0f} / +{dU[1]:.1f} eV{ok}")
    Pr(f"    -> If E_R sits within ~1 eV of threshold, K(thermal) ~ {rows[-1][2]:.0e} (= +{np.log(rows[-1][2]):.0f} e-folds)."
       " If E_R is a few keV above threshold, K(thermal) ~ 1.2-1.4. The data do not fix E_R.")
    Pr("    Even K=2e8 raises the D2-molecule rate only from 3e-64 to ~6e-56 /s: the resonance changes WHAT to"
       " detect, not WHETHER cold D-D is detectable (confirms R3 Table D).")
    Pr("    Products if real: Gamma_ee/Gamma_tot ~ 10/11 -> per fusion 0.91 e+e- pair (E_e+ + E_e- = 22.8 MeV,"
       " continuum), 2 x 511 keV back-to-back after e+ stops, bremsstrahlung up to ~20 MeV; p and n branches"
       " each ~4.5%. A neutron-only or p-only detector under-counts by x11; a 511-511 coincidence pair plus"
       " a >3 MeV gamma/electron window is needed.")
    Es = np.logspace(-2, 4, 400)
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    for i, E_R in enumerate((3000, 1000, 300, 0.0)):
        C = (Kcal - 1) * ((2500 - E_R) ** 2 + Gout ** 2 / 4)
        ax.loglog(Es, resonance_K(Es, E_R, C), color=COL[i], lw=2, label=f"E$_R$ = {E_R:.0f} eV")
    ax.axvline(0.04, color="grey", lw=0.8, ls=":")
    ax.text(0.045, 3e7, "thermal pair", fontsize=7, color="#52514e")
    ax.axvline(2500, color="grey", lw=0.8, ls=":")
    ax.text(2600, 3e7, "calibration\n(2.5 keV)", fontsize=7, color="#52514e")
    ax.set_xlabel("pair relative energy E$_{cm}$ (eV)")
    ax.set_ylabel("S-factor multiplier K(E)")
    ax.set_title("M1-B6: threshold-resonance multiplier, calibrated to Szczecin 2.5 keV excess", fontsize=9)
    ax.legend(fontsize=7, frameon=False)
    ax.grid(color="#e6e6e6", lw=0.6)
    fig.tight_layout()
    fig.savefig(os.path.join(P.OUT, "m1_resonance.png"), dpi=140)
    return rows


# ------------------------------------------------------------------ B7
def part_B7():
    Pr("\nB7. Thermal-spike 'plateau' (UC Davis/LBNL 2026; Szczecin 2024-26) - not available beam-free")
    Yp = P.thick_target_yield(2500, "ZrD2", 340) + 0 * P.thick_target_yield(2500, "PdD0.9", 340)
    Pr(f"    Plateau yield per incident D taken = thick-target yield at E_d=2.5 keV with U_e=340 eV (ZrD2):"
       f" {Yp:.1e} fusions/D")
    need = 2.7e-3 / Yp
    Pr(f"    To reach the C3 proton threshold (2.7e-3 /s) a device needs {need:.1e} keV deuterons/s"
       f" = {need*e_C*1e6:.1f} uA of keV D+ -> an ion beam. Cold sources of keV deuterons:")
    rec = sum(phi * 6.1e22 * s_el * 1e-24 for _, phi, _, s_el, _, _ in NBANDS)
    Pr(f"      cosmic knock-on recoils in 1 cm3 PdD: {rec:.0e} /s (B5); LiD fracture at 1e9 V/m: ~5e13 D+ per cm2"
       " once per fresh face (B1). Neither gives a sustained 1e14 /s.")
    Pr("    A cold drive (current pulse, laser <=1 W) raises local T by <=1e2 K; a thermal spike needs"
       " ~keV into ~nm^3 (T ~ 1e4 K). Equivalent spike rate from cold means: 0.")
    Pr(f"    Conclusion: the plateau mechanism contributes <= {1e-3*Yp:.0e} fusions/s in a beam-free device.")


def main():
    part_B0()
    part_B1()
    part_B2()
    part_B3()
    part_B4()
    b5 = part_B5()
    part_B6()
    part_B7()
    with open(os.path.join(P.OUT, "m1_nonthermal.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    return b5


if __name__ == "__main__":
    main()
