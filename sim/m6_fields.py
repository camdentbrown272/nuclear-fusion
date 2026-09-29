"""M6 Q3: local field enhancement |E_loc/E_0| at tips (spheroidal asperities) and nanogaps
(sphere dimers) for Pd, PdD and Au; relevance to (a) the d-d Coulomb barrier and
(b) surface chemistry / loading.

Models (quasi-static, valid for feature sizes << lambda/(2 pi n) ~ 90 nm):
 * prolate spheroid (asperity on a plane = half spheroid by image symmetry): field just outside the
   tip E_tip/E_0 = (eps/eps_d) / (1 + L (eps/eps_d - 1)), L = depolarisation factor.
 * sphere dimer, field along the axis: exact axisymmetric multipole solution (Legendre series with
   axial translation theorem), truncated at l_max, convergence shown.
Run: python3 sim/m6_fields.py -> docs/models/figs/m6_fields.png, m6_fields.txt
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import gammaln, eval_legendre
from m6_common import eps, savefig, Tee, e, eps0, c0, amu, me, kB, hbar

out = Tee("m6_fields.txt")


def e1(m, l):
    return eps(m, np.array([l]))[0]


def depol_prolate(aspect):
    """Depolarisation factor along the long axis of a prolate spheroid, aspect = a/b >= 1."""
    aspect = np.asarray(aspect, float)
    ecc = np.sqrt(np.clip(1 - 1 / aspect ** 2, 1e-12, 1))
    L = (1 - ecc ** 2) / ecc ** 2 * (np.log((1 + ecc) / (1 - ecc)) / (2 * ecc) - 1)
    return np.where(aspect < 1.0001, 1 / 3, L)


def tip_enh(em, ed, aspect):
    L = depol_prolate(aspect)
    ein = 1 / (1 + L * (em / ed - 1))
    return np.abs(em / ed * ein)


def dimer_gap_field(em, ed, R, g, lmax=150):
    """Two identical spheres radius R, surface gap g, uniform applied field E0 along the axis.
    Returns |E|/E0 at the gap midpoint and just outside the sphere surface on the axis."""
    d = 2 * R + g                    # centre separation
    ls = np.arange(1, lmax + 1)
    n = len(ls)
    # Potential outside sphere j: sum_l B_l^j r_j^-(l+1) P_l(cos th_j) (R=1 units).
    # Axial translation: for a point near centre A, with centre B at +D z_hat from A,
    #   P_l(cos th_B)/r_B^(l+1) = sum_n (-1)^l C(l+n,n) r_A^n P_n(cos th_A) / D^(l+n+1)
    # and with B at -D z_hat: (-1)^n C(l+n,n) ... (verified numerically below).
    Dn = d / R
    lg = lambda a: gammaln(a + 1)
    T_plus = np.zeros((n, n))   # sphere 2 (at +D) multipole l -> regular n about sphere 1
    T_minus = np.zeros((n, n))  # sphere 1 (at -D rel. to sphere 2) -> regular n about sphere 2
    for i, l in enumerate(ls):
        for k, m in enumerate(ls):
            c = np.exp(lg(l + m) - lg(l) - lg(m) - (l + m + 1) * np.log(Dn))
            T_plus[k, i] = (-1) ** l * c
            T_minus[k, i] = (-1) ** m * c
    # n=0 (monopole of regular expansion) does not affect the field on a neutral sphere: skip.
    alpha = np.array([m * (ed - em) / (m * em + (m + 1) * ed) for m in ls])   # B_m = alpha_m C_m R^(2m+1)
    # external regular coefficients: applied field -E0 z = -E0 r P1  -> C_1 = -1
    C0 = np.zeros(n)
    C0[0] = -1.0
    # unknowns B1 (sphere 1), B2 (sphere 2):  B1 = alpha*(C0 + T_plus B2),  B2 = alpha*(C0 + T_minus B1)
    I = np.eye(n)
    Aalpha = np.diag(alpha)
    M = np.block([[I, -Aalpha @ T_plus], [-Aalpha @ T_minus, I]])
    rhs = np.concatenate([alpha * C0, alpha * C0])
    B = np.linalg.solve(M.astype(complex), rhs.astype(complex))
    B1, B2 = B[:n], B[n:]

    def phi(z):   # potential on axis (R = 1 units), sphere 1 centre at 0, sphere 2 at Dn
        r1, r2 = abs(z), abs(z - Dn)
        c1 = np.sign(z) if z != 0 else 1.0
        c2 = np.sign(z - Dn)
        val = -z
        val += np.sum(B1 * c1 ** ls / r1 ** (ls + 1))
        val += np.sum(B2 * c2 ** ls / r2 ** (ls + 1))
        return val

    def Ez(z, hstep=1e-6):
        return -(phi(z + hstep) - phi(z - hstep)) / (2 * hstep)
    zmid = 1 + g / (2 * R)
    zsurf = 1 + 1e-4 * min(g / R, 0.05)
    return abs(Ez(zmid)), abs(Ez(zsurf))


# ---------------- verification ----------------
out("VERIFICATION")
# translation theorem check at a random point
rng = np.random.default_rng(1)
D = 3.0
for trial in range(2):
    x, zz = rng.uniform(-0.4, 0.4, 2)
    rA = np.hypot(x, zz)
    cA = zz / rA
    for l in (1, 3):
        rB = np.hypot(x, zz - D)
        cB = (zz - D) / rB
        lhs = eval_legendre(l, cB) / rB ** (l + 1)
        rhs = sum((-1) ** l * np.exp(gammaln(l + m + 1) - gammaln(l + 1) - gammaln(m + 1)) * rA ** m * eval_legendre(m, cA) / D ** (l + m + 1) for m in range(0, 80))
        rB2 = np.hypot(x, zz + D)
        cB2 = (zz + D) / rB2
        lhs2 = eval_legendre(l, cB2) / rB2 ** (l + 1)
        rhs2 = sum((-1) ** m * np.exp(gammaln(l + m + 1) - gammaln(l + 1) - gammaln(m + 1)) * rA ** m * eval_legendre(m, cA) / D ** (l + m + 1) for m in range(0, 80))
        out(f"  translation theorem l={l}: +D {lhs:.8f} vs {rhs:.8f};  -D {lhs2:.8f} vs {rhs2:.8f}")
# isolated sphere limit: huge gap -> surface field 3 eps/(eps+2 eps_d) relative
em, ed = e1("Au", 0.53), e1("vac", 0.53)
mid, surf = dimer_gap_field(em, ed, 20e-9, 2000e-9, lmax=10)
out(f"  dimer, gap = 100 R: surface field {surf:.3f} vs single sphere |3 eps/(eps+2)| = {abs(3*em/(em+2*ed)):.3f}")
out(f"  spheroid aspect 1 -> sphere: {tip_enh(em, ed, 1.0):.3f} (same formula)")
for lm in (50, 100, 200, 300):
    out(f"  convergence Au dimer R=20 nm, g=1 nm, 633 nm: lmax={lm}: |E_gap|/E0 = {dimer_gap_field(e1('Au',0.633), 1.0, 20e-9, 1e-9, lm)[0]:.2f}")

# ---------------- tip enhancement ----------------
out("\nTIP (prolate-spheroid asperity) enhancement |E_tip/E_0| (quasi-static; add the flat-surface factor")
out("|E_surf/E_inc| <= 2 for the incident field; retardation/radiation damping lower large-aspect Au values).")
asp = np.array([1, 2, 3, 5, 10, 20])
for m, med in (("Pd", "D2O"), ("PdD", "D2O"), ("PdD", "vac"), ("Au", "D2O"), ("Au", "vac")):
    for l in (0.633, 0.785, 1.064):
        v = tip_enh(e1(m, l), e1(med, l), asp)
        out(f"  {m:3s}/{med:3s} {l*1e3:5.0f} nm: " + "  ".join(f"a/b={a:>2}: {x:6.1f}" for a, x in zip(asp, v)))

# ---------------- gap enhancement ----------------
out("\nNANOGAP (sphere dimer, R = 20 nm) |E_gap/E_0| at the gap centre, field along dimer axis:")
gaps = np.array([0.5, 1, 2, 5, 10, 20])
gap_tab = {}
for m, med in (("PdD", "D2O"), ("PdD", "vac"), ("Au", "D2O"), ("Au", "vac")):
    for l in (0.633, 0.785):
        vals = [dimer_gap_field(e1(m, l), e1(med, l), 20e-9, g * 1e-9, 250)[0] for g in gaps]
        gap_tab[(m, med, l)] = vals
        out(f"  {m:3s}/{med:3s} {l*1e3:4.0f} nm: " + "  ".join(f"g={g:>4} nm: {x:6.1f}" for g, x in zip(gaps, vals)))
out("  Gaps < ~0.5 nm: electron tunnelling quenches the enhancement (Savage et al., Nature 491, 574 (2012);")
out("  Zhu et al., Nat. Commun. 7, 11495 (2016)) -> realistic ceiling |E/E0| ~ 100-200 for Au, ~10-30 for Pd/PdD.")

# ---------------- (a) relevance to the nuclear barrier ----------------
out("\n(a) RELEVANCE TO THE d-d COULOMB BARRIER")
k_e = 1 / (4 * np.pi * eps0)
E_coul5fm = k_e * e / (5e-15) ** 2
E_coul05A = k_e * e / (0.5e-10) ** 2
out(f"  Coulomb field of a deuteron at 5 fm: {E_coul5fm:.2e} V/m; at 0.5 A: {E_coul05A:.2e} V/m")
cases = [("CW 30 mW on 1 mm^2 (Letts-type)", 0.03 / 1e-6),
         ("CW 1 W on 0.1 mm^2 (upper CW, ~1e7 W/m^2)", 1.0 / 1e-7),
         ("ps pulse at Pd ablation threshold ~0.1 J/cm^2 / 1 ps (out of 'cold' scope)", 0.1e4 / 1e-12)]
m_D = 2.01410 * amu
for name, I in cases:
    E0 = np.sqrt(2 * I / (c0 * eps0))
    for enh in (10, 100):
        El = E0 * enh
        w = 2 * np.pi * c0 / 0.785e-6
        Up_D = e ** 2 * El ** 2 / (4 * m_D * w ** 2) / e
        Up_e = e ** 2 * El ** 2 / (4 * me * w ** 2) / e
        dV_pair = El * 0.74e-10 * (0.74e-10 / 10e-9)     # differential potential across D2 at a 10 nm tip (field gradient)
        out(f"  {name}: E0 = {E0:.2e} V/m, x{enh}: E_loc = {El:.2e} V/m = {El/E_coul5fm:.1e} of the 5 fm field; "
            f"eV across 0.74 A: {El*0.74e-10:.1e}; relative (gradient) shift on a D2 pair at 10 nm tip: {dV_pair:.1e} eV; "
            f"ponderomotive U_p(D) = {Up_D:.1e} eV, U_p(e) = {Up_e:.1e} eV")
# effect on penetration: X = sqrt(E_G/U); dlnP = 0.5 sqrt(E_G/U) dU/U
EG = 985.8e3
for U in (34.0, 100.0):
    for dU in (1e-6, 1e-3, 1.0):
        out(f"  Using M0's U_eff = {U:.0f} eV: adding dU = {dU:g} eV multiplies the rate by exp({0.5*np.sqrt(EG/U)*dU/U:.2e})")
out("  => even the most extreme optical near-field (CW) shifts pair energies by <1e-4 eV; the rate changes")
out("     by a factor 1 + O(1e-5). Negligible, as expected (R4 sec. 2.8).")

# ---------------- (b) relevance to surface chemistry / loading ----------------
out("\n(b) RELEVANCE TO SURFACE CHEMISTRY AND LOADING")
lam = 0.785
Eph = 1.23984 / lam * e
for P, area_cm2 in ((0.03, 0.01), (0.03, 0.5), (0.3, 0.5)):
    flux_ph = P / Eph / area_cm2
    out(f"  {P*1e3:.0f} mW over {area_cm2} cm^2: photon flux {flux_ph:.2e} cm^-2 s^-1 vs electrolysis at 0.1 A/cm^2 = {0.1/e:.2e} e cm^-2 s^-1")
out("  Hot-carrier (plasmon) quantum yield for H2 dissociation / desorption on Au/Pd: 1e-4 - 1e-2 [BK:")
out("  Mukherjee et al., Nano Lett. 13, 240 (2013); Zhou et al., Science 362, 69 (2018)].")
for P, area in ((0.03, 0.5), (0.3, 0.5)):
    for qy in (1e-4, 1e-2):
        rate = P / Eph / area * qy
        out(f"  {P*1e3:.0f} mW / {area} cm^2, QY {qy:g}: photo-driven D events {rate:.1e} cm^-2 s^-1 = "
            f"{rate/(1e16):.1e} of a 1e16 D cm^-2 s^-1 permeation flux")
# photothermal modulation of the exit-face recombination (C3)
kPd = 72.0   # W/m/K, Pd (CRC); PdD lower ~ 30-50 [BK]
for P, A_abs, w in ((0.03, 0.3, 0.5e-3), (0.3, 0.3, 0.5e-3), (0.3, 0.9, 0.5e-3)):
    for kap in (kPd, 30.0):
        dT = P * A_abs / (2 * np.sqrt(np.pi) * kap * w)
        out(f"  photothermal: {P*1e3:.0f} mW, A={A_abs}, spot w={w*1e3:.1f} mm, kappa={kap:.0f} W/m/K (semi-infinite): dT = {dT:.3f} K")
Ea = 0.5  # eV, D2 recombinative desorption from Pd, 0.4-0.9 eV [BK: Behm, Christmann & Ertl, Surf. Sci. 99, 320 (1980)]
T = 300.0
out(f"  Arrhenius sensitivity of recombination at 300 K: dln k_r/dT = Ea/kT^2 = {Ea*e/(kB*T**2)*100:.1f} %/K (Ea=0.5 eV), {0.9*e/(kB*T**2)*100:.1f} %/K (Ea=0.9 eV)")
out("  => a chopped 0.1-1 W laser is a usable, contact-free MODULATOR of exit-face recombination/flux")
out("     (few % to tens of %), acting thermally; plasmonic structuring raises A from ~0.3 to ~0.9 (x3).")

# ---------------- figure ----------------
fig, ax = plt.subplots(1, 2, figsize=(12, 4.3))
aa = np.linspace(1, 20, 100)
for m, med, l in (("PdD", "D2O", 0.785), ("PdD", "vac", 0.785), ("Au", "D2O", 0.785), ("Au", "vac", 0.633)):
    ax[0].plot(aa, tip_enh(e1(m, l), e1(med, l), aa), label=f"{m}/{med}, {l*1e3:.0f} nm")
ax[0].set(xlabel="asperity aspect ratio a/b", ylabel="|E_tip / E_0|", yscale="log", title="Lightning-rod (spheroid) enhancement")
for (m, med, l), v in gap_tab.items():
    ax[1].loglog(gaps, v, "o-", label=f"{m}/{med}, {l*1e3:.0f} nm")
ax[1].axvspan(0.1, 0.5, color="grey", alpha=0.2, label="tunnelling-quenched")
ax[1].set(xlabel="gap (nm)", ylabel="|E_gap / E_0|", title="Sphere-dimer gap field (R = 20 nm)")
for a in ax:
    a.legend(fontsize=7)
    a.grid(alpha=0.3, which="both")
savefig(fig, "m6_fields.png")
out.save()
