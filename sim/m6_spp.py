"""M6 Q1: surface-plasmon-polariton (SPP) dispersion for Pd, PdD, Au, Au-on-Pd at vacuum and
D2O interfaces; which roughness wavevectors couple which photon energies; reconstruction and
consistency test of the ENEA roughness-PSD band; light-source requirements.

Run: python3 sim/m6_spp.py   -> docs/models/figs/m6_spp_*.png, m6_spp.txt
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import newton
from m6_common import eps, eV_to_um, um_to_eV, savefig, Tee, kB, e, hbar, c0, F_HYD

out = Tee("m6_spp.txt")

E = np.linspace(0.66, 3.4, 400)       # photon energy, eV (0.36-1.88 um; J&C range)
lam = eV_to_um(E)


def kspp_2(em, ed):
    """Normalised SPP wavevector (k/k0) of a single metal/dielectric interface."""
    return np.sqrt(em * ed / (em + ed))


def kz(epsx, kx):
    q = np.sqrt(epsx - kx ** 2 + 0j)
    return np.where(q.imag < 0, -q, q)


def kspp_3(ed, ef, es, t_k0, guess=None):
    """TM bound mode of dielectric / film (thickness t, normalised by k0) / substrate.
    Numerator of the pole of r_123 (H-field Fresnel coefficients):
      (p1+p2)(p2+p3) + (p1-p2)(p2-p3) exp(2 i kz2 t) = 0, p_i = kz_i/eps_i.
    Solved by continuation in thickness from a thick film (single Au/dielectric SPP)."""
    def F(kx, t):
        k1, k2, k3 = kz(ed, kx), kz(ef, kx), kz(es, kx)
        p1, p2, p3 = k1 / ed, k2 / ef, k3 / es
        return (p1 + p2) * (p2 + p3) + (p1 - p2) * (p2 - p3) * np.exp(2j * k2 * t)
    kx = kspp_2(ef, ed)
    for t in np.geomspace(max(t_k0, 20.0), t_k0, 25) if t_k0 < 20 else [t_k0]:
        kx = newton(lambda k: F(k, t), kx, tol=1e-13, maxiter=500)
    return kx


interfaces = [
    ("Pd / vacuum", "Pd", "vac"),
    ("PdD / vacuum", "PdD", "vac"),
    ("Pd / D2O", "Pd", "D2O"),
    ("PdD / D2O", "PdD", "D2O"),
    ("Au / vacuum", "Au", "vac"),
    ("Au / D2O", "Au", "D2O"),
    ("Ni / vacuum", "Ni", "vac"),
]
res = {}
for name, m, d in interfaces:
    em, ed = eps(m, lam), eps(d, lam)
    ks = kspp_2(em, ed)
    res[name] = (ks, ed)

# Au film on PdD under D2O (Letts-type Au overlay), thickness 10, 20, 50 nm
for t_nm in (10, 20, 50):
    ks = np.empty_like(lam, dtype=complex)
    for i, l in enumerate(lam):
        ed, ef, es = eps("D2O", np.array([l]))[0], eps("Au", np.array([l]))[0], eps("PdD", np.array([l]))[0]
        ks[i] = kspp_3(ed, ef, es, 2 * np.pi * t_nm * 1e-3 / l)
    res[f"Au {t_nm} nm on PdD / D2O"] = (ks, eps("D2O", lam))

# ---------------- verification: thick film -> single interface ----------------
l = 0.785
ed, ef, es = eps("D2O", np.array([l]))[0], eps("Au", np.array([l]))[0], eps("PdD", np.array([l]))[0]
k_thick = kspp_3(ed, ef, es, 2 * np.pi * 0.4 / l, kspp_2(ef, ed))
out("VERIFY 3-layer solver, Au 400 nm on PdD vs Au/D2O single interface at 785 nm:",
    np.round(k_thick, 6), np.round(kspp_2(ef, ed), 6))
t0 = kspp_3(ed, es * 1.0, es, 2 * np.pi * 0.02 / l, kspp_2(es, ed))
out("VERIFY film = substrate material (PdD/PdD) reproduces PdD/D2O:", np.round(t0, 6),
    np.round(kspp_2(es, ed), 6))

# ---------------- table at common laser lines ----------------
lasers = [0.405, 0.532, 0.633, 0.670, 0.785, 0.830, 1.064, 1.55]
out("\nSPP properties at laser lines. n_eff = Re k_spp/k0; L = 1/(2 Im k_spp) propagation length;")
out("d_d = 1/Im kz in dielectric (field decay length); Lam_0 = grating period for 1st-order coupling at")
out("normal incidence = lambda/n_eff; f_0 = 1/Lam_0 (cycles/m); Q = Re k/(2 Im k).")
hdr = f"{'interface':28s} {'lam nm':>6s} {'E eV':>5s} {'n_eff':>7s} {'L um':>7s} {'d_d nm':>7s} {'Lam_0 nm':>8s} {'f_0 1/m':>9s} {'Q':>6s}"
out(hdr)
table = {}
for name in res:
    for l in lasers:
        if "on PdD" in name:
            t_nm = int(name.split()[1])
            ed, ef, es = eps("D2O", np.array([l]))[0], eps("Au", np.array([l]))[0], eps("PdD", np.array([l]))[0]
            ks = kspp_3(ed, ef, es, 2 * np.pi * t_nm * 1e-3 / l, kspp_2(ef, ed))
            edd = ed
        else:
            m, d = [x for x in interfaces if x[0] == name][0][1:]
            em, edd = eps(m, np.array([l]))[0], eps(d, np.array([l]))[0]
            ks = kspp_2(em, edd)
        k0 = 2 * np.pi / l
        L = 1 / (2 * ks.imag * k0)
        dd = 1 / (kz(edd, ks).imag * k0)
        Lam0 = l / ks.real
        Q = ks.real / (2 * ks.imag)
        table[(name, l)] = dict(neff=ks.real, L=L, dd=dd, Lam0=Lam0, Q=Q)
        out(f"{name:28s} {l*1e3:6.0f} {1.23984/l:5.2f} {ks.real:7.4f} {L:7.2f} {dd*1e3:7.0f} {Lam0*1e3:8.0f} {1/(Lam0*1e-6):9.3e} {Q:6.1f}")

# ---------------- ENEA band reconstruction ----------------
out("\nENEA roughness-PSD band (R1 G4: 'micron to sub-micron', ~1e5-1e7 m^-1, band NOT verified).")
out("Two readings: (A) spatial frequency f = 1/Lambda in cycles/m -> q = 2 pi f; (B) q itself in rad/m.")
for label, qlo, qhi in [("A: f=1e5-1e7 /m", 2 * np.pi * 1e5, 2 * np.pi * 1e7), ("B: q=1e5-1e7 rad/m", 1e5, 1e7)]:
    for name in ("Pd / D2O", "PdD / D2O", "Au / D2O", "PdD / vacuum"):
        ks, ed = res[name]
        k0 = 2 * np.pi / (lam * 1e-6)
        kre = ks.real * k0
        nd = np.sqrt(ed.real)
        # first-order grating coupling possible at some angle if G in [kre - nd k0, kre + nd k0]
        Gmin, Gmax = kre - nd * k0, kre + nd * k0
        ok_any = (Gmax >= qlo) & (Gmin <= qhi)
        ok_norm = (kre >= qlo) & (kre <= qhi)
        Eany = E[ok_any]
        En = E[ok_norm]
        out(f"  {label:20s} {name:14s}: normal-incidence 1st order couples E = "
            f"{(f'{En.min():.2f}-{En.max():.2f} eV' if En.size else 'none in 0.66-3.4 eV')}; "
            f"some angle: {(f'{Eany.min():.2f}-{Eany.max():.2f} eV' if Eany.size else 'none')}")
ks, ed = res["PdD / D2O"]
k0 = 2 * np.pi / (lam * 1e-6)
fmax = (ks.real * k0).max() / (2 * np.pi)
out(f"  Highest f that couples at normal incidence for E <= 3.4 eV (PdD/D2O): {fmax:.2e} /m (Lambda = {1e9/fmax:.0f} nm).")
out("  => under reading A, f > ~3.5e6 /m (Lambda < ~285 nm) couples no visible/NIR photon in first order;")
out("     above ~5 eV Re eps(Pd) > -eps_d and no bound SPP exists. The top half-decade is not plasmonic.")
d = (ks.real - np.sqrt(ed.real)) * k0 / (2 * np.pi)
for Ev in (1.17, 1.58, 1.96, 2.33):
    i = np.argmin(abs(E - Ev))
    out(f"  Minimum f for 1st-order coupling at grazing incidence, PdD/D2O, E={Ev} eV: (n_eff-n_d)/lambda = {d[i]:.2e} /m"
        f"; SPP linewidth in f: Im k/pi = {ks.imag[i]*k0[i]/np.pi:.2e} /m")
out("  => because Pd/PdD SPPs lie within ~4 % of the light line and have Q ~ 10-30, EVERY spatial frequency")
out("     from ~1e5 to ~3e6 /m couples some visible/NIR photon at some angle: the band is not selective.")
out("  f = 1e5 /m at normal incidence would need lambda ~ 13 um in D2O, where D2O absorbs (alpha ~ 1e3 cm^-1)")
out("  and the Pd 'SPP' is a Zenneck-like wave extending ~lambda into the medium (see mid-IR check below).")

# mid-IR check of how loosely bound Pd SPPs are at 3 and 10 um (vacuum side, C3)
for l in (3.0, 10.0):
    em = eps("Pd", np.array([l]))[0]
    ks = kspp_2(em, 1.0)
    k0 = 2 * np.pi / l
    out(f"  Pd/vacuum at {l} um: n_eff-1 = {ks.real-1:.2e}, field decay into vacuum = {1/(kz(1.0, ks).imag*k0):.1f} um, L = {1/(2*ks.imag*k0)/1e3:.2f} mm")

# ---------------- thermal SPP population (dark cell) ----------------
T = 300.0
out("\nThermal (dark-cell) SPP occupation n_BE = 1/(exp(E/kT)-1) at 300 K:")
for Ev in (0.1, 0.5, 1.0, 1.55, 1.96, 2.5):
    n = 1 / np.expm1(Ev * e / (kB * T))
    out(f"  E = {Ev:4.2f} eV: n = {n:.2e}")
out("  => visible/NIR SPPs are not populated without external illumination (n < 1e-16 at 1 eV).")
out("  A PSD correlation in dark electrolytic cells cannot be a plasmon-coupling effect; if real it is")
out("  a proxy for metallurgy (etch pits, grain/texture, dislocation outcrops) that also sets the PSD.")

# ---------------- dual-laser pairs for Letts beat frequencies ----------------
out("\nDual-laser pairs for the Letts beat frequencies (Delta lambda = lambda^2 Delta f / c):")
for lc in (0.670, 0.785, 0.830):
    for df in (8.3e12, 15.3e12, 20.4e12):
        f1 = c0 / (lc * 1e-6)
        l2 = c0 / (f1 - df) * 1e6
        out(f"  centre {lc*1e3:.0f} nm, df = {df/1e12:4.1f} THz: lambda1 = {lc*1e3:.1f} nm, lambda2 = {l2*1e3:.1f} nm (Delta = {(l2-lc)*1e3:.1f} nm)")

# ---------------- figures ----------------
fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
cols = plt.cm.tab10(np.arange(10))
for i, (name, (ks, ed)) in enumerate(res.items()):
    ax[0].plot(E, ks.real, color=cols[i], label=name)
    k0 = 2 * np.pi / lam
    ax[1].semilogy(E, 1 / (2 * ks.imag * k0), color=cols[i], label=name)
    ax[2].semilogy(E, ks.real * k0 * 1e6 / (2 * np.pi), color=cols[i], label=name)
ax[0].plot(E, np.sqrt(eps("D2O", lam).real), "k--", lw=1, label="light line in D2O")
ax[0].plot(E, np.ones_like(E), "k:", lw=1, label="light line in vacuum")
ax[0].set(xlabel="photon energy (eV)", ylabel="Re k_spp / k0", title="SPP effective index")
ax[1].set(xlabel="photon energy (eV)", ylabel="propagation length L (um)", title="SPP propagation length 1/(2 Im k)")
ax[2].axhspan(1e5, 1e7, color="orange", alpha=0.15, label="ENEA band, reading A (f=1e5-1e7 /m)")
ax[2].axhspan(1e5 / (2 * np.pi), 1e7 / (2 * np.pi), color="purple", alpha=0.08, label="reading B (q=1e5-1e7 rad/m)")
ax[2].set(xlabel="photon energy (eV)", ylabel="f_spp = Re k_spp / 2pi  (cycles / m)",
          title="Roughness spatial frequency that couples at normal incidence", ylim=(1e4, 2e7))
ax[0].legend(fontsize=6.5)
ax[2].legend(fontsize=6.5, loc="lower right")
for a in ax:
    a.grid(alpha=0.3)
savefig(fig, "m6_spp_dispersion.png")

# dielectric functions
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
lp = np.linspace(0.3, 1.9, 300)
for m in ("Pd", "Pd_Palm", "PdD", "Au", "Ni", "Cu"):
    lq = lp[lp < 1.68] if m == "Pd_Palm" else lp
    ev = eps(m, lq)
    ax[0].plot(lq, ev.real, label=m)
    ax[1].plot(lq, ev.imag, label=m)
ax[0].set(xlabel="wavelength (um)", ylabel="Re eps", title="Real permittivity")
ax[1].set(xlabel="wavelength (um)", ylabel="Im eps", title=f"Imag permittivity (PdD = {F_HYD} x Pd, [BK])")
for a in ax:
    a.legend(fontsize=7)
    a.grid(alpha=0.3)
savefig(fig, "m6_eps.png")
out.save()
