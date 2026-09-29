"""M6 Q1-Q3: deterministic gratings on PdD, Pd and Au-over-PdD.  RCWA (TM) absorptance maps,
optimum period/depth at common laser lines, tolerances, and near-surface field enhancement.

Run: python3 sim/m6_gratings.py  -> docs/models/figs/m6_grating_*.png, m6_gratings.txt
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
from m6_common import eps, savefig, Tee
from m6_rcwa import rcwa_tm, sinusoid_layers

out = Tee("m6_gratings.txt")
NH = 30        # harmonics -NH..NH (convergence shown below)
NS = 16        # staircase slices


def e1(m, l):
    return eps(m, np.array([l]))[0]


def system(name, l, period, depth, theta=0.0, want_field=False, N=NH, ns=NS):
    """Return RCWA result for a named surface at wavelength l (um)."""
    if name == "PdD/D2O":
        es, ed, under = e1("PdD", l), e1("D2O", l), None
        sub = es
        gm = es
    elif name == "Pd/D2O":
        gm, ed, under = e1("Pd", l), e1("D2O", l), None
        sub = gm
    elif name == "PdD/vac":
        gm, ed, under = e1("PdD", l), 1.0 + 0j, None
        sub = gm
    elif name == "Au(20nm)/PdD, D2O":
        # sinusoidal grating cut in the top of a 20 nm (mean) Au overlayer on PdD:
        # crest-to-trough in Au, then a uniform 20 nm Au film, then PdD
        gm, ed = e1("Au", l), e1("D2O", l)
        under = [{"d": 0.020, "eps": gm}]
        sub = e1("PdD", l)
    elif name == "Au(thick)/D2O":
        gm, ed, under = e1("Au", l), e1("D2O", l), None
        sub = gm
    else:
        raise KeyError(name)
    lay = sinusoid_layers(depth, gm, ed, ns, under)
    return rcwa_tm(l, period, theta, ed, lay, sub, N=N, want_field=want_field)


# ------------------------------ verification ------------------------------
out("VERIFICATION")
l, P, hgt = 0.785, 0.575, 0.040
ref = None
for N in (10, 20, 30, 40, 60):
    A = system("PdD/D2O", l, P, hgt, N=N)["A"]
    out(f"  PdD/D2O 785 nm, P=575 nm, h=40 nm: N={N:3d} harmonics  A = {A:.4f}")
for ns in (8, 16, 32):
    A = system("PdD/D2O", l, P, hgt, ns=ns)["A"]
    out(f"  slices={ns:3d}  A = {A:.4f}")
for N in (20, 30, 40, 60):
    A = system("Au(20nm)/PdD, D2O", l, 0.57, 0.030, N=N)["A"]
    out(f"  Au(20nm)/PdD 785 nm, P=570 nm, h=30 nm: N={N:3d}  A = {A:.4f}")
# flat limit vs Fresnel
em, ed = e1("PdD", l), e1("D2O", l)
kz1, kz2 = np.sqrt(ed), np.sqrt(em)
rp = (em * kz1 - ed * kz2) / (em * kz1 + ed * kz2)
out(f"  flat PdD/D2O: Fresnel A = {1-abs(rp)**2:.5f}; RCWA(h=1e-6) A = {system('PdD/D2O', l, P, 1e-6)['A']:.5f}")

# ------------------------------ maps ------------------------------
lams = np.linspace(0.50, 1.10, 49)
pers = np.linspace(0.25, 0.90, 53)
maps = {}
depth_for = {"PdD/D2O": 0.040, "Au(20nm)/PdD, D2O": 0.030, "PdD/vac": 0.040}
for nm in depth_for:
    Amap = np.zeros((len(pers), len(lams)))
    for j, Pp in enumerate(pers):
        for i, ll in enumerate(lams):
            Amap[j, i] = system(nm, ll, Pp, depth_for[nm], N=20, ns=12)["A"]
    maps[nm] = Amap
fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
for a, nm in zip(ax, maps):
    flat = np.array([system(nm, ll, 0.5, 1e-6, N=5, ns=2)["A"] for ll in lams])
    im = a.pcolormesh(lams * 1e3, pers * 1e3, maps[nm] / flat[None, :], shading="auto", cmap="magma")
    plt.colorbar(im, ax=a, label="A_grating / A_flat")
    a.set(xlabel="wavelength (nm)", ylabel="period (nm)", title=f"{nm}, h = {depth_for[nm]*1e3:.0f} nm, normal incidence, TM")
savefig(fig, "m6_grating_maps.png")

# ------------------------------ optimum at laser lines ------------------------------
out("\nOPTIMUM SINUSOIDAL GRATINGS (normal incidence, TM/p-polarised; grating vector in plane of incidence)")
out(f"{'surface':20s} {'lam':>5s} {'A_flat':>7s} {'P_opt':>6s} {'h_opt':>6s} {'A_max':>6s} {'P tol(A>=half gain)':>20s} {'h range(A>=0.9Amax)':>21s} {'|E|^2 crest/flat':>16s}")
rows = []
for nm in ("PdD/D2O", "Pd/D2O", "PdD/vac", "Au(20nm)/PdD, D2O", "Au(thick)/D2O"):
    for ll in (0.633, 0.785, 1.064):
        A_flat = system(nm, ll, 0.5, 1e-6, N=5, ns=2)["A"]
        best = (0, None, None)
        for hh in (0.01, 0.02, 0.03, 0.04, 0.06, 0.08, 0.10):
            r = minimize_scalar(lambda Pp: -system(nm, ll, Pp, hh)["A"],
                                bounds=(0.6 * ll / 1.45, 1.05 * ll), method="bounded", options={"xatol": 1e-3})
            if -r.fun > best[0]:
                best = (-r.fun, r.x, hh)
        Amax, Popt, hopt = best
        # refine depth
        r = minimize_scalar(lambda hh: -system(nm, ll, Popt, hh)["A"], bounds=(0.3 * hopt, 2 * hopt), method="bounded")
        hopt = r.x
        Amax = -r.fun
        Ps = np.linspace(Popt - 0.2 * ll, Popt + 0.2 * ll, 81)
        Aps = np.array([system(nm, ll, Pp, hopt)["A"] for Pp in Ps])
        half = A_flat + 0.5 * (Amax - A_flat)
        okP = Ps[Aps >= half]
        hs = np.linspace(0.2 * hopt, 2.5 * hopt, 47)
        Ahs = np.array([system(nm, ll, Popt, hh)["A"] for hh in hs])
        okh = hs[Ahs >= 0.9 * Amax]
        fr = system(nm, ll, Popt, hopt, want_field=True)
        ff = system(nm, ll, Popt, 1e-6, want_field=True)
        gainE = fr["E2norm"].max() / ff["E2norm"].max()
        rows.append((nm, ll, A_flat, Popt, hopt, Amax, okP.min(), okP.max(), okh.min(), okh.max(), gainE, fr["E2norm"].max()))
        out(f"{nm:20s} {ll*1e3:5.0f} {A_flat:7.3f} {Popt*1e3:6.0f} {hopt*1e3:6.1f} {Amax:6.3f} "
            f"{okP.min()*1e3:8.0f}-{okP.max()*1e3:<8.0f}nm   {okh.min()*1e3:6.1f}-{okh.max()*1e3:<6.1f}nm   {gainE:8.2f} (abs {fr['E2norm'].max():.2f})")

out("\nNotes: A = absorptance; |E|^2 normalised to incident |E|^2 in the superstrate, evaluated on the plane")
out("touching the crests (lower bound of the true local maximum, which sits at the crest surface).")

# ------------------------------ angular scan (C3 vacuum face / C1 via window) ------------------------------
out("\nANGLE TOLERANCE at the 785 nm optimum (A vs incidence angle):")
for nm in ("PdD/D2O", "Au(20nm)/PdD, D2O"):
    Popt, hopt = [(r[3], r[4]) for r in rows if r[0] == nm and abs(r[1] - 0.785) < 1e-6][0]
    ths = np.linspace(0, 10, 21)
    As = [system(nm, 0.785, Popt, hopt, theta=t)["A"] for t in ths]
    Am = max(As)
    half_ok = ths[np.array(As) >= 0.5 * (Am + system(nm, 0.785, 0.5, 1e-6, N=5, ns=2)["A"])]
    out(f"  {nm:20s}: A(0) = {As[0]:.3f}; A >= half-gain for theta <= {half_ok.max():.1f} deg")

# ------------------------------ figure: spectra at optimum ------------------------------
fig, ax = plt.subplots(1, 2, figsize=(12, 4.3))
lsp = np.linspace(0.55, 1.1, 111)
for nm in ("PdD/D2O", "Au(20nm)/PdD, D2O", "PdD/vac"):
    Popt, hopt = [(r[3], r[4]) for r in rows if r[0] == nm and abs(r[1] - 0.785) < 1e-6][0]
    ax[0].plot(lsp * 1e3, [system(nm, ll, Popt, hopt)["A"] for ll in lsp], label=f"{nm}: P={Popt*1e3:.0f}, h={hopt*1e3:.0f} nm")
    ax[0].plot(lsp * 1e3, [system(nm, ll, 0.5, 1e-6, N=5, ns=2)["A"] for ll in lsp], "--", color=ax[0].lines[-1].get_color(), lw=1)
ax[0].set(xlabel="wavelength (nm)", ylabel="absorptance", title="Grating optimised for 785 nm (solid) vs flat (dashed)")
ax[0].legend(fontsize=7)
for nm in ("PdD/D2O", "Au(20nm)/PdD, D2O"):
    Popt, hopt = [(r[3], r[4]) for r in rows if r[0] == nm and abs(r[1] - 0.785) < 1e-6][0]
    fr = system(nm, 0.785, Popt, hopt, want_field=True)
    ax[1].plot(fr["x"], fr["E2norm"], label=nm)
ax[1].set(xlabel="x / period (crest at 0.5)", ylabel="|E|^2 / |E_inc|^2 at crest plane", title="Near-field intensity at resonance, 785 nm")
ax[1].legend(fontsize=7)
for a in ax:
    a.grid(alpha=0.3)
savefig(fig, "m6_grating_spectra.png")
np.savetxt(__import__("os").path.join(__import__("m6_common").FIGS, "m6_grating_optima.csv"),
           np.array([[r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8], r[9], r[10], r[11]] for r in rows]),
           delimiter=",", header="rows in order of m6_gratings.txt table: lam_um,A_flat,P_opt_um,h_opt_um,A_max,P_lo,P_hi,h_lo,h_hi,E2_gain,E2_abs",
           fmt="%.4g")
out.save()
