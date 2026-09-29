"""M6 Q6: segmented active face of the detector-facing membrane (C3).

 1. D flux partition between sectors with different exit skins: 2-D steady diffusion across the
    membrane (x across sectors, z through thickness), electrochemical (Robin) entry condition,
    detailed-balance desorption exit condition with sector-specific effective sticking s.
    Linear Fickian diffusion in x (no explicit alpha/beta interface) -> indicative only; hand-off to M3.
 2. Collimated Si detector per sector: acceptance map vs source position (deterministic quadrature),
    cross-talk between neighbouring sectors, geometric efficiency, minimum sector size / gap.
Run: python3 sim/m6_sectors.py -> docs/models/figs/m6_sector_*.png, m6_sectors.txt
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from m6_common import savefig, Tee, amu, kB, NA, e

out = Tee("m6_sectors.txt")

# ------------------------------------------------------------------------------------
# isotherm (same [BK] anchors as m6_skins.py)
xs_iso = np.array([0.015, 0.60, 0.65, 0.70, 0.80, 0.90, 0.95, 1.00])
lnf_iso = np.log(np.array([0.04, 0.045, 1.0, 10.0, 3e2, 1e4, 1e5, 1e6]))
Ks = 0.015 / np.sqrt(0.04)


def fug_and_deriv(x):
    x = np.clip(x, 1e-9, 0.999)
    f = np.where(x < 0.015, (x / Ks) ** 2, np.exp(np.interp(x, xs_iso, lnf_iso)))
    slope = np.interp(x, 0.5 * (xs_iso[1:] + xs_iso[:-1]), np.diff(lnf_iso) / np.diff(xs_iso))
    # piecewise-linear ln f: use the segment slope
    seg = np.clip(np.searchsorted(xs_iso, x) - 1, 0, len(xs_iso) - 2)
    slope = (np.diff(lnf_iso) / np.diff(xs_iso))[seg]
    df = np.where(x < 0.015, 2 * x / Ks ** 2, f * slope)
    return f, df


T = 300.0
mD2 = 4.028 * amu
Z1 = 1e5 / np.sqrt(2 * np.pi * mD2 * kB * T) * 1e-4        # cm^-2 s^-1 bar^-1
Dd = 5e-7                                                  # cm^2/s  (R1)
nPd = 12.02 / 106.42 * NA                                  # cm^-3


def solve_membrane(sectors, w_mm=3.0, g_mm=0.5, h_um=50.0, j_abs=0.05, x_eq=0.95, gap_s=0.1,
                   dx_um=25.0, nz=20, back_factor=None):
    """sectors: list of (name, s). Returns grid and solution x(xpos, z) (z=0 entry, z=h exit)."""
    n = len(sectors)
    W = n * w_mm + (n + 1) * g_mm
    nx = int(round(W * 1e3 / dx_um))
    dx = W / 10 / nx            # cm
    h = h_um * 1e-4
    dz = h / nz
    xc = (np.arange(nx) + 0.5) * dx * 10    # mm
    s_exit = np.full(nx, gap_s)
    fac_in = np.ones(nx)
    labels = np.full(nx, -1)
    for i, (nm, s) in enumerate(sectors):
        x0 = g_mm + i * (w_mm + g_mm)
        m = (xc >= x0) & (xc < x0 + w_mm)
        s_exit[m] = s
        labels[m] = i
        if back_factor and nm in back_factor:
            fac_in[m] = back_factor[nm]
    J0 = j_abs / e * fac_in                   # D cm^-2 s^-1 at x = 0 (entry)
    N = nx * nz
    idx = lambda i, k: k * nx + i
    # linear diffusion operator (D n_Pd) Laplacian with no-flux sides
    rows, cols, vals = [], [], []
    cx, cz = Dd * nPd / dx ** 2, Dd * nPd / dz ** 2
    for k in range(nz):
        for i in range(nx):
            p = idx(i, k)
            diag = 0.0
            if i > 0:
                rows.append(p); cols.append(idx(i - 1, k)); vals.append(cx); diag -= cx
            if i < nx - 1:
                rows.append(p); cols.append(idx(i + 1, k)); vals.append(cx); diag -= cx
            if k > 0:
                rows.append(p); cols.append(idx(i, k - 1)); vals.append(cz); diag -= cz
            if k < nz - 1:
                rows.append(p); cols.append(idx(i, k + 1)); vals.append(cz); diag -= cz
            rows.append(p); cols.append(p); vals.append(diag)
    Lop = sp.csr_matrix((vals, (rows, cols)), shape=(N, N))
    bot = np.array([idx(i, 0) for i in range(nx)])
    top = np.array([idx(i, nz - 1) for i in range(nx)])

    def resid(u):
        r = Lop @ u
        # entry: J_in = J0 (1 - x/x_eq) into the first cell (per volume: / dz)
        r[bot] += J0 * (1 - u[bot] / x_eq) / dz
        # exit: -2 s Z1 f(x) out of the last cell (approximating the surface value by the cell value)
        f, df = fug_and_deriv(u[top])
        r[top] -= 2 * s_exit * Z1 * f / dz
        return r

    def jac(u):
        f, df = fug_and_deriv(u[top])
        d = np.zeros(N)
        d[bot] += -J0 / x_eq / dz
        d[top] += -2 * s_exit * Z1 * df / dz
        return Lop + sp.diags(d)
    u = np.full(N, 0.3)
    for it in range(200):
        r = resid(u)
        du = spla.spsolve(jac(u).tocsc(), -r)
        lam = 1.0
        base = np.linalg.norm(r)
        while lam > 1e-4:
            un = np.clip(u + lam * du, 1e-9, 0.999)
            if np.linalg.norm(resid(un)) < base * (1 - 1e-4 * lam) or lam < 2e-4:
                break
            lam *= 0.5
        u = un
        if np.max(np.abs(lam * du)) < 1e-10:
            break
    U = u.reshape(nz, nx)
    f, _ = fug_and_deriv(U[-1])
    Jout = 2 * s_exit * Z1 * f
    Jin = J0 * (1 - U[0] / x_eq)
    return dict(xc=xc, U=U, Jout=Jout, Jin=Jin, labels=labels, iters=it, resnorm=np.linalg.norm(resid(u)),
                h=h, dx=dx, W=W)


# ---------------- verification: 1-D limits ----------------
out("1. SECTOR FLUX PARTITION (2-D steady diffusion, h = 50 um, entry j_abs = 0.05 A/cm^2 absorbed, x_eq = 0.95)")
out("VERIFY uniform skin -> 1-D solution (analytic: J = J0(1 - x_in/x_eq) = D n (x_in - x_out)/h = 2 s Z1 f(x_out))")
from scipy.optimize import brentq
for s in (0.1, 1e-7):
    r = solve_membrane([("u", s)] * 2, g_mm=0.0, gap_s=s, dx_um=100)
    J0 = 0.05 / e
    h = 50e-4

    def F(xo):
        fo = fug_and_deriv(np.array([xo]))[0][0]
        J = 2 * s * Z1 * fo
        xin = xo + J * h / (Dd * nPd)
        return J - J0 * (1 - xin / 0.95)
    xo = brentq(F, 1e-9, 0.95)
    out(f"  s = {s:g}: 2-D x_out = {r['U'][-1].mean():.5f}, x_in = {r['U'][0].mean():.5f}, J = {r['Jout'].mean():.3e};"
        f" 1-D x_out = {xo:.5f} (Newton iters {r['iters']}, |res| {r['resnorm']:.1e})")

SECT = [("(a) bare Pd", 0.1), ("(b) PdO 10 nm", 1e-6), ("(c) Au|Pd|PdO", 1e-6), ("(d) Pd/CaO", 0.1),
        ("(e) Ni/Cu", 1e-2), ("(f) Au 20 nm", 1e-7)]
cases = {
    "gaps bare Pd": dict(gap_s=0.1),
    "gaps Au-capped (s=1e-7)": dict(gap_s=1e-7),
}
results = {}
for cname, kw in cases.items():
    r = solve_membrane(SECT, back_factor={"(c) Au|Pd|PdO": 0.1}, **kw)
    results[cname] = r
    out(f"\n  Case: {cname}  (Newton iters {r['iters']}, |res| {r['resnorm']:.1e}); sector (c) entry reduced x0.1 by back-face Au [assumption]")
    out(f"  {'sector':16s} {'s':>7s} {'J_out mean':>11s} {'x_out mean':>10s} {'x_in mean':>9s} {'edge influence width (mm)':>30s}")
    for i, (nm, s) in enumerate(SECT):
        m = r["labels"] == i
        xo = r["U"][-1][m]
        # transition width at the sector edge: distance from the edge over which x_out reaches 90 % of its interior value
        interior = np.median(xo)
        xpos = r["xc"][m] - r["xc"][m][0]
        if interior > 0.05:
            edge = f"{xpos[np.argmax(np.abs(xo - interior) < 0.1 * interior)]:.2f}"
        else:
            edge = "n/a (drained)"
        out(f"  {nm:16s} {s:7.0e} {r['Jout'][m].mean():11.2e} {xo.mean():10.3f} {r['U'][0][m].mean():9.3f} {edge:>30s}")
    tot = r["Jout"].sum()
    out("  Share of total exit flux by sector: " + ", ".join(f"{nm.split()[0]} {r['Jout'][r['labels']==i].sum()/tot*100:.1f}%" for i, (nm, s) in enumerate(SECT)) +
        f", gaps {r['Jout'][r['labels']==-1].sum()/tot*100:.1f}%")
out("\n  => sectors do NOT share the same D flux or loading: high-s skins (bare Pd, Pd/CaO, Ni/Cu) drain their")
out("     footprint to x < 0.05 at the entry and ~1e-4 at the exit and carry 3-4x the flux of the low-s skins")
out("     (PdO, Au), which sit on the plateau/beta boundary (x ~ 0.55-0.65). Bare inter-sector gaps drain too.")
out("     Lateral coupling is confined to <~0.3 mm (a few membrane thicknesses) from each sector edge, so sectors of >= 2 mm")
out("     are independent -- but a sector comparison confounds skin chemistry with (x, J), which differ by orders")
out("     of magnitude. Equalise with a D2 backfill (m6_skins D) or accept and log per-sector (x_out, J) estimates.")

# 1-D sweep of j_abs for the bare-Pd sector to show the drained regime persists
out("\n  Bare-Pd exit, 1-D: absorbed current needed to reach x_in at the entry face (h = 50 um):")
for xin in (0.6, 0.8, 0.9):
    J = Dd * nPd * xin / 50e-4      # x_out ~ 0
    out(f"   x_in = {xin}: J = {J:.2e} D cm^-2 s^-1 = {J*e:.2f} A/cm^2 absorbed (plus J0-saturation: j_abs >= that/(1-x_in/0.95))")

fig, ax = plt.subplots(2, 1, figsize=(11, 6.5), sharex=True)
for cname, r in results.items():
    ax[0].plot(r["xc"], r["U"][-1], label=f"exit face, {cname}")
    ax[0].plot(r["xc"], r["U"][0], "--", label=f"entry face, {cname}")
    ax[1].semilogy(r["xc"], np.maximum(r["Jout"], 1e8), label=cname)
r = results["gaps bare Pd"]
for i, (nm, s) in enumerate(SECT):
    x0 = 0.5 + i * 3.5
    for a in ax:
        a.axvspan(x0, x0 + 3, color=f"C{i}", alpha=0.08)
    ax[0].text(x0 + 1.5, 0.97, nm, ha="center", fontsize=7, rotation=0)
ax[0].set(ylabel="loading x = D/Pd", ylim=(0, 1.02), title="Six-sector membrane, h = 50 um, j_abs = 0.05 A/cm^2")
ax[1].set(xlabel="position across membrane (mm)", ylabel="exit flux (D cm^-2 s^-1)")
ax[0].legend(fontsize=6, loc="center right")
ax[1].legend(fontsize=7)
savefig(fig, "m6_sector_flux.png")

# =====================================================================================
out("\n2. COLLIMATED Si DETECTOR PER SECTOR: acceptance, cross-talk, efficiency")


def acceptance(rho, D, r_det, r_ap, z_ap, n=160):
    """Fraction of 4pi from a point source at radial offset rho (membrane plane z=0) that reaches a
    detector disc (radius r_det at height D, centred on the axis) through a circular aperture
    (radius r_ap at height z_ap, centred). Tube walls between aperture and detector are absorbing;
    since both discs are coaxial and the tube radius >= both, the straight path stays inside."""
    # polar quadrature on the detector disc
    rr = (np.arange(n) + 0.5) / n * r_det
    th = (np.arange(2 * n) + 0.5) / (2 * n) * 2 * np.pi
    R, TH = np.meshgrid(rr, th, indexing="ij")
    dA = (r_det / n) * (2 * np.pi / (2 * n)) * R
    px, py = R * np.cos(TH), R * np.sin(TH)
    out_ = []
    for r0 in np.atleast_1d(rho):
        vx, vy = px - r0, py
        dist2 = vx ** 2 + vy ** 2 + D ** 2
        t = z_ap / D
        ix, iy = r0 + t * vx, t * vy
        ok = ix ** 2 + iy ** 2 <= r_ap ** 2
        dOmega = D / np.sqrt(dist2) * dA / dist2
        out_.append(np.sum(dOmega * ok) / (4 * np.pi))
    return np.array(out_)


# verification: no aperture restriction, on-axis -> 0.5(1 - D/sqrt(D^2 + r^2))
D, rd = 20.0, 3.0
a0 = acceptance(0.0, D, rd, 1e9, 1.0)[0]
out(f"VERIFY on-axis disc solid angle: quadrature {a0:.6f} vs analytic {0.5*(1-D/np.hypot(D,rd)):.6f}")

configs = []
for D in (10.0, 20.0, 30.0):
    for r_det in (2.0, 3.0, 4.0):               # detector radius mm (12.6, 28, 50 mm^2)
        for z_ap in (0.5, 1.0, 2.0):            # aperture (knife-edge) height above membrane, mm
            r_ap = r_det                        # aperture radius = detector radius (straight tube)
            configs.append((D, r_det, z_ap, r_ap))

out("  Layout: circular sectors (patches) of radius r_s on a hexagonal pitch p = 2 r_s + g; each detector coaxial with")
out("  its patch; aperture radius = detector radius; cross-talk = counts from the 6 nearest patches / own-patch counts")
out("  for equal emission per unit area. Penumbra (geometric) radius r_p = r_ap + 2 r_ap z_ap/(D - z_ap).")
rows = []
rq = np.linspace(0, 12, 481)
for D, r_det, z_ap, r_ap in configs:
    A = acceptance(rq, D, r_det, r_ap, z_ap, n=60)
    r_s = r_ap + 0.0                             # patch radius = aperture radius
    own = np.trapezoid(A * 2 * np.pi * rq * (rq <= r_s), rq)
    r_p = r_ap + 2 * r_ap * z_ap / (D - z_ap)
    for g in (0.25, 0.5, 1.0, 2.0):
        p = 2 * r_s + g
        # counts from a neighbouring patch centred at distance p: integrate A over that disc (2-D quadrature)
        xs = np.linspace(p - r_s, p + r_s, 121)
        ys = np.linspace(-r_s, r_s, 121)
        X, Y = np.meshgrid(xs, ys)
        inside = (X - p) ** 2 + Y ** 2 <= r_s ** 2
        Ar = np.interp(np.hypot(X, Y), rq, A)
        nb = np.sum(Ar * inside) * (xs[1] - xs[0]) * (ys[1] - ys[0])
        eff = own / (np.pi * r_s ** 2)
        rows.append((D, r_det, z_ap, g, 6 * nb / own, eff, r_p))
out(f"  {'D mm':>5s} {'r_det':>5s} {'z_ap':>5s} {'gap':>5s} {'cross-talk (6 nb)':>18s} {'eff (Omega/4pi, own patch)':>27s} {'penumbra r_p':>12s}")
for rw in rows:
    if rw[3] in (0.5, 1.0) and rw[2] in (0.5, 2.0):
        out(f"  {rw[0]:5.0f} {rw[1]:5.1f} {rw[2]:5.1f} {rw[3]:5.2f} {rw[4]*100:17.2f}% {rw[5]*100:26.3f}% {rw[6]:12.2f}")
# recommended
out("\n  Minimum gap for cross-talk <= 1 % (all 6 neighbours combined):")
for D in (10.0, 20.0, 30.0):
    for r_det in (3.0,):
        for z_ap in (0.5, 1.0, 2.0):
            ok = [rw for rw in rows if rw[0] == D and rw[1] == r_det and rw[2] == z_ap and rw[4] <= 0.01]
            gmin = min([rw[3] for rw in ok]) if ok else np.nan
            eff = [rw[5] for rw in rows if rw[0] == D and rw[1] == r_det and rw[2] == z_ap][0]
            out(f"   D = {D:4.0f} mm, r_det = r_patch = {r_det} mm, aperture {z_ap} mm above membrane: g_min = {gmin} mm; efficiency {eff*100:.2f} %")
# uncollimated comparison
for D, rdet in ((5.0, 12.5), (10.0, 12.5)):
    out(f"  Reference: one uncollimated Ø{2*rdet:.0f} mm detector at {D:.0f} mm, on-axis source: Omega/4pi = {0.5*(1-D/np.hypot(D,rdet))*100:.1f} %")
out("  Minimum sector area: a Ø6 mm patch (0.28 cm^2) with Ø6 mm (28 mm^2) detector at 10 mm gives ~2 % efficiency;")
out("  six such patches + 1 mm gaps fit a Ø22 mm active face (hex: 1 centre + 5, or 2x3 array).")

# figure: acceptance profiles
fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
for D, r_det, z_ap in ((10.0, 3.0, 0.5), (10.0, 3.0, 2.0), (20.0, 3.0, 0.5), (30.0, 3.0, 2.0)):
    A = acceptance(rq, D, r_det, r_det, z_ap, n=60)
    ax[0].plot(rq, A * 100, label=f"D={D:.0f}, r_det={r_det}, aperture z={z_ap} mm")
ax[0].axvline(3.0, color="k", lw=0.7)
ax[0].set(xlabel="source offset from axis (mm)", ylabel="acceptance Omega/4pi (%)", title="Collimated-detector acceptance", xlim=(0, 6))
ax[0].legend(fontsize=7)
for D in (10.0, 20.0, 30.0):
    for z_ap, ls in ((0.5, "-"), (2.0, "--")):
        gg = [rw[3] for rw in rows if rw[0] == D and rw[1] == 3.0 and rw[2] == z_ap]
        ct = [rw[4] * 100 for rw in rows if rw[0] == D and rw[1] == 3.0 and rw[2] == z_ap]
        ax[1].semilogy(gg, ct, ls, marker="o", label=f"D={D:.0f} mm, aperture z={z_ap} mm")
ax[1].axhline(1, color="k", lw=0.7)
ax[1].set(xlabel="inter-sector gap (mm)", ylabel="cross-talk from 6 neighbours (%)", title="r_patch = r_det = 3 mm")
ax[1].legend(fontsize=7)
for a in ax:
    a.grid(alpha=0.3, which="both")
savefig(fig, "m6_sector_collimation.png")
out.save()
