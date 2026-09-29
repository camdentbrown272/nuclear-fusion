"""
M2 — current distribution and loading maps for the four cell families.

Run:  python3 sim/m2_current.py
Writes docs/models/figs/m2_current.txt and m2_*.png (current/loading maps,
DFM uniformity, area-fraction table).

Primary distribution: Laplace with Dirichlet electrodes.
Secondary: nonlinear kinetic BC from m2_common.Kinetics (Volmer-Tafel-Heyrovsky
on Pd in 0.1 M LiOD) -> local i -> f -> x via the beta-PdD isotherm.
Optional bubble screening: Bruggeman kappa_eff = kappa (1-eps)^1.5 with a
drift-flux void fraction in the gas plume / gap (see bubble_* functions).
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from m2_common import (FIGS, SURFACES, Kinetics, Polarization, kappa_LiOD, F, Vm_gas,
                       x_of_f)
from m2_fv import EL, CA, AN, INS
from m2_geom import (coax_wire, planar_foil, dfm_disk, eccentric_coax_ratio,
                     mesh_ripple, sphere_resistance)

OUT = []


def log(s=""):
    print(s)
    OUT.append(s)


KAPPA = kappa_LiOD(0.1, 25.0)          # 0.0121 S/cm
POL = {k: Polarization(Kinetics(**v)) for k, v in SURFACES.items()}
UB = 1.0                               # cm/s effective bubble-swarm rise velocity (0.3-3)


def stats(cell, i):
    A = cell.c_area
    im = np.sum(A * i) / A.sum()
    within = np.sum(A[np.abs(i / im - 1) <= 0.05]) / A.sum()
    return dict(mean=im, max=i.max() / im, min=i.min() / im, w5=within,
                cv=np.sqrt(np.sum(A * (i / im - 1) ** 2) / A.sum()))


def area_frac(cell, i, pol, thr):
    x = pol.x_of_i(np.maximum(i, 1e-12))
    return np.sum(cell.c_area[x >= thr]) / cell.c_area.sum()


# ------------------------------------------------------------------ bubbles
def apply_bubbles_coax(cell, meta, i_face, ub=UB, d0=0.05, spread=0.1, kappa0=KAPPA):
    """D2 plume around a vertical wire: gas flow accumulated from the tip upward,
    plume thickness d(z) = d0 + spread (z - ztip); void eps = jg/(jg+ub)."""
    a, zt = meta["a"], meta["ztip"]
    lat = cell.c_axis == "r"
    zf, Af, i = cell.c_z[lat], cell.c_area[lat], i_face[lat]
    o = np.argsort(zf)
    Itip = np.sum((cell.c_area * i_face)[~lat])
    Qz = (Itip + np.cumsum((Af * i)[o])) * Vm_gas * 1e6 / (2 * F)     # cm^3/s
    zs = zf[o]
    kap = np.full(cell.mat.shape, kappa0)
    R, Z = np.meshgrid(cell.rc, cell.zc, indexing="ij")
    d = d0 + spread * np.clip(Z - zt, 0, None)
    Q = np.interp(Z, zs, Qz, left=0.0)
    Ap = np.pi * ((a + d) ** 2 - a**2)
    jg = Q / Ap
    eps = np.clip(jg / (jg + ub), 0, 0.6)
    inside = (R > a) & (R < a + d) & (Z > zt)
    kap[inside] = kappa0 * (1 - eps[inside]) ** 1.5
    return kap, eps * inside


def apply_bubbles_dfm(cell, meta, I, ub=UB, kappa0=KAPPA):
    """D2 from the upward-facing disk rises as a column of radius Rd through the
    gap to the (porous) anode. Uniform void in r<Rd, 0<z<g."""
    Rd, g = meta["Rd"], meta["g"]
    jg = I / (np.pi * Rd**2) * Vm_gas * 1e6 / (2 * F)
    eps = min(jg / (jg + ub), 0.6)
    kap = np.full(cell.mat.shape, kappa0)
    R, Z = np.meshgrid(cell.rc, cell.zc, indexing="ij")
    kap[(R < Rd) & (Z > 0) & (Z < g)] = kappa0 * (1 - eps) ** 1.5
    return kap, eps


def rebuild(cell, kap):
    from m2_fv import FVCell
    return FVCell(cell.re, cell.ze, cell.mat, kap, cell.geom)


# ================================================================ VERIFICATION
def verification():
    log("=" * 78)
    log("VERIFICATION")
    log("=" * 78)
    # V1 concentric cylinders (plates) : R = ln(Ra/a)/(2 pi kappa L), uniform i
    a, Ra, L = 0.05, 1.0, 3.0
    Ran = np.log(Ra / a) / (2 * np.pi * KAPPA * L)
    for hf in [1.0, 0.5, 0.25]:
        c, m = coax_wire(a=a, L=L, Ra=Ra, ends="plates", hfac=hf)
        r = c.run(1.0)
        s = stats(c, r["i"])
        log(f"V1 coax between plates hfac={hf}: R_num={r['Va']:.4f} ohm, R_analytic={Ran:.4f}"
            f" (err {100*(r['Va']/Ran-1):+.2f}%), max/mean={s['max']:.5f}")
    # V2 Newman disk: i/i_avg = 0.5/sqrt(1-(r/R)^2) for disk in infinite insulating plane
    for hf in [1.0, 0.5]:
        c, m = dfm_disk(Rd=1.0, Rc=8.0, g=8.0, H=8.5, hfac=hf)
        r = c.run(1.0)
        im = np.sum(c.c_area * r["i"]) / c.c_area.sum()
        rr = c.c_r
        errs = []
        for rq in [0.0, 0.5, 0.8, 0.9]:
            k = np.argmin(np.abs(rr - rq))
            errs.append((rr[k], r["i"][k] / im, 0.5 / np.sqrt(1 - rr[k] ** 2)))
        log(f"V2 Newman disk (Rc=8Rd, g=8Rd) hfac={hf}: " +
            "; ".join(f"r={a_:.2f}: {b_:.3f} vs {c_:.3f}" for a_, b_, c_ in errs))
        # Newman resistance of a disk to infinity 1/(4 kappa Rd)
        log(f"    R_num={r['Va']:.3f} ohm vs 1/(4 kappa Rd)={1/(4*KAPPA*1.0):.3f} (finite cell adds ~g/(kappa pi Rc^2)"
            f"={8/(KAPPA*np.pi*64):.3f}) ")
    # V3 sphere-in-sphere analytic
    log(f"V3 sphere a=0.1 cm in b=1.0 cm: R = {sphere_resistance(0.1, 1.0, KAPPA):.2f} ohm, uniform by symmetry")
    # V4 secondary 1-D planar (tube cell) vs analytic series: Va = I R_ohm - eta(i)
    c, m = dfm_disk(Rd=1.0, Rc=1.0, g=1.0)
    for iav in [0.01, 0.1, 0.5]:
        I = iav * np.pi
        r = c.run(I, POL["typical"])
        Va_an = iav * 1.0 / KAPPA - POL["typical"].eta_of_i(iav)
        log(f"V4 1-D secondary (tube, g=1cm) i={iav}: Va_num={r['Va']:.4f} V, analytic={Va_an:.4f} V")
    # V5 eccentric cylinders small-e limit: max/mean ~ 1 + 2 a e/(Ra^2)
    for e in [0.02, 0.05, 0.1, 0.2]:
        mx, mn = eccentric_coax_ratio(0.05, 1.0, e)
        log(f"V5 eccentric coax a=0.5mm Ra=10mm e={10*e:.1f}mm: max/mean={mx:.4f} min/mean={mn:.4f}"
            f"  (small-e estimate 1+2ae/Ra^2={1+2*0.05*e:.4f})")


# ================================================================ SCENARIOS
IAVG = [0.01, 0.05, 0.1, 0.3, 0.5]


def run_case(name, builder, kw, bubbles=None, currents=IAVG, store=None):
    cell, meta = builder(**kw)
    rows = []
    rp = cell.run(1.0)
    sp_ = stats(cell, rp["i"])
    for iav in currents:
        I = iav * cell.area_cath
        c2 = cell
        r = cell.run(I, POL["typical"])
        eps = 0.0
        if bubbles is not None:
            for _ in range(4):           # fixed-point on bubble conductivity
                if bubbles == "coax":
                    kap, epsf = apply_bubbles_coax(cell, meta, r["i"])
                    eps = float(np.nanmax(epsf))
                else:
                    kap, eps = apply_bubbles_dfm(cell, meta, I)
                c2 = rebuild(cell, kap)
                r = c2.run(I, POL["typical"])
        s = stats(c2, r["i"])
        P = c2.dissipation(r["phi"], r["Va"], r["i"])
        fr = {k: (area_frac(c2, r["i"], POL[k], 0.90), area_frac(c2, r["i"], POL[k], 0.95))
              for k in ["typical", "good"]}
        rows.append(dict(iav=iav, I=I, Va=r["Va"], P=P, eps=eps, **s, fr=fr,
                         xmin={k: float(POL[k].x_of_i(r['i'].min())) for k in POL},
                         xmax={k: float(POL[k].x_of_i(r['i'].max())) for k in POL}))
        if store is not None:
            store[(name, iav)] = (c2, r, meta)
    return cell, meta, sp_, rows


def scenarios():
    log("")
    log("=" * 78)
    log(f"SCENARIOS  (0.1 M LiOD, 25 C, kappa={KAPPA*1e3:.2f} mS/cm; secondary kinetics = 'typical'"
        f" surface; bubble u_b={UB} cm/s where noted)")
    log("=" * 78)
    cases = [
        # name, builder, kwargs, bubbles
        ("A1 coax free tip a=0.5mm L=3cm Ra=10mm", coax_wire, dict(a=0.05, L=3, Ra=1.0, ends="free"), None),
        ("A1b  + bubble plume", coax_wire, dict(a=0.05, L=3, Ra=1.0, ends="free"), "coax"),
        ("A2 coax between PTFE end plates", coax_wire, dict(a=0.05, L=3, Ra=1.0, ends="plates"), None),
        ("A2b  + bubble plume", coax_wire, dict(a=0.05, L=3, Ra=1.0, ends="plates"), "coax"),
        ("A3 coax tip in PTFE sleeve", coax_wire, dict(a=0.05, L=3, Ra=1.0, ends="sleeve"), None),
        ("A4 coax free tip a=0.125mm", coax_wire, dict(a=0.0125, L=3, Ra=1.0, ends="free"), None),
        ("A5 coax free tip a=1mm Ra=20mm", coax_wire, dict(a=0.1, L=3, Ra=2.0, ends="free"), None),
        ("A6 coax plates a=1mm L=5 Ra=5mm", coax_wire, dict(a=0.1, L=5, Ra=0.5, ends="plates"), None),
        ("A7 coax free, anode half length", coax_wire, dict(a=0.05, L=3, Ra=1.0, ends="free", anode_len=1.5), None),
        ("B1 foil 100um W=10mm g=5mm open", planar_foil, dict(t=0.01, W=1.0, g=0.5), None),
        ("B2 foil 100um W=10mm g=5mm framed", planar_foil, dict(t=0.01, W=1.0, g=0.5, framed=True), None),
        ("C1 DFM flush, wide floor Rc=3Rd, g=10mm", dfm_disk, dict(Rd=1.0, Rc=3.0, g=1.0), None),
        ("C2 DFM collar well 5mm, Rc=3Rd, g=10mm", dfm_disk, dict(Rd=1.0, Rc=3.0, g=1.0, well=0.5), None),
        ("C3 DFM collar well 10mm, g=15mm", dfm_disk, dict(Rd=1.0, Rc=3.0, g=1.5, well=1.0), None),
        ("C4 DFM tube + mesh anode, g=10mm", dfm_disk, dict(Rd=1.0, Rc=1.0, g=1.0), None),
        ("C4b  + bubble column", dfm_disk, dict(Rd=1.0, Rc=1.0, g=1.0), "dfm"),
        ("C5 DFM tube + mesh anode, g=3mm", dfm_disk, dict(Rd=1.0, Rc=1.0, g=0.3), None),
        ("C6 DFM tube + ring anode, g=10mm", dfm_disk, dict(Rd=1.0, Rc=1.0, g=1.0, anode="ring"), None),
        ("C7 DFM tube + ring anode, g=20mm", dfm_disk, dict(Rd=1.0, Rc=1.0, g=2.0, anode="ring"), None),
        ("C8 DFM tube, small disk anode Ran=Rd/2", dfm_disk, dict(Rd=1.0, Rc=1.0, g=1.0, Ran=0.5), None),
        ("C9 DFM tube, seal overhang 0.5mm (crevice 0.2mm)", dfm_disk, dict(Rd=1.0, g=1.0, step=0.05, crevice=0.02), None),
        ("C10 DFM tube, seal recess 0.5mm (flush annulus)", dfm_disk, dict(Rd=1.0, g=1.0, step=-0.05), None),
    ]
    store = {}
    summary = []
    for name, b, kw, bub in cases:
        cell, meta, sp_, rows = run_case(name, b, kw, bub, store=store)
        summary.append((name, sp_, rows))
        log("")
        log(f"{name}:  area={cell.area_cath:.3f} cm^2 ({'per cm depth' if cell.geom=='cart' else 'total'});"
            f" PRIMARY max/mean={sp_['max']:.3f} min/mean={sp_['min']:.3f} area within +-5%={100*sp_['w5']:.0f}%")
        log("   i_avg   I(A)   Va(V)  P_ohm(W) eps   max/mean min/mean  +-5%  | x>=0.90 typ/good/exc | x>=0.95 typ/good/exc | x range typical")
        for r in rows:
            log(f"   {r['iav']:5.2f} {r['I']:6.3f} {r['Va']:7.2f} {r['P']:8.3f} {r['eps']:.2f}  {r['max']:7.3f} {r['min']:7.3f}  {100*r['w5']:4.0f}% |"
                f" {100*r['fr']['typical'][0]:4.0f}%/{100*r['fr']['good'][0]:4.0f}%/{100*r['fr']['exceptional'][0]:4.0f}% | {100*r['fr']['typical'][1]:4.0f}%/{100*r['fr']['good'][1]:4.0f}%/{100*r['fr']['exceptional'][1]:4.0f}% |"
                f" {r['xmin']['typical']:.3f}-{r['xmax']['typical']:.3f}")
    return summary, store


# ================================================================ PARAMETER SWEEPS
def dfm_sweeps():
    log("")
    log("=" * 78)
    log("DFM SWEEPS (primary distribution; secondary only improves it, Wa<<1 at loading currents)")
    log("=" * 78)
    res = {}
    gaps = [0.2, 0.3, 0.5, 1.0, 1.5, 2.0, 3.0]
    for Rd in [0.5, 1.0, 1.25]:
        for label, kw in [("flush wide", dict(Rc=3 * Rd)), ("ring in tube", dict(Rc=Rd, anode="ring")),
                          ("mesh in tube", dict(Rc=Rd))]:
            out = []
            for g in gaps:
                c, m = dfm_disk(Rd=Rd, g=g, **kw)
                r = c.run(1.0)
                s = stats(c, r["i"])
                out.append((g, s["max"], s["min"], s["w5"]))
            res[(Rd, label)] = out
            log(f"Rd={10*Rd:.1f}mm {label:14s}: " + " ".join(f"g={10*g:.0f}mm:[{mx:.3f},{mn:.3f}]" for g, mx, mn, w in out))
    # collar depth sweep
    wells = [0.0, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0]
    for Rd in [0.5, 1.0, 1.25]:
        out = []
        for w in wells:
            c, m = dfm_disk(Rd=Rd, Rc=3 * Rd, g=w + 0.5, well=w)
            r = c.run(1.0)
            s = stats(c, r["i"])
            out.append((w, s["max"], s["min"], s["w5"]))
        res[(Rd, "well")] = out
        log(f"Rd={10*Rd:.1f}mm collar depth h (anode h+5mm): " + " ".join(f"h={10*w:.1f}mm:[{mx:.3f},{mn:.3f}]" for w, mx, mn, _ in out))
    # seal-step tolerance
    steps = [-0.1, -0.05, -0.02, -0.01, 0.0, 0.01, 0.02, 0.05, 0.1]
    out = []
    for st in steps:
        c, m = dfm_disk(Rd=1.0, g=1.0, step=st, crevice=0.02) if st != 0 else dfm_disk(Rd=1.0, g=1.0)
        r = c.run(1.0)
        s = stats(c, r["i"])
        # fraction of disk area within +-5% and with i >= 0.5 mean
        A = c.c_area
        im = np.sum(A * r["i"]) / A.sum()
        low = np.sum(A[r["i"] < 0.5 * im]) / A.sum()
        out.append((st, s["max"], s["min"], s["w5"], low))
    res["step"] = out
    log("Seal step (Rd=10mm tube, g=10mm; +: wall overhangs disk by s with 0.2mm crevice, -: flush insulating annulus of width s)")
    for st, mx, mn, w5, low in out:
        log(f"   s={10*st:+.1f}mm: max/mean={mx:.3f} min/mean={mn:.3f} area within +-5%={100*w5:.1f}%  area with i<0.5 mean={100*low:.2f}%")
    # ring anode: required gap for +-5% vs Rd
    log("Ring anode in tube: smallest g (mm) giving +-5% over the whole disk")
    for Rd in [0.5, 0.75, 1.0, 1.25]:
        gg = None
        for g in np.arange(0.2, 4.01, 0.1):
            c, m = dfm_disk(Rd=Rd, Rc=Rd, g=g, anode="ring")
            r = c.run(1.0)
            s = stats(c, r["i"])
            if s["max"] <= 1.05 and s["min"] >= 0.95:
                gg = g
                break
        log(f"   Rd={10*Rd:.1f}mm: g_min={'%.0f' % (10*gg) if gg else '>40'} mm  (g/Rd={gg/Rd if gg else float('nan'):.2f})")
    # mesh ripple
    log("Mesh/helix anode ripple on the cathode (primary, leading Fourier mode 2exp(-2pi d/p)):")
    for p in [0.1, 0.2, 0.5, 1.0]:
        log("   pitch %.0f mm: " % (10 * p) + " ".join(f"d={10*d:.0f}mm:{100*mesh_ripple(p, d):.2g}%" for d in [0.2, 0.3, 0.5, 1.0]))
    # eccentricity tolerance for coax
    log("Coax eccentricity tolerance (primary, infinite coax): e for +-5%:")
    for a, Ra in [(0.0125, 0.5), (0.05, 0.5), (0.05, 1.0), (0.1, 1.0), (0.1, 2.0)]:
        es = np.linspace(0.001, Ra - a - 1e-3, 2000)
        ok = [e for e in es if eccentric_coax_ratio(a, Ra, e)[0] <= 1.05]
        log(f"   a={10*a:.2f}mm Ra={10*Ra:.0f}mm : e_max={10*ok[-1]:.2f} mm  (e/Ra={ok[-1]/Ra:.2f})")
    return res


def coax_sweeps():
    log("")
    log("=" * 78)
    log("COAX SWEEPS (free tip): primary tip peak and fraction of wire length within +-5%")
    log("=" * 78)
    for a in [0.0125, 0.025, 0.05, 0.1]:
        for Ra in [0.5, 1.0, 2.0]:
            c, m = coax_wire(a=a, L=3.0, Ra=Ra, ends="free")
            r = c.run(1.0)
            s = stats(c, r["i"])
            r2 = c.run(0.1 * c.area_cath, POL["typical"])
            s2 = stats(c, r2["i"])
            log(f"   a={10*a:.3f}mm Ra={10*Ra:.0f}mm: primary max/mean={s['max']:.2f} within+-5%={100*s['w5']:.0f}% |"
                f" secondary@0.1A/cm2 max/mean={s2['max']:.2f} within+-5%={100*s2['w5']:.0f}%")


def thresholds():
    log("")
    log("=" * 78)
    log("LOADING THRESHOLDS: current density needed for x>=0.90 / 0.95 (uniform, no permeation drain)")
    log("=" * 78)
    for k, pol in POL.items():
        out = []
        for thr in [0.85, 0.90, 0.95]:
            ok = pol.x >= thr
            out.append(f"x>={thr:.2f}: " + (f"{pol.i[ok][0]*1e3:8.1f} mA/cm2" if ok.any() else "   not reached"))
        log(f"   {k:12s} " + " | ".join(out) + f" | x at 1 A/cm2 = {float(pol.x_of_i(1.0)):.3f}")


def wagner():
    log("")
    log("=" * 78)
    log("WAGNER NUMBER Wa = kappa (d eta/d i) / L   (typical kinetics, 0.1 M LiOD)")
    log("=" * 78)
    pol = POL["typical"]
    for iav in [0.001, 0.01, 0.1, 0.5]:
        de = pol.eta_of_i(iav * 1.01) - pol.eta_of_i(iav)
        Rct = -de / (0.01 * iav)                 # ohm cm^2
        ell = KAPPA * Rct
        log(f"   i={iav*1e3:6.1f} mA/cm2: R_ct={Rct:8.3f} ohm cm2, kappa R_ct={10*ell:.3f} mm ->"
            f" Wa(L=a=0.5mm)={ell/0.05:.3f}, Wa(L=Rd=10mm)={ell/1.0:.4f}")


# ================================================================ FIGURES
def figures(summary, store, sweeps):
    pol = POL["typical"]
    # --- coax: i and x along wire
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    for nm, iav, ls in [("A1 coax free tip a=0.5mm L=3cm Ra=10mm", 0.1, "-"),
                        ("A1b  + bubble plume", 0.1, "--"),
                        ("A2 coax between PTFE end plates", 0.1, "-"),
                        ("A2b  + bubble plume", 0.1, "--"),
                        ("A3 coax tip in PTFE sleeve", 0.1, ":")]:
        c, r, m = store[(nm, iav)]
        lat = c.c_axis == "r"
        z = c.c_z[lat] - m["ztip"]
        o = np.argsort(z)
        im = np.sum(c.c_area * r["i"]) / c.c_area.sum()
        ax[0].plot(z[o] * 10, r["i"][lat][o] / im, ls, label=nm.split(" ", 1)[1][:34])
        for k, col in [("typical", None)]:
            ax[1].plot(z[o] * 10, POL[k].x_of_i(r["i"][lat][o]), ls, label=nm.split(" ", 1)[1][:34])
    ax[0].set_xlabel("height above wire tip (mm)"); ax[0].set_ylabel("i / i_mean (lateral surface)")
    ax[0].set_title("Coax wire, i_avg = 100 mA/cm$^2$ (secondary)"); ax[0].legend(fontsize=7); ax[0].grid(alpha=.3)
    ax[1].set_xlabel("height above wire tip (mm)"); ax[1].set_ylabel("x = D/Pd ('typical' surface)")
    ax[1].set_title("Local loading along the wire"); ax[1].grid(alpha=.3)
    c, r, m = store[("A1 coax free tip a=0.5mm L=3cm Ra=10mm", 0.1)]
    phi = c.field(r["phi"])
    Rg, Zg = np.meshgrid(c.rc, c.zc, indexing="ij")
    cs = ax[2].contourf(Rg * 10, Zg * 10, phi, 30, cmap="viridis")
    ax[2].contourf(Rg * 10, Zg * 10, np.where(c.mat == CA, 1, np.nan), colors="k")
    ax[2].contourf(Rg * 10, Zg * 10, np.where(c.mat == AN, 1, np.nan), colors="r")
    plt.colorbar(cs, ax=ax[2], label="electrolyte potential (V)")
    ax[2].set_xlabel("r (mm)"); ax[2].set_ylabel("z (mm)"); ax[2].set_title("A1 potential map (wire black, anode red)")
    fig.tight_layout(); fig.savefig(os.path.join(FIGS, "m2_coax_maps.png"), dpi=130); plt.close(fig)

    # --- DFM radial profiles
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    for nm in ["C1 DFM flush, wide floor Rc=3Rd, g=10mm", "C2 DFM collar well 5mm, Rc=3Rd, g=10mm",
               "C3 DFM collar well 10mm, g=15mm", "C4 DFM tube + mesh anode, g=10mm",
               "C6 DFM tube + ring anode, g=10mm", "C7 DFM tube + ring anode, g=20mm",
               "C8 DFM tube, small disk anode Ran=Rd/2", "C9 DFM tube, seal overhang 0.5mm (crevice 0.2mm)",
               "C10 DFM tube, seal recess 0.5mm (flush annulus)"]:
        c, r, m = store[(nm, 0.1)]
        o = np.argsort(c.c_r)
        im = np.sum(c.c_area * r["i"]) / c.c_area.sum()
        ax[0].plot(c.c_r[o] * 10, r["i"][o] / im, label=nm.split(" ", 1)[1][:40])
        ax[1].plot(c.c_r[o] * 10, POL["good"].x_of_i(r["i"][o]), label=nm.split(" ", 1)[1][:40])
    ax[0].axhspan(0.95, 1.05, color="g", alpha=0.15)
    ax[0].set_ylim(0, 3); ax[0].set_xlabel("r (mm)"); ax[0].set_ylabel("i / i_mean")
    ax[0].set_title("DFM disk Rd=10 mm, i_avg=100 mA/cm$^2$ (secondary); band = ±5%"); ax[0].legend(fontsize=6.5); ax[0].grid(alpha=.3)
    ax[1].set_xlabel("r (mm)"); ax[1].set_ylabel("x = D/Pd ('good' surface)"); ax[1].set_title("Local entry-side loading (no permeation drain)")
    ax[1].axhline(0.9, color="k", lw=.8); ax[1].grid(alpha=.3)
    fig.tight_layout(); fig.savefig(os.path.join(FIGS, "m2_dfm_profiles.png"), dpi=130); plt.close(fig)

    # --- DFM sweeps
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    for key, st in [((1.0, "flush wide"), "-"), ((1.0, "ring in tube"), "--"), ((1.0, "mesh in tube"), ":"),
                    ((0.5, "ring in tube"), "--"), ((1.25, "ring in tube"), "--")]:
        d = np.array(sweeps[key])
        ax[0].semilogy(d[:, 0] * 10, d[:, 1], st, marker="o", label=f"Rd={10*key[0]:.1f}mm {key[1]} max/mean")
    ax[0].axhline(1.05, color="k", lw=.8); ax[0].set_xlabel("anode gap g (mm)"); ax[0].set_ylabel("i_max / i_mean (primary)")
    ax[0].legend(fontsize=7); ax[0].grid(alpha=.3); ax[0].set_title("DFM: edge peak vs gap")
    for Rd in [0.5, 1.0, 1.25]:
        d = np.array(sweeps[(Rd, "well")])
        ax[1].plot(d[:, 0] * 10, d[:, 1], marker="o", label=f"Rd={10*Rd:.1f} mm max/mean")
        ax[1].plot(d[:, 0] * 10, d[:, 2], marker="s", ls="--", label=f"Rd={10*Rd:.1f} mm min/mean")
    ax[1].axhspan(0.95, 1.05, color="g", alpha=.15); ax[1].set_ylim(0.4, 2.5)
    ax[1].set_xlabel("PTFE collar (well) depth h (mm)"); ax[1].set_ylabel("i / i_mean"); ax[1].legend(fontsize=7)
    ax[1].grid(alpha=.3); ax[1].set_title("DFM in wide cell: collar depth")
    d = np.array(sweeps["step"])
    ax[2].plot(d[:, 0] * 10, d[:, 1], "o-", label="max/mean")
    ax[2].plot(d[:, 0] * 10, d[:, 2], "s-", label="min/mean")
    ax[2].plot(d[:, 0] * 10, d[:, 3], "^-", label="area fraction within ±5%")
    ax[2].set_xlabel("seal step s (mm)  (+ overhang/crevice, − flush annulus)"); ax[2].legend(fontsize=7); ax[2].grid(alpha=.3)
    ax[2].set_title("Tube DFM: seal-step tolerance")
    fig.tight_layout(); fig.savefig(os.path.join(FIGS, "m2_dfm_sweeps.png"), dpi=130); plt.close(fig)

    # --- area-fraction bar chart
    fig, ax = plt.subplots(figsize=(12, 5.5))
    names = [s[0] for s in summary]
    yv = np.arange(len(names))
    for k, (iav, col) in enumerate([(0.1, "#4c72b0"), (0.3, "#dd8452"), (0.5, "#55a868")]):
        vals = [next(r for r in s[2] if r["iav"] == iav)["fr"]["good"][0] for s in summary]
        ax.barh(yv + (k - 1) * 0.27, np.array(vals) * 100, 0.27, color=col, label=f"i_avg={int(iav*1e3)} mA/cm$^2$")
    ax.set_yticks(yv); ax.set_yticklabels([n[:48] for n in names], fontsize=7); ax.invert_yaxis()
    ax.set_xlabel("% of cathode area with x ≥ 0.90 ('good' surface)"); ax.legend(fontsize=8); ax.grid(alpha=.3, axis="x")
    fig.tight_layout(); fig.savefig(os.path.join(FIGS, "m2_area_fraction.png"), dpi=130); plt.close(fig)

    # --- foil profile
    fig, ax = plt.subplots(figsize=(6.5, 4))
    for nm in ["B1 foil 100um W=10mm g=5mm open", "B2 foil 100um W=10mm g=5mm framed"]:
        c, r, m = store[(nm, 0.1)]
        face = c.c_axis == "z"
        o = np.argsort(c.c_r[face])
        im = np.sum(c.c_area * r["i"]) / c.c_area.sum()
        ax.plot(c.c_r[face][o] * 10, r["i"][face][o] / im, label=nm[:34])
    ax.set_xlabel("x from foil centre (mm)"); ax.set_ylabel("i / i_mean (broad face)"); ax.legend(fontsize=8); ax.grid(alpha=.3)
    ax.set_title("Planar foil, i_avg=100 mA/cm$^2$")
    fig.tight_layout(); fig.savefig(os.path.join(FIGS, "m2_foil.png"), dpi=130); plt.close(fig)


if __name__ == "__main__":
    verification()
    thresholds()
    wagner()
    summary, store = scenarios()
    sw = dfm_sweeps()
    coax_sweeps()
    figures(summary, store, sw)
    with open(os.path.join(FIGS, "m2_current.txt"), "w") as fh:
        fh.write("\n".join(OUT) + "\n")
