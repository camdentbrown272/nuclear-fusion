"""
M2 geometry builders for the four cell families of the brief (cm units).

  (a) coax_wire   : Pd wire (radius a, length L) inside cylindrical Pt mesh anode
  (b) planar_foil : Pd foil (thickness t, width W) between two parallel anode meshes
                    (planar 2-D cross-section per unit depth, symmetric quarter)
  (c) dfm_disk    : detector-facing Pd disk (radius Rd) forming the cell floor,
                    anode disk / mesh / ring at gap g, optional insulating well
                    (collar), optional seal-step defect
  (d) sphere      : analytic reference
Each builder returns an FVCell plus a dict of metadata.
"""
import numpy as np
from m2_fv import FVCell, make_edges, EL, CA, AN, INS


def coax_wire(a=0.05, L=3.0, Ra=1.0, kappa=0.012, ends="free", anode_ext=0.5,
              bottom_gap=1.0, top_len=0.3, hfac=1.0, anode_len=None):
    """ends: 'free'   -> wire tip hangs free at z = bottom_gap above cell floor;
                         anode spans [bottom_gap - anode_ext, top] (or anode_len);
             'plates' -> wire spans between two insulating plates at z=0 and z=L,
                         anode spans the same height (the ideal 1-D coax);
             'sleeve' -> free tip covered by a PTFE sleeve of 3a outer radius, 2 mm long.
    Wire top: goes through the free surface (no flux) at z = H."""
    t_an = 0.05
    if ends == "plates":
        H = L
        zk = [0.0, L]
        ztip = 0.0
    else:
        ztip = bottom_gap
        H = ztip + L + top_len * 0      # wire emerges through free surface at z=H
        zk = [0.0, ztip, H]
    za0 = max(0.0, ztip - anode_ext) if ends != "plates" else 0.0
    za1 = H
    if anode_len is not None and ends != "plates":
        za0 = max(0.0, ztip + L / 2 - anode_len / 2)
        za1 = min(H, ztip + L / 2 + anode_len / 2)
    zk = sorted(set(zk + [za0, za1]))
    if ends == "sleeve":
        zk = sorted(set(zk + [ztip - 0.2]))
    hmin = a / 6 * hfac
    rkeys = [0.0, a, Ra, Ra + t_an]
    rhs = [a / 3 * hfac, hmin, 0.02 * hfac, 0.02 * hfac]
    if ends == "sleeve":
        rkeys = [0.0, a, 3 * a, Ra, Ra + t_an]
        rhs = [a / 3 * hfac, hmin, hmin, 0.02 * hfac, 0.02 * hfac]
    re_ = make_edges(rkeys, rhs, 0.08 * hfac)
    zhs = []
    for z in zk:
        zhs.append(hmin if (abs(z - ztip) < 1e-12 and ends != "plates") else 0.03 * hfac)
    ze = make_edges(zk, zhs, 0.1 * hfac)
    rc, zc = 0.5 * (re_[1:] + re_[:-1]), 0.5 * (ze[1:] + ze[:-1])
    R, Z = np.meshgrid(rc, zc, indexing="ij")
    mat = np.full(R.shape, EL)
    mat[(R < a) & (Z > ztip)] = CA
    mat[(R > Ra)] = INS
    mat[(R > Ra) & (Z > za0) & (Z < za1)] = AN
    if ends == "sleeve":
        mat[(R < 3 * a) & (Z > ztip - 0.2) & (Z < ztip)] = INS
    cell = FVCell(re_, ze, mat, kappa, "cyl")
    return cell, dict(kind="coax", a=a, L=L, Ra=Ra, ends=ends, ztip=ztip, H=H)


def planar_foil(t=0.01, W=1.0, g=0.5, kappa=0.012, framed=False, Wcell=None, Wa=None, hfac=1.0):
    """Quarter-domain: x in [0, Wcell/2], y in [0, g+ta]. Foil occupies x<W/2, y<t/2.
    framed=True: foil edges are embedded in an insulating frame whose walls continue
    to the anode (Wcell = W) -> ideal 1-D channel."""
    ta = 0.05
    if Wcell is None:
        Wcell = W if framed else W + 2 * g + 2.0
    if Wa is None:
        Wa = Wcell
    hmin = min(t / 4, 0.005) * hfac
    xk = sorted(set([0.0, W / 2, Wa / 2, Wcell / 2]))
    xhs = [0.05 * hfac if x < W / 2 - 1e-9 else hmin for x in xk]
    xhs[0] = 0.05 * hfac
    xe = make_edges(xk, xhs, 0.05 * hfac)
    yk = [0.0, t / 2, g, g + ta]
    ye = make_edges(yk, [t / 6, hmin, 0.02 * hfac, 0.02 * hfac], 0.05 * hfac)
    xc, yc = 0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1])
    X, Y = np.meshgrid(xc, yc, indexing="ij")
    mat = np.full(X.shape, EL)
    mat[(X < W / 2) & (Y < t / 2)] = CA
    if framed:
        mat[(X > W / 2) & (Y < g)] = INS
    mat[Y > g] = INS
    mat[(Y > g) & (X < Wa / 2)] = AN
    cell = FVCell(xe, ye, mat, kappa, "cart")
    return cell, dict(kind="foil", t=t, W=W, g=g, framed=framed)


def dfm_disk(Rd=1.0, g=1.0, kappa=0.012, Rc=None, anode="disk", Ran=None, well=0.0,
             H=None, step=0.0, crevice=0.02, hfac=1.0, sag=0.0):
    """Pd disk (radius Rd) at the cell floor z=0, wetted face up, back face in vacuum.
    Rc     : cell (electrolyte) radius; Rc=Rd -> 'tube' cell (insulating wall meets disk at 90 deg)
    anode  : 'disk'  -> thin plate/mesh of radius Ran at z=g (default Ran=Rc)
             'ring'  -> ring of 1 mm square section, centred at radius Ran (default Rc-0.1) and z=g
    well   : depth of an insulating collar of inner radius Rd above the disk (cm);
             the cell widens to Rc above it. (well=0 & Rc>Rd -> flush disk in wide floor)
    sag    : dished mesh anode, anode surface at z = g - sag (1 - r^2/Rc^2) (cm) (parallelism test)
    step   : seal-step defect for the tube cell (cm). step>0: insulating wall bore is
             Rd-step and the disk continues under it in a crevice of height `crevice`
             (gasket gap). step<0: bore is Rd+|step|, leaving a coplanar insulating annulus.
    """
    if Rc is None:
        Rc = Rd
    if H is None:
        H = g + 0.5
    if Ran is None:
        Ran = Rc if anode == "disk" else Rc - 0.1
    tA = 0.05
    Rbore = Rd - step if step > 0 else Rd + abs(step)
    Rcell = max(Rc, Rbore) if step != 0 else Rc
    hmin = 0.004 * hfac
    rk = {0.0, Rd, Rcell}
    if step != 0:
        rk |= {Rbore}
    if anode == "disk":
        rk |= {min(Ran, Rcell)}
    else:
        rk |= {Ran - 0.05, Ran + 0.05}
    rk = sorted(x for x in rk if x <= Rcell + 1e-12)
    rhs = [0.05 * hfac if x == 0.0 else (hmin if abs(x - Rd) < 1e-9 or abs(x - Rbore) < 1e-9 else 0.02 * hfac) for x in rk]
    re_ = make_edges(rk, rhs, 0.06 * hfac)
    zk = {0.0, g, g + tA, H}
    if well > 0:
        zk.add(min(well, g))
    if step > 0:
        zk.add(crevice)
    if anode == "ring":
        zk |= {g - 0.05, g + 0.05}
        zk.discard(g + tA)
    if sag != 0:
        zk |= {g - abs(sag) - 0.01}
    zk = sorted(z for z in zk if z <= H + 1e-12)
    zhs = [hmin if z == 0.0 else 0.02 * hfac for z in zk]
    if sag != 0:
        zhs = [hmin if (z == 0.0 or g - abs(sag) - 0.011 <= z <= g + tA + 1e-9) else 0.02 * hfac for z in zk]
    ze0 = make_edges(zk, zhs, 0.06 * hfac)
    # prepend one layer of solid below z=0 for the disk / floor
    ze = np.concatenate([[-0.01], ze0])
    rc, zc = 0.5 * (re_[1:] + re_[:-1]), 0.5 * (ze[1:] + ze[:-1])
    Rg, Zg = np.meshgrid(rc, zc, indexing="ij")
    mat = np.full(Rg.shape, EL)
    below = Zg < 0
    mat[below] = INS
    mat[below & (Rg < Rd)] = CA
    if well > 0:
        mat[(Rg > Rd) & (Zg > 0) & (Zg < well)] = INS
    if step > 0:        # wall bore Rbore < Rd; crevice of height `crevice` under the wall
        mat[(Rg > Rbore) & (Zg > crevice)] = INS
        mat[(Rg > Rd) & (Zg > 0) & (Zg < crevice)] = INS
    if anode == "disk":
        zlow = g - sag * (1 - (Rg / Rcell) ** 2)
        mat[(Zg > zlow) & (Zg < g + tA) & (Rg < Ran)] = AN
        mat[Zg > g + tA] = INS          # nothing above the (porous) anode carries current
    else:
        mat[(np.abs(Zg - g) < 0.05) & (np.abs(Rg - Ran) < 0.05)] = AN
    cell = FVCell(re_, ze, mat, kappa, "cyl")
    return cell, dict(kind="dfm", Rd=Rd, g=g, Rc=Rc, anode=anode, well=well, step=step)


# --------------------------------------------------------------- analytic helpers
def eccentric_coax_ratio(a, Ra, e):
    """Primary current density on the inner cylinder of two eccentric cylinders
    (radii a < Ra, centre offset e). Exact via bipolar coordinates:
    Phi = C ln|z-c|/|z+c|, poles at +-c. Returns (i_max/i_mean, i_min/i_mean)."""
    if e == 0:
        return 1.0, 1.0
    # centres x1 (inner), x2 (outer) on real axis: x1^2 - a^2 = c^2 = x2^2 - Ra^2, x1 - x2 = e
    # -> x1 = (Ra^2 - a^2 + e^2) / (2e)  ... solve directly
    x1 = (Ra**2 - a**2 - e**2) / (2 * e)
    x2 = x1 + e
    c = np.sqrt(x1**2 - a**2)
    th = np.linspace(0, 2 * np.pi, 4001)
    z = x1 + a * np.exp(1j * th)
    E = np.abs(1 / (z - c) - 1 / (z + c))
    return E.max() / E.mean(), E.min() / E.mean()


def mesh_ripple(pitch, dist):
    """Relative peak-to-mean ripple on a planar electrode facing an array of line
    electrodes (period `pitch`) at distance `dist` (primary distribution, leading
    Fourier mode): 2 exp(-2 pi dist / pitch)."""
    return 2 * np.exp(-2 * np.pi * dist / pitch)


def sphere_resistance(a, b, kappa):
    return (1 / a - 1 / b) / (4 * np.pi * kappa)
