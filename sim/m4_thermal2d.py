"""
M4 — axisymmetric (r,z) finite-volume thermal model of a closed electrolytic
cell inside four calorimeter types.

  ISO   : isoperibolic "metal Dewar" (F&P style): cell in a blackened vacuum gap
          to a bath, sensor = thermistor in the electrolyte, lid leaks via neck+leads.
  FLOW  : mass-flow: cell sleeve in contact with a water-cooled jacket (sides+bottom),
          lid insulated with foam; calibration constant = captured fraction.
  SEEB0 : Seebeck envelope, no heat spreader: cell sits in air inside a box whose
          walls are thermoelectric modules (Storms-type).
  SEEB1 : Seebeck envelope with an isothermal Cu shell clamped to the cell by a
          gap pad, 4-pi module coverage, all leads thermally anchored to the shell.
  SEEB2 : as SEEB1 but leads run straight from the lid through the wall (not anchored).

For every design and for each heat-source location (cathode, anode, electrolyte
Joule heating, recombiner, electrolyte surface, calibration heaters) the model
returns the calorimeter's calibration constant, so that the position dependence
(the physical basis of Shanahan's CCS critique) can be quantified.

Run:  python3 sim/m4_thermal2d.py
"""
import os
import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, os.path.dirname(__file__))
from m4_params import MAT, TEC_ALPHA, TEC_K, TEC_A, TEC_T, SIGMA_SB  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), "..", "docs", "models", "figs")
os.makedirs(OUT, exist_ok=True)
mm = 1e-3


# ============================================================ generic FV solver
def make_edges(breaks, h, hmin=None):
    """Cell edges containing every breakpoint, spacing <= h."""
    b = np.unique(np.round(np.asarray(breaks, float), 9))
    e = [b[0]]
    for a, c in zip(b[:-1], b[1:]):
        n = max(1, int(np.ceil((c - a) / h - 1e-9)))
        e.extend(list(np.linspace(a, c, n + 1)[1:]))
    return np.array(e)


class Axi:
    """Axisymmetric FV conduction on a tensor grid with per-cell (kr, kz, rho*c)."""

    def __init__(self, re, ze):
        self.re, self.ze = re, ze
        self.rc, self.zc = 0.5 * (re[1:] + re[:-1]), 0.5 * (ze[1:] + ze[:-1])
        self.nr, self.nz = len(self.rc), len(self.zc)
        self.dr, self.dz = np.diff(re), np.diff(ze)
        n = (self.nr, self.nz)
        self.kr = np.full(n, np.nan)
        self.kz = np.full(n, np.nan)
        self.rhoc = np.full(n, np.nan)
        self.label = np.full(n, "", dtype=object)
        self.vol = np.pi * (re[1:] ** 2 - re[:-1] ** 2)[:, None] * self.dz[None, :]
        self.dirichlet = np.zeros(n, bool)      # cell held at T = 0 (sink)
        self.links = []                         # (mask, G_total) conductance to T=0

    def paint(self, r0, r1, z0, z1, mat=None, k=None, rhoc=None, label=None, kr=None, kz=None):
        R, Z = np.meshgrid(self.rc, self.zc, indexing="ij")
        m = (R > r0) & (R < r1) & (Z > z0) & (Z < z1)
        if mat is not None:
            k, rhoc = MAT[mat]
        if k is not None:
            self.kr[m] = k
            self.kz[m] = k
        if kr is not None:
            self.kr[m] = kr
        if kz is not None:
            self.kz[m] = kz
        if rhoc is not None:
            self.rhoc[m] = rhoc
        if label is not None:
            self.label[m] = label
        return m

    def mask(self, r0, r1, z0, z1):
        R, Z = np.meshgrid(self.rc, self.zc, indexing="ij")
        return (R > r0) & (R < r1) & (Z > z0) & (Z < z1)

    def idx(self, i, j):
        return i * self.nz + j

    def assemble(self):
        assert not np.isnan(self.kr).any(), "unpainted cells"
        nr, nz = self.nr, self.nz
        N = nr * nz
        rows, cols, vals = [], [], []
        diag = np.zeros(N)
        # radial faces
        self.Gr = np.zeros((nr - 1, nz))
        for i in range(nr - 1):
            rf = self.re[i + 1]
            A = 2 * np.pi * rf * self.dz
            R1 = np.log(rf / self.rc[i]) / (2 * np.pi * self.kr[i, :] * self.dz) if self.rc[i] > 0 else 0
            R2 = np.log(self.rc[i + 1] / rf) / (2 * np.pi * self.kr[i + 1, :] * self.dz)
            self.Gr[i, :] = 1.0 / (R1 + R2)
        self.Gz = np.zeros((nr, nz - 1))
        A = self.vol / self.dz[None, :]
        for j in range(nz - 1):
            R1 = 0.5 * self.dz[j] / (self.kz[:, j] * A[:, j])
            R2 = 0.5 * self.dz[j + 1] / (self.kz[:, j + 1] * A[:, j + 1])
            self.Gz[:, j] = 1.0 / (R1 + R2)
        I = np.arange(N).reshape(nr, nz)

        def add(a, b, G):
            rows.extend([a, b, a, b]); cols.extend([b, a, a, b]); vals.extend([-G, -G, G, G])

        for i in range(nr - 1):
            for j in range(nz):
                add(I[i, j], I[i + 1, j], self.Gr[i, j])
        for i in range(nr):
            for j in range(nz - 1):
                add(I[i, j], I[i, j + 1], self.Gz[i, j])
        K = sp.csr_matrix((vals, (rows, cols)), shape=(N, N))
        # links to ground
        self.Glink = np.zeros(N)
        for m, Gtot in self.links:
            w = self.vol[m] / self.vol[m].sum()
            self.Glink[I[m]] += Gtot * w
        K = K + sp.diags(self.Glink)
        # Dirichlet cells: replace rows
        d = self.dirichlet.ravel()
        self.Kfull = K.tocsr()
        Kd = K.tolil()
        for n in np.where(d)[0]:
            Kd.rows[n] = [n]
            Kd.data[n] = [1.0]
        self.K = Kd.tocsr()
        self.C = self.rhoc.ravel() * self.vol.ravel()
        self._lu = spla.splu(self.K.tocsc())

    def solve(self, q):
        b = q.ravel().copy()
        b[self.dirichlet.ravel()] = 0
        return self._lu.solve(b).reshape(self.nr, self.nz)

    def sink_flux(self, T, cells):
        """Heat flowing into the Dirichlet cells selected by boolean mask `cells` (W)."""
        r = self.Kfull @ T.ravel()
        # For a Dirichlet cell at T=0, -(K T)_n is the heat delivered into it.
        return -(r[cells.ravel()]).sum()

    def link_flux(self, T, m=None):
        if m is None:
            return (self.Glink * T.ravel()).sum()
        return (self.Glink * T.ravel())[m.ravel()].sum()

    def transient(self, q, t_end, dt, probe, n_out=400):
        """Implicit Euler from T=0 with constant source q; returns t and probe(T) history."""
        C = self.C.copy()
        d = self.dirichlet.ravel()
        C[d] = 0.0
        A = (sp.diags(C / dt) + self.K).tocsc()
        lu = spla.splu(A)
        T = np.zeros(self.nr * self.nz)
        b0 = q.ravel().copy(); b0[d] = 0
        ts, ys = [0.0], [probe(T.reshape(self.nr, self.nz))]
        nsteps = int(round(t_end / dt))
        every = max(1, nsteps // n_out)
        for s in range(1, nsteps + 1):
            T = lu.solve(C / dt * T + b0)
            if s % every == 0:
                ts.append(s * dt); ys.append(probe(T.reshape(self.nr, self.nz)))
        return np.array(ts), np.array(ys)


# ============================================================ geometry
G = dict(
    Ri=20 * mm,          # electrolyte (liner inner) radius -> 60 mL over 48 mm
    t_liner=1.0 * mm,    # PTFE liner
    t_wall=2.0 * mm,     # 316L wall
    z_bot=4 * mm,        # 316L bottom thickness
    h_el=48 * mm,        # electrolyte height
    h_head=30 * mm,      # headspace
    t_lid=13 * mm,       # 316L lid (CF-flange class)
    r_cath=0.5 * mm, z_c0=15 * mm, z_c1=45 * mm,   # 1 mm dia Pd cathode, 30 mm long
    r_an=14 * mm,        # Pt anode helix radius
    r_rec=10 * mm, h_rec=15 * mm,                  # recombiner basket radius / height
    r_heat=3 * mm,       # sheathed calibration heater 3 mm off axis
)


def cell_z(g):
    z0 = 0.0
    z_liner = z0 + g["z_bot"]
    z_el0 = z_liner + g["t_liner"]
    z_el1 = z_el0 + g["h_el"]
    z_lid0 = z_el1 + g["h_head"]
    z_lid1 = z_lid0 + g["t_lid"]
    return z0, z_liner, z_el0, z_el1, z_lid0, z_lid1


def paint_cell(m, g, k_el, k_head):
    Ri, tl, tw = g["Ri"], g["t_liner"], g["t_wall"]
    Rc = Ri + tl + tw
    z0, zl, ze0, ze1, zlid0, zlid1 = cell_z(g)
    m.paint(0, Rc, z0, zlid1, mat="ss316", label="wall")
    m.paint(0, Ri + tl, zl, zlid0, mat="ptfe", label="liner")
    m.paint(0, Ri, ze0, ze1, k=k_el, rhoc=MAT["d2o"][1], label="el")
    m.paint(0, Ri, ze1, zlid0, k=k_head, rhoc=MAT["gas"][1], label="head")
    # recombiner basket hanging from lid on a 4 mm SS rod
    zr1 = zlid0 - 5 * mm
    zr0 = zr1 - g["h_rec"]
    m.paint(0, g["r_rec"], zr0, zr1, mat="cat", label="rec")
    m.paint(0, 2 * mm, zr1, zlid0, mat="ss316", label="rod")
    m.paint(0, g["r_cath"], ze0 + g["z_c0"] - 5 * mm, ze0 + g["z_c1"] - 5 * mm, mat="pd", label="cath")
    m.paint(g["r_an"], g["r_an"] + 0.5 * mm, ze0 + g["z_c0"] - 5 * mm, ze0 + g["z_c1"] - 5 * mm,
            mat="pt", label="an")
    return Rc, zlid1


def sources(m, g):
    """Unit-power (1 W) source maps for each location."""
    z0, zl, ze0, ze1, zlid0, zlid1 = cell_z(g)
    za, zb = ze0 + g["z_c0"] - 5 * mm, ze0 + g["z_c1"] - 5 * mm
    R, Z = np.meshgrid(m.rc, m.zc, indexing="ij")
    S = {}

    def norm(w):
        w = np.where(np.isfinite(w), w, 0.0)
        return w / w.sum()

    S["cathode"] = norm(m.vol * (m.label == "cath"))
    S["anode"] = norm(m.vol * (m.label == "an"))
    jm = (m.label == "el") & (R > g["r_cath"]) & (R < g["r_an"]) & (Z > za) & (Z < zb)
    S["joule"] = norm(m.vol * jm / np.maximum(R, 1e-6) ** 2)   # J^2/sigma ~ 1/r^2
    S["recomb"] = norm(m.vol * (m.label == "rec"))
    sm = (m.label == "el") & (Z > ze1 - 2 * mm)
    S["surface"] = norm(m.vol * sm)
    hm = (R > g["r_heat"] - 0.5 * mm) & (R < g["r_heat"] + 0.5 * mm) & (Z > za) & (Z < zb)
    S["heater_cath"] = norm(m.vol * hm)
    S["heater_rec"] = S["recomb"].copy()    # heater wound inside the basket
    S["uniform_el"] = norm(m.vol * (m.label == "el"))
    return S


def build(design, g=G, h=1.0 * mm, k_el=10.0, k_head=0.5, shell="cu", t_shell=6 * mm,
          f_side=0.85, f_top=0.60, f_bot=0.85, G_leads=0.012, G_neck=0.030, h_rad=5.0,
          k_pad=3.0, tec_K=TEC_K):
    """Return (model, info) for one design."""
    Ri, tl, tw = g["Ri"], g["t_liner"], g["t_wall"]
    Rc = Ri + tl + tw
    z0, zl, ze0, ze1, zlid0, zlid1 = cell_z(g)
    zr1 = zlid0 - 5 * mm; zr0 = zr1 - g["h_rec"]
    zc0, zc1 = ze0 + g["z_c0"] - 5 * mm, ze0 + g["z_c1"] - 5 * mm
    rb = [0, g["r_cath"], g["r_heat"] - 0.5 * mm, g["r_heat"] + 0.5 * mm, 2 * mm, g["r_rec"], g["r_an"],
          g["r_an"] + 0.5 * mm, Ri, Ri + tl, Rc]
    zb = [z0, zl, ze0, ze1, ze1 - 2 * mm, zlid0, zlid1, zr0, zr1, zc0, zc1]
    info = dict(design=design)

    if design == "ISO":
        gap = 10 * mm
        Rout, zbot, ztop = Rc + gap, z0 - gap, zlid1 + h
        re = make_edges(rb + [Rout], h); ze = make_edges(zb + [zbot, ztop], h)
        m = Axi(re, ze)
        # vacuum gap with blackened surfaces: conductance h_rad per area of the cell surface
        k_side = h_rad * Rc * np.log(Rout / Rc)
        m.paint(0, 1, -1, 1, k=h_rad * gap, rhoc=1.0, label="vac")
        m.paint(Rc, Rout, -1, 1, k=k_side, rhoc=1.0, label="vac")
        paint_cell(m, g, k_el, k_head)
        m.paint(0, 1, zlid1, 1, k=1e-6, rhoc=1.0, label="void")
        # bath: outermost radial column and bottom row held at 0
        m.dirichlet[-1, :] = True
        m.dirichlet[:, 0] = True
        m.dirichlet[:, -1] = True
        lid_top = m.mask(0, Rc, zlid1 - h, zlid1)
        m.links.append((lid_top, G_neck + G_leads))
        info["probe_r"], info["probe_z"] = 10 * mm, ze0 + 24 * mm
    elif design == "FLOW":
        foam = 30 * mm
        Rj = Rc + 0.5 * mm + 2 * mm
        Rout, zbot, ztop = Rj + foam, z0 - 2.5 * mm - foam, zlid1 + foam
        re = make_edges(rb + [Rc + 0.5 * mm, Rj, Rout], h)
        ze = make_edges(zb + [z0 - 0.5 * mm, z0 - 2.5 * mm, zbot, ztop], h)
        m = Axi(re, ze)
        m.paint(0, 1, -1, 1, mat="foam", label="foam")
        m.paint(0, Rc + 0.5 * mm, z0 - 0.5 * mm, zlid0, mat="pad", label="pad")
        jac = m.paint(0, Rj, z0 - 2.5 * mm, zlid0 - 5 * mm, mat="cu", label="jacket")
        m.paint(0, Rc + 0.5 * mm, z0 - 0.5 * mm, zlid0 - 5 * mm, mat="pad", label="pad")
        paint_cell(m, g, k_el, k_head)
        m.paint(0, Rc, zlid0, zlid1, mat="ss316", label="wall")
        m.paint(Rc, Rc + 0.5 * mm, zlid0 - 5 * mm, zlid1, mat="foam", label="foam")
        m.dirichlet[jac & (m.label == "jacket")] = True
        amb = np.zeros_like(m.dirichlet); amb[-1, :] = True; amb[:, 0] = True; amb[:, -1] = True
        m.dirichlet |= amb
        info["jacket"] = (m.label == "jacket") & m.dirichlet
        info["ambient"] = amb & ~info["jacket"]
        lid_top = m.mask(0, Rc, zlid1 - h, zlid1)
        m.links.append((lid_top, G_leads))
        info["probe_r"], info["probe_z"] = 10 * mm, ze0 + 24 * mm
    elif design in ("SEEB0", "SEEB1", "SEEB2"):
        t_pad = 0.5 * mm
        top_gap = 20 * mm
        if design == "SEEB0":
            Rs0 = Rc + 12 * mm      # air gap all round, no spreader
            zs_bot, zs_top = z0 - 12 * mm, zlid1 + top_gap
            Rt0 = Rs0
        else:
            Rs0 = Rc + t_pad
            zs_bot, zs_top = z0 - t_pad, zlid1 + top_gap
            Rt0 = Rs0 + t_shell
        Rt1 = Rt0 + TEC_T
        zt_b1 = zs_bot - (t_shell if design != "SEEB0" else 0)   # outer face of shell bottom
        zt_t0 = zs_top + (t_shell if design != "SEEB0" else 0)
        zt_b0, zt_t1 = zt_b1 - TEC_T, zt_t0 + TEC_T
        hs = 1.0 * mm   # one row of isothermal sink (held at 0) outside the modules
        re = make_edges(rb + [Rs0, Rt0, Rt1, Rt1 + hs], h)
        ze = make_edges(zb + [zs_bot, zs_top, zt_b1, zt_t0, zt_b0, zt_t1, zt_b0 - hs, zt_t1 + hs], h)
        m = Axi(re, ze)
        k_mod = tec_K * TEC_T / TEC_A
        kf = MAT["foam"][0]
        m.paint(0, 1, -1, 1, mat="foam", label="corner")
        # isothermal sink (water-jacketed Cu box) outside the modules
        m.paint(Rt1, 1, -1, 1, mat="cu", label="sink")
        m.paint(0, 1, -1, zt_b0, mat="cu", label="sink")
        m.paint(0, 1, zt_t1, 1, mat="cu", label="sink")
        # TEC layers (anisotropic: modules conduct across, not along)
        m.paint(Rt0, Rt1, zt_b1, zt_t0, kr=f_side * k_mod + (1 - f_side) * kf, kz=0.05,
                rhoc=1.5e6, label="tec_side")
        m.paint(0, Rt0, zt_b0, zt_b1, kz=f_bot * k_mod + (1 - f_bot) * kf, kr=0.05,
                rhoc=1.5e6, label="tec_bot")
        m.paint(0, Rt0, zt_t0, zt_t1, kz=f_top * k_mod + (1 - f_top) * kf, kr=0.05,
                rhoc=1.5e6, label="tec_top")
        if design == "SEEB0":
            m.paint(0, Rt0, zs_bot, zs_top, k=0.05, rhoc=1.2e3, label="air")  # air + radiation
            # cell stands on 3 PEEK feet (lumped as a 5 mm thick 10% fill disc)
            m.paint(0, Rc, zs_bot, z0, k=0.1 * 0.25 + 0.9 * 0.05, rhoc=1.3e5, label="feet")
        else:
            m.paint(0, Rt0, zt_b1, zt_t0, mat=shell, label="shell")
            m.paint(0, Rs0, zs_bot, zs_top, k=0.05, rhoc=1.2e3, label="air")
            m.paint(0, Rs0, zs_bot, z0, k=k_pad, rhoc=2.5e6, label="pad")
            m.paint(Rc, Rs0, z0, zlid1, k=k_pad, rhoc=2.5e6, label="pad")
        paint_cell(m, g, k_el, k_head)
        m.dirichlet[-1, :] = True
        m.dirichlet[:, 0] = True
        m.dirichlet[:, -1] = True
        info.update(f=dict(tec_side=f_side, tec_bot=f_bot, tec_top=f_top), k_mod=k_mod)
        if design == "SEEB1":
            anchor = m.mask(0, Rt0, zs_top, zt_t0) & (m.label == "shell")
            m.links.append((anchor, G_leads))
        else:
            lid_top = m.mask(0, Rc, zlid1 - h, zlid1)
            m.links.append((lid_top, G_leads))
        info["probe_r"], info["probe_z"] = 10 * mm, ze0 + 24 * mm
    else:
        raise ValueError(design)
    m.assemble()
    info["Rc"] = Rc
    return m, info


def signal(m, info, T):
    """Design-specific raw calorimeter output for temperature field T (sink at 0)."""
    d = info["design"]
    if d == "ISO":
        i = np.argmin(abs(m.rc - info["probe_r"])); j = np.argmin(abs(m.zc - info["probe_z"]))
        return T[i, j]                                 # K
    if d == "FLOW":
        return m.sink_flux(T, info["jacket"])          # W captured by coolant
    # Seebeck: V = alpha/A_m * sum f * dT_across * dA  (per unit alpha/A_m -> K m^2)
    V = 0.0
    for lab, f in info["f"].items():
        mk = m.label == lab
        if lab == "tec_side":
            for j in np.where(mk.any(axis=0))[0]:
                ii = np.where(mk[:, j])[0]
                i0 = ii[0] - 1
                q = m.Gr[i0, j] * (T[i0, j] - T[ii[0], j])          # W into layer at this z
                kr = m.kr[ii[0], j]
                V += f * q / (kr / TEC_T)                           # = f * dT * A
        else:
            for i in np.where(mk.any(axis=1))[0]:
                jj = np.where(mk[i, :])[0]
                if lab == "tec_bot":
                    jin = jj[-1] + 1
                    q = m.Gz[i, jj[-1]] * (T[i, jin] - T[i, jj[-1]])
                else:
                    jin = jj[0] - 1
                    q = m.Gz[i, jin] * (T[i, jin] - T[i, jj[0]])
                kz = m.kz[i, jj[0]]
                V += f * q / (kz / TEC_T)
    return V * TEC_ALPHA / TEC_A        # volts


# ============================================================ analysis
POS = ["cathode", "anode", "joule", "surface", "recomb", "heater_cath", "heater_rec"]


def cal_constants(design, **kw):
    m, info = build(design, **kw)
    S = sources(m, G)
    res = {}
    for p in POS + ["uniform_el"]:
        T = m.solve(S[p])
        res[p] = signal(m, info, T)
        if p == "cathode":
            i = np.argmin(abs(m.rc - info["probe_r"])); j = np.argmin(abs(m.zc - info["probe_z"]))
            res["_Tel_cathode"] = T[i, j]
            res["_Tmax"] = T.max()
    return res, m, info


def rel(res, ref="heater_cath"):
    return {p: res[p] / res[ref] - 1 for p in POS}


def verify_line_source(txt):
    """Analytic check: line source on axis of a long cylinder, Dirichlet outer wall.
    T(r) = q'/(2 pi k) ln(R/r).  Also energy conservation."""
    R, L, k = 20 * mm, 400 * mm, 1.0
    out = []
    for h in [2 * mm, 1 * mm, 0.5 * mm]:
        re = make_edges([0, 0.5 * mm, R], h); ze = make_edges([0, L], 4 * mm)
        m = Axi(re, ze)
        m.paint(0, 1, -1, 1, k=k, rhoc=1.0)
        m.dirichlet[-1, :] = True
        m.assemble()
        q = np.zeros((m.nr, m.nz)); q[0, :] = 1.0 / m.nz     # 1 W over L -> q' = 2.5 W/m
        T = m.solve(q)
        j = m.nz // 2
        r_test = m.rc[m.nr // 2]
        # adiabatic ends -> 1D radial; last cell centre held at 0 so reference R_eff = rc[-1]
        Tan = (1 / L) / (2 * np.pi * k) * np.log(m.rc[-1] / r_test)
        err = T[m.nr // 2, j] / Tan - 1
        cons = m.sink_flux(T, m.dirichlet) - 1.0
        out.append((h / mm, err, cons))
        txt.append(f"  line-source test h={h/mm:.1f} mm: T(r={r_test/mm:.1f} mm) rel.err={err:+.2e}, "
                   f"energy balance err={cons:+.1e} W")
    return out


def verify_finite_cylinder(txt):
    """2D analytic check: uniform source q in a finite cylinder (radius a, height L), T=0 on all faces.
    T(r,z) = sum_{m odd} 4/(m pi) * q L^2/(k m^2 pi^2) * sin(m pi z/L) * [1 - I0(m pi r/L)/I0(m pi a/L)]"""
    from scipy.special import i0e
    k = 1.0
    for h in [2 * mm, 1 * mm, 0.5 * mm]:
        re = make_edges([0, 30 * mm], h); ze = make_edges([0, 60 * mm], h)
        m = Axi(re, ze)
        m.paint(0, 1, -1, 1, k=k, rhoc=1.0)
        m.dirichlet[-1, :] = True; m.dirichlet[:, 0] = True; m.dirichlet[:, -1] = True
        m.assemble()
        q = m.vol * (~m.dirichlet); q = q / q.sum()          # 1 W uniform in the interior
        qv = 1.0 / m.vol[~m.dirichlet].sum()
        T = m.solve(q)
        a = m.rc[-1]; z0, z1 = m.zc[0], m.zc[-1]; L = z1 - z0
        i, j = 0, m.nz // 2
        r, z = m.rc[i], m.zc[j] - z0
        Tan = 0.0
        for mm_ in range(1, 400, 2):
            b = mm_ * np.pi / L
            ratio = i0e(b * r) / i0e(b * a) * np.exp(b * (r - a))
            Tan += 4 / (mm_ * np.pi) * qv / (k * b * b) * np.sin(b * z) * (1 - ratio)
        txt.append(f"  finite-cylinder 2D test h={h/mm:.1f} mm: T(axis, mid-height) FV={T[i, j]:.5f} K, "
                   f"analytic={Tan:.5f} K, rel.err={T[i, j]/Tan-1:+.2e}; energy balance "
                   f"{m.sink_flux(T, m.dirichlet)-1:+.1e} W")


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    txt = ["M4 thermal 2D (axisymmetric FV) results", "=" * 60, "", "Verification"]
    verify_line_source(txt)
    verify_finite_cylinder(txt)

    # ---- grid convergence for the headline numbers
    txt.append("")
    txt.append("Grid convergence: recombiner-vs-cathode-heater calibration offset (%)")
    conv = {}
    for d in ["ISO", "FLOW", "SEEB1"]:
        row = []
        for h in [2 * mm, 1 * mm, 0.5 * mm]:
            r, _, _ = cal_constants(d, h=h)
            row.append(100 * (r["recomb"] / r["heater_cath"] - 1))
        conv[d] = row
        txt.append(f"  {d:6s} h=2,1,0.5 mm: " + ", ".join(f"{x:+.4f}" for x in row))

    # ---- baseline table
    designs = ["ISO", "FLOW", "SEEB0", "SEEB2", "SEEB1"]
    base = {}
    txt += ["", "Baseline calibration constants (k_el=10 W/mK, k_head=0.5 W/mK)",
            "  units: ISO K/W (probe), FLOW captured fraction, SEEB V/W"]
    hdr = "  design  " + " ".join(f"{p:>11s}" for p in POS) + "   Tel(K/W)  Tmax(K/W)"
    txt.append(hdr)
    for d in designs:
        r, m, info = cal_constants(d)
        base[d] = r
        # energy check
        txt.append(f"  {d:6s}  " + " ".join(f"{r[p]:11.5g}" for p in POS) +
                   f"   {r['_Tel_cathode']:.3g}   {r['_Tmax']:.3g}")
    txt += ["", "Position dependence relative to the cathode-position calibration heater (%)"]
    txt.append("  design  " + " ".join(f"{p:>11s}" for p in POS))
    for d in designs:
        rr = rel(base[d])
        txt.append(f"  {d:6s}  " + " ".join(f"{100*rr[p]:+11.4f}" for p in POS))

    # ---- sensitivity sweeps
    txt += ["", "Sweep: electrolyte effective conductivity (stirring) k_el -> recomb/cathode offset (%) "
            "and cathode-cal-constant change vs k_el=10 (%)"]
    kels = [1, 3, 10, 30, 100]
    sweep_k = {d: [] for d in designs}
    for d in designs:
        row = []
        for ke in kels:
            r, _, _ = cal_constants(d, k_el=ke)
            row.append((100 * (r["recomb"] / r["heater_cath"] - 1),
                        100 * (r["cathode"] / base[d]["cathode"] - 1)))
        sweep_k[d] = row
        txt.append(f"  {d:6s} " + "  ".join(f"k={k:>3}: {a:+.3f}/{b:+.3f}" for k, (a, b) in zip(kels, row)))

    txt += ["", "Sweep: headspace effective conductivity k_head (condensation heat-pipe) -> recomb offset (%)"]
    khs = [0.1, 0.5, 2, 10]
    for d in designs:
        row = []
        for kh in khs:
            r, _, _ = cal_constants(d, k_head=kh)
            row.append(100 * (r["recomb"] / r["heater_cath"] - 1))
        txt.append(f"  {d:6s} " + "  ".join(f"kh={k}: {a:+.3f}" for k, a in zip(khs, row)))

    txt += ["", "SEEB1 design sweeps (recomb and joule offsets vs cathode heater, %)"]
    sweeps = []
    for lab, kw in [("baseline Cu 6 mm, f_top 0.60", {}),
                    ("Al 6 mm shell", dict(shell="al6061")),
                    ("Cu 3 mm shell", dict(t_shell=3 * mm)),
                    ("Cu 10 mm shell", dict(t_shell=10 * mm)),
                    ("uniform coverage f=0.85 all faces", dict(f_top=0.85)),
                    ("top coverage f=0.40", dict(f_top=0.40)),
                    ("uniform coverage + Cu 10 mm", dict(f_top=0.85, t_shell=10 * mm)),
                    ("pad k=1 W/mK", dict(k_pad=1.0)),
                    ("leads 5x (G=0.06 W/K)", dict(G_leads=0.06)),
                    ("low-K modules (0.25 W/K)", dict(tec_K=0.25)),
                    ("low-K modules + uniform + Cu10", dict(tec_K=0.25, f_top=0.85, t_shell=10 * mm))]:
        r, m, info = cal_constants("SEEB1", **kw)
        rr = rel(r)
        sweeps.append((lab, r["heater_cath"], rr["recomb"], rr["joule"], rr["surface"], r["_Tel_cathode"]))
        txt.append(f"  {lab:38s} S={r['heater_cath']*1e3:7.2f} mV/W  recomb {100*rr['recomb']:+.4f}  "
                   f"joule {100*rr['joule']:+.4f}  surface {100*rr['surface']:+.4f}  "
                   f"Tel {r['_Tel_cathode']:.2f} K/W")
    txt += ["", "SEEB2 (leads not anchored) sweep of lead conductance"]
    for Gl in [0.003, 0.012, 0.03, 0.06]:
        r, _, _ = cal_constants("SEEB2", G_leads=Gl)
        rr = rel(r)
        txt.append(f"  G_leads={Gl:.3f} W/K: recomb {100*rr['recomb']:+.4f} %  joule {100*rr['joule']:+.4f} %")

    # ---- transients (step response of signal) for the cathode heater and recombiner
    txt += ["", "Step-response time constants (s): t63, t99, t99.9 of final signal"]
    tr = {}
    for d in ["ISO", "FLOW", "SEEB1", "SEEB0"]:
        m, info = build(d)
        S = sources(m, G)
        tr[d] = {}
        for p in ["heater_cath", "recomb"]:
            final = signal(m, info, m.solve(S[p]))
            t, y = m.transient(S[p], t_end=40000 if d in ("ISO", "SEEB0") else 12000, dt=10.0,
                               probe=lambda T: signal(m, info, T))
            y = y / final
            def cross(level):
                k = np.argmax(y >= level)
                return t[k] if y[k] >= level else np.nan
            t63, t99, t999 = cross(1 - np.exp(-1)), cross(0.99), cross(0.999)
            tr[d][p] = (t, y, t63, t99, t999)
            txt.append(f"  {d:6s} {p:12s} t63={t63:7.0f}  t99={t99:7.0f}  t99.9={t999:7.0f}")
        C_tot = np.nansum(m.C[~m.dirichlet.ravel()])
        txt.append(f"  {d:6s} total heat capacity inside sink = {C_tot:.0f} J/K")

    with open(os.path.join(OUT, "m4_thermal2d.txt"), "w") as fh:
        fh.write("\n".join(txt) + "\n")
    print("\n".join(txt))

    # ---- figures
    fig, ax = plt.subplots(figsize=(8, 4.2))
    x = np.arange(len(POS))
    w = 0.16
    for n, d in enumerate(designs):
        rr = rel(base[d])
        ax.bar(x + (n - 2) * w, [100 * rr[p] for p in POS], w, label=d)
    ax.axhspan(-0.1, 0.1, color="0.85", zorder=0, label="±0.1 % target")
    ax.set_xticks(x); ax.set_xticklabels(POS, rotation=20)
    ax.set_ylabel("calibration constant vs cathode heater (%)")
    ax.set_yscale("symlog", linthresh=0.1)
    ax.legend(fontsize=8, ncol=3)
    ax.set_title("Position dependence of calibration constant (k_el=10, k_head=0.5 W/mK)")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "m4_position_dependence.png"), dpi=130)

    fig, ax = plt.subplots(figsize=(6.5, 4))
    for d in designs:
        ax.plot(kels, [abs(a) for a, b in sweep_k[d]], "o-", label=d)
    ax.axhline(0.1, color="k", ls=":", label="0.1 %")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("effective electrolyte conductivity k_el (W/m/K)  [stirring]")
    ax.set_ylabel("|recombiner − cathode| cal. offset (%)")
    ax.legend(fontsize=8); ax.set_title("Recombiner-vs-cathode offset vs stirring")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "m4_offset_vs_stirring.png"), dpi=130)

    fig, ax = plt.subplots(figsize=(6.5, 4))
    for d in tr:
        t, y = tr[d]["heater_cath"][:2]
        ax.plot(t / 3600, y, label=f"{d} cathode heater")
        t, y = tr[d]["recomb"][:2]
        ax.plot(t / 3600, y, "--", label=f"{d} recombiner")
    ax.set_xlabel("time (h)"); ax.set_ylabel("signal / final"); ax.set_xlim(0, 6)
    ax.legend(fontsize=7); ax.set_title("Step response")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "m4_step_response.png"), dpi=130)

    # temperature map for SEEB1 with recombiner source
    m, info = build("SEEB1")
    S = sources(m, G)
    T = m.solve(S["recomb"] * 1.5 + S["joule"] * 1.0 + S["cathode"] * 0.5 + S["anode"] * 0.5)
    fig, ax = plt.subplots(figsize=(4.2, 6))
    pc = ax.pcolormesh(m.re / mm, m.ze / mm, T.T, shading="flat", cmap="inferno")
    fig.colorbar(pc, label="T − T_sink (K)")
    ax.set_xlabel("r (mm)"); ax.set_ylabel("z (mm)"); ax.set_aspect("equal")
    ax.set_title("SEEB1, 3.5 W (1 A, 3.5 V): T field")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "m4_seeb1_Tmap.png"), dpi=130)
    return base, sweeps, tr


if __name__ == "__main__":
    main()
