"""
M2 finite-volume solver for primary / secondary current distribution.

Laplace equation div(kappa grad Phi) = 0 in the electrolyte, on a tensor-product
(r,z) grid (axisymmetric, geom='cyl') or (x,y) grid (planar per unit depth,
geom='cart'). Each cell has a material code:
    0 electrolyte, 1 cathode (Pd), 2 anode (Pt), 3 insulator (PTFE/glass/Kel-F)
Faces between electrolyte and:
    electrolyte : conductance A / (d1/k1 + d2/k2)
    insulator / domain edge : zero flux
    anode       : Dirichlet Phi = Va (electrolyte potential at anode surface)
    cathode     : primary  -> Dirichlet Phi = 0
                  secondary-> nonlinear Butler-Volmer-type BC  h(Phi_P - Phi_s) = i(eta),
                               eta = -Phi_s  (cathode metal at 0 V)
Newton iteration on Phi; total current is imposed by root-finding on Va.
Units: cm, S/cm, A, V.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import brentq

EL, CA, AN, INS = 0, 1, 2, 3


def make_edges(keys, hs, hmax, growth=1.12):
    """1-D cell edges passing exactly through each key coordinate.
    hs[k] is the target cell size at keys[k]; sizes grow geometrically away
    from each key up to hmax."""
    keys = np.asarray(keys, float)
    edges = [keys[0]]
    for k in range(len(keys) - 1):
        a, b = keys[k], keys[k + 1]
        L = b - a
        ha, hb = hs[k], hs[k + 1]
        # local target size h(s) = min(ha g^(n), hb g^(m), hmax) built by marching
        pts = [a]
        s = a
        while True:
            da, db = s - a, b - s
            h = min(ha + (growth - 1) * da, hb + (growth - 1) * db, hmax)
            if s + h >= b - 0.5 * min(h, hb):
                break
            s += h
            pts.append(s)
        pts.append(b)
        pts = np.array(pts)
        # rescale interior spacing so the segment closes exactly
        d = np.diff(pts)
        d *= L / d.sum()
        edges.extend(list(a + np.cumsum(d)))
    return np.array(edges)


class FVCell:
    def __init__(self, re, ze, mat, kappa, geom="cyl", Va_mask=None):
        self.re, self.ze = np.asarray(re), np.asarray(ze)
        self.nr, self.nz = len(re) - 1, len(ze) - 1
        self.rc = 0.5 * (self.re[1:] + self.re[:-1])
        self.zc = 0.5 * (self.ze[1:] + self.ze[:-1])
        self.dr = np.diff(self.re)
        self.dz = np.diff(self.ze)
        self.mat = mat.copy()
        self.kap = np.broadcast_to(kappa, mat.shape).astype(float).copy() if np.ndim(kappa) else np.full(mat.shape, float(kappa))
        self.geom = geom
        self._build()

    # ---------------------------------------------------------------- geometry
    def _area_r(self, i_edge, j):
        """area of radial face at edge index i_edge (between cells i_edge-1, i_edge)."""
        if self.geom == "cyl":
            return 2 * np.pi * self.re[i_edge] * self.dz[j]
        return self.dz[j]

    def _area_z(self, i, j_edge):
        if self.geom == "cyl":
            return np.pi * (self.re[i + 1] ** 2 - self.re[i] ** 2)
        return self.dr[i]

    def _build(self):
        nr, nz, mat, kap = self.nr, self.nz, self.mat, self.kap
        idx = -np.ones((nr, nz), int)
        el = mat == EL
        idx[el] = np.arange(el.sum())
        self.idx, self.n = idx, el.sum()
        rows, cols, vals = [], [], []
        b = np.zeros(self.n)            # coefficient multiplying Va
        cath = []                       # (cell, area, h, r, z, orientation)
        self.efaces = []                # (p, q, G) electrolyte-electrolyte for dissipation
        diag = np.zeros(self.n)

        def link(p_ij, q_ij, A, d1, d2, axis):
            (i, j), (k, l) = p_ij, q_ij
            p = idx[i, j]
            mq = mat[k, l]
            if mq == EL:
                q = idx[k, l]
                if q < p:
                    return
                G = A / (d1 / kap[i, j] + d2 / kap[k, l])
                rows.extend([p, p, q, q]); cols.extend([p, q, q, p]); vals.extend([G, -G, G, -G])
                self.efaces.append((p, q, G))
            elif mq == AN:
                G = A * kap[i, j] / d1
                diag[p] += G
                b[p] += G
                self.anode_faces.append((p, G))
            elif mq == CA:
                h = kap[i, j] / d1
                # face position
                if axis == "r":
                    rf = self.re[max(i, k)]
                    zf = self.zc[j]
                else:
                    rf = self.rc[i]
                    zf = self.ze[max(j, l)]
                cath.append((p, A, h, rf, zf, axis))

        self.anode_faces = []
        for i in range(nr):
            for j in range(nz):
                if mat[i, j] != EL:
                    continue
                if i + 1 < nr:
                    link((i, j), (i + 1, j), self._area_r(i + 1, j), self.dr[i] / 2, self.dr[i + 1] / 2, "r")
                if i - 1 >= 0 and mat[i - 1, j] != EL:
                    link((i, j), (i - 1, j), self._area_r(i, j), self.dr[i] / 2, self.dr[i - 1] / 2, "r")
                if j + 1 < nz:
                    link((i, j), (i, j + 1), self._area_z(i, j + 1), self.dz[j] / 2, self.dz[j + 1] / 2, "z")
                if j - 1 >= 0 and mat[i, j - 1] != EL:
                    link((i, j), (i, j - 1), self._area_z(i, j), self.dz[j] / 2, self.dz[j - 1] / 2, "z")
        A = sp.coo_matrix((vals, (rows, cols)), shape=(self.n, self.n)).tocsr()
        self.L = (A + sp.diags(diag)).tocsr()
        self.bVa = b
        c = np.array([(p, a, h, r, z) for (p, a, h, r, z, ax) in cath]) if cath else np.zeros((0, 5))
        self.c_cell = c[:, 0].astype(int)
        self.c_area = c[:, 1]
        self.c_h = c[:, 2]
        self.c_r = c[:, 3]
        self.c_z = c[:, 4]
        self.c_axis = np.array([ax for (*_, ax) in cath])
        self.area_cath = self.c_area.sum()

    # ---------------------------------------------------------------- solvers
    def solve_primary(self, Va=1.0):
        h, A, p = self.c_h, self.c_area, self.c_cell
        M = self.L + sp.coo_matrix((A * h, (p, p)), shape=(self.n, self.n)).tocsr()
        phi = spla.spsolve(M.tocsc(), self.bVa * Va)
        i_face = h * phi[p]              # A/cm^2
        return phi, i_face

    def _face_eta(self, phiP, pol):
        """Solve h(phiP + eta) = i(eta) per face by vectorised bisection (eta<=0)."""
        h = self.c_h
        lo = np.full_like(phiP, pol.eta[-1])      # most negative
        hi = np.zeros_like(phiP)
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            i_m, _ = pol.i_of_eta(mid)
            g = h * (phiP + mid) - i_m            # increasing in eta
            hi = np.where(g > 0, mid, hi)
            lo = np.where(g > 0, lo, mid)
        eta = 0.5 * (lo + hi)
        return eta

    def solve_secondary(self, Va, pol, phi0=None, tol=1e-9, maxit=60):
        n, p, A, h = self.n, self.c_cell, self.c_area, self.c_h
        if phi0 is None:
            phi0, _ = self.solve_primary(Va)
        phi = phi0.copy()
        for it in range(maxit):
            eta = self._face_eta(phi[p], pol)
            i, didE = pol.i_of_eta(eta)
            gp = -didE                               # >0
            flux = A * h * (phi[p] + eta)            # = A i
            dflux = A * h * gp / (h + gp)
            res = self.L @ phi - self.bVa * Va
            np.add.at(res, p, flux)
            J = self.L + sp.coo_matrix((dflux, (p, p)), shape=(n, n)).tocsr()
            dphi = spla.spsolve(J.tocsc(), -res)
            phi += dphi
            if np.max(np.abs(dphi)) < tol * max(1.0, Va):
                break
        eta = self._face_eta(phi[p], pol)
        i, _ = pol.i_of_eta(eta)
        return phi, i, eta

    def run(self, I_target, pol=None):
        """Impose total cathode current I_target (A). pol=None -> primary."""
        if pol is None:
            phi1, i1 = self.solve_primary(1.0)
            I1 = np.sum(self.c_area * i1)
            s = I_target / I1
            return dict(Va=s, phi=phi1 * s, i=i1 * s, eta=np.zeros_like(i1))
        cache = {}

        def Itot(Va):
            phi0 = cache.get("phi")
            phi, i, eta = self.solve_secondary(Va, pol, phi0=phi0 if phi0 is not None else None)
            cache["phi"] = phi
            cache["res"] = (phi, i, eta)
            return np.sum(self.c_area * i)

        # bracket
        eta_guess = pol.eta_of_i(I_target / self.area_cath)
        phi1, i1 = self.solve_primary(1.0)
        Rprim = 1.0 / np.sum(self.c_area * i1)          # ohm (primary resistance)
        Va_lo = -eta_guess * 0.5
        Va_hi = -eta_guess + I_target * Rprim * 1.5 + 0.2
        cache["phi"] = phi1 * Va_hi
        f = lambda V: np.log(Itot(V) / I_target)
        while f(Va_hi) < 0:
            Va_hi *= 2
        Va = brentq(f, Va_lo, Va_hi, xtol=1e-7, rtol=1e-7)
        Itot(Va)
        phi, i, eta = cache["res"]
        return dict(Va=Va, phi=phi, i=i, eta=eta, Rprim=Rprim)

    def dissipation(self, phi, Va, i_face=None):
        """Ohmic power in the electrolyte (W): sum over faces G dPhi^2 (+anode half cells)."""
        P = 0.0
        for (p, q, G) in self.efaces:
            P += G * (phi[p] - phi[q]) ** 2
        for (p, G) in self.anode_faces:
            P += G * (phi[p] - Va) ** 2
        if i_face is not None:                      # cathode half-cells
            P += np.sum(self.c_area * i_face**2 / self.c_h)
        return P

    def field(self, phi):
        out = np.full((self.nr, self.nz), np.nan)
        out[self.mat == EL] = phi
        return out
