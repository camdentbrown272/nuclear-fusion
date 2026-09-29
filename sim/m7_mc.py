"""M7 simplified coupled e+/e-/photon Monte Carlo (numpy-vectorised).

Physics (all cross-sections from m7_common):
  charged:  condensed history, continuous collision loss (ICRU-37/ESTAR Bethe
            formula) + soft radiative loss (k < kc); discrete bremsstrahlung
            photons for k > kc sampled from the screened Bethe-Heitler
            spectrum; multiple scattering by the Highland formula applied
            incrementally in cumulative path length (PDG eq. 34.15); positron
            annihilation in flight (Heitler) and at rest (2 x 511 keV,
            back-to-back, isotropic; 3-gamma o-Ps decay neglected).
            No delta-ray transport and no energy-loss straggling
            (collisional); range straggling therefore comes only from
            bremsstrahlung and scattering.
  photons:  Woodcock (delta) tracking; photoabsorption (local deposit, no
            fluorescence), Klein-Nishina Compton (recoil electron
            transported), pair production (energy shared uniformly, both
            leptons transported). Coherent scattering neglected.
Geometry: axis-aligned boxes and tubes (annular cylinders along x, y or z);
the first solid in the list that contains a point defines its region.
Units: cm, MeV.
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from m7_common import MATERIALS, ME  # noqa: E402
from m7_ipc import rotate  # noqa: E402

EPS_PUSH = 1e-6      # cm, push across boundaries


# --------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------
class Box:
    def __init__(self, c, h):
        self.c = np.asarray(c, float)
        self.h = np.asarray(h, float)

    def inside(self, P):
        return np.all(np.abs(P - self.c) <= self.h, axis=1)

    def ray(self, P, D):
        t_best = np.full(len(P), np.inf)
        rel = P - self.c
        for k in range(3):
            dk = D[:, k]
            with np.errstate(divide='ignore', invalid='ignore'):
                for sgn in (-1, 1):
                    t = (sgn * self.h[k] - rel[:, k]) / dk
                    ok = (t > 1e-12) & np.isfinite(t)
                    Q = rel + t[:, None] * D
                    o = [j for j in range(3) if j != k]
                    ok &= (np.abs(Q[:, o[0]]) <= self.h[o[0]] + 1e-9) & (np.abs(Q[:, o[1]]) <= self.h[o[1]] + 1e-9)
                    t_best = np.where(ok & (t < t_best), t, t_best)
        return t_best


class Tube:
    def __init__(self, axis, c, rmin, rmax, hl):
        self.a = axis
        self.o = [j for j in range(3) if j != axis]
        self.c = np.asarray(c, float)
        self.rmin, self.rmax, self.hl = rmin, rmax, hl

    def inside(self, P):
        rel = P - self.c
        u = rel[:, self.a]
        r2 = rel[:, self.o[0]] ** 2 + rel[:, self.o[1]] ** 2
        return (np.abs(u) <= self.hl) & (r2 <= self.rmax ** 2) & (r2 >= self.rmin ** 2)

    def ray(self, P, D):
        rel = P - self.c
        u, du = rel[:, self.a], D[:, self.a]
        wx, wy = rel[:, self.o[0]], rel[:, self.o[1]]
        dx, dy = D[:, self.o[0]], D[:, self.o[1]]
        A = dx * dx + dy * dy
        B = wx * dx + wy * dy
        W2 = wx * wx + wy * wy
        t_best = np.full(len(P), np.inf)
        with np.errstate(divide='ignore', invalid='ignore'):
            for R in (self.rmax, self.rmin):
                if R <= 0:
                    continue
                C = W2 - R * R
                disc = B * B - A * C
                okd = (disc >= 0) & (A > 1e-14)
                sq = np.sqrt(np.where(okd, disc, 0))
                for t in ((-B - sq) / A, (-B + sq) / A):
                    ok = okd & (t > 1e-12) & (np.abs(u + t * du) <= self.hl + 1e-9)
                    t_best = np.where(ok & (t < t_best), t, t_best)
            for sgn in (-1, 1):
                t = (sgn * self.hl - u) / du
                r2 = (wx + t * dx) ** 2 + (wy + t * dy) ** 2
                ok = (t > 1e-12) & np.isfinite(t) & (r2 <= self.rmax ** 2 + 1e-9) & (r2 >= self.rmin ** 2 - 1e-9)
                t_best = np.where(ok & (t < t_best), t, t_best)
        return t_best


class Geometry:
    def __init__(self, world_half, fill='AIR'):
        self.world = Box([0, 0, 0], world_half)
        self.fill = fill
        self.solids, self.mats, self.names = [], [], []

    def add(self, solid, mat, name):
        self.solids.append(solid)
        self.mats.append(mat)
        self.names.append(name)

    def finalize(self):
        # region 0 = fill, region i+1 = solid i
        self.reg_names = ['fill'] + self.names
        mats = [self.fill] + self.mats
        self.mat_list = sorted(set(mats), key=lambda s: (s != 'VAC', s))
        self.reg_mat = np.array([self.mat_list.index(m) for m in mats])
        return self

    def locate(self, P):
        reg = np.zeros(len(P), int)
        for i in range(len(self.solids) - 1, -1, -1):
            reg[self.solids[i].inside(P)] = i + 1
        reg[~self.world.inside(P)] = -1
        return reg

    def ray(self, P, D):
        t = self.world.ray(P, D)
        for s in self.solids:
            t = np.minimum(t, s.ray(P, D))
        return t


# --------------------------------------------------------------------------
# Cross-section tables
# --------------------------------------------------------------------------
class Tables:
    def __init__(self, mat_names, kc=0.02, Tmax=40.0):
        self.kc = kc
        self.lT = np.linspace(np.log(1e-3), np.log(Tmax), 300)
        T = np.exp(self.lT)
        self.lE = np.linspace(np.log(1e-3), np.log(Tmax), 300)
        E = np.exp(self.lE)
        n = len(mat_names)
        z = lambda: np.zeros((n, len(T)))  # noqa: E731
        self.rho, self.X0 = np.zeros(n), np.full(n, np.inf)
        self.scol_e, self.scol_p, self.srad_soft, self.sbrem, self.sann, self.range = z(), z(), z(), z(), z(), z()
        self.mu_ph, self.mu_inc, self.mu_pair = z(), z(), z()
        self.mats = []
        for i, nm in enumerate(mat_names):
            m = MATERIALS[nm]
            self.mats.append(m)
            if m is None:
                self.range[i] = np.inf
                continue
            rho = m.rho
            self.rho[i] = rho
            self.X0[i] = m.X0 / rho
            self.scol_e[i] = m.S_col(T) * rho
            self.scol_p[i] = m.S_col(T, positron=True) * rho
            soft, hard = np.zeros(len(T)), np.zeros(len(T))
            for j, t in enumerate(T):
                kk = np.geomspace(1e-5, min(kc, t) * 0.999999, 120)
                soft[j] = np.trapezoid(kk * m.brem_dsdk(t, kk), kk)
                if t > kc:
                    kh = np.geomspace(kc, t * 0.999999, 160)
                    hard[j] = np.trapezoid(m.brem_dsdk(t, kh), kh)
            self.srad_soft[i] = soft * rho
            self.sbrem[i] = hard * rho
            self.sann[i] = m.ann_inflight(T) * rho
            stot = self.scol_e[i] + m.S_rad(T) * rho
            rr = np.concatenate([[T[0] / stot[0]], T[0] / stot[0] + np.cumsum(0.5 * (1 / stot[1:] + 1 / stot[:-1]) * np.diff(T))])
            self.range[i] = rr
            self.mu_ph[i] = m.mu_photo(E) * rho
            self.mu_inc[i] = m.mu_incoh(E) * rho
            self.mu_pair[i] = m.mu_pair(E) * rho
        self.mu_tot = self.mu_ph + self.mu_inc + self.mu_pair
        self.mu_max = self.mu_tot.max(axis=0) * 1.0001 + 1e-12

    def idx(self, lx):
        f = (lx - self.lT[0]) / (self.lT[1] - self.lT[0])
        f = np.clip(f, 0, len(self.lT) - 1.000001)
        i = f.astype(int)
        return i, f - i

    def get(self, tab, m, T):
        i, w = self.idx(np.log(np.clip(T, 1e-3, None)))
        return tab[m, i] * (1 - w) + tab[m, i + 1] * w


# --------------------------------------------------------------------------
# Monte Carlo
# --------------------------------------------------------------------------
class MC:
    def __init__(self, geom, sensitive, rng, kc=0.02, tcut=0.02, gcut=0.01, Tmax=40.0):
        self.g = geom.finalize()
        self.tab = Tables(self.g.mat_list, kc=kc, Tmax=Tmax)
        self.rng = rng
        self.tcut, self.gcut = tcut, gcut
        self.sensitive = list(sensitive)
        self.det_of_reg = np.full(len(self.g.reg_names), -1)
        for d, nm in enumerate(self.sensitive):
            self.det_of_reg[self.g.reg_names.index(nm)] = d

    # --- tallies
    def reset(self, nev):
        self.edep = np.zeros((nev, len(self.sensitive)))
        self.ann = []        # (ev, reg, inflight)
        self.ann_pos = []
        self.brem = []       # (ev, reg, k)
        self.esc_g = []      # (ev, E)
        self.esc_q = []      # (ev, q, T) charged leaving world

    def _dep(self, ev, reg, dE):
        d = self.det_of_reg[np.maximum(reg, 0)]
        ok = (d >= 0) & (reg >= 0) & (dE > 0)
        if ok.any():
            np.add.at(self.edep, (ev[ok], d[ok]), dE[ok])

    def run(self, nev, charged=None, photons=None, max_rounds=50):
        """charged: dict ev,q,P,D,T ; photons: dict ev,P,D,E"""
        self.reset(nev)
        qst = [charged] if charged is not None else []
        gst = [photons] if photons is not None else []
        for _ in range(max_rounds):
            if not qst and not gst:
                break
            newg, newq = [], []
            for c in qst:
                if len(c['ev']):
                    newg += self._charged(**c)
            for p in gst:
                if len(p['ev']):
                    newq += self._photons(**p)
            qst, gst = newq, newg
        self.ann = np.array(self.ann, dtype=float).reshape(-1, 3) if len(self.ann) else np.zeros((0, 3))
        self.ann_pos = np.concatenate(self.ann_pos) if len(self.ann_pos) else np.zeros((0, 3))
        return self

    # ------------------------------------------------------------- charged
    def _charged(self, ev, q, P, D, T):
        g, tb, rng = self.g, self.tab, self.rng
        ev, q, P, D, T = ev.copy(), q.copy(), P.copy(), D.copy(), T.copy()
        xc = np.full(len(ev), 1e-7)       # cumulative path in X0
        out_g = []
        ph = {'ev': [], 'P': [], 'D': [], 'E': []}

        def emit(e, p, d, E):
            ph['ev'].append(e); ph['P'].append(p); ph['D'].append(d); ph['E'].append(E)

        def annihilate_rest(sel):
            if not sel.any():
                return
            e, p, r = ev[sel], P[sel], reg[sel]
            self.ann += [(a, b, 0) for a, b in zip(e, r)]
            self.ann_pos.append(p.copy())
            n = sel.sum()
            c = rng.uniform(-1, 1, n)
            f = rng.uniform(0, 2 * np.pi, n)
            s = np.sqrt(1 - c * c)
            d1 = np.stack([s * np.cos(f), s * np.sin(f), c], 1)
            emit(e, p, d1, np.full(n, ME))
            emit(e.copy(), p.copy(), -d1, np.full(n, ME))

        while len(ev):
            reg = g.locate(P)
            out = reg < 0
            if out.any():
                self.esc_q += list(zip(ev[out], q[out], T[out], D[out, 2]))
                keep = ~out
                ev, q, P, D, T, xc, reg = ev[keep], q[keep], P[keep], D[keep], T[keep], xc[keep], reg[keep]
                if not len(ev):
                    break
            m = g.reg_mat[reg]
            rho = tb.rho[m]
            pos = q > 0
            S = np.where(pos, tb.get(tb.scol_p, m, T), tb.get(tb.scol_e, m, T)) + tb.get(tb.srad_soft, m, T)
            R = tb.get(tb.range, m, T)
            sb = tb.get(tb.sbrem, m, T)
            sa = np.where(pos, tb.get(tb.sann, m, T), 0.0)
            s_geo = g.ray(P, D) + EPS_PUSH
            s_phys = np.where(rho > 0, np.maximum(0.1 * R, 2e-6), np.inf)
            with np.errstate(divide='ignore'):
                s_b = np.where(sb > 0, -np.log(rng.uniform(size=len(ev))) / sb, np.inf)
                s_a = np.where(sa > 0, -np.log(rng.uniform(size=len(ev))) / sa, np.inf)
            s = np.minimum.reduce([s_geo, s_phys, s_b, s_a])
            s = np.where(np.isfinite(s), s, 1.0)
            dT = S * s
            Tm = np.clip(T - 0.5 * dT, 1e-3, None)
            S2 = np.where(pos, tb.get(tb.scol_p, m, Tm), tb.get(tb.scol_e, m, Tm)) + tb.get(tb.srad_soft, m, Tm)
            dT = S2 * s
            stop = dT >= T - self.tcut
            frac = np.where(stop & (dT > 0), np.clip(T / np.maximum(dT, 1e-30), 0, 1), 1.0)
            P = P + D * (s * frac)[:, None]
            dep = np.where(stop, T, dT)
            self._dep(ev, reg, dep)
            T = np.where(stop, 0.0, T - dT)
            # multiple scattering
            ms = (rho > 0) & ~stop
            if ms.any():
                x0 = xc[ms]
                x1 = x0 + s[ms] / tb.X0[m[ms]]
                H = lambda x: x * (1 + 0.038 * np.log(x)) ** 2  # noqa: E731
                var = np.clip(H(x1) - H(x0), 0, None)
                xc[ms] = x1
                Ti = T[ms]
                p = np.sqrt(Ti * (Ti + 2 * ME))
                beta = p / (Ti + ME)
                th0 = 13.6 / (beta * p) * np.sqrt(var)
                th = th0 * np.sqrt(-2 * np.log(rng.uniform(size=ms.sum())))
                cth = np.where(th0 > 1.0, rng.uniform(-1, 1, ms.sum()), np.cos(np.minimum(th, np.pi)))
                D[ms] = rotate(D[ms], cth, rng.uniform(0, 2 * np.pi, ms.sum()))
            # discrete bremsstrahlung
            isb = (~stop) & (s == s_b) & (T > tb.kc)
            if isb.any():
                idx = np.where(isb)[0]
                k = self._sample_brem(m[idx], T[idx])
                E = T[idx] + ME
                u = rng.uniform(size=len(idx))
                thg = np.minimum(np.sqrt(u / (1 - u)) * ME / E, np.pi)
                dg = rotate(D[idx], np.cos(thg), rng.uniform(0, 2 * np.pi, len(idx)))
                emit(ev[idx], P[idx].copy(), dg, k)
                self.brem += list(zip(ev[idx], reg[idx], k))
                T[idx] -= k
            # annihilation in flight
            isa = (~stop) & (s == s_a) & pos & ~isb
            killed = np.zeros(len(ev), bool)
            if isa.any():
                idx = np.where(isa)[0]
                self._ann_flight(idx, ev, P, D, T, reg, emit)
                killed[idx] = True
            # energy cut
            low = (~killed) & (T < self.tcut)
            if low.any():
                self._dep(ev[low], reg[low], T[low])
                T[low] = 0
            done = killed | (T <= 0)
            annihilate_rest(done & ~killed & pos)
            keep = ~done
            ev, q, P, D, T, xc = ev[keep], q[keep], P[keep], D[keep], T[keep], xc[keep]
        if ph['ev']:
            out_g.append({'ev': np.concatenate(ph['ev']), 'P': np.concatenate(ph['P']),
                          'D': np.concatenate(ph['D']), 'E': np.concatenate(ph['E'])})
        return out_g

    def _sample_brem(self, m, T):
        tb, rng = self.tab, self.rng
        k = np.zeros(len(T))
        todo = np.arange(len(T))
        mats = tb.mats
        for _ in range(200):
            if not len(todo):
                break
            mm, TT = m[todo], T[todo]
            kk = tb.kc * (TT / tb.kc) ** rng.uniform(size=len(todo))
            f = np.zeros(len(todo))
            fe = np.zeros(len(todo))
            for mi in np.unique(mm):
                s = mm == mi
                f[s] = kk[s] * mats[mi].brem_dsdk(TT[s], kk[s])
                fe[s] = 1.15 * tb.kc * mats[mi].brem_dsdk(TT[s], np.full(s.sum(), tb.kc))
            acc = rng.uniform(size=len(todo)) * fe <= f
            k[todo[acc]] = kk[acc]
            todo = todo[~acc]
        k[todo] = tb.kc
        return k

    def _ann_flight(self, idx, ev, P, D, T, reg, emit):
        rng = self.rng
        g = T[idx] / ME + 1
        e0 = 1 / (g + 1 + np.sqrt(g * g - 1))
        eps = np.zeros(len(idx))
        todo = np.arange(len(idx))
        for _ in range(200):
            if not len(todo):
                break
            gg, ee0 = g[todo], e0[todo]
            e = ee0 * np.exp(rng.uniform(size=len(todo)) * np.log((1 - ee0) / ee0))
            gf = 1 - e + (2 * gg * e - 1) / (e * (gg + 1) ** 2)
            acc = rng.uniform(size=len(todo)) <= gf
            eps[todo[acc]] = e[acc]
            todo = todo[~acc]
        eps[todo] = 0.5
        Etot = (g + 1) * ME
        k1, k2 = eps * Etot, (1 - eps) * Etot
        sg = np.sqrt(g * g - 1)
        c1 = np.clip((g + 1 - 1 / eps) / sg, -1, 1)
        c2 = np.clip((g + 1 - 1 / (1 - eps)) / sg, -1, 1)
        phi = rng.uniform(0, 2 * np.pi, len(idx))
        d1 = rotate(D[idx], c1, phi)
        d2 = rotate(D[idx], c2, phi + np.pi)
        emit(ev[idx], P[idx].copy(), d1, k1)
        emit(ev[idx].copy(), P[idx].copy(), d2, k2)
        self.ann += [(a, b, 1) for a, b in zip(ev[idx], reg[idx])]
        self.ann_pos.append(P[idx].copy())

    # ------------------------------------------------------------- photons
    def _photons(self, ev, P, D, E):
        g, tb, rng = self.g, self.tab, self.rng
        ev, P, D, E = ev.copy(), P.copy(), D.copy(), E.copy()
        q = {'ev': [], 'q': [], 'P': [], 'D': [], 'T': []}

        def spawn(e, qq, p, d, t):
            ok = t > self.tcut
            if (~ok).any():
                r = g.locate(p[~ok])
                self._dep(e[~ok], r, t[~ok])
                # low-energy positrons: annihilate at rest immediately
                lp = (~ok) & (qq > 0)
                if lp.any():
                    q['ev'].append(e[lp]); q['q'].append(qq[lp]); q['P'].append(p[lp]); q['D'].append(d[lp])
                    q['T'].append(np.full(lp.sum(), 1e-6))
            if ok.any():
                q['ev'].append(e[ok]); q['q'].append(qq[ok]); q['P'].append(p[ok]); q['D'].append(d[ok]); q['T'].append(t[ok])

        while len(ev):
            lE = np.log(np.clip(E, 1e-3, None))
            i, w = tb.idx(lE)
            mumax = tb.mu_max[i] * (1 - w) + tb.mu_max[i + 1] * w
            s = -np.log(rng.uniform(size=len(ev))) / mumax
            P = P + D * s[:, None]
            reg = g.locate(P)
            out = reg < 0
            if out.any():
                self.esc_g += list(zip(ev[out], E[out], D[out, 2]))
            keep = ~out
            ev, P, D, E, reg, i, w = ev[keep], P[keep], D[keep], E[keep], reg[keep], i[keep], w[keep]
            if not len(ev):
                break
            m = g.reg_mat[reg]
            mph = tb.mu_ph[m, i] * (1 - w) + tb.mu_ph[m, i + 1] * w
            minc = tb.mu_inc[m, i] * (1 - w) + tb.mu_inc[m, i + 1] * w
            mpr = tb.mu_pair[m, i] * (1 - w) + tb.mu_pair[m, i + 1] * w
            mt = mph + minc + mpr
            mmx = tb.mu_max[i] * (1 - w) + tb.mu_max[i + 1] * w
            u = rng.uniform(size=len(ev)) * mmx
            real = u < mt
            photo = real & (u < mph)
            comp = real & (u >= mph) & (u < mph + minc)
            pair = real & (u >= mph + minc)
            dead = np.zeros(len(ev), bool)
            if photo.any():
                self._dep(ev[photo], reg[photo], E[photo])
                dead |= photo
            if comp.any():
                idx = np.where(comp)[0]
                eps, cth = self._kn(E[idx])
                phi = rng.uniform(0, 2 * np.pi, len(idx))
                Dn = rotate(D[idx], cth, phi)
                Ee = E[idx] * (1 - eps)
                pe = E[idx][:, None] * D[idx] - (eps * E[idx])[:, None] * Dn
                pe /= np.maximum(np.linalg.norm(pe, axis=1), 1e-30)[:, None]
                spawn(ev[idx], np.full(len(idx), -1), P[idx].copy(), pe, Ee)
                E[idx] = eps * E[idx]
                D[idx] = Dn
            if pair.any():
                idx = np.where(pair)[0]
                avail = E[idx] - 2 * ME
                tp = avail * rng.uniform(size=len(idx))
                for qq, tt in ((1, tp), (-1, avail - tp)):
                    th = np.minimum(ME / (tt + ME) * np.sqrt(-2 * np.log(rng.uniform(size=len(idx)))), np.pi)
                    dd = rotate(D[idx], np.cos(th), rng.uniform(0, 2 * np.pi, len(idx)))
                    spawn(ev[idx], np.full(len(idx), qq), P[idx].copy(), dd, tt)
                dead |= pair
            low = (~dead) & (E < self.gcut)
            if low.any():
                self._dep(ev[low], reg[low], E[low])
                dead |= low
            keep = ~dead
            ev, P, D, E = ev[keep], P[keep], D[keep], E[keep]
        if q['ev']:
            return [{k: np.concatenate(v) for k, v in q.items()}]
        return []

    def _kn(self, E):
        rng = self.rng
        k = E / ME
        eps = np.zeros(len(E))
        cth = np.zeros(len(E))
        todo = np.arange(len(E))
        for _ in range(500):
            if not len(todo):
                break
            kk = k[todo]
            e0 = 1 / (1 + 2 * kk)
            a1 = -np.log(e0)
            a2 = 0.5 * (1 - e0 * e0)
            u1, u2, u3 = rng.uniform(size=(3, len(todo)))
            e = np.where(u1 < a1 / (a1 + a2), np.exp(-a1 * u2), np.sqrt(e0 * e0 + (1 - e0 * e0) * u2))
            t = (1 - e) / (kk * e)
            s2 = t * (2 - t)
            gfun = 1 - e * s2 / (1 + e * e)
            acc = u3 <= gfun
            eps[todo[acc]] = e[acc]
            cth[todo[acc]] = 1 - t[acc]
            todo = todo[~acc]
        return eps, np.clip(cth, -1, 1)
