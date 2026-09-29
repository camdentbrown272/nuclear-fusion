#!/usr/bin/env python3
"""M5 neutron detection for a tabletop cell: moderated 3He bank (MC), EJ-309 (PSD), bubble detectors.

Neutron transport MC (numpy, vectorised):
  * Epithermal/fast: elastic scattering, isotropic in the CM, on free H (Gammel n-p cross-section)
    and C / O (coarse ENDF/B-VIII digitisations, +/-15 %); inelastic channels neglected.
  * Thermal (E < 0.2 eV): one-group model. Scattering isotropic with an effective transport cross-section
    per H calibrated to Lamarsh's water thermal diffusion coefficient D = 0.16 cm; absorption 1/v with
    Maxwellian averaging (x 0.886) at 293 K.
  * 3He tubes: 1/v absorption, sigma(0.0253 eV) = 5333 b; tube walls neglected.
Validation: Fermi age of 2 MeV neutrons in water to 1.46 eV (Lamarsh: 27 cm2 for fission neutrons)
and thermal diffusion length in water (2.85 cm) and polyethylene (literature ~2.1-2.3 cm).

    python3 sim/m5_neutron.py      # -> figs/m5_neutron_*.png, figs/m5_neutron.txt
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m5_stopping as st  # noqa: E402

FIGS = st.FIGS
RNG = np.random.default_rng(1234)
NA = 6.02214076e23
E_TH = 0.2e-6  # MeV: below this, neutron joins the thermal group
MAXW = 0.886  # <sigma_a> for 1/v absorber over a Maxwellian flux at 293 K relative to sigma(0.0253 eV)

# --- cross sections (barns), energy in MeV
_EC = np.array([1e-9, 1e-2, 0.1, 0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 10.0, 14.0, 20.0, 30.0])
_SC = np.array([4.74, 4.72, 4.5, 3.9, 3.4, 2.6, 2.2, 1.72, 1.6, 1.9, 1.5, 1.2, 1.3, 1.1, 0.85, 0.7, 0.6])
_SO = np.array([3.8, 3.8, 3.7, 3.5, 3.3, 5.5, 1.8, 1.3, 1.5, 1.2, 1.3, 1.0, 1.2, 1.1, 0.9, 0.7, 0.6])


def sig_H(E):
    E = np.maximum(E, 1e-9)
    return (3 * np.pi / (1.206 * E + (-1.86 + 0.09415 * E + 0.0001306 * E ** 2) ** 2)
            + np.pi / (1.206 * E + (0.4223 + 0.13 * E) ** 2))


def sig_C(E):
    return np.interp(np.log(np.maximum(E, 1e-9)), np.log(_EC), _SC)


def sig_O(E):
    return np.interp(np.log(np.maximum(E, 1e-9)), np.log(_EC), _SO)


# thermal-group constants
SIG_TR_H_TH = 28.6   # b, effective isotropic transport xs per bound H, calibrated to water D=0.16 cm (see validate)
SIG_A_H = 0.332 * MAXW
SIG_A_C = 0.0035 * MAXW
SIG_A_O = 0.00019 * MAXW
SIG_A_HE3 = 5333.0 * MAXW


class Medium:
    def __init__(self, nH, nC=0.0, nO=0.0):
        self.nH, self.nC, self.nO = nH, nC, nO

    def sigma_t(self, E):
        """macroscopic total (1/cm) per neutron energy array (fast/epithermal and thermal)."""
        th = E < E_TH
        fast = (self.nH * sig_H(E) + self.nC * sig_C(E) + self.nO * sig_O(E)) * 1e-24
        therm = (self.nH * (SIG_TR_H_TH + SIG_A_H) + self.nC * (4.7 + SIG_A_C) + self.nO * (3.8 + SIG_A_O)) * 1e-24
        return np.where(th, therm, fast)

    def collide(self, E, d):
        """Process a collision: returns new E, new direction, absorbed flag."""
        n = len(E)
        th = E < E_TH
        absorbed = np.zeros(n, bool)
        # thermal: absorption probability
        sa = (self.nH * SIG_A_H + self.nC * SIG_A_C + self.nO * SIG_A_O)
        st_ = (self.nH * (SIG_TR_H_TH + SIG_A_H) + self.nC * (4.7 + SIG_A_C) + self.nO * (3.8 + SIG_A_O))
        absorbed[th] = RNG.uniform(0, 1, th.sum()) < sa / st_
        # fast: choose nucleus
        sH = self.nH * sig_H(E)
        sC = self.nC * sig_C(E)
        sO = self.nO * sig_O(E)
        u = RNG.uniform(0, 1, n) * (sH + sC + sO)
        A = np.where(u < sH, 1.0, np.where(u < sH + sC, 12.0, 16.0))
        # isotropic CM scattering -> lab energy and direction
        muc = RNG.uniform(-1, 1, n)
        Enew = E * (A ** 2 + 2 * A * muc + 1) / (A + 1) ** 2
        mulab = (1 + A * muc) / np.sqrt(A ** 2 + 2 * A * muc + 1)
        dnew = rotate(d, mulab)
        # thermal-group scattering: isotropic, energy stays thermal
        iso = random_dirs(n)
        dnew = np.where(th[:, None], iso, dnew)
        Enew = np.where(th, 0.0253e-6, np.where(Enew < E_TH, 0.0253e-6, Enew))
        return Enew, dnew, absorbed


def random_dirs(n):
    mu = RNG.uniform(-1, 1, n)
    ph = RNG.uniform(0, 2 * np.pi, n)
    s = np.sqrt(1 - mu ** 2)
    return np.stack([s * np.cos(ph), s * np.sin(ph), mu], axis=1)


def rotate(d, mu):
    """Rotate unit vectors d by polar angle arccos(mu) with random azimuth."""
    n = len(mu)
    ph = RNG.uniform(0, 2 * np.pi, n)
    s = np.sqrt(np.maximum(1 - mu ** 2, 0))
    # orthonormal basis
    a = np.where(np.abs(d[:, 2:3]) < 0.9, np.array([[0, 0, 1.0]]), np.array([[1.0, 0, 0]]))
    u = np.cross(d, a)
    u /= np.linalg.norm(u, axis=1)[:, None]
    v = np.cross(d, u)
    out = mu[:, None] * d + s[:, None] * (np.cos(ph)[:, None] * u + np.sin(ph)[:, None] * v)
    return out / np.linalg.norm(out, axis=1)[:, None]


def hdpe(rho=0.95):
    n = rho / 14.027 * NA
    return Medium(2 * n, n)


def water():
    n = 0.998 / 18.015 * NA
    return Medium(2 * n, 0.0, n)


# ------------------------------------------------------------------------------ validation (infinite medium)
def infinite_medium(med, E0, n=40000, stop_E=None):
    """Track in infinite medium; return r^2 at slowing-down past stop_E (MeV) or at thermal absorption."""
    x = np.zeros((n, 3))
    d = random_dirs(n)
    E = np.full(n, E0)
    alive = np.ones(n, bool)
    r2 = np.zeros(n)
    for _ in range(20000):
        idx = np.nonzero(alive)[0]
        if len(idx) == 0:
            break
        s = -np.log(RNG.uniform(0, 1, len(idx))) / med.sigma_t(E[idx])
        x[idx] += d[idx] * s[:, None]
        En, dn, ab = med.collide(E[idx], d[idx])
        if stop_E is not None:
            done = En < stop_E
        else:
            done = ab
        r2[idx[done]] = np.sum(x[idx[done]] ** 2, axis=1)
        alive[idx[done]] = False
        E[idx], d[idx] = En, dn
    return r2


def thermal_diffusion_length(med, n=40000):
    """Thermal neutrons born at the origin; L^2 = <r^2>/6 at absorption."""
    x = np.zeros((n, 3))
    d = random_dirs(n)
    E = np.full(n, 0.0253e-6)
    alive = np.ones(n, bool)
    r2 = np.zeros(n)
    for _ in range(100000):
        idx = np.nonzero(alive)[0]
        if len(idx) == 0:
            break
        s = -np.log(RNG.uniform(0, 1, len(idx))) / med.sigma_t(E[idx])
        x[idx] += d[idx] * s[:, None]
        En, dn, ab = med.collide(E[idx], d[idx])
        r2[idx[ab]] = np.sum(x[idx[ab]] ** 2, axis=1)
        alive[idx[ab]] = False
        d[idx] = dn
    return np.sqrt(r2.mean() / 6)


# ------------------------------------------------------------------------------ 3He bank MC
def he3_bank(R_c=8.0, t_front=4.0, t_back=6.0, n_tubes=None, tube_r=1.27, p_atm=4.0, H=50.0, L_active=40.0,
             n=20000, E0=2.45, source="center", ext_spectrum=None, pitch=4.0, t_shield=0.0):
    """Cylindrical HDPE moderator (inner radius R_c = cavity for the cell, thickness t_front up to the tube axes,
    t_back behind), N tubes parallel to z on a ring. Returns detection efficiency (fraction of source neutrons
    absorbed in 3He). source='center' (isotropic point) or 'external' (isotropic inward flux on outer surface)."""
    med = hdpe()
    r_ring = R_c + t_front
    R_m = r_ring + t_back           # moderator outer radius (a Cd/boron layer here absorbs thermal neutrons)
    R_o = R_m + t_shield            # outer radius incl. borated-HDPE shield (thermals absorbed on contact)
    if n_tubes is None:
        n_tubes = max(int(2 * np.pi * r_ring / pitch), 4)
    ang = 2 * np.pi * np.arange(n_tubes) / n_tubes
    tx, ty = r_ring * np.cos(ang), r_ring * np.sin(ang)
    nHe = p_atm * 101325 / (1.380649e-23 * 293) * 1e-6  # /cm3
    if source == "center":
        x = np.zeros((n, 3))
        d = random_dirs(n)
        E = np.full(n, E0)
    else:
        # isotropic flux entering the outer cylindrical surface (cosine-weighted inward)
        phi = RNG.uniform(0, 2 * np.pi, n)
        z = RNG.uniform(-H / 2, H / 2, n)
        x = np.stack([R_o * np.cos(phi) * 0.99999, R_o * np.sin(phi) * 0.99999, z], axis=1)
        nrm = -np.stack([np.cos(phi), np.sin(phi), np.zeros(n)], axis=1)
        mu = np.sqrt(RNG.uniform(0, 1, n))
        d = rotate(nrm, mu)
        E = ext_spectrum(n)
    alive = np.ones(n, bool)
    det = np.zeros(n, bool)
    for _ in range(6000):
        idx = np.nonzero(alive)[0]
        if len(idx) == 0:
            break
        xi, di, Ei = x[idx], d[idx], E[idx]
        r = np.hypot(xi[:, 0], xi[:, 1])
        in_cav = r < R_c - 1e-9
        # distance to cylinder boundaries along direction
        a = di[:, 0] ** 2 + di[:, 1] ** 2
        b = xi[:, 0] * di[:, 0] + xi[:, 1] * di[:, 1]

        def dist_cyl(R):
            c = r ** 2 - R ** 2
            disc = b ** 2 - a * c
            with np.errstate(invalid="ignore", divide="ignore"):
                s1 = (-b - np.sqrt(np.maximum(disc, 0))) / a
                s2 = (-b + np.sqrt(np.maximum(disc, 0))) / a
            s = np.where(s1 > 1e-7, s1, np.where(s2 > 1e-7, s2, np.inf))
            return np.where((disc > 0) & (a > 1e-12), s, np.inf)

        s_c = dist_cyl(R_c)
        s_o = dist_cyl(R_o)
        if t_shield > 0:
            s_o = np.minimum(s_o, dist_cyl(R_m))
        with np.errstate(divide="ignore"):
            s_z = np.where(di[:, 2] > 0, (H / 2 - xi[:, 2]) / di[:, 2], np.where(di[:, 2] < 0, (-H / 2 - xi[:, 2]) / di[:, 2], np.inf))
        s_bound = np.minimum(np.minimum(s_c, s_o), s_z)
        s_col = np.where(in_cav, np.inf, -np.log(RNG.uniform(0, 1, len(idx))) / med.sigma_t(Ei))
        # tube intersections (only in HDPE region)
        s_tube = np.full(len(idx), np.inf)
        chord = np.zeros(len(idx))
        for k in range(n_tubes):
            px, py = xi[:, 0] - tx[k], xi[:, 1] - ty[k]
            bb = px * di[:, 0] + py * di[:, 1]
            cc = px ** 2 + py ** 2 - tube_r ** 2
            disc = bb ** 2 - a * cc
            ok = (disc > 0) & (a > 1e-12)
            with np.errstate(invalid="ignore", divide="ignore"):
                sq = np.sqrt(np.maximum(disc, 0))
                s1 = (-bb - sq) / a
                s2 = (-bb + sq) / a
            hit = ok & (s1 > 1e-7)
            better = hit & (s1 < s_tube)
            s_tube = np.where(better, s1, s_tube)
            chord = np.where(better, s2 - s1, chord)
        # z-extent of tube active length
        ztube = xi[:, 2] + di[:, 2] * s_tube
        tube_ok = np.abs(ztube) < L_active / 2
        s_tube = np.where(tube_ok & ~in_cav, s_tube, np.inf)
        s_min = np.minimum(np.minimum(s_col, s_bound), s_tube)
        newx = xi + di * s_min[:, None]
        ev_tube = s_tube <= np.minimum(s_col, s_bound)
        ev_col = (s_col < s_bound) & ~ev_tube
        ev_bnd = ~ev_tube & ~ev_col
        # tube crossing
        Et = Ei[ev_tube]
        sigHe = SIG_A_HE3 * np.sqrt(0.0253e-6 / np.maximum(np.where(Et < E_TH, 0.0253e-6, Et), 1e-12)) / MAXW
        sigHe = np.where(Et < E_TH, SIG_A_HE3, sigHe)
        pabs = 1 - np.exp(-nHe * sigHe * 1e-24 * chord[ev_tube])
        absd = RNG.uniform(0, 1, ev_tube.sum()) < pabs
        ti = idx[ev_tube]
        det[ti[absd]] = True
        alive[ti[absd]] = False
        # survivors: exit tube (move by chord)
        keep = ti[~absd]
        x[keep] = newx[ev_tube][~absd] + di[ev_tube][~absd] * (chord[ev_tube][~absd] + 1e-6)[:, None]
        # collisions
        ci = idx[ev_col]
        x[ci] = newx[ev_col]
        En, dn, ab = med.collide(Ei[ev_col], di[ev_col])
        E[ci], d[ci] = En, dn
        alive[ci[ab]] = False
        if t_shield > 0:  # thermal neutrons in the borated shield are absorbed
            rr = np.hypot(x[ci, 0], x[ci, 1])
            alive[ci[(rr > R_m) & (En < E_TH)]] = False
        # boundaries
        bi = idx[ev_bnd]
        xb = newx[ev_bnd]
        rb = np.hypot(xb[:, 0], xb[:, 1])
        out = (rb > R_o - 1e-6) | (np.abs(xb[:, 2]) > H / 2 - 1e-6)
        if t_shield > 0:  # Cd/boron sheet at R_m kills thermal neutrons crossing it
            out |= (np.abs(rb - R_m) < 1e-5) & (Ei[ev_bnd] < E_TH)
        alive[bi[out]] = False
        x[bi] = xb + di[ev_bnd] * 1e-6
    return det.mean()


# ------------------------------------------------------------------------------ cosmic spectrum (external)
def cosmic_spectrum_sampler(which="fast"):
    """Sea-level cosmic neutron flux model (NYC, Gordon et al. IEEE TNS 51 (2004) 3427; JESD89A).
    Components (n cm^-2 h^-1): thermal <0.4 eV ~4, epithermal 0.4 eV-1 MeV ~10, 1-10 MeV ~9, >10 MeV 13.
    Only fast/epithermal reach the 3He behind a Cd/borated liner. Assumption uncertainty +/-50 %."""
    comps = {"epi": (10.0, 0.4e-6, 1.0), "fast": (9.0, 1.0, 10.0), "hi": (13.0, 10.0, 30.0)}

    def sample(n):
        w = np.array([comps[k][0] for k in comps])
        k = RNG.choice(len(w), size=n, p=w / w.sum())
        lo = np.array([comps[c][1] for c in comps])[k]
        hi = np.array([comps[c][2] for c in comps])[k]
        return np.exp(RNG.uniform(np.log(lo), np.log(hi)))

    return sample, sum(c[0] for c in comps.values()) / 3600.0  # flux cm-2 s-1


# ------------------------------------------------------------------------------ EJ-309
def ej309(n=200000, E0=2.45, dist=10.0, radius=6.35, length=12.7, thr_MeVee=0.10, E_arr=None, iso=False):
    """Single-scatter + multiple-scatter approx for a cylindrical EJ-309 cell facing a point source.
    EJ-309: rho 0.959 g/cm3, H 5.43e22, C 4.35e22 /cm3 (Eljen data sheet,
    https://eljentechnology.com/products/liquid-scintillators/ej-301-ej-309). Proton light (MeVee):
    L = 0.817 Ep - 2.63 (1 - exp(-0.297 Ep^0.9)) (Enqvist et al. NIM A 715 (2013) 79, EJ-309 fit)."""
    nH, nC = 5.43e22, 4.35e22
    geo = 0.5 * (1 - dist / np.sqrt(dist ** 2 + radius ** 2))
    # track through cylinder along ~axis (source on axis); chord ~ length for small angles
    mu = RNG.uniform(np.sqrt(1 - (radius / np.hypot(dist, radius)) ** 2), 1, n)
    L = length / mu
    if iso:  # isotropic external field: mean chord of a convex body = 4V/S
        L = RNG.exponential(4 * np.pi * radius ** 2 * length / (2 * np.pi * radius ** 2 + 2 * np.pi * radius * length), n)
    light = np.zeros(n)
    E = np.full(n, E0) if E_arr is None else E_arr.copy()
    depth = np.zeros(n)
    alive = np.ones(n, bool)
    for _ in range(6):  # up to 6 scatters (continuing roughly forward)
        sH = nH * sig_H(E) * 1e-24
        sC = nC * sig_C(E) * 1e-24
        s = -np.log(RNG.uniform(0, 1, n)) / (sH + sC)
        depth = depth + s
        inside = alive & (depth < L)
        onH = RNG.uniform(0, 1, n) < sH / (sH + sC)
        frac = RNG.uniform(0, 1, n)
        Ep = np.where(onH, E * frac, 0.0)
        Lp = np.maximum(0.817 * Ep - 2.63 * (1 - np.exp(-0.297 * Ep ** 0.9)), 0)
        light += np.where(inside, Lp, 0)
        E = np.where(inside & onH, E * (1 - frac), np.where(inside, E * (1 - 0.28 * RNG.uniform(0, 1, n)), E))
        alive = inside & (E > 0.1)
    intrinsic = np.mean(light > thr_MeVee)
    return geo, intrinsic, light


def ej309_bkg(radius=6.35, length=12.7, window=(0.10, 0.80), gamma_cps=150.0, psd_misid=1e-3):
    """Cosmic fast-neutron background in the 2.45 MeV-neutron light window (0.1-0.8 MeVee; the 2.45 MeV
    proton-recoil edge is at ~0.73 MeVee) of one 5"x5" EJ-309 cell, plus gamma leakage through PSD.
    gamma singles ~150 cps at sea level unshielded (assumption), PSD misID 1e-3 at 100 keVee (typical for EJ-309)."""
    samp = lambda n: np.exp(RNG.uniform(np.log(1.0), np.log(30.0), n))
    flux = 22.0 / 3600  # 1-30 MeV (Gordon 2004 model above)
    S = 2 * np.pi * radius ** 2 + 2 * np.pi * radius * length
    En = samp(200000)
    _, _, light = ej309(n=200000, E_arr=En, iso=True, radius=radius, length=length)
    frac = np.mean((light > window[0]) & (light < window[1]))
    n_rate = flux * S / 4 * frac
    return dict(neutron=n_rate, gamma=gamma_cps * psd_misid, total=n_rate + gamma_cps * psd_misid, frac=frac)


# ------------------------------------------------------------------------------ channel summary for m5_stats
_BANK = dict(R_c=8.0, t_front=4.0, t_back=6.0)
T_SHIELD = 20.0  # cm borated HDPE outside a 1 mm Cd sheet


def channel_summary():
    """Baseline numbers used by m5_stats (recomputed quickly)."""
    eff = he3_bank(**_BANK, t_shield=T_SHIELD, n=8000)
    samp, flux = cosmic_spectrum_sampler()
    ext = he3_bank(**_BANK, n=8000, source="external", ext_spectrum=samp, t_shield=T_SHIELD)
    R_o = _BANK["R_c"] + _BANK["t_front"] + _BANK["t_back"] + T_SHIELD
    area = 2 * np.pi * R_o * 50.0  # outer lateral area (cm2)
    # isotropic field of fluence rate phi: inward current phi/4 per cm2; x1.5 for top/bottom (not modelled)
    bkg_cosmic = flux / 4 * area * ext * 1.5
    # muon-induced neutrons in the shield and nearby building materials: ~10 % of the hadronic
    # component at sea level for low-Z shielding (assumption); a muon veto removes most of them
    bkg_mu = 0.1 * bkg_cosmic
    n_tubes = int(2 * np.pi * (_BANK["R_c"] + _BANK["t_front"]) / 4.0)
    bkg_int = n_tubes * 1.0 / 3600  # 1 count/h/tube intrinsic alpha + electronic after PSA (assumption)
    geo, intr, _ = ej309()
    ej_eff = 2 * geo * intr
    ej_bkg = 2 * ej309_bkg()["total"]
    return dict(eff=eff, bkg=bkg_cosmic + bkg_mu + bkg_int, bkg_cosmic=bkg_cosmic, bkg_mu=bkg_mu, bkg_int=bkg_int,
                ext_resp=ext, ej_eff=ej_eff, ej_bkg=ej_bkg, n_tubes=n_tubes)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = open(os.path.join(FIGS, "m5_neutron.txt"), "w")

    def P(*a):
        print(*a)
        print(*a, file=out)

    P("== Validation (infinite media)")
    w = water()
    Lw = thermal_diffusion_length(w, 20000)
    Lp = thermal_diffusion_length(hdpe(), 20000)
    P(f"   thermal diffusion length: water {Lw:.2f} cm (Lamarsh 2.85 cm; calibration input D=0.16 cm), "
      f"HDPE {Lp:.2f} cm (literature ~2.1-2.3 cm; prediction)")
    for E0, lab in [(2.0, "2.0 MeV"), (2.45, "2.45 MeV")]:
        r2 = infinite_medium(w, E0, 20000, stop_E=1.46e-6)
        P(f"   Fermi age {lab} -> 1.46 eV in water: {r2.mean() / 6:.1f} cm2 (Lamarsh fission-spectrum value 27 cm2)")
    r2 = infinite_medium(hdpe(), 2.45, 20000, stop_E=1.46e-6)
    P(f"   Fermi age 2.45 MeV -> 1.46 eV in HDPE: {r2.mean() / 6:.1f} cm2")

    P("\n== 3He bank efficiency vs moderator thickness in front of tubes (2.45 MeV point source at cavity centre)")
    P("   cavity radius 8 cm, tubes 1\" x 40 cm active, 4 atm, pitch 4 cm on ring, 6 cm HDPE behind")
    tf = [1, 2, 3, 4, 5, 6, 8, 10]
    effs = []
    for t in tf:
        e = he3_bank(R_c=8.0, t_front=t, t_back=6.0, n=6000)
        effs.append(e)
        P(f"   t_front={t:4.1f} cm: eff = {e:.3f}")
    tb = [2, 4, 6, 10]
    effb = [he3_bank(R_c=8.0, t_front=4.0, t_back=t, n=6000) for t in tb]
    P("   back-moderator scan (t_front=4): " + ", ".join(f"{t} cm: {e:.3f}" for t, e in zip(tb, effb)))
    pr = [(4.0, 4.0), (10.0, 4.0), (4.0, 3.0)]
    for p_atm, pitch in pr:
        e = he3_bank(R_c=8.0, t_front=4.0, t_back=6.0, p_atm=p_atm, pitch=pitch, n=6000)
        P(f"   p={p_atm} atm, pitch {pitch} cm: eff = {e:.3f}")
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(tf, effs, "o-")
    ax.set_xlabel("HDPE thickness between cavity and tube axes (cm)")
    ax.set_ylabel("detection efficiency (2.45 MeV, 4π)")
    ax.set_title("Moderated ³He ring, cavity r = 8 cm")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "m5_neutron_he3.png"), dpi=130)
    plt.close(fig)

    P("   cavity-radius scan (t_front=4, t_back=6, 20 cm shield): " + ", ".join(
        f"R_c={rc} cm: {he3_bank(R_c=rc, t_front=4, t_back=6, t_shield=20, n=6000):.3f}" for rc in [5, 8, 12, 16, 20]))
    samp, flux = cosmic_spectrum_sampler()
    P("   outer borated-HDPE shield scan (response per neutron entering the OUTER surface; rate scales with its area):")
    for tsh in [0, 10, 20, 30]:
        R_o = 8 + 4 + 6 + tsh
        ext = he3_bank(R_c=8.0, t_front=4.0, t_back=6.0, t_shield=tsh, source="external", ext_spectrum=samp, n=8000)
        rate = flux / 4 * (2 * np.pi * R_o * 50 * 1.5) * ext
        P(f"     shield {tsh:3d} cm: response {ext:.4f} -> cosmic bkg {rate:.3f} cps")
    cs = channel_summary()
    P("\n== Baseline bank (t_front=4, t_back=6, cavity r=8 cm, 1 mm Cd + 20 cm borated HDPE shield)")
    P(f"   tubes: {cs['n_tubes']}; eff(2.45 MeV, centre) = {cs['eff']:.3f}")
    P(f"   response to external cosmic spectrum (per entering neutron) = {cs['ext_resp']:.3f}")
    P(f"   background: cosmic {cs['bkg_cosmic']:.3f} + muon-induced {cs['bkg_mu']:.3f} + intrinsic {cs['bkg_int']:.3f}"
      f" = {cs['bkg']:.3f} cps (unshielded bank: ~1.8 cps)")
    eb = ej309_bkg()
    P(f"   EJ-309 5x5 cell background in 0.1-0.8 MeVee n-band: cosmic n {eb['neutron']:.3f} + gamma leak {eb['gamma']:.3f} cps")
    geo, intr, light = ej309()
    P(f"\n== EJ-309 5\"x5\" at 10 cm: geometric {geo:.3f}, intrinsic (>0.1 MeVee) {intr:.3f}, per cell {geo * intr:.4f};"
      f" 2 cells {2 * geo * intr:.4f}")
    # bubble detectors
    sens = 3.0          # bubbles / uSv (BD-PND, 1.9-3.7 b/uSv)
    h10 = 416e-12 * 1e6  # uSv cm2 at 2.5 MeV (ICRP 74 H*(10) conversion)
    for r in [3.0, 5.0]:
        per_n = sens * h10 / (4 * np.pi * r ** 2)
        P(f"   BD-PND at {r} cm: {per_n:.2e} bubbles per emitted neutron")
    bkg_b = sens * 0.009 * 24  # sea-level neutron H*(10) ~9 nSv/h (UNSCEAR 2000)
    P(f"   BD-PND background ~{bkg_b:.2f} bubbles/day (cosmic neutron dose ~9 nSv/h)")

    P("\n== Cosmic background systematics for the 3He bank")
    B = cs["bkg_cosmic"]
    for days in [1, 14, 30]:
        N = B * days * 86400
        P(f"   {days:3d} d: N_cosmic = {N:.3g}; statistical 1 sigma = {100 / np.sqrt(N):.3f} %;"
          f" a 1 hPa pressure change shifts it by 0.72 % = {0.0072 * np.sqrt(N):.1f} sigma")
    out.close()
    return cs


if __name__ == "__main__":
    main()
