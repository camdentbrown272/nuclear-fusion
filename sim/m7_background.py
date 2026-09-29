"""M7 / Q2-Q3: background model and minimum detectable pair-production rate
for the 511-511 coincidence channel and the 12-30 MeV calorimetric channel.

Method
------
1. Cosmic-muon knock-on (delta-ray) showers in the shield: MC of delta rays
   (d2N/dT dx = 0.1535 (Z/A) F/T^2 MeV^-1 (g/cm2)^-1, PDG eq. 34.8, T = 1.2 MeV
   .. 1.1 GeV for a 4 GeV muon) in a laterally infinite slab (5 cm HDPE inner
   liner + 10 cm Pb), tallying photons/positrons that cross into the cavity.
   Cavity fluence rate phi(E) = 4 J_in(E) (isotropic enclosure).
2. Detector response to an isotropic photon field inside the cavity
   (sphere source, cosine-law inward emission; fluence = N / (pi R^2)):
   R_c(E) = 511-511 coincidences per unit fluence, R_hi(E) = sum in window.
3. Terrestrial gamma lines leaking through the shield (concrete room model),
   radon, cosmogenic beta+ activation, Michel positrons from stopped mu+,
   direct muons through the crystals (chord MC), hadronic component
   (ratio to muonic, from literature), accidentals.
4. MDA: 30 d total, on/off modulation (15 d on / 15 d off), Asimov 5 sigma
   (Li & Ma eq. 17 with alpha = 1; Cowan et al. EPJC 71, 1554 (2011)).
Outputs: docs/models/figs/m7_background.txt, m7_background.png, m7_mda.png
"""
import os
import sys
import time
import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(__file__))
from m7_mc import MC, Geometry, Box  # noqa: E402
from m7_common import MATERIALS, ME, PDG_DEDX_MIN_BK  # noqa: E402
from m7_detectors import smear, win511, coinc_mask, TAU2  # noqa: E402
from m7_transport import ipc_source  # noqa: E402
import m7_geom  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), '..', 'docs', 'models', 'figs')

# ------------------------------------------------------------------ inputs
# Cosmic muons at sea level (PDG Cosmic-ray review 2024, sec. 30.3 [BK]):
J1 = 1.0 / 60.0                        # cm^-2 s^-1 through a horizontal surface
I_V = J1 / (np.pi / 2)                 # vertical intensity with cos^2 law
PHI_MU = 2 * np.pi / 3 * I_V           # omnidirectional fluence rate (track length / cm3 / s)
T_DMIN, T_DMAX = 1.2, 1100.0           # delta-ray energies able to make pairs; T_max for 4 GeV mu
STOP_RATE = 0.012                      # stopped muons / kg / s at sea level [est., range 0.005-0.02]
MU_PLUS = 0.56                         # mu+/(mu+ + mu-) (charge ratio 1.27, PDG [BK])
MU_TAU = 2.197e-6
# Hadronic (neutron/proton) contribution to cosmogenic positron production at the
# surface relative to muonic. Literature [BK: Heusser, Annu. Rev. Nucl. Part. Sci.
# 45, 543 (1995); Semkow et al., NIM A 489, 519 (2002)]: an active muon veto
# reduces shielded-surface Ge backgrounds (incl. the 511 keV line) by only
# ~2-5x, i.e. the unvetoable (mostly hadronic) part is 20-50 % of the total.
F_HAD = {'central': 0.5, 'low': 0.25, 'high': 1.0}
VETO_INEFF = {'central': 3e-3, 'low': 1e-3, 'high': 1e-2}   # prompt, hermetic 5 cm plastic
VETO_WINDOW = 20e-6                    # s, extended veto after every veto hit (Michel suppression)
VETO_RATE_PER_M2 = 170.0               # muons/m2/s through a ~horizontal+vertical veto surface [~1 cm^-2 min^-1]
# Sites (muon-flux factor, hadron factor) [BK: Heusser 1995 fig. 3; Mei & Hime PRD 73, 053004 (2006)]
SITES = {'sea level': (1.0, 1.0), 'shallow 30 m w.e.': (0.15, 1e-3), 'deep >1 km w.e.': (1e-5, 0.0)}
# Terrestrial gamma field from a concrete room (UNSCEAR 2000 Annex B typical
# building concrete: 40K 400, 226Ra 40, 232Th 30 Bq/kg [BK]); uncollided fluence
# at the centre of a void in an infinite source medium = S_v / mu; the room
# is taken as 0.5 of that (floor + walls + ceiling with openings).
CONC_RHO = 2.35
CONC = {'K': 400.0, 'U': 40.0, 'Th': 30.0}
LINES = [(2.6145, 'Th', 0.359), (1.7645, 'U', 0.153), (2.2042, 'U', 0.049), (1.4608, 'K', 0.1066)]
ROOM_FACTOR = 0.5
PB_T, HDPE_T = 10.0, 5.0              # shield: 10 cm Pb outside 5 cm (borated) HDPE
BUILDUP_PAIR = 1.3                    # fraction of scattered photons still >1.02 MeV ~ small-angle buildup

T_ON = 15 * 86400.0                   # 30 d, 50 % duty
DAY = 86400.0


def mu_concrete(E):
    # ordinary concrete ~ (O 0.53, Si 0.34, Ca 0.044, Al 0.034, Na 0.029, Fe 0.014, K 0.013, H 0.01)
    from m7_common import Material
    m = Material('CONC', CONC_RHO, {'O': 0.53, 'Si': 0.34, 'Al': 0.078, 'Na': 0.029, 'Fe': 0.014, 'H': 0.01}, 'mass')
    return float(m.mu_total(E)[0]) * CONC_RHO


def shield_transmission(E, pb=PB_T, pe=HDPE_T):
    return float(np.exp(-MATERIALS['PB'].mu_total(E)[0] * MATERIALS['PB'].rho * pb
                        - MATERIALS['HDPE'].mu_total(E)[0] * MATERIALS['HDPE'].rho * pe))


# ------------------------------------------------------------------ 1. shield showers
def slab_showers(n=3000, seed=5, pb=PB_T, pe=HDPE_T):
    rng = np.random.default_rng(seed)
    H = (pb + pe) / 2
    g = Geometry([200, 200, H + 1e-3], 'VAC')
    g.add(Box([0, 0, -H + pe / 2], [200, 200, pe / 2]), 'HDPE', 'hdpe')
    g.add(Box([0, 0, -H + pe + pb / 2], [200, 200, pb / 2]), 'PB', 'pb')
    mc = MC(g, [], rng, kc=0.1, tcut=0.1, gcut=0.1, Tmax=1200.0)
    lay = [('HDPE', pe, -H, -H + pe), ('PB', pb, -H + pe, H)]
    areal = np.array([MATERIALS[m].ZA * MATERIALS[m].rho * t for m, t, _, _ in lay])
    li = rng.choice(2, size=n, p=areal / areal.sum())
    z = np.array([rng.uniform(lay[i][2], lay[i][3]) for i in li])
    Pp = np.stack([rng.uniform(-50, 50, n), rng.uniform(-50, 50, n), z], 1)
    T = T_DMIN * (T_DMAX / T_DMIN) ** rng.uniform(size=n)
    w = np.log(T_DMAX / T_DMIN) / (T * (1 / T_DMIN - 1 / T_DMAX))
    c = rng.uniform(-1, 1, n)
    f = rng.uniform(0, 2 * np.pi, n)
    D = np.stack([np.sqrt(1 - c * c) * np.cos(f), np.sqrt(1 - c * c) * np.sin(f), c], 1)
    mc.run(n, charged={'ev': np.arange(n), 'q': -np.ones(n, int), 'P': Pp, 'D': D, 'T': T})
    eg = np.array(mc.esc_g).reshape(-1, 3)
    inward = eg[:, 2] < 0
    ev = eg[inward, 0].astype(int)
    E = eg[inward, 1]
    wt = w[ev]
    R_delta = PHI_MU * 0.1535 * areal.sum() * (1 / T_DMIN - 1 / T_DMAX)     # delta rays / cm2 / s
    eq = np.array(mc.esc_q).reshape(-1, 4)
    pos_in = (eq[:, 1] > 0) & (eq[:, 3] < 0) if len(eq) else np.zeros(0, bool)
    J_pos = R_delta * np.sum(w[eq[pos_in, 0].astype(int)]) / n if len(eq) else 0.0
    return E, wt * R_delta / n, R_delta, J_pos     # photon energies and current weights [cm^-2 s^-1]


# ------------------------------------------------------------------ 2. responses
def sphere_source(n, R, E, rng):
    c = rng.uniform(-1, 1, n)
    f = rng.uniform(0, 2 * np.pi, n)
    s = np.sqrt(1 - c * c)
    nrm = np.stack([s * np.cos(f), s * np.sin(f), c], 1)
    P = R * nrm
    # cosine-law inward direction about -nrm
    ct = np.sqrt(rng.uniform(size=n))
    from m7_ipc import rotate
    D = rotate(-nrm, ct, rng.uniform(0, 2 * np.pi, n))
    return P, D


def responses(layout, energies, n=30000, seed=9, n_hi=200000):
    dk, kw, pairs = layout
    mat = m7_geom.DETECTORS[dk][0]
    rng = np.random.default_rng(seed)
    Rs = 16.0
    out = {}
    for E0 in energies:
        g, sens = m7_geom.build(L_um=100, det=dk, encl=None, world=30.0, **kw)
        mc = MC(g, sens, rng, kc=0.05, tcut=0.05, gcut=0.02, Tmax=max(40.0, 1.1 * E0))
        nn = n if E0 < 1.05 else (n_hi if E0 <= 30 else n_hi // 4)
        P, D = sphere_source(nn, Rs, E0, rng)
        mc.run(nn, photons={'ev': np.arange(nn), 'P': P, 'D': D, 'E': np.full(nn, E0)})
        Es = smear(mc.edep[:, 2:], mat, rng)
        fl = nn / (np.pi * Rs ** 2)
        lo, hi = win511(mat)
        S = Es.sum(1)
        out[E0] = {'c': coinc_mask(Es, mat, pairs).sum() / fl,
                   's': ((Es > lo) & (Es < hi)).sum() / fl,          # singles, summed over detectors
                   'h5': ((S > 5) & (S < 30)).sum() / fl,
                   'h12': ((S > 12) & (S < 30)).sum() / fl,
                   'nc': coinc_mask(Es, mat, pairs).sum(), 'n': nn}
    return out


def fold(E, w, resp, key):
    Es = np.array(sorted(resp))
    R = np.array([resp[e][key] for e in Es])
    Ri = np.interp(np.log(E), np.log(Es), R)
    Ri[E < Es[0] * 0.97] = 0
    return float(np.sum(w * Ri))


# ------------------------------------------------------------------ 3. other terms
def central_points(g, names, n, rng):
    """Uniform points (by mass) in named regions; returns positions and total mass [g]."""
    g.finalize()
    ids = [g.reg_names.index(nm) for nm in names]
    box = 16.0
    m = 400000
    Pt = rng.uniform(-box, box, size=(m, 3))
    reg = g.locate(Pt)
    sel = np.isin(reg, ids)
    rho = np.array([MATERIALS[g.mat_list[g.reg_mat[r]]].rho if g.mat_list[g.reg_mat[r]] != 'VAC' else 0 for r in reg[sel]])
    vol = (2 * box) ** 3 / m
    mass = rho.sum() * vol
    idx = rng.choice(np.where(sel)[0], size=n, p=rho / rho.sum())
    return Pt[idx] + rng.uniform(-0.5, 0.5, (n, 3)) * 0.0, mass


CENTRAL = ['foil', 'flange', 'chamber', 'topplate', 'electrolyte', 'electrolyte2', 'cellclamp', 'cellwall',
           'cellbottom', 'si_dE', 'si_E', 'si_pcb']


def michel(layout, n=20000, seed=13):
    dk, kw, pairs = layout
    mat = m7_geom.DETECTORS[dk][0]
    rng = np.random.default_rng(seed)
    g, sens = m7_geom.build(L_um=100, det=dk, **kw)
    Pp, mass = central_points(g, CENTRAL, n, rng)
    g, sens = m7_geom.build(L_um=100, det=dk, **kw)
    mc = MC(g, sens, rng, Tmax=60.0)
    x = np.linspace(0, 1, 2001)
    cdf = np.cumsum(x ** 2 * (3 - 2 * x))
    cdf /= cdf[-1]
    T = np.interp(rng.uniform(size=n), cdf, x) * 52.83 - ME
    T = np.clip(T, 0.01, None)
    c = rng.uniform(-1, 1, n)
    f = rng.uniform(0, 2 * np.pi, n)
    D = np.stack([np.sqrt(1 - c * c) * np.cos(f), np.sqrt(1 - c * c) * np.sin(f), c], 1)
    mc.run(n, charged={'ev': np.arange(n), 'q': np.ones(n, int), 'P': Pp, 'D': D, 'T': T})
    Es = smear(mc.edep[:, 2:], mat, rng)
    S = Es.sum(1)
    return mass, coinc_mask(Es, mat, pairs).mean(), ((S > 12) & (S < 30)).mean(), ((S > 5) & (S < 30)).mean()


def central_bplus(layout, n=40000, seed=17):
    """Coincidence acceptance for low-energy beta+ (cosmogenic 58Co/56Co-like,
    <T> ~ 0.3-0.6 MeV) uniformly distributed in the central hardware."""
    dk, kw, pairs = layout
    mat = m7_geom.DETECTORS[dk][0]
    rng = np.random.default_rng(seed)
    g, sens = m7_geom.build(L_um=100, det=dk, **kw)
    Pp, mass = central_points(g, CENTRAL, n, rng)
    g, sens = m7_geom.build(L_um=100, det=dk, **kw)
    mc = MC(g, sens, rng)
    T = rng.uniform(0.05, 1.0, n)
    c = rng.uniform(-1, 1, n)
    f = rng.uniform(0, 2 * np.pi, n)
    D = np.stack([np.sqrt(1 - c * c) * np.cos(f), np.sqrt(1 - c * c) * np.sin(f), c], 1)
    mc.run(n, charged={'ev': np.arange(n), 'q': np.ones(n, int), 'P': Pp, 'D': D, 'T': T})
    Es = smear(mc.edep[:, 2:], mat, rng)
    return mass, coinc_mask(Es, mat, pairs).mean()


def muon_chords(layout, n=400000, seed=19):
    """Direct cosmic muons through the crystals: rate in sum-energy windows."""
    dk, kw, pairs = layout
    mat = m7_geom.DETECTORS[dk][0]
    rng = np.random.default_rng(seed)
    g, sens = m7_geom.build(L_um=100, det=dk, **kw)
    g.finalize()
    dets = [s for s, nm in zip(g.solids, g.names) if nm.startswith('det')]
    A = 60.0                               # generate on z = +40 plane, 120 x 120 cm
    # directions: pdf ∝ cos^3 on the downward hemisphere
    u = rng.uniform(size=n)
    ct = (1 - u) ** 0.25                   # cdf of cos^3 over cos in [0,1]: 1 - c^4
    ph = rng.uniform(0, 2 * np.pi, n)
    st = np.sqrt(1 - ct ** 2)
    D = np.stack([st * np.cos(ph), st * np.sin(ph), -ct], 1)
    P0 = np.stack([rng.uniform(-A, A, n), rng.uniform(-A, A, n), np.full(n, 40.0)], 1)
    rate = J1 * (2 * A) ** 2 / n
    dedx = PDG_DEDX_MIN_BK[mat if mat != 'GE' else 'GE'] * MATERIALS[mat].rho * 1.10
    Edep = np.zeros(n)
    for s in dets:
        t1 = s.ray(P0, D)
        hit = np.isfinite(t1)
        P1 = P0 + D * (np.where(hit, t1, 0) + 1e-6)[:, None]
        t2 = np.where(hit, s.ray(P1, D), 0)
        Edep += t2 * dedx
    Es = Edep + np.random.default_rng(seed + 1).normal(size=n) * 0.03 * Edep
    return {'any': rate * (Edep > 0).sum(), 'h12': rate * ((Es > 12) & (Es < 30)).sum(),
            'h5': rate * ((Es > 5) & (Es < 30)).sum(), 'mean_dep': Edep[Edep > 0].mean()}


def asimov_onoff(s, b):
    """Median significance for on/off counting with equal exposure (Li & Ma 17, alpha=1)."""
    non, noff = s + b, b
    if s <= 0:
        return 0.0
    t1 = non * np.log(2 * non / (non + noff)) if non > 0 else 0
    t2 = noff * np.log(2 * noff / (non + noff)) if noff > 0 else 0
    return np.sqrt(2 * (t1 + t2))


def mda(B, eff, t_on=T_ON, z=5.0):
    """Minimum source rate [pairs/s] for median 5 sigma, on/off with t_on = t_off."""
    b = B * t_on
    s = brentq(lambda s: asimov_onoff(s, b) - z, 1e-6, 1e9)
    # also require >= 5 expected signal counts (for b -> 0 the Asimov formula alone
    # would allow a discovery on ~2 counts; we keep the stricter)
    s = max(s, 5.0)
    return s / (eff * t_on), s, b


def main(layout_name=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from m7_detectors import LAYOUTS
    t0 = time.time()
    lines = []
    P = lines.append
    layout_name = layout_name or 'BGO 3x3 quad (+-x,+-z)'
    lay = LAYOUTS[layout_name]
    dk = lay[0]
    mat = m7_geom.DETECTORS[dk][0]
    P(f"Background model for layout: {layout_name}")
    P(f"muon inputs: J1 = {J1*60:.2f} /cm2/min, I_v = {I_V*1e4:.0f} /m2/s/sr, omni fluence rate {PHI_MU:.4f} /cm2/s")

    # 1. shield showers
    E_g, w_g, R_delta, J_pos = slab_showers()
    P(f"\n[1] delta rays in 5 cm HDPE + 10 cm Pb slab: {R_delta:.3f} /cm2/s (unvetoed)")
    J_tot = w_g.sum()
    J_gt1 = w_g[E_g > 1.022].sum()
    J_511 = w_g[np.abs(E_g - ME) < 0.003].sum()
    P(f"    inward photon current into cavity: all >0.1 MeV {J_tot:.2e}, >1.022 MeV {J_gt1:.2e}, 511 line {J_511:.2e} /cm2/s; positrons {J_pos:.2e} /cm2/s")
    phi_E, phi_w = E_g, 4 * w_g                  # cavity fluence rate contributions
    # 2. responses
    energies = [0.2, 0.3, 0.4, 0.48, 0.511, 0.55, 0.7, 0.9, 1.1, 1.3, 1.46, 1.76, 2.2, 2.61, 3.5, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0]
    resp = responses(lay, energies)
    P("\n[2] response to isotropic cavity photons (per unit fluence, cm2):")
    P("    E[MeV]  R_coinc      R_singles511  R_12-30MeV  (raw coinc counts / photons simulated)")
    for e in energies:
        r = resp[e]
        P(f"    {e:6.3f}  {r['c']:.3e}   {r['s']:.3e}     {r['h12']:.3e}   {r['nc']}/{r['n']}")
    Bmu_c = fold(phi_E, phi_w, resp, 'c')
    Bmu_s = fold(phi_E, phi_w, resp, 's')
    Bmu_h12 = fold(phi_E, phi_w, resp, 'h12')
    Bmu_h5 = fold(phi_E, phi_w, resp, 'h5')
    P(f"\n    muon-shower cavity photons (unvetoed, sea level): coinc {Bmu_c:.2e} /s, singles-511 {Bmu_s:.2e} /s, "
      f"sum 12-30 MeV {Bmu_h12:.2e} /s, 5-30 MeV {Bmu_h5:.2e} /s")
    # Michel positrons from mu+ stopped in the shield: same showers, delayed.
    # Energy ratio (Michel e+ energy / delta-ray energy deposited in shield):
    m_pb = MATERIALS['PB'].rho * 1e3 * 1e-3      # kg per (cm3*1e3) placeholder
    E_delta = PHI_MU * 0.1535 * (MATERIALS['PB'].ZA * MATERIALS['PB'].rho * PB_T + MATERIALS['HDPE'].ZA * MATERIALS['HDPE'].rho * HDPE_T) * np.log(T_DMAX / T_DMIN)  # MeV/cm2/s
    E_mich = STOP_RATE * MU_PLUS * (MATERIALS['PB'].rho * PB_T + MATERIALS['HDPE'].rho * HDPE_T) * 1e-3 * 36.9   # MeV/cm2/s (per cm2 of wall)
    r_mich = E_mich / E_delta
    P(f"    Michel-e+ energy injected in shield / delta-ray energy = {r_mich:.3f} (delayed, tau = 2.2 us)")
    del m_pb

    # 3. terrestrial lines
    P("\n[3] terrestrial gamma lines leaking through 10 cm Pb + 5 cm HDPE (room: concrete K/U/Th = 400/40/30 Bq/kg)")
    Bter = {'c': 0, 's': 0, 'h12': 0}
    for E, ser, yld in LINES:
        S_v = CONC[ser] * CONC_RHO * 1e-3 * yld
        phi_lab = ROOM_FACTOR * S_v / mu_concrete(E)
        tr = shield_transmission(E)
        phi_c = phi_lab * tr * BUILDUP_PAIR
        rc = np.interp(np.log(E), np.log(energies), [resp[e]['c'] for e in energies])
        rs = np.interp(np.log(E), np.log(energies), [resp[e]['s'] for e in energies])
        Bter['c'] += phi_c * rc
        Bter['s'] += phi_c * rs
        P(f"    {E:.4f} MeV ({ser}): lab fluence {phi_lab:.3f} /cm2/s, shield transmission {tr:.2e}, cavity {phi_c:.2e} /cm2/s -> coinc {phi_c*rc:.2e} /s, singles {phi_c*rs:.2e} /s")
    # radon in cavity (flushed): 214Bi lines; singles only matter
    P("    radon: cavity flushed with LN2 boil-off N2 to <1 Bq/m3 (22 L cavity -> <0.02 Bq 214Bi): <1e-4 /s singles, coincidences negligible")

    # 4. Michel positrons from mu+ stopping in the central hardware
    mass_c, eps_mich_c, eps_mich_h12, eps_mich_h5 = michel(lay)
    R_mich_c = STOP_RATE * MU_PLUS * mass_c * 1e-3
    P(f"\n[4] central hardware mass {mass_c:.0f} g; stopped mu+ there {R_mich_c:.2e} /s; Michel e+ -> coinc eff {eps_mich_c:.4f}, 12-30 MeV {eps_mich_h12:.4f}")
    # 5. cosmogenic beta+ in central hardware (SS/Cu: 58Co, 56Co, 48V, 52Mn; saturation)
    mass_b, eps_bp = central_bplus(lay)
    R_bplus = 1.5e-4 * mass_c * 1e-3            # e+/s: ~50 atoms/kg/d x ~0.25 beta+ branch (SS, [BK] Cebrian 2017)
    P(f"[5] cosmogenic beta+ in central hardware: {R_bplus:.1e} e+/s x coinc acceptance {eps_bp:.4f}")
    # 6. muons through crystals
    mu = muon_chords(lay)
    P(f"[6] cosmic muons through crystals: {mu['any']:.2f} /s (mean deposit {mu['mean_dep']:.0f} MeV); in 12-30 MeV {mu['h12']:.3f} /s; 5-30 MeV {mu['h5']:.3f} /s")
    # 7. hadrons in crystals (surface): sea-level neutron flux >20 MeV ~2.5e-3 /cm2/s outdoors
    #    (Gordon et al., IEEE TNS 51, 3427 (2004), JESD89A: 13 n/cm2/h above 10 MeV [BK]),
    #    x0.6 for one concrete floor above; non-elastic Sigma ~0.05/cm (BGO); mean chord 4V/S
    R_c, L_c = m7_geom.DETECTORS[dk][1], m7_geom.DETECTORS[dk][2]
    ndet = 2 * len(lay[2])
    V = np.pi * R_c ** 2 * L_c
    Sarea = 2 * np.pi * R_c * L_c + 2 * np.pi * R_c ** 2
    chord = 4 * V / Sarea
    sig = {'BGO': 0.053, 'NAI': 0.035, 'LABR3': 0.045, 'GE': 0.040}[mat]
    had_int = 2.5e-3 * 0.6 * (Sarea / 4) * (1 - np.exp(-sig * chord)) * ndet
    f_vis = {'h12': (0.3, 0.15, 0.5), 'h5': (0.5, 0.3, 0.7)}
    P(f"[7] hadron interactions in crystals (surface): {had_int:.3f} /s; visible fraction 12-30 MeV assumed 0.3 (0.15-0.5)")

    # ------------------------------------------------------------------ budgets
    def budget(site, case='central', veto=True, win=VETO_WINDOW):
        fm, fh = SITES[site]
        vi = VETO_INEFF[case] if veto else 1.0
        fhad = F_HAD[case] * fh
        k = {'central': 1.0, 'low': 1 / 3, 'high': 3.0}[case]      # generic MC/normalisation uncertainty on muonic terms
        delay = np.exp(-win / MU_TAU) if veto else 1.0
        acc2 = 2 * TAU2.get(mat, 20e-9) * 0.5 * 0.5                   # singles in window <=0.5 /s each
        c = {
            'mu showers (prompt)': k * fm * Bmu_c * vi,
            'mu+ Michel in shield (delayed)': k * fm * Bmu_c * r_mich * max(delay, vi),
            'hadronic (n,p) showers': k * Bmu_c * fhad,
            'Michel e+ in central hardware': fm * R_mich_c * eps_mich_c * max(delay, vi),
            'terrestrial gamma pair conversion': Bter['c'] * (1 if case == 'central' else (0.5 if case == 'low' else 2.0)),
            'cosmogenic beta+ (central)': R_bplus * eps_bp * (fm + fh) / 2,
            'accidentals (2 tau R1 R2)': acc2,
        }
        h = {
            'direct muons in crystals': fm * mu['h12'] * vi,
            'mu showers (prompt)': k * fm * Bmu_h12 * vi,
            'mu+ Michel (shield + crystals + central)': fm * (k * Bmu_h12 * r_mich + STOP_RATE * MU_PLUS * V * ndet * MATERIALS[mat].rho * 1e-3 * 0.4
                                                             + R_mich_c * eps_mich_h12) * max(delay, vi),
            'hadrons in crystals': fh * had_int * f_vis['h12'][{'central': 0, 'low': 1, 'high': 2}[case]],
            'hadronic showers in shield': k * Bmu_h12 * fhad,
        }
        return c, h

    P("\n==== BACKGROUND BUDGETS (counts/s) ====")
    table = {}
    for site in SITES:
        for veto in (True, False):
            if site != 'sea level' and not veto:
                continue
            for case in ('central', 'low', 'high'):
                c, h = budget(site, case, veto)
                table[(site, veto, case)] = (sum(c.values()), sum(h.values()), c, h)
            c, h = table[(site, veto, 'central')][2:]
            tag = f"{site}, {'veto (20 us ext.)' if veto else 'NO veto'}"
            P(f"\n-- {tag}")
            P("   511-511 coincidence:")
            for k_, v in c.items():
                P(f"      {k_:38s} {v:.2e} /s  ({v*DAY:.2f} /day)")
            P(f"      TOTAL {sum(c.values()):.2e} /s = {sum(c.values())*DAY:.2f} /day  (range {table[(site,veto,'low')][0]*DAY:.2f}-{table[(site,veto,'high')][0]*DAY:.2f} /day)")
            P("   calorimeter sum 12-30 MeV:")
            for k_, v in h.items():
                P(f"      {k_:38s} {v:.2e} /s  ({v*DAY:.1f} /day)")
            P(f"      TOTAL {sum(h.values()):.2e} /s = {sum(h.values())*DAY:.1f} /day  (range {table[(site,veto,'low')][1]*DAY:.1f}-{table[(site,veto,'high')][1]*DAY:.1f} /day)")
    # short veto window (1 us) to show why the extended window matters
    c1, h1 = budget('sea level', 'central', True, win=1e-6)
    P(f"\n   sea level with a 1 us veto window instead of 20 us: coinc {sum(c1.values())*DAY:.1f} /day, 12-30 MeV {sum(h1.values())*DAY:.0f} /day")
    P(f"   veto dead time: ~{VETO_RATE_PER_M2*3.0*VETO_WINDOW*100:.1f} % for ~3 m2 of veto at {VETO_RATE_PER_M2:.0f} /m2/s")

    # terrestrial term vs shield
    P("\n   terrestrial 511-511 term vs shield (per day):")
    for pb, cu in ((5, 0), (10, 0), (15, 0), (10, 5)):
        tot = 0.0
        for E, ser, yld in LINES:
            S_v = CONC[ser] * CONC_RHO * 1e-3 * yld
            phi_lab = ROOM_FACTOR * S_v / mu_concrete(E)
            tr = shield_transmission(E, pb=pb) * np.exp(-MATERIALS['CU'].mu_total(E)[0] * MATERIALS['CU'].rho * cu)
            rc = np.interp(np.log(E), np.log(energies), [resp[e]['c'] for e in energies])
            tot += phi_lab * tr * BUILDUP_PAIR * rc
        mass_pb = 11.35e-3 * ((2 * (15 + HDPE_T + pb + cu)) ** 3 - (2 * (15 + HDPE_T + cu)) ** 3)
        P(f"     Pb {pb:2d} cm + Cu {cu} cm: {tot*DAY:.2f} /day   (Pb mass for a 30 cm cavity + 5 cm HDPE: {mass_pb:.0f} kg)")

    # signal efficiencies for this layout and MDA
    rng = np.random.default_rng(77)
    g, sens = m7_geom.build(L_um=100, det=dk, **lay[1])
    mcs = MC(g, sens, rng)
    nsig = 8000
    mcs.run(nsig, charged=ipc_source(nsig, 100, rng))
    Es = smear(mcs.edep[:, 2:], mat, rng)
    Ssum = Es.sum(1)
    eff = {'511-511 coincidence': coinc_mask(Es, mat, lay[2]).mean(),
           'calorimeter 12-30 MeV': ((Ssum > 12) & (Ssum < 30)).mean()}
    both = (coinc_mask(Es, mat, lay[2]) | ((Ssum > 12) & (Ssum < 30))).mean()
    P(f"\n[8] IPC signal efficiency (100 um foil, uniform depth, {nsig} events): " + ", ".join(f"{k} {v:.4f}" for k, v in eff.items())
      + f"; OR of both {both:.4f}")
    P("\n==== MDA: e+e- pairs produced per second in the foil, 5 sigma median, 30 d (15 d on / 15 d off) ====")
    P("site / veto | channel | B [/day] (range) | eps | MDA [pairs/s] (range: low-B .. high-B)")
    mdas = {}
    for site in SITES:
        for veto in (True, False):
            if site != 'sea level' and not veto:
                continue
            for ci, ch in enumerate(eff):
                Bc = table[(site, veto, 'central')][ci]
                Bl = table[(site, veto, 'low')][ci]
                Bh = table[(site, veto, 'high')][ci]
                m_c = mda(Bc, eff[ch])[0]
                m_l = mda(Bl, eff[ch])[0]
                m_h = mda(Bh, eff[ch])[0]
                mdas[(site, veto, ch)] = (m_c, m_l, m_h, Bc)
                P(f"{site}, {'veto' if veto else 'no veto'} | {ch} | {Bc*DAY:.3g} ({Bl*DAY:.3g}-{Bh*DAY:.3g}) | {eff[ch]:.4f} | "
                  f"{m_c:.2e} ({m_l:.2e}-{m_h:.2e})")
    # combined (Stouffer-like: independent channels; approx MDA from summed Fisher information)
    P("\ncombined H_Cz test (both channels, same R): MDA where Z_A(511)^2 + Z_A(hi)^2 = 25")
    for site in SITES:
        chs = list(eff)
        B = [table[(site, True, 'central')][i] for i in range(2)]
        f = lambda R: sum(asimov_onoff(eff[c] * R * T_ON, B[i] * T_ON) ** 2 for i, c in enumerate(chs)) - 25  # noqa: E731
        R0 = brentq(f, 1e-7, 1e3)
        P(f"   {site}: {R0:.2e} pairs/s")
        mdas[(site, 'comb')] = R0

    # figures
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))
    bins = np.logspace(-1, 3, 61)
    h, _ = np.histogram(phi_E, bins, weights=phi_w)
    ax[0].step(np.sqrt(bins[1:] * bins[:-1]), h / np.diff(bins), where='mid')
    ax[0].set_xscale('log')
    ax[0].set_yscale('log')
    ax[0].set_xlabel('photon energy [MeV]')
    ax[0].set_ylabel('cavity fluence rate [cm$^{-2}$ s$^{-1}$ MeV$^{-1}$]')
    ax[0].set_title('muon delta-ray showers, 10 cm Pb + 5 cm HDPE (unvetoed)')
    Ees = np.array(energies)
    ax[1].loglog(Ees, [max(resp[e]['c'], 1e-4) for e in energies], 'o-', label='511-511 coinc.')
    ax[1].loglog(Ees, [max(resp[e]['h12'], 1e-4) for e in energies], 's-', label='sum 12-30 MeV')
    ax[1].loglog(Ees, [max(resp[e]['s'], 1e-4) for e in energies], '^-', label='singles in 511 window')
    ax[1].set_xlabel('photon energy [MeV]')
    ax[1].set_ylabel('response [counts per unit fluence, cm$^2$]')
    ax[1].set_title(f'{layout_name}: response to isotropic cavity photons')
    ax[1].legend(fontsize=8)
    labels, vals = [], []
    for site in SITES:
        for veto in (True, False):
            if site != 'sea level' and not veto:
                continue
            c, h_ = budget(site, 'central', veto)
            labels.append(f"{site}\n{'veto' if veto else 'no veto'}")
            vals.append(c)
    keys = list(vals[0])
    bottom = np.zeros(len(vals))
    for k_ in keys:
        v = np.array([x[k_] * DAY for x in vals])
        ax[2].bar(range(len(vals)), v, bottom=bottom, label=k_)
        bottom += v
    ax[2].set_xticks(range(len(vals)))
    ax[2].set_xticklabels(labels, fontsize=8)
    ax[2].set_ylabel('511-511 background [counts/day]')
    ax[2].set_title('coincidence background budget (central)')
    ax[2].legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'm7_background.png'), dpi=130)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    days = np.linspace(2, 120, 60)
    for (site, veto, ch), st in [(('sea level', True, '511-511 coincidence'), '-'), (('sea level', True, 'calorimeter 12-30 MeV'), '--'),
                                 (('shallow 30 m w.e.', True, '511-511 coincidence'), '-'), (('shallow 30 m w.e.', True, 'calorimeter 12-30 MeV'), '--')]:
        B = mdas[(site, veto, ch)][3]
        ax.loglog(days, [mda(B, eff[ch], t_on=d / 2 * DAY)[0] for d in days], st, label=f'{site}: {ch}')
    ax.set_xlabel('total run time [days] (50 % on)')
    ax.set_ylabel('MDA [e+e- pairs / s in foil], 5 sigma median')
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(fontsize=8)
    ax.set_title(f'{layout_name}')
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'm7_mda.png'), dpi=130)
    plt.close(fig)
    P(f"\nrun time {time.time()-t0:.0f} s")
    return lines, table, resp, (phi_E, phi_w), mu


if __name__ == '__main__':
    lines, *_ = main(sys.argv[1] if len(sys.argv) > 1 else None)
    txt = '\n'.join(lines)
    open(os.path.join(OUT, 'm7_background.txt'), 'w').write(txt + '\n')
    print(txt)
