"""M7 / Q2-Q3: detection efficiencies of candidate gamma-detector layouts
around the reference foil assembly (m7_geom), from the coupled MC.

  * 511-511 back-to-back coincidence efficiency for
      - a pure annihilation point source at the foil centre ("annihilation
        at the foil", as the brief asks),
      - a 22Na calibration source at the foil (what is actually measured),
      - E0 IPC pairs born uniformly in a 100 um foil (the physics signal).
  * high-energy (5-30 MeV, 12-30 MeV) calorimetric efficiency for IPC pairs.
Outputs: docs/models/figs/m7_detectors.txt, m7_eff_layouts.png, m7_sum_spectra.png
"""
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from m7_mc import MC  # noqa: E402
from m7_ipc import fermi  # noqa: E402
from m7_common import ME  # noqa: E402
import m7_geom  # noqa: E402
from m7_transport import ipc_source  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), '..', 'docs', 'models', 'figs')

# Energy resolution FWHM/E at 662 keV [BK: Saint-Gobain/Scionix/Ortec data sheets]
RES662 = {'BGO': 0.105, 'NAI': 0.070, 'LABR3': 0.029}
RESFLOOR = {'BGO': 0.02, 'NAI': 0.015, 'LABR3': 0.006}
# coincidence resolving time 2*tau [s] (digitiser CFD, [BK])
TAU2 = {'BGO': 20e-9, 'NAI': 10e-9, 'LABR3': 2e-9, 'GE': 50e-9}
# Temperature coefficient of light output [%/K] [BK: Saint-Gobain; Moszynski]
TEMPCO = {'BGO': -1.2, 'NAI': -0.3, 'LABR3': -0.01, 'GE': 0.0}


def fwhm(mat, E):
    E = np.asarray(E, float)
    if mat == 'GE':      # HPGe: 1.0 keV electronic (+) 1.9 keV statistical at 1332 keV -> 1.55 keV at 511
        return np.sqrt(0.0010 ** 2 + 0.0019 ** 2 * np.clip(E, 0, None) / 1.332)
    r = np.maximum(RES662[mat] * np.sqrt(0.662 / np.clip(E, 1e-6, None)), RESFLOOR[mat])
    return r * E


def win511(mat):
    if mat == 'GE':
        return 0.511 - 0.004, 0.511 + 0.004     # includes ~2.5 keV Doppler FWHM
    w = float(fwhm(mat, 0.511))
    return 0.511 - w, 0.511 + w                   # +-1 FWHM = 98 % of a Gaussian peak


def smear(ed, mat, rng):
    s = fwhm(mat, ed) / 2.3548
    return np.where(ed > 0, ed + rng.normal(size=ed.shape) * s, 0.0)


def coinc_mask(E, mat, pairs):
    lo, hi = win511(mat)
    inw = (E > lo) & (E < hi)
    m = np.zeros(len(E), bool)
    for a, b in pairs:
        m |= inw[:, a] & inw[:, b]
    return m


def point_511(g, sens, n, rng, pos=(0, 0, 0)):
    mc = MC(g, sens, rng)
    c = rng.uniform(-1, 1, n)
    f = rng.uniform(0, 2 * np.pi, n)
    s = np.sqrt(1 - c * c)
    d = np.stack([s * np.cos(f), s * np.sin(f), c], 1)
    P = np.tile(pos, (n, 1)).astype(float)
    ev = np.arange(n)
    mc.run(n, photons={'ev': np.concatenate([ev, ev]), 'P': np.concatenate([P, P]),
                       'D': np.concatenate([d, -d]), 'E': np.full(2 * n, ME)})
    return mc


def na22(g, sens, n, rng, L_um=100.0):
    """n 22Na decays at the foil centre: beta+ (90.3 %) + 1274.5 keV gamma (99.9 %)."""
    mc = MC(g, sens, rng)
    Q = 0.5459
    T = np.linspace(1e-4, Q - 1e-4, 2000)
    p = np.sqrt(T * (T + 2 * ME))
    E = T + ME
    w = p * E * (Q - T) ** 2 * fermi(-10, E / ME)
    cdf = np.cumsum(w)
    cdf /= cdf[-1]
    isb = rng.uniform(size=n) < 0.903
    nb = isb.sum()
    Tb = np.interp(rng.uniform(size=nb), cdf, T)
    iso = lambda k: (lambda c, f: np.stack([np.sqrt(1 - c * c) * np.cos(f), np.sqrt(1 - c * c) * np.sin(f), c], 1))(
        rng.uniform(-1, 1, k), rng.uniform(0, 2 * np.pi, k))  # noqa: E731
    z0 = -L_um * 1e-4 / 2
    ev = np.arange(n)
    mc.run(n, charged={'ev': ev[isb], 'q': np.ones(nb, int), 'P': np.tile([0, 0, z0], (nb, 1)),
                       'D': iso(nb), 'T': Tb},
           photons={'ev': ev, 'P': np.tile([0, 0, z0], (n, 1)), 'D': iso(n), 'E': np.full(n, 1.2745)})
    # note: run() transports charged first, then photons, so both are included
    return mc


LAYOUTS = {
    # name: (detector key, build kwargs, coincidence pairs)
    'BGO 3x3 pair':       ('BGO3x3', {}, [(0, 1)]),
    'NaI 3x3 pair':       ('NAI3x3', {}, [(0, 1)]),
    'LaBr3 3x3 pair':     ('LABR3x3', {}, [(0, 1)]),
    'HPGe 100% pair':     ('HPGE100', {}, [(0, 1)]),
    'BGO 4x4 pair':       ('BGO4x4', {}, [(0, 1)]),
    'NaI 5x5 pair':       ('NAI5x5', {}, [(0, 1)]),
    'BGO 3x3 quad (+-x,+-y)': ('BGO3x3', {'four': True}, [(0, 1), (2, 3)]),
    'BGO 3x3 z-pair':     ('BGO3x3', {'axes': 'zpair'}, [(0, 1)]),
    'NaI 5x5 z-pair':     ('NAI5x5', {'axes': 'zpair'}, [(0, 1)]),
    'BGO 3x3 quad (+-x,+-z)': ('BGO3x3', {'axes': 'xz'}, [(0, 1), (2, 3)]),
}


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    t0 = time.time()
    lines = []
    P = lines.append
    rng = np.random.default_rng(2024)
    res = {}
    spectra = {}
    n_pt, n_ipc = 40000, 6000
    for name, (dk, kw, pairs) in LAYOUTS.items():
        mat = m7_geom.DETECTORS[dk][0]
        g, sens = m7_geom.build(L_um=100, det=dk, **kw)
        nd = len(sens) - 2
        # (1) pure 511 pair at the foil centre
        mc = point_511(g, sens, n_pt, rng)
        E = smear(mc.edep[:, 2:], mat, rng)
        e_pt = coinc_mask(E, mat, pairs).mean()
        lo, hi = win511(mat)
        e_single = ((E > lo) & (E < hi)).any(1).mean()
        # (2) 22Na at the foil
        g, sens = m7_geom.build(L_um=100, det=dk, **kw)
        mc = na22(g, sens, n_pt, rng)
        E = smear(mc.edep[:, 2:], mat, rng)
        e_na = coinc_mask(E, mat, pairs).mean()
        # (3) IPC events
        g, sens = m7_geom.build(L_um=100, det=dk, **kw)
        mc = MC(g, sens, rng)
        mc.run(n_ipc, charged=ipc_source(n_ipc, 100, rng))
        E = smear(mc.edep[:, 2:], mat, rng)
        e_ipc = coinc_mask(E, mat, pairs).mean()
        S = E.sum(1)
        e5 = ((S > 5) & (S < 30)).mean()
        e12 = ((S > 12) & (S < 30)).mean()
        e_any5 = (E.max(1) > 5).mean()
        spectra[name] = S
        res[name] = (e_pt, e_na, e_ipc, e5, e12, e_any5, e_single)
        vol = np.pi * m7_geom.DETECTORS[dk][1] ** 2 * m7_geom.DETECTORS[dk][2] * nd
        P(f"{name:26s} ({nd} x {dk}, {vol:5.0f} cm3 total): 511-511 coinc eff: point@foil {e_pt:.4f} | 22Na@foil {e_na:.4f} "
          f"| IPC {e_ipc:.4f} ;  IPC sum-E in 5-30 MeV {e5:.4f}, 12-30 MeV {e12:.4f}; any det >5 MeV {e_any5:.4f}; "
          f"window {lo*1e3:.0f}-{hi*1e3:.0f} keV")
    # (4) z-offset scan of the BGO pair
    P("\nBGO 3x3 pair: dependence on detector-axis height zc (foil at z=0, electrolyte below)")
    for zc in [1.0, -0.5, -2.0]:
        g, sens = m7_geom.build(L_um=100, det='BGO3x3', zc=zc)
        mc = MC(g, sens, rng)
        mc.run(n_ipc, charged=ipc_source(n_ipc, 100, rng))
        E = smear(mc.edep[:, 2:], 'BGO', rng)
        P(f"  zc = {zc:+.1f} cm: IPC 511-511 {coinc_mask(E, 'BGO', [(0,1)]).mean():.4f}; 12-30 MeV {((E.sum(1)>12)&(E.sum(1)<30)).mean():.4f}")
    # (5) foil thickness dependence for the BGO pair
    P("\nBGO 3x3 pair: dependence on foil thickness (uniform source depth)")
    for L in [25, 300]:
        g, sens = m7_geom.build(L_um=L, det='BGO3x3')
        mc = MC(g, sens, rng)
        mc.run(n_ipc, charged=ipc_source(n_ipc, L, rng))
        E = smear(mc.edep[:, 2:], 'BGO', rng)
        P(f"  L = {L} um: IPC 511-511 {coinc_mask(E, 'BGO', [(0,1)]).mean():.4f}; 12-30 MeV {((E.sum(1)>12)&(E.sum(1)<30)).mean():.4f}")
    # (6) aluminium chamber instead of SS (less e+ stopping near foil?)
    g, sens = m7_geom.build(L_um=100, det='BGO3x3', chamber_mat='AL', chamber_wall=0.3)
    mc = MC(g, sens, rng)
    mc.run(n_ipc, charged=ipc_source(n_ipc, 100, rng))
    E = smear(mc.edep[:, 2:], 'BGO', rng)
    P(f"\nBGO 3x3 pair with 3 mm Al chamber wall instead of 1.5 mm SS: IPC 511-511 {coinc_mask(E,'BGO',[(0,1)]).mean():.4f}; "
      f"12-30 MeV {((E.sum(1)>12)&(E.sum(1)<30)).mean():.4f}")
    P(f"\nstatistical 1-sigma: point/22Na ~{np.sqrt(0.05/n_pt):.4f}, IPC ~{np.sqrt(0.05/n_ipc):.4f} (at eff 0.05)")
    P(f"run time {time.time()-t0:.0f} s")

    # figures
    fig, ax = plt.subplots(figsize=(10, 4.5))
    names = list(res)
    x = np.arange(len(names))
    for i, (lab, j) in enumerate([('511-511, point at foil', 0), ('511-511, 22Na at foil', 1), ('511-511, IPC pairs', 2),
                                  ('IPC sum 5-30 MeV', 3), ('IPC sum 12-30 MeV', 4)]):
        ax.bar(x + (i - 2) * 0.16, [res[k][j] for k in names], 0.16, label=lab)
    ax.set_xticks(x)
    ax.set_xticklabels([k.replace(' pair', '\npair').replace(' quad', '\nquad') for k in names], fontsize=7)
    ax.set_ylabel('efficiency per event')
    ax.set_title('Gamma-detector layouts around the foil (MC)')
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'm7_eff_layouts.png'), dpi=130)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for k in ['BGO 3x3 pair', 'NaI 3x3 pair', 'BGO 4x4 pair', 'BGO 3x3 quad (+-x,+-y)']:
        S = spectra[k]
        ax.hist(S[S > 0.05], bins=120, range=(0, 30), histtype='step', label=k, density=False)
    ax.set_yscale('log')
    ax.set_xlabel('summed deposited energy (smeared) [MeV]')
    ax.set_ylabel(f'events / 0.25 MeV ({n_ipc} IPC events)')
    ax.set_title('IPC signal in the calorimeter sum')
    ax.axvline(12, color='k', ls=':')
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'm7_sum_spectra.png'), dpi=130)
    plt.close(fig)
    txt = '\n'.join(lines)
    open(os.path.join(OUT, 'm7_detectors.txt'), 'w').write(txt + '\n')
    print(txt)
    return res


if __name__ == '__main__':
    main()
