"""M7 / Q3 option: dedicated near-4pi NaI(Tl) well calorimeter around a small
permeation cell (no Si telescope, no UHV chamber). Tests how much of the
23.85 MeV IPC energy a hermetic calorimeter recovers, and the resulting MDA.

Geometry (cm): NaI(Tl) cylinder 8" x 8" (R 10.16, half-length 10.16) with a
5.2 cm diameter well, 13.2 cm deep from the top, 0.5 mm Al liner. Inside, a
mini cell: Pd foil (r = 1 cm) at the calorimeter centre; D2O below (r 1.5 cm,
2 cm deep) in PTFE (3 mm); above the foil a 2 cm high evacuated space closed
by a 1 mm SS cap (D2 permeate pumped through a 6 mm tube, not modelled).
Output: docs/models/figs/m7_wellcal.txt, m7_wellcal.png
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from m7_mc import MC, Geometry, Tube  # noqa: E402
from m7_detectors import smear  # noqa: E402
from m7_transport import ipc_source  # noqa: E402
from m7_background import mda, J1, VETO_INEFF, SITES, T_ON, DAY  # noqa: E402
from m7_common import MATERIALS, PDG_DEDX_MIN_BK  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), '..', 'docs', 'models', 'figs')
R_NAI, HL_NAI, R_WELL, ZW = 10.16, 10.16, 2.6, -3.0
S_QUAD = 2 * (2 * np.pi * 5.08 * 10.16 + 2 * np.pi * 5.08 ** 2)  # surface of the reference BGO 4x4 pair [cm2]


def _bmu_h12():
    """Unvetoed muon-shower rate in 12-30 MeV for the BGO 4x4 pair, from m7_background.txt."""
    try:
        for ln in open(os.path.join(OUT, 'm7_background.txt')):
            if 'muon-shower cavity photons' in ln:
                return float(ln.split('sum 12-30 MeV')[1].split('/s')[0])
    except (OSError, IndexError, ValueError):
        pass
    return 0.165


BMU_H12 = _bmu_h12()


def build(L_um=100.0):
    L = L_um * 1e-4
    g = Geometry([30, 30, 30], 'AIR')
    g.add(Tube(2, [0, 0, -L / 2], 0, 1.2, L / 2), 'PD', 'foil')
    g.add(Tube(2, [0, 0, 1.0], 0, 1.2, 1.0), 'VAC', 'vac')
    g.add(Tube(2, [0, 0, 1.0], 1.2, 1.3, 1.0), 'SS316', 'capwall')
    g.add(Tube(2, [0, 0, 2.05], 0, 1.3, 0.05), 'SS316', 'cap')
    g.add(Tube(2, [0, 0, -L - 1.0], 0, 1.5, 1.0), 'D2O', 'electrolyte')
    g.add(Tube(2, [0, 0, -L - 1.0], 1.5, 1.8, 1.0), 'PTFE', 'cellwall')
    g.add(Tube(2, [0, 0, -L - 2.15], 0, 1.8, 0.15), 'PTFE', 'cellbottom')
    g.add(Tube(2, [0, 0, (HL_NAI + ZW) / 2], R_WELL, R_WELL + 0.05, (HL_NAI - ZW) / 2), 'AL', 'liner')
    g.add(Tube(2, [0, 0, ZW - 0.025], 0, R_WELL + 0.05, 0.025), 'AL', 'linerbottom')
    g.add(Tube(2, [0, 0, (HL_NAI + ZW) / 2], 0, R_WELL + 0.05, (HL_NAI - ZW) / 2), 'AIR', 'well')
    g.add(Tube(2, [0, 0, 0], 0, R_NAI, HL_NAI), 'NAI', 'nai')
    g.add(Tube(2, [0, 0, 0], 0, R_NAI + 0.1, HL_NAI + 0.1), 'AL', 'can')
    return g


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rng = np.random.default_rng(31)
    lines = []
    P = lines.append
    n = 5000
    g = build()
    mc = MC(g, ['nai'], rng)
    mc.run(n, charged=ipc_source(n, 100, rng))
    E = smear(mc.edep[:, 0], 'NAI', rng)
    e12 = ((E > 12) & (E < 30)).mean()
    e5 = ((E > 5) & (E < 30)).mean()
    pk = ((E > 23.85 * 0.93) & (E < 23.85 * 1.07)).mean()
    P(f"8\"x8\" NaI well calorimeter, IPC pairs in 100 um foil: eff 12-30 MeV {e12:.3f}; 5-30 MeV {e5:.3f}; "
      f"full-energy peak 23.85 MeV +-7 % {pk:.3f}; mean deposit {mc.edep.mean():.2f} MeV")
    # backgrounds (same physics as m7_background, scaled to this crystal)
    V = np.pi * R_NAI ** 2 * 2 * HL_NAI - np.pi * R_WELL ** 2 * (HL_NAI - ZW)
    S = 2 * np.pi * R_NAI * 2 * HL_NAI + 2 * np.pi * R_NAI ** 2
    chord = 4 * V / S
    # direct muons: chord MC through the cylinder
    m = 400000
    A = 40.0
    u = rng.uniform(size=m)
    ct = (1 - u) ** 0.25
    ph = rng.uniform(0, 2 * np.pi, m)
    st = np.sqrt(1 - ct ** 2)
    D = np.stack([st * np.cos(ph), st * np.sin(ph), -ct], 1)
    P0 = np.stack([rng.uniform(-A, A, m), rng.uniform(-A, A, m), np.full(m, 30.0)], 1)
    T = Tube(2, [0, 0, 0], 0, R_NAI, HL_NAI)
    t1 = T.ray(P0, D)
    hit = np.isfinite(t1)
    t2 = np.where(hit, T.ray(P0 + D * (np.where(hit, t1, 0) + 1e-6)[:, None], D), 0)
    dep = t2 * PDG_DEDX_MIN_BK['NAI'] * MATERIALS['NAI'].rho * 1.1
    rate = J1 * (2 * A) ** 2 / m
    mu_any, mu12 = rate * hit.sum(), rate * ((dep > 12) & (dep < 30)).sum()
    had = 2.5e-3 * 0.6 * (S / 4) * (1 - np.exp(-0.035 * chord))
    stop = 0.012 * 0.56 * V * MATERIALS['NAI'].rho * 1e-3
    P(f"muons through crystal {mu_any:.2f} /s, in 12-30 MeV window {mu12:.3f} /s; hadron interactions {had:.3f} /s; stopped mu+ {stop:.3f} /s")
    # shield-shower photons: scale BGO-quad result by crystal surface (sphere-source response ~ area x efficiency)
    P("\nMDA (pairs/s in foil, 5 sigma median, 30 d on/off), 12-30 MeV window:")
    res = {}
    for site, (fm, fh) in SITES.items():
        for case in ('central', 'low', 'high'):
            vi = VETO_INEFF[case]
            fv = {'central': 0.3, 'low': 0.15, 'high': 0.5}[case]
            k = {'central': 1.0, 'low': 1 / 3, 'high': 3.0}[case]
            fhad = {'central': 0.5, 'low': 0.25, 'high': 1.0}[case]
            B = fm * mu12 * vi + fh * had * fv + fm * stop * 0.4 * max(np.exp(-20 / 2.197), vi) \
                + k * BMU_H12 * (S / S_QUAD) * (fm * vi + fhad * fh)
            res[(site, case)] = B
        Bc = res[(site, 'central')]
        m_c = mda(Bc, e12)[0]
        m_l = mda(res[(site, 'low')], e12)[0]
        m_h = mda(res[(site, 'high')], e12)[0]
        P(f"  {site:20s}: B = {Bc*DAY:9.1f} /day ({res[(site,'low')]*DAY:.1f}-{res[(site,'high')]*DAY:.1f}); MDA {m_c:.2e} ({m_l:.2e}-{m_h:.2e}) pairs/s")
    fig, ax = plt.subplots(figsize=(7, 4.3))
    ax.hist(E[E > 0.05], bins=130, range=(0, 26), histtype='step')
    ax.axvspan(12, 26, color='g', alpha=0.08)
    ax.set_yscale('log')
    ax.set_xlabel('deposited energy (NaI resolution) [MeV]')
    ax.set_ylabel(f'events / 0.2 MeV ({n} IPC events)')
    ax.set_title('8"x8" NaI well calorimeter: IPC signal (23.85 MeV = full absorption)')
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'm7_wellcal.png'), dpi=130)
    plt.close(fig)
    txt = '\n'.join(lines)
    open(os.path.join(OUT, 'm7_wellcal.txt'), 'w').write(txt + '\n')
    print(txt)


if __name__ == '__main__':
    main()
