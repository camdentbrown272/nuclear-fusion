"""M7 / Q1(a-c): transport of E0 IPC e+e- pairs born in a Pd foil
(electrolyte on one side, UHV on the other) in the reference geometry
(m7_geom.build).

Outputs: docs/models/figs/m7_transport.txt, m7_annihilation.png,
m7_brems.png, m7_edep.png
Run time ~2-4 min on 4 CPUs (single-threaded numpy).
"""
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from m7_mc import MC  # noqa: E402
from m7_ipc import sample_pairs, TSUM  # noqa: E402
from m7_common import MATERIALS  # noqa: E402
import m7_geom  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), '..', 'docs', 'models', 'figs')


def ipc_source(n, L_um, rng, depth='uniform', r_src=1.0):
    """Primary e+ and e- for n IPC events. depth: 'uniform' | 'front' | 'back'."""
    L = L_um * 1e-4
    Tp, Tm, dp, dm = sample_pairs(n, rng)
    r = r_src * np.sqrt(rng.uniform(size=n))
    ph = rng.uniform(0, 2 * np.pi, n)
    if depth == 'uniform':
        z = -L * rng.uniform(size=n)
    elif depth == 'front':      # within 0.5 um of the vacuum face
        z = -0.5e-4 * rng.uniform(size=n)
    else:                       # within 0.5 um of the electrolyte face
        z = -L + 0.5e-4 * rng.uniform(size=n)
    P = np.stack([r * np.cos(ph), r * np.sin(ph), z], 1)
    ev = np.arange(n)
    return {'ev': np.concatenate([ev, ev]), 'q': np.concatenate([np.ones(n, int), -np.ones(n, int)]),
            'P': np.concatenate([P, P]), 'D': np.concatenate([dp, dm]), 'T': np.concatenate([Tp, Tm])}


def group_of(names):
    lut = {}
    for grp, regs in m7_geom.REGION_GROUPS.items():
        for r in regs:
            lut[r] = grp
    return [lut.get(n, 'other') for n in names]


def run_config(L_um, depth='uniform', n=20000, seed=11, det='BGO3x3', **kw):
    rng = np.random.default_rng(seed)
    g, sens = m7_geom.build(L_um=L_um, det=det, **kw)
    mc = MC(g, sens, rng)
    src = ipc_source(n, L_um, rng, depth)
    mc.run(n, charged=src)
    return mc, sens


def summarize_ann(mc, n):
    names = mc.g.reg_names
    grp = group_of(names)
    ann = mc.ann
    res = {}
    for gname in list(m7_geom.REGION_GROUPS) + ['escaped world (>30 cm)']:
        res[gname] = 0
    for ev, reg, fl in ann:
        res[grp[int(reg)]] += 1
    esc_pos = sum(1 for e in mc.esc_q if e[1] > 0)
    res['escaped world (>30 cm)'] = esc_pos
    inflight = int(np.sum(ann[:, 2])) if len(ann) else 0
    # annihilations from the primary positron only are not separable from
    # pair-produced secondary positrons; report per IPC event
    return {k: v / n for k, v in res.items()}, inflight / n, len(ann) / n


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    lines = []
    P = lines.append
    t0 = time.time()
    P("M7 Q1 transport of E0 IPC pairs (T+ + T- = %.3f MeV) in reference geometry" % TSUM)
    P("Per IPC event. 'annihilations' include positrons from secondary pair production (so sums can exceed 1).")
    # ------------------------------------------------ analytic orientation
    pd = MATERIALS['PD']
    for L in [25, 50, 100, 300]:
        x = L * 1e-4 * pd.rho
        P(f"L = {L:3d} um Pd = {x:.3f} g/cm2 = {x/pd.X0:.4f} X0; CSDA range of 1 MeV e+ in Pd = "
          f"{pd.csda_range(1.0, True)[0]/pd.rho*1e4:.0f} um, of 11.4 MeV = {pd.csda_range(11.4, True)[0]/pd.rho*1e4:.0f} um")
    for m in ['D2O', 'PTFE', 'SS316', 'SI', 'BGO']:
        mm = MATERIALS[m]
        P(f"CSDA range of 11.4 MeV e- in {m:6s}: {mm.csda_range(11.4)[0]/mm.rho:.2f} cm ; X0 = {mm.X0/mm.rho:.2f} cm")
    # ------------------------------------------------ main scan
    configs = [(25, 'uniform'), (50, 'uniform'), (100, 'uniform'), (300, 'uniform'),
               (100, 'front'), (100, 'back')]
    n = 10000
    table = {}
    keep = {}
    for L, dep in configs:
        mc, sens = run_config(L, dep, n=n, seed=100 + L + len(dep))
        fr, infl, nann = summarize_ann(mc, n)
        table[(L, dep)] = (fr, infl, nann)
        keep[(L, dep)] = mc
        P(f"\n-- L = {L} um, source depth = {dep}: annihilations/event = {nann:.3f}, in flight = {infl:.4f}")
        for k, v in fr.items():
            P(f"   {k:32s} {v:7.4f}")
        # brems produced in foil vs all
        br = np.array(mc.brem).reshape(-1, 3)
        foil_id = mc.g.reg_names.index('foil')
        kf = br[br[:, 1] == foil_id, 2]
        P(f"   bremsstrahlung photons (k>20 keV) created per event: total {len(br)/n:.2f}, in foil {len(kf)/n:.3f}; "
          f"energy radiated: total {br[:,2].sum()/n:.3f} MeV, in foil {kf.sum()/n:.4f} MeV")
        # deposits
        ed = mc.edep
        si = ed[:, 0] + ed[:, 1]
        P(f"   Si telescope: P(E-det dep > 50 keV) = {(ed[:,1]>0.05).mean():.4f}; P(dE-det > 5 keV) = {(ed[:,0]>0.005).mean():.4f}; "
          f"median E-det dep (hits) = {np.median(ed[ed[:,1]>0.05,1]) if (ed[:,1]>0.05).any() else 0:.3f} MeV; "
          f"P(E-det dep in 2.6-3.1 MeV p window) = {((ed[:,1]>2.6)&(ed[:,1]<3.1)).mean():.2e}")
        bg = ed[:, 2:]
        P(f"   BGO 3x3 pair at +-x: P(any det > 5 MeV) = {(bg.max(1)>5).mean():.4f}; P(sum > 12 MeV) = {(bg.sum(1)>12).mean():.4f}; "
          f"<sum dep> = {bg.sum(1).mean():.2f} MeV")
    # ------------------------------------------------ summary table for doc
    P("\nSUMMARY: fraction of IPC events whose positron annihilates in each region (uniform depth)")
    hdr = ['foil', 'electrolyte', 'PTFE cell', 'SS flange/chamber', 'Si telescope + PCB', 'gamma detectors (crystal+can)', 'air / outside', 'enclosure (HDPE liner/moderator)', 'escaped world (>30 cm)']
    P("L_um | " + " | ".join(hdr) + " | in-flight")
    for L in [25, 50, 100, 300]:
        fr, infl, _ = table[(L, 'uniform')]
        P(f"{L} | " + " | ".join(f"{fr[h]:.3f}" for h in hdr) + f" | {infl:.3f}")
    # ------------------------------------------------ figures
    mc = keep[(100, 'uniform')]
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    labels = []
    for L in [25, 100, 300]:
        fr, infl, _ = table[(L, 'uniform')]
        labels.append(L)
    x = np.arange(len(hdr))
    w = 0.25
    for i, L in enumerate([25, 100, 300]):
        fr, infl, _ = table[(L, 'uniform')]
        ax[0].bar(x + (i - 1) * w, [fr[h] for h in hdr], w, label=f'L = {L} um')
    ax[0].set_xticks(x)
    ax[0].set_xticklabels([h.replace(' (crystal+can)', '').replace(' ', '\n', 1) for h in hdr], fontsize=7)
    ax[0].set_ylabel('fraction of IPC events')
    ax[0].set_title('(a) where the positron annihilates')
    ax[0].legend(fontsize=8)
    ap = mc.ann_pos
    ax[1].hist(ap[:, 2], bins=120, range=(-6, 8), histtype='step', lw=1.3)
    ax[1].set_xlabel('z of annihilation vertex [cm] (foil at z=0, electrolyte z<0)')
    ax[1].set_ylabel('annihilations / bin (10k events)')
    ax[1].set_yscale('log')
    ax[1].set_title('annihilation vertex z (L = 100 um, uniform)')
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'm7_annihilation.png'), dpi=130)
    plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    for L in [25, 300]:
        m = keep[(L, 'uniform')]
        br = np.array(m.brem).reshape(-1, 3)
        foil_id = m.g.reg_names.index('foil')
        bins = np.linspace(0, 23, 47)
        h_all, _ = np.histogram(br[:, 2], bins)
        h_f, _ = np.histogram(br[br[:, 1] == foil_id, 2], bins)
        c = 0.5 * (bins[1:] + bins[:-1])
        ax[0].step(c, h_all / n / np.diff(bins), where='mid', label=f'all materials, L={L} um')
        ax[0].step(c, h_f / n / np.diff(bins), where='mid', ls='--', label=f'in foil only, L={L} um')
    ax[0].set_yscale('log')
    ax[0].set_xlabel('bremsstrahlung photon energy k [MeV]')
    ax[0].set_ylabel('photons / IPC event / MeV')
    ax[0].set_title('(b) bremsstrahlung produced')
    ax[0].legend(fontsize=7)
    m = keep[(100, 'uniform')]
    eg = np.array(m.esc_g).reshape(-1, 3)
    bins = np.linspace(0, 24, 97)
    h, _ = np.histogram(eg[:, 1], bins)
    c = 0.5 * (bins[1:] + bins[:-1])
    ax[1].step(c, h / n / np.diff(bins), where='mid', label='photons leaving 50 cm world box (L=100 um)')
    ax[1].set_yscale('log')
    ax[1].set_xlabel('photon energy [MeV]')
    ax[1].set_ylabel('photons / IPC event / MeV')
    ax[1].set_title('photon spectrum escaping the assembly (incl. 511 keV)')
    ax[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'm7_brems.png'), dpi=130)
    plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    ed = m.edep
    ax[0].hist(ed[ed[:, 1] > 0.005, 1] * 1e3, bins=100, range=(0, 2000), histtype='step', label='Si E (1000 um)')
    ax[0].hist(ed[ed[:, 0] > 0.0005, 0] * 1e3, bins=100, range=(0, 2000), histtype='step', label='Si dE (25 um)')
    ax[0].axvspan(2600, 3100, color='r', alpha=0.1)
    ax[0].set_xlabel('deposited energy [keV]')
    ax[0].set_ylabel('events / 20 keV (10k IPC events)')
    ax[0].set_yscale('log')
    ax[0].set_title('(c) Si telescope deposits (p window 2.6-3.1 MeV is off-scale right)')
    ax[0].legend(fontsize=8)
    s = ed[:, 2:].sum(1)
    ax[1].hist(ed[:, 2][ed[:, 2] > 0.01], bins=120, range=(0, 24), histtype='step', label='one BGO 3"x3"')
    ax[1].hist(s[s > 0.01], bins=120, range=(0, 24), histtype='step', label='sum of both BGO')
    ax[1].set_yscale('log')
    ax[1].set_xlabel('deposited energy [MeV]')
    ax[1].set_ylabel('events / 0.2 MeV')
    ax[1].set_title('(c) BGO deposits (no resolution smearing)')
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'm7_edep.png'), dpi=130)
    plt.close(fig)
    P(f"\nrun time {time.time()-t0:.0f} s")
    txt = '\n'.join(lines)
    open(os.path.join(OUT, 'm7_transport.txt'), 'w').write(txt + '\n')
    print(txt)


if __name__ == '__main__':
    main()
