"""M7 / Q5: calibration with exempt-quantity sources.

Exempt quantities: 10 CFR 30.71 Schedule B (US NRC,
https://www.nrc.gov/reading-rm/doc-collections/cfr/part030/part030-0071.html)
[BK values]; sources sold as exempt check sources under a distributor licence
(10 CFR 32.18/32.19) need no user licence. Nuclear data: NNDC/ENSDF [BK].

Computes, for the recommended layout, the full-energy-peak efficiency of a
point source at the foil position and the time to reach 1 % statistics.
Outputs: docs/models/figs/m7_calibration.txt
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from m7_mc import MC  # noqa: E402
from m7_detectors import smear, fwhm, coinc_mask, na22, LAYOUTS  # noqa: E402
import m7_geom  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), '..', 'docs', 'models', 'figs')
UCI = 3.7e4   # Bq per uCi

# nuclide: (Schedule-B exempt quantity [uCi], lines [(E MeV, yield)], what it calibrates)
SOURCES = {
    'Na-22':  (10.0, [(0.511, 1.806), (1.2745, 0.9994)], '511-511 coincidence efficiency at the foil position; 1275 keV energy point; timing (CRT)'),
    'Cs-137': (10.0, [(0.6617, 0.851)], 'energy scale and resolution at 662 keV; gain-stabilisation reference (removed during runs)'),
    'Co-60':  (1.0, [(1.1732, 0.9985), (1.3325, 0.9998)], 'energy scale 1.17/1.33 MeV; 2.505 MeV sum peak tests summing/coincidence logic'),
    'Am-241': (0.01, [(0.0595, 0.359)], 'Si telescope: 5.486 MeV alpha energy scale and dE-E thresholds; low-energy gamma threshold'),
    'Ba-133': (10.0, [(0.356, 0.62), (0.081, 0.33)], 'low-energy scale / threshold, 80-400 keV'),
    'Ge-68':  (0.1, [(0.511, 1.78)], '(not in Schedule B -> 0.1 uCi default) harder e+ (1.9 MeV end-point) for e+ transport/annihilation-location check'),
    'Sr-90':  (0.1, [], 'beta- to 2.28 MeV (Y-90): electron response of Si telescope and scintillator entrance windows'),
}


def peak_eff(g, sens, mat, E, n, rng):
    mc = MC(g, sens, rng)
    c = rng.uniform(-1, 1, n)
    f = rng.uniform(0, 2 * np.pi, n)
    s = np.sqrt(1 - c * c)
    D = np.stack([s * np.cos(f), s * np.sin(f), c], 1)
    mc.run(n, photons={'ev': np.arange(n), 'P': np.zeros((n, 3)), 'D': D, 'E': np.full(n, E)})
    Es = smear(mc.edep[:, 2:], mat, rng)
    w = fwhm(mat, E)
    return float(((Es > E - w) & (Es < E + w)).sum() / n)


def main(layout='BGO 3x3 quad (+-x,+-z)'):
    rng = np.random.default_rng(55)
    dk, kw, pairs = LAYOUTS[layout]
    mat = m7_geom.DETECTORS[dk][0]
    lines = [f"Calibration plan, layout: {layout}"]
    P = lines.append
    g, sens = m7_geom.build(L_um=100, det=dk, **kw)
    mcn = na22(g, sens, 40000, rng)
    eps_na = coinc_mask(smear(mcn.edep[:, 2:], mat, rng), mat, pairs).mean()
    effs = {}
    for E in (0.0595, 0.356, 0.6617, 1.1732, 1.3325, 1.4608, 2.6145):
        g, sens = m7_geom.build(L_um=100, det=dk, **kw)
        effs[E] = peak_eff(g, sens, mat, E, 40000, rng)
    P(f"22Na at foil: 511-511 coincidences per decay = {eps_na:.4f}")
    P("full-energy-peak efficiency (sum over detectors, source at foil): " + ", ".join(f"{E*1e3:.0f} keV {e:.4f}" for E, e in effs.items()))
    P("\nnuclide | exempt qty | activity used | key rate | time to 1e4 counts | calibrates")
    for nuc, (q, lines_, what) in SOURCES.items():
        A = q * UCI
        use = A
        if nuc == 'Na-22':
            use = 0.1 * UCI                     # keep singles < ~5 kcps in BGO
            rate = use * 0.903 * eps_na
            key = f"{rate:.0f} coinc/s"
        elif nuc == 'Am-241':
            om = (1 - 2.0 / np.sqrt(2.0 ** 2 + 1.2 ** 2)) / 2
            rate = use * om
            key = f"{rate:.1f} alpha/s in Si at 2 cm"
        elif nuc == 'Sr-90':
            om = (1 - 2.0 / np.sqrt(2.0 ** 2 + 1.2 ** 2)) / 2
            rate = use * 2 * om
            key = f"{rate:.0f} e-/s into Si at 2 cm"
        elif nuc == 'Ge-68':
            rate = use * 0.89 * eps_na / 0.903 * 0.9
            key = f"~{rate:.0f} coinc/s"
        else:
            E, y = lines_[0]
            e = effs.get(E, np.interp(E, sorted(effs), [effs[k] for k in sorted(effs)]))
            rate = use * y * e
            key = f"{rate:.0f} peak counts/s ({E*1e3:.0f} keV)"
        P(f"{nuc} | {q} uCi | {use/UCI:.2g} uCi = {use:.0f} Bq | {key} | {1e4/rate/60:.1f} min | {what}")
    P("\nNatural (no licence): thoriated W welding rods (2 % ThO2, 10 CFR 40.13(c)(1)(iii)) for 2614 keV; "
      "KCl (1 kg = 16 kBq 40K) for 1461 keV; cosmic muons (MIP peak ~9 MeV/cm in BGO) and delayed Michel "
      "electrons (52.8 MeV end-point) for the 5-60 MeV calorimeter scale.")
    for E, e in effs.items():
        if E == 2.6145:
            A_rod = 170.0                       # Bq Th-232 in a 150 mm x 2.4 mm 2 % thoriated rod (48 mg ThO2, 4.06 kBq/g Th)
            r = A_rod * 0.359 * e * 5          # 5 rods taped at the foil position
            P(f"5 thoriated rods at the foil (~{5*A_rod:.0f} Bq Th-232 in equilibrium): 2614 keV peak {r:.2f} /s -> 1e4 counts in {1e4/r/3600:.1f} h")
        if E == 1.4608:
            r = 16.3e3 * 0.1066 * e * 0.3       # 1 kg KCl bag around the chamber, ~30 % self/geometry factor
            P(f"1 kg KCl around the chamber: 1461 keV peak ~{r:.1f} /s")
    txt = '\n'.join(lines)
    open(os.path.join(OUT, 'm7_calibration.txt'), 'w').write(txt + '\n')
    print(txt)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'BGO 3x3 quad (+-x,+-z)')
