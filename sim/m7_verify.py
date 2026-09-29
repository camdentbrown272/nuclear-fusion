"""M7 verification suite: analytic limits and reference comparisons for the
physics library and the MC. Output: docs/models/figs/m7_verify.txt"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from m7_common import MATERIALS, ME, self_test  # noqa: E402
from m7_mc import MC, Geometry, Tube, Box  # noqa: E402
from m7_ipc import sample_pairs  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), '..', 'docs', 'models', 'figs')


def iso(n, rng):
    c = rng.uniform(-1, 1, n)
    f = rng.uniform(0, 2 * np.pi, n)
    s = np.sqrt(1 - c * c)
    return np.stack([s * np.cos(f), s * np.sin(f), c], 1)


def main():
    lines = []
    P = lines.append
    self_test(P)
    rng = np.random.default_rng(99)

    # 1. NaI 3"x3" 662 keV point source at 10 cm: MC vs analytic ray-trace
    P("\n[V1] 3\"x3\" NaI, 662 keV point source 10 cm from the face")
    g = Geometry([30, 30, 30], 'VAC')
    T = Tube(2, [0, 0, 10 + 3.81], 0, 3.81, 3.81)
    g.add(T, 'NAI', 'nai')
    n = 200000
    D = iso(n, rng)
    mc = MC(g, ['nai'], rng)
    mc.run(n, photons={'ev': np.arange(n), 'P': np.zeros((n, 3)), 'D': D, 'E': np.full(n, 0.662)})
    e = mc.edep[:, 0]
    t1 = T.ray(np.zeros((n, 3)), D)
    hit = np.isfinite(t1)
    t2 = np.where(hit, T.ray(D * (np.where(hit, t1, 0) + 1e-7)[:, None], D), 0)
    mu = MATERIALS['NAI'].mu_total(0.662)[0] * MATERIALS['NAI'].rho
    P(f"     interaction probability: MC {(e>0.005).mean():.4f} vs analytic {(1-np.exp(-mu*t2)).mean():.4f}; "
      f"peak/total {(e>0.655).sum()/(e>0.005).sum():.3f} (Heath NaI catalogue ~0.5-0.55 at 10 cm [BK])")

    # 2. electron depth-dose in water: R50 and practical range
    P("\n[V2] broad-beam electrons in water: R50 and Rp vs clinical relations (AAPM TG-25: E0 = 2.33 R50; Rp = 0.52 E0 - 0.3) [BK]")
    for E0 in (10.0, 20.0):
        g = Geometry([60, 60, 60], 'VAC')
        dz = 0.25
        names = []
        for i in range(60):
            g.add(Box([0, 0, (i + 0.5) * dz], [50, 50, dz / 2]), 'H2O', f'L{i}')
            names.append(f'L{i}')
        n = 3000
        mc = MC(g, names, np.random.default_rng(3))
        mc.run(n, charged={'ev': np.arange(n), 'q': -np.ones(n, int), 'P': np.tile([0, 0, 1e-5], (n, 1)),
                           'D': np.tile([0, 0, 1.0], (n, 1)), 'T': np.full(n, E0)})
        dose = mc.edep.sum(0) / n / dz
        z = (np.arange(60) + 0.5) * dz
        imax = np.argmax(dose)
        R50 = z[imax:][np.argmax(dose[imax:] < 0.5 * dose.max())]
        d = np.gradient(dose, z)
        i = np.argmin(d)
        Rp = z[i] - dose[i] / d[i]
        P(f"     E0 = {E0:4.1f} MeV: R50 MC {R50:.2f} cm vs {E0/2.33:.2f}; Rp MC {Rp:.2f} vs {0.52*E0-0.3:.2f} cm; CSDA {MATERIALS['H2O'].csda_range(E0)[0]/0.997:.2f} cm")

    # 3. Highland multiple scattering through 300 um Pd
    P("\n[V3] 10 MeV e- through 300 um Pd: projected theta0")
    L = 0.03
    g = Geometry([5, 5, L / 2 + 1e-4], 'VAC')
    g.add(Box([0, 0, 0], [5, 5, L / 2]), 'PD', 'foil')
    n = 20000
    mc = MC(g, ['foil'], np.random.default_rng(4))
    mc.run(n, charged={'ev': np.arange(n), 'q': -np.ones(n, int), 'P': np.tile([0, 0, -L / 2 + 1e-6], (n, 1)),
                       'D': np.tile([0, 0, 1.0], (n, 1)), 'T': np.full(n, 10.0)})
    q = np.array(mc.esc_q)
    fw = q[:, 3] > 0
    th = np.arccos(np.clip(q[fw, 3], -1, 1))
    x = L * 12.02 / 9.20
    p = np.sqrt(10 * 10 + 2 * ME * 10)
    hl = 13.6 / (p / (10 + ME) * p) * np.sqrt(x) * (1 + 0.038 * np.log(x))
    th0 = np.sqrt(np.mean(th[th < 3 * 1.41 * hl] ** 2) / 2)
    P(f"     MC {th0*1e3:.0f} mrad vs Highland {hl*1e3:.0f} mrad (PDG accuracy +-11 %); mean exit energy {q[fw,2].mean():.2f} MeV")

    # 4. energy conservation for IPC pairs in a closed BGO block
    P("\n[V4] energy conservation: IPC pairs at the centre of a 80 cm BGO cube")
    g = Geometry([45, 45, 45], 'VAC')
    g.add(Box([0, 0, 0], [40, 40, 40]), 'BGO', 'bgo')
    n = 1000
    Tp, Tm, dp, dm = sample_pairs(n, rng)
    ev = np.arange(n)
    mc = MC(g, ['bgo'], np.random.default_rng(5))
    mc.run(n, charged={'ev': np.concatenate([ev, ev]), 'q': np.concatenate([np.ones(n, int), -np.ones(n, int)]),
                       'P': np.zeros((2 * n, 3)), 'D': np.concatenate([dp, dm]), 'T': np.concatenate([Tp, Tm])})
    esc = np.array(mc.esc_g).reshape(-1, 3)[:, 1].sum() + (np.array(mc.esc_q)[:, 2].sum() if len(mc.esc_q) else 0)
    tot = mc.edep.sum() + esc
    exp = (Tp + Tm).sum() + 2 * ME * n
    P(f"     deposited+escaped {tot/n:.4f} MeV/event vs expected {exp/n:.4f} (T+ + T- + 2 m_e); ratio {tot/exp:.5f}")
    P(f"     containment: <E_dep> = {mc.edep.mean():.3f} MeV = {mc.edep.sum()/exp*100:.2f} % (23.85 MeV line in a hermetic calorimeter)")

    # 5. positron annihilation in flight in water (fraction of positrons slowing from T0)
    P("\n[V5] in-flight annihilation probability, e+ fully stopped in water")
    for T0 in (1.0, 10.0):
        g = Geometry([60, 60, 60], 'VAC')
        g.add(Box([0, 0, 0], [55, 55, 55]), 'H2O', 'w')
        n = 4000
        mc = MC(g, [], np.random.default_rng(6))
        mc.run(n, charged={'ev': np.arange(n), 'q': np.ones(n, int), 'P': np.zeros((n, 3)), 'D': iso(n, rng), 'T': np.full(n, T0)})
        a = np.array(mc.ann)
        prim = a[a[:, 0] < n]
        # analytic: P = 1 - exp(-int n_e sigma dx) over the CSDA path
        m = MATERIALS['H2O']
        Tg = np.geomspace(0.02, T0, 400)
        integ = np.trapezoid(m.ann_inflight(Tg) / (m.S_col(Tg, True) + m.S_rad(Tg)), Tg)
        P(f"     T0 = {T0:4.1f} MeV: MC {prim[:,2].sum()/n:.4f} vs CSDA integral {1-np.exp(-integ):.4f}")
    txt = '\n'.join(lines)
    open(os.path.join(OUT, 'm7_verify.txt'), 'w').write(txt + '\n')
    print(txt)


if __name__ == '__main__':
    main()
