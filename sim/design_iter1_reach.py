"""T1 reach for iteration-1 rev B (pre-registration v1 §9).

Median-5sigma D+D fusion rate per membrane for the pooled T1 test and for a single membrane.
alpha = 1.15e-7 (T1 weight 0.4 of global 2.87e-7; ADR-006), >=10 net events floor, P2 = 42 d x 0.8 live.
eps = PID efficiency (exit 0.22, entry 0.162 at 12 um; ADR-005 §2.2) x grid 0.80 x septum 0.9;
proton branch 0.5. Background per cell per day includes electrolyte recoils and septum/grid (n,xp).
"""
import numpy as np
from scipy.stats import poisson

ALPHA = 1.15e-7
T_LIVE = 42 * 0.8 * 86400.0
EPS = {"exit": 0.22 * 0.8 * 0.9, "entry": 0.162 * 0.8 * 0.9}


def ncrit(B):
    n = 0
    while poisson.sf(n - 1, B) > ALPHA:
        n += 1
    return n


def s_median(B):
    nc = max(ncrit(B), int(np.ceil(B)) + 10)
    s = 0.0
    while poisson.median(B + s) < nc:
        s += 0.05
    return s, nc


def reach(bkg_per_day, eps, n_mem):
    B = n_mem * bkg_per_day * T_LIVE / 86400.0
    s, nc = s_median(B)
    return 2 * s / (n_mem * eps * T_LIVE), B, nc


if __name__ == "__main__":
    for bkg in (0.03, 0.05, 0.09):
        for face, eps in EPS.items():
            r1, B1, n1 = reach(bkg, eps, 1)
            r8, B8, n8 = reach(bkg, eps, 8)
            print(f"bkg {bkg:.2f}/d {face:5s}: single {r1:.1e} (B={B1:.1f}, ncrit={n1}); "
                  f"pooled-8 per membrane {r8:.1e} (B={B8:.1f}, ncrit={n8}) fusions/s")
