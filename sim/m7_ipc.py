"""M7 / Q1: kinematics of E0 internal pair creation (IPC) for the proposed
0+ (d+d threshold) -> 0+ (g.s.) transition in 4He.

Transition energy
    Q = 2 m_d c^2 - m_alpha c^2 = 23.8465 MeV   (AME2020 masses [BK];
    m_d = 1875.612942 MeV, m_alpha = 3727.379378 MeV, CODATA 2018)
A 0+ -> 0+ transition cannot emit a single photon. It proceeds by internal
conversion (impossible for bare/He nuclei in a metal to any useful degree)
or internal pair creation, so the full energy appears as an e+e- pair:
    T+ + T- = Q - 2 m_e = 22.824 MeV   (recoil <= 0.08 MeV neglected).

Distribution (first Born approximation; Oppenheimer & Schwinger, Phys. Rev.
56, 1066 (1939); Thomas, Phys. Rev. 58, 714 (1940); Church & Weneser,
Phys. Rev. 103, 1035 (1956); review: Schlueter, Soff & Greiner, Phys. Rep. 75,
327 (1981); Wilkinson, Nucl. Phys. A133, 1 (1969)).
In units m_e = c = 1, with total energies E+, E- and momenta p+, p-,
E+ + E- = W0 = Q/m_e, theta = e+e- opening angle:

    d^2 W / (dE+ dcos(theta)) ∝ p+ p- (E+ E- - 1 + p+ p- cos(theta))      (1)

Derivation check: the E0 amplitude is (nuclear monopole matrix element) x
(time component of the lepton current), |u_bar gamma^0 v|^2 summed over spins
= 4 (E+E- + p+.p- - m^2); phase space d^3p+ d^3p-/(E+E-) delta(...) gives
p+ p- dE+ dOmega+ dOmega-. Integrating over cos(theta):

    dW/dE+ ∝ p+ p- (E+ E- - 1)                                          (2)

and the angular correlation for fixed energies is
    W(theta) ∝ 1 + [p+p-/(E+E- - 1)] cos(theta)  ≈ 1 + beta+ beta- cos(theta).
Coulomb (Fermi-function) corrections for Z = 2 change the spectrum by <~1 %
(checked below) and are neglected in the transport.
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from m7_common import ME, ALPHA  # noqa: E402

M_D = 1875.612942
M_ALPHA = 3727.379378
Q_DD = 2 * M_D - M_ALPHA            # 23.8465 MeV
W0 = Q_DD / ME                       # total pair energy in m_e units
TSUM = Q_DD - 2 * ME                 # kinetic sum 22.824 MeV


def marginal_Eplus(Ep, W=W0):
    """dW/dE+ (unnormalised), E+ total energy in m_e units."""
    Em = W - Ep
    pp = np.sqrt(np.clip(Ep ** 2 - 1, 0, None))
    pm = np.sqrt(np.clip(Em ** 2 - 1, 0, None))
    return pp * pm * (Ep * Em - 1)


def fermi(Z, E):
    """Non-relativistic-form Fermi function 2 pi eta / (1 - exp(-2 pi eta)),
    eta = Z alpha / beta (Z>0 electron, Z<0 positron)."""
    p = np.sqrt(np.clip(E ** 2 - 1, 1e-12, None))
    eta = Z * ALPHA * E / p
    x = np.clip(2 * np.pi * eta, -50, 50)
    return -x / np.expm1(-x)


def sample_pairs(n, rng, W=W0):
    """Return T+ [MeV], T- [MeV], unit vectors d+ (n,3), d- (n,3)."""
    # E+ by rejection from eq. (2)
    grid = np.linspace(1.0, W - 1.0, 4001)
    fmax = marginal_Eplus(grid, W).max() * 1.001
    Ep = np.empty(0)
    while Ep.size < n:
        x = rng.uniform(1.0, W - 1.0, 2 * n)
        acc = rng.uniform(0, fmax, 2 * n) < marginal_Eplus(x, W)
        Ep = np.concatenate([Ep, x[acc]])
    Ep = Ep[:n]
    Em = W - Ep
    pp = np.sqrt(Ep ** 2 - 1)
    pm = np.sqrt(Em ** 2 - 1)
    a = Ep * Em - 1.0
    b = pp * pm
    # cos(theta) from density (a + b c) on [-1,1]: inverse CDF of linear pdf
    u = rng.uniform(size=n)
    # CDF: [a(c+1) + b(c^2-1)/2] / (2a)  -> solve b/2 c^2 + a c + (a - b/2 - 2 a u) = 0
    A = b / 2
    B = a
    C = a - b / 2 - 2 * a * u
    disc = np.sqrt(np.clip(B ** 2 - 4 * A * C, 0, None))
    c = np.where(A > 1e-12, (-B + disc) / (2 * A), 2 * u - 1)
    c = np.clip(c, -1, 1)
    # isotropic e+ direction
    dp = isotropic(n, rng)
    dm = rotate(dp, c, rng.uniform(0, 2 * np.pi, n))
    return (Ep - 1) * ME, (Em - 1) * ME, dp, dm


def isotropic(n, rng):
    c = rng.uniform(-1, 1, n)
    ph = rng.uniform(0, 2 * np.pi, n)
    s = np.sqrt(1 - c ** 2)
    return np.stack([s * np.cos(ph), s * np.sin(ph), c], axis=1)


def rotate(d, cth, phi):
    """Rotate unit vectors d by polar angle acos(cth), azimuth phi."""
    sth = np.sqrt(np.clip(1 - cth ** 2, 0, None))
    ux, uy, uz = d[:, 0], d[:, 1], d[:, 2]
    perp = np.sqrt(np.clip(1 - uz ** 2, 0, None))
    small = perp < 1e-8
    perp_s = np.where(small, 1.0, perp)
    nx = ux * cth + sth * (ux * uz * np.cos(phi) - uy * np.sin(phi)) / perp_s
    ny = uy * cth + sth * (uy * uz * np.cos(phi) + ux * np.sin(phi)) / perp_s
    nz = uz * cth - perp * sth * np.cos(phi)
    sg = np.sign(uz + 1e-300)
    nx = np.where(small, sth * np.cos(phi), nx)
    ny = np.where(small, sth * np.sin(phi), ny)
    nz = np.where(small, sg * cth, nz)
    out = np.stack([nx, ny, nz], axis=1)
    return out / np.linalg.norm(out, axis=1)[:, None]


def main(outdir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rng = np.random.default_rng(7001)
    lines = []
    P = lines.append
    P(f"Q(d+d -> 4He g.s.) = {Q_DD:.4f} MeV ; kinetic sum T+ + T- = {TSUM:.4f} MeV ; W0 = {W0:.3f} m_e")
    recoil_max = (np.sqrt((Q_DD) ** 2) ** 2) / (2 * M_ALPHA)
    P(f"max 4He recoil energy (pair momenta parallel) = {recoil_max*1e3:.0f} keV (neglected)")

    # analytic marginal
    Ep = np.linspace(1, W0 - 1, 20001)
    f = marginal_Eplus(Ep)
    f /= np.trapezoid(f, (Ep - 1) * ME)
    T = (Ep - 1) * ME
    cdf = np.concatenate([[0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(T))])
    for thr in [0.1, 0.3, 0.5, 1.0, 2.0, 5.0]:
        P(f"P(T+ < {thr:4.1f} MeV) = {np.interp(thr, T, cdf):.4f}")
    # <cos theta> analytic: int (b/3)*2 / int 2a  weighted by marginal-norm
    pp = np.sqrt(Ep ** 2 - 1)
    pm = np.sqrt((W0 - Ep) ** 2 - 1)
    a = Ep * (W0 - Ep) - 1
    b = pp * pm
    num = np.trapezoid(pp * pm * (2 * b / 3), Ep)
    den = np.trapezoid(pp * pm * 2 * a, Ep)
    P(f"<cos theta_+-> analytic = {num/den:.4f}  (ultra-relativistic limit 1/3)")
    # Fermi-function effect
    fc = f * fermi(-2, Ep) * fermi(2, W0 - Ep)
    fc /= np.trapezoid(fc, T)
    P(f"Coulomb (Z=2 Fermi fn) change of <T+>: {np.trapezoid(T*fc,T)-np.trapezoid(T*f,T):+.3f} MeV "
      f"(of {np.trapezoid(T*f,T):.3f}); max local shape change {np.max(np.abs(fc[(T>0.5)&(T<TSUM-0.5)]/f[(T>0.5)&(T<TSUM-0.5)]-1))*100:.2f} % for 0.5<T+<{TSUM-0.5:.1f} MeV")

    # MC sample
    n = 400000
    Tp, Tm, dp, dm = sample_pairs(n, rng)
    cth = np.sum(dp * dm, axis=1)
    th = np.degrees(np.arccos(np.clip(cth, -1, 1)))
    P(f"MC n={n}: <T+> = {Tp.mean():.3f} MeV, <T-> = {Tm.mean():.3f} MeV, sum check max|T+ + T- - Tsum| = {np.max(np.abs(Tp+Tm-TSUM)):.2e}")
    P(f"MC <cos theta> = {cth.mean():.4f} ; median opening angle = {np.median(th):.1f} deg ; P(theta>90) = {(th>90).mean():.3f}; P(theta>150) = {(th>150).mean():.3f}; P(theta<30) = {(th<30).mean():.3f}")
    minv = np.sqrt(2 * ME ** 2 + 2 * ((Tp + ME) * (Tm + ME) - np.sqrt(Tp * (Tp + 2 * ME)) * np.sqrt(Tm * (Tm + 2 * ME)) * cth))
    P(f"e+e- invariant mass: mean {minv.mean():.2f} MeV, 5-95 % range {np.percentile(minv,5):.2f}-{np.percentile(minv,95):.2f} MeV (continuum, no peak)")
    asym = np.abs(Tp - Tm) / (Tp + Tm)
    P(f"energy asymmetry |T+-T-|/sum: mean {asym.mean():.3f}; P(both > 5 MeV) = {((Tp>5)&(Tm>5)).mean():.3f}; P(min(T+,T-) < 1 MeV) = {(np.minimum(Tp,Tm)<1).mean():.3f}")

    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    ax[0].hist(Tp, bins=115, range=(0, 23), density=True, histtype='step', lw=1.5, label='MC e$^+$')
    ax[0].plot(T, f, 'k--', lw=1, label='eq. (2) analytic')
    ax[0].plot(T, fc, 'r:', lw=1, label='with Z=2 Fermi fn.')
    ax[0].set_xlabel('positron kinetic energy T$_+$ [MeV]')
    ax[0].set_ylabel('dN/dT$_+$ [1/MeV]')
    ax[0].set_title('E0 IPC energy sharing (T$_+$+T$_-$=22.82 MeV)')
    ax[0].legend(fontsize=8)
    ax[1].hist(th, bins=90, range=(0, 180), density=True, histtype='step', lw=1.5, label='E0, all energies (MC)')
    tt = np.linspace(0, 180, 200)
    ur = (1 + np.cos(np.radians(tt))) * np.sin(np.radians(tt))
    ur /= np.trapezoid(ur, tt)
    ax[1].plot(tt, ur, 'k--', lw=1, label=r'$(1+\cos\theta)\sin\theta$ (UR limit)')
    ax[1].set_xlabel(r'e$^+$e$^-$ opening angle $\theta$ [deg]')
    ax[1].set_ylabel(r'dN/d$\theta$ [1/deg]')
    ax[1].set_title('Angular correlation')
    ax[1].legend(fontsize=8)
    h = ax[2].hist2d(Tp, th, bins=[46, 36], range=[[0, 23], [0, 180]], cmap='viridis')
    ax[2].set_xlabel('T$_+$ [MeV]')
    ax[2].set_ylabel(r'$\theta$ [deg]')
    ax[2].set_title('joint distribution (MC)')
    fig.colorbar(h[3], ax=ax[2])
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'm7_ipc_kinematics.png'), dpi=130)
    plt.close(fig)
    txt = '\n'.join(lines)
    open(os.path.join(outdir, 'm7_ipc_kinematics.txt'), 'w').write(txt + '\n')
    print(txt)


if __name__ == '__main__':
    out = os.path.join(os.path.dirname(__file__), '..', 'docs', 'models', 'figs')
    os.makedirs(out, exist_ok=True)
    main(out)
