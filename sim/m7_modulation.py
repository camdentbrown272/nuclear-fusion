"""M7 / Q4: flux-modulation (lock-in / on-off) protocol.

1. Diffusion time of D through Pd foils of 25-300 um; permeation transfer
   function (entry-face concentration modulated, exit face at c = 0 into UHV):
       H(w) = J_out(w)/J_out(0) = q / sinh(q),  q = sqrt(i w L^2 / D)
   (Crank, The Mathematics of Diffusion, 2nd ed. 1975, ch. 4 [BK]).
2. Efficiency of a square-wave (50 % duty) flux modulation analysed with a
   template regression, eta(P) = Var[r(t)] / 0.25, where r(t) is the exit flux
   normalised to its steady-state ON value.
3. Noise budget vs period with slow background drifts (barometric pressure,
   temperature) -> recommended period window.
4. End-to-end toy study: 30 d of 10-min bins with Poisson counts, barometric
   drift, a temperature-coupled gain artifact synchronous with the drive, and
   an injected signal; compares (i) naive on/off, (ii) regression with
   covariates, (iii) 2x2 factorial (D2O/H2O x on/off) with current steering.
5. Look-elsewhere control for the channel combination.
Outputs: docs/models/figs/m7_modulation.txt, m7_modulation.png, m7_toy.png
"""
import os
import sys
import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
OUT = os.path.join(os.path.dirname(__file__), '..', 'docs', 'models', 'figs')
KB = 8.617333e-5

# D in Pd, alpha phase (Voelkl & Alefeld, in "Hydrogen in Metals I", Springer
# 1978 [BK]); H for comparison. D0 [cm2/s], Ea [eV]
DIFF = {'D': (1.7e-3, 0.206), 'H': (2.9e-3, 0.230)}


def Dcoef(iso='D', TC=25.0):
    D0, Ea = DIFF[iso]
    return D0 * np.exp(-Ea / (KB * (TC + 273.15)))


def H_perm(f, L, D):
    """Exit-flux transfer function for a slab driven at the entry face."""
    w = 2 * np.pi * np.asarray(f, float)
    q = np.sqrt(1j * w * L ** 2 / D)
    with np.errstate(over='ignore', invalid='ignore'):
        h = np.where(np.abs(q) < 1e-6, 1.0 + 0j, q / np.sinh(q))
    return np.where(np.isfinite(h), h, 0)


def step_response(t, L, D, nterm=200):
    """Exit flux after a step in entry concentration (J/J_inf)."""
    t = np.asarray(t, float)
    n = np.arange(1, nterm + 1)
    s = 1 + 2 * np.sum(((-1.0) ** n)[None, :] * np.exp(-(n[None, :] * np.pi) ** 2 * D * t[:, None] / L ** 2), axis=1)
    return np.clip(s, 0, 1)


def square_response(P, L, D, extra_tau=0.0, nper=6, npts=4000):
    """Periodic steady-state exit flux for a 50 % square-wave entry drive,
    via Fourier series of the drive times H(f) (and an optional first-order
    surface/electrochemical lag extra_tau). Returns t, r(t) over one period."""
    t = np.linspace(0, P, npts, endpoint=False)
    r = np.full(npts, 0.5)
    for k in range(1, 400, 2):   # odd harmonics of a 0/1 square wave: (2/(pi k)) sin
        f = k / P
        h = H_perm(f, L, D) / (1 + 1j * 2 * np.pi * f * extra_tau)
        r += (2 / (np.pi * k)) * np.imag(h * np.exp(1j * 2 * np.pi * f * t))
    return t, r


def eta(P, L, D, extra_tau=0.0):
    t, r = square_response(P, L, D, extra_tau)
    return np.var(r) / 0.25


def toy(seed=3, days=30, P_h=2.0, S_rate=0.0, B_rate=2000.0, art=0.05, steer=True, factorial=True):
    """Toy experiment in one counting channel (rates per day). Returns z-scores
    for three analyses. art = fractional background change when the cell heats
    (no current steering) -> gain shift -> window-rate change."""
    rng = np.random.default_rng(seed)
    dt = 10 / 1440.0                        # 10-min bins [d]
    n = int(days / dt)
    t = np.arange(n) * dt
    # pseudo-random block schedule: each block of 2 periods randomly ON-OFF or OFF-ON
    per = P_h / 24.0
    nb = int(np.ceil(days / (2 * per)))
    order = rng.integers(0, 2, nb)
    phase = ((t / per).astype(int)) % 2
    blk = (t / (2 * per)).astype(int)
    on = np.where(order[blk] == 0, phase == 0, phase == 1).astype(float)
    # isotope schedule for the factorial design: alternate D2O / H2O every 5 days
    iso_D = ((t // 5) % 2 == 0).astype(float) if factorial else np.ones(n)
    # barometric drift (Ornstein-Uhlenbeck, sigma 8 hPa, corr 2 d) with -0.7 %/hPa (hadronic)
    p = np.zeros(n)
    a = np.exp(-dt / 2.0)
    for i in range(1, n):
        p[i] = a * p[i - 1] + np.sqrt(1 - a * a) * 8.0 * rng.normal()
    baro = 1 - 0.007 * p
    heat = on if not steer else np.zeros(n)          # with steering the heat load is constant
    lam = B_rate * dt * baro * (1 + art * heat) + S_rate * dt * on * iso_D
    k = rng.poisson(lam)
    # (i) naive on/off (D2O periods only)
    selD = iso_D > 0
    non, noff = k[selD & (on > 0)].sum(), k[selD & (on == 0)].sum()
    ton, toff = (selD & (on > 0)).sum(), (selD & (on == 0)).sum()
    z_naive = (non / ton - noff / toff) / np.sqrt(non / ton ** 2 + noff / toff ** 2)
    # (ii) Poisson GLM with covariates (pressure measured, heat-proxy = cell power measured)
    cols = [np.ones(n), p] + ([heat] if heat.any() else []) + [on * iso_D]
    X = np.column_stack(cols)
    beta = glm_poisson(k, X, offset=np.log(dt))
    cov = np.linalg.inv(X.T @ (np.exp(X @ beta + np.log(dt))[:, None] * X))
    z_glm = beta[-1] / np.sqrt(cov[-1, -1])
    # (iii) factorial: interaction (D on - D off) - (H on - H off), pressure covariate
    if factorial:
        X3 = np.column_stack([np.ones(n), p, on, iso_D, on * iso_D])
        b3 = glm_poisson(k, X3, offset=np.log(dt))
        c3 = np.linalg.inv(X3.T @ (np.exp(X3 @ b3 + np.log(dt))[:, None] * X3))
        z_fac = b3[4] / np.sqrt(c3[4, 4])
    else:
        z_fac = np.nan
    return z_naive, z_glm, z_fac


def glm_poisson(k, X, offset, it=50):
    beta = np.zeros(X.shape[1])
    beta[0] = np.log(max(k.mean(), 1e-9)) - offset if np.isscalar(offset) else np.log(max(k.mean(), 1e-9)) - offset.mean()
    for _ in range(it):
        eta_ = X @ beta + offset
        mu = np.exp(eta_)
        g = X.T @ (k - mu)
        Hm = X.T @ (mu[:, None] * X)
        step = np.linalg.solve(Hm, g)
        beta += step
        if np.max(np.abs(step)) < 1e-9:
            break
    return beta


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    lines = []
    P = lines.append
    P("M7 Q4 flux-modulation protocol")
    for TC in (25.0, 60.0):
        P(f"D(D in Pd, alpha phase, {TC:.0f} C) = {Dcoef('D', TC):.2e} cm2/s ; D(H) = {Dcoef('H', TC):.2e} cm2/s ; ratio D/H = {Dcoef('D', TC)/Dcoef('H', TC):.2f}")
    Ls = [25, 50, 100, 200, 300]
    P("\nL[um] | t_lag=L2/6D [s] | tau1=L2/(pi2 D) [s] | t90 step [s] | f_3dB [mHz] | P for eta>=0.9 [min] (D x1 / x0.3 / x3)")
    fs = np.logspace(-6, 0, 3000)
    res = {}
    for L_um in Ls:
        L = L_um * 1e-4
        D = Dcoef('D', 25.0)
        tt = np.logspace(-2, 6, 4000)
        sr = step_response(tt, L, D)
        t90 = tt[np.argmax(sr >= 0.9)]
        h = np.abs(H_perm(fs, L, D))
        f3 = fs[np.argmax(h < 1 / np.sqrt(2))]
        Pmin = []
        for fac in (1.0, 0.3, 3.0):
            Pg = np.logspace(0, 5, 120)
            et = np.array([eta(p, L, D * fac) for p in Pg])
            Pmin.append(Pg[np.argmax(et >= 0.9)] / 60)
        res[L_um] = (L ** 2 / (6 * D), L ** 2 / (np.pi ** 2 * D), t90, f3, Pmin)
        P(f"{L_um:5d} | {L**2/(6*D):9.1f} | {L**2/(np.pi**2*D):9.1f} | {t90:9.1f} | {f3*1e3:8.2f} | "
          f"{Pmin[0]:6.1f} / {Pmin[1]:6.1f} / {Pmin[2]:6.1f}")
    P("\nWith an additional first-order electrochemical/surface lag tau_s (to be measured by the RGA step test):")
    for ts in (60.0, 600.0, 1800.0):
        Pg = np.logspace(1, 5.5, 150)
        et = np.array([eta(p, 100e-4, Dcoef('D'), ts) for p in Pg])
        P(f"  tau_s = {ts/60:5.1f} min, L = 100 um: P(eta>=0.9) = {Pg[np.argmax(et>=0.9)]/60:6.1f} min; eta(P=2h) = {eta(7200,100e-4,Dcoef('D'),ts):.3f}; eta(P=6h) = {eta(21600,100e-4,Dcoef('D'),ts):.3f}")

    # noise vs period: drift spectrum of barometric effect on a hadron-dominated channel
    P("\nDrift budget (lock-in at frequency 1/P, background B, 30 d):")
    for B in (5.0, 50.0, 1000.0):     # counts/day
        for Ph in (0.5, 2.0, 6.0, 24.0):
            # OU pressure noise sigma 8 hPa, corr time 2 d; coefficient -0.7 %/hPa on the whole background
            per = Ph / 24.0
            f = 1 / per
            tc = 2.0
            Sp = 2 * 8.0 ** 2 * tc / (1 + (2 * np.pi * f * tc) ** 2)      # hPa^2 / (1/d)
            # square-wave demodulation of drift: variance of (on - off) mean rate difference over T=30 d
            var_drift = (0.007 * B) ** 2 * Sp * (8 / np.pi ** 2) / 30.0 * 2
            var_pois = 4 * B / 30.0
            P(f"  B = {B:6.0f}/d, P = {Ph:4.1f} h: sigma(on-off rate) Poisson {np.sqrt(var_pois):.3f}/d, uncorrected barometric drift {np.sqrt(var_drift):.3f}/d")

    # toy studies
    P("\nToy experiments (30 d, 10-min bins, P = 2 h randomized blocks, D2O/H2O alternating every 5 d):")
    cases = [('no signal, no steering (5 % heat artifact)', dict(S_rate=0, steer=False)),
             ('no signal, current steering', dict(S_rate=0, steer=True)),
             ('signal 100/d on B 2000/d, steering', dict(S_rate=100, steer=True)),
             ('signal 100/d, no steering', dict(S_rate=100, steer=False))]
    toyres = {}
    for lab, kw in cases:
        zs = np.array([toy(seed=s, **kw) for s in range(60)])
        toyres[lab] = zs
        P(f"  {lab:48s}: mean z  naive {np.nanmean(zs[:,0]):+.2f} (sd {np.nanstd(zs[:,0]):.2f}) | GLM w/ covariates {np.nanmean(zs[:,1]):+.2f} "
          f"(sd {np.nanstd(zs[:,1]):.2f}) | factorial interaction {np.nanmean(zs[:,2]):+.2f} (sd {np.nanstd(zs[:,2]):.2f})")

    # look-elsewhere
    P("\nLook-elsewhere: one-sided global 5 sigma (p = 2.87e-7) with Bonferroni over N pre-registered tests")
    for N in (1, 2, 4, 8, 12, 24):
        P(f"  N = {N:2d}: local threshold z = {stats.norm.isf(2.87e-7 / N):.2f}")

    # figures
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    for L_um in Ls:
        L = L_um * 1e-4
        ax[0].semilogx(fs * 3600, np.abs(H_perm(fs, L, Dcoef('D'))), label=f'{L_um} um')
    ax[0].set_xlabel('modulation frequency [1/h]')
    ax[0].set_ylabel('|J_out(f)/J_out(0)|')
    ax[0].set_title('permeation transfer function (D in Pd, 25 C)')
    ax[0].legend(fontsize=8)
    ax[0].set_xlim(1e-2, 1e3)
    Pg = np.logspace(0.5, 5, 150)
    for L_um in Ls:
        ax[1].semilogx(Pg / 60, [eta(p, L_um * 1e-4, Dcoef('D')) for p in Pg], label=f'{L_um} um')
    ax[1].semilogx(Pg / 60, [eta(p, 100e-4, Dcoef('D'), 1800) for p in Pg], 'k--', label='100 um + 30 min surface lag')
    ax[1].axhline(0.9, color='gray', ls=':')
    ax[1].set_xlabel('square-wave period P [min]')
    ax[1].set_ylabel('efficiency eta = Var r / 0.25')
    ax[1].set_title('lock-in efficiency vs period')
    ax[1].legend(fontsize=7)
    t, r = square_response(7200, 300e-4, Dcoef('D'))
    t2, r2 = square_response(7200, 100e-4, Dcoef('D'), 1800)
    ax[2].plot(t / 60, (t % 7200 < 3600).astype(float), 'k:', label='drive (current to foil)')
    ax[2].plot(t / 60, r, label='exit flux, 300 um')
    ax[2].plot(t2 / 60, r2, label='exit flux, 100 um + 30 min lag')
    ax[2].set_xlabel('time in period [min]')
    ax[2].set_ylabel('normalised flux')
    ax[2].set_title('P = 2 h, periodic steady state')
    ax[2].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'm7_modulation.png'), dpi=130)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 4.3))
    labs = list(toyres)
    for j, nm in enumerate(['naive on/off', 'GLM + covariates', 'factorial D/H x on/off']):
        ax.errorbar(np.arange(len(labs)) + (j - 1) * 0.2, [np.nanmean(toyres[l][:, j]) for l in labs],
                    yerr=[np.nanstd(toyres[l][:, j]) for l in labs], fmt='o', label=nm)
    ax.axhline(0, color='k', lw=0.5)
    ax.set_xticks(range(len(labs)))
    ax.set_xticklabels([l.replace(', ', '\n') for l in labs], fontsize=7)
    ax.set_ylabel('z-score (mean +- sd over 60 toys)')
    ax.set_title('Artifact rejection: 5 % heat-coupled rate artifact, B = 2000/d')
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'm7_toy.png'), dpi=130)
    plt.close(fig)
    txt = '\n'.join(lines)
    open(os.path.join(OUT, 'm7_modulation.txt'), 'w').write(txt + '\n')
    print(txt)


if __name__ == '__main__':
    main()
