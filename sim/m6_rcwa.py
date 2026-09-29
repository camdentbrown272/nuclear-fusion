"""M6: minimal 1D-grating RCWA for TM (p) polarisation, with Li's factorisation rules.

Geometry: superstrate (incidence medium, z<0) / grating layers / uniform layers / semi-infinite
substrate. Fields e^{-i w t}; u = eta0*H_y; lengths normalised by k0 (z' = k0 z).

Layer equations (derived in docs/models/M6-surface-microgeometry.md, section 2.3):
   dU/dz' = i A^{-1} S_x,  dS_x/dz' = i (I - K E^{-1} K) U,
   d2U/dz'2 = A^{-1}(K E^{-1} K - I) U,
with E = Toeplitz[eps], A = Toeplitz[1/eps], K = diag(kx_n/k0) (Lalanne & Morris 1996;
Li 1996 inverse rule).  A stable top-down admittance recursion (S_x = Y U) gives the
specular/diffracted reflection; the substrate is absorbing, so absorptance = 1 - sum R_n.
"""
import numpy as np


def _toeplitz(coef, N):
    """Toeplitz matrix T[m,n] = c_{m-n} for harmonics -N..N from Fourier coefficients
    array indexed c[k + 2N] for k = -2N..2N."""
    idx = np.arange(-N, N + 1)
    return coef[(idx[:, None] - idx[None, :]) + 2 * N]


def fourier_step(eps_in, eps_out, f, N, center=0.5):
    """Fourier coefficients (k=-2N..2N) of a period-1 function equal to eps_in on an
    interval of fill fraction f centred at x=center, eps_out elsewhere."""
    k = np.arange(-2 * N, 2 * N + 1)
    c = (eps_in - eps_out) * f * np.sinc(k * f) * np.exp(-2j * np.pi * k * center)
    c = c.astype(complex)
    c[2 * N] += eps_out
    return c


def _kz(epsm, kx):
    """Normalised gamma for uniform medium: fields ~ exp(-gamma z'); decays/propagates to +z."""
    kz = np.sqrt(epsm - kx ** 2 + 0j)
    kz = np.where(kz.imag < 0, -kz, kz)
    return -1j * kz, kz


def _layer_modes(eps_coef_or_scalar, K, N, uniform):
    n = 2 * N + 1
    if uniform:
        epsm = eps_coef_or_scalar
        gam, kz = _kz(epsm, np.diag(K))
        W = np.eye(n, dtype=complex)
        V = (1j / epsm) * np.diag(gam)
        return W, V, gam
    ce, ci = eps_coef_or_scalar
    E = _toeplitz(ce, N)
    A = _toeplitz(ci, N)
    Ainv = np.linalg.inv(A)
    M = Ainv @ (K @ np.linalg.solve(E, K) - np.eye(n))
    lam, W = np.linalg.eig(M)
    gam = np.sqrt(lam + 0j)
    gam = np.where(gam.real < 0, -gam, gam)
    # purely oscillatory modes: choose the branch that propagates toward +z (Im gamma < 0)
    osc = np.abs(gam.real) < 1e-12 * np.maximum(1, np.abs(gam))
    gam = np.where(osc & (gam.imag > 0), -gam, gam)
    V = 1j * A @ W @ np.diag(gam)
    return W, V, gam


def rcwa_tm(lam_um, period_um, theta_deg, eps_sup, layers, eps_sub, N=30, want_field=False):
    """layers: list of dicts {'d': thickness_um, 'eps_in':, 'eps_out':, 'f': fill, 'center':}
    (grating slices) or {'d':, 'eps':} (uniform). Returns dict with R orders, A, fields.
    """
    k0 = 2 * np.pi / lam_um
    nsup = np.sqrt(eps_sup + 0j)
    kx0 = (nsup * np.sin(np.radians(theta_deg))).real
    orders = np.arange(-N, N + 1)
    kxn = kx0 + orders * lam_um / period_um
    K = np.diag(kxn).astype(complex)
    n = 2 * N + 1
    # substrate admittance
    Ws, Vs, _ = _layer_modes(eps_sub, K, N, True)
    Y = Vs @ np.linalg.inv(Ws)
    for L in reversed(layers):
        d = L["d"] * k0
        if "eps" in L:
            W, V, g = _layer_modes(L["eps"], K, N, True)
        else:
            ce = fourier_step(L["eps_in"], L["eps_out"], L["f"], N, L.get("center", 0.5))
            ci = fourier_step(1 / L["eps_in"], 1 / L["eps_out"], L["f"], N, L.get("center", 0.5))
            W, V, g = _layer_modes((ce, ci), K, N, False)
        X = np.diag(np.exp(-g * d))
        YW = Y @ W
        R = np.linalg.solve(V + YW, V - YW)
        XRX = X @ R @ X
        I = np.eye(n)
        Y = V @ (I - XRX) @ np.linalg.inv(W @ (I + XRX))
    gam0, kz0 = _kz(eps_sup, kxn)
    V0 = (1j / eps_sup) * np.diag(gam0)
    ainc = np.zeros(n, complex)
    ainc[N] = 1.0
    r = np.linalg.solve(V0 + Y, (V0 - Y) @ ainc)
    flux = (kz0 / eps_sup).real
    Rn = flux * np.abs(r) ** 2 / flux[N]
    Rn[flux <= 0] = 0.0
    out = {"R": Rn, "orders": orders, "A": 1 - Rn.sum(), "R0": Rn[N], "r": r, "kx": kxn}
    if want_field:
        # field just above the z=0 plane (crest tops) sampled across one period
        x = np.linspace(0, 1, 401)
        ph = np.exp(2j * np.pi * np.outer(x, orders) * 1.0) * np.exp(2j * np.pi * x[:, None] * kx0 * period_um / lam_um)
        U = ainc + r
        Sx = V0 @ (ainc - r) / 1.0
        # S_x = -i A dU/dz'; uniform medium: E_x amplitudes = S_x ; E_z = -K U / eps
        Ex = ph @ Sx
        Ez = ph @ (-(kxn * U) / eps_sup)
        E2 = np.abs(Ex) ** 2 + np.abs(Ez) ** 2
        Einc2 = 1.0 / abs(eps_sup)   # |E_inc|^2 = |u|^2/|n|^2
        out["x"] = x
        out["E2norm"] = E2 / Einc2
    return out


def sinusoid_layers(depth_um, eps_metal, eps_diel, nslice=20, under=None):
    """Staircase approximation to z_s(x) = (h/2)(1+cos 2pi x/L): metal below the surface.
    Crest (metal top) at x = 0.5 period, z=0; trough at z=h."""
    layers = []
    for j in range(nslice):
        zmid = (j + 0.5) / nslice
        # metal where cos(2 pi x) < 2 z/h - 1, centred at x = 0.5
        f = 1 - np.arccos(np.clip(2 * zmid - 1, -1, 1)) / np.pi
        layers.append({"d": depth_um / nslice, "eps_in": eps_metal, "eps_out": eps_diel, "f": f, "center": 0.5})
    if under:
        layers += under
    return layers
