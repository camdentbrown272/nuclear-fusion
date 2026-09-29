"""M6 Q1-Q2: roughness power spectral density (PSD) targets.

 * Builds the recommended 'ring' PSD (isotropic, band-limited around the SPP-coupling spatial frequency
   f0 = 1/P_opt from m6_gratings.py) and a synthetic surface realising it; verifies rms and PSD.
 * Compares with a generic etched-metal PSD (k-correlation / ABC model, [BK] typical parameters) and the
   reconstructed ENEA band, to show what an AFM acceptance test should measure.
PSD convention: 2-D isotropic PSD C(q) with rms^2 = integral C(q) d^2q/(2 pi)^2, q = 2 pi f.
Run: python3 sim/m6_psd.py -> docs/models/figs/m6_psd.png, m6_psd.txt
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from m6_common import savefig, Tee

out = Tee("m6_psd.txt")
rng = np.random.default_rng(6)

L = 40e-6        # patch side (m)
n = 1024
dx = L / n
fx = np.fft.fftfreq(n, dx)
FX, FY = np.meshgrid(fx, fx)
F = np.hypot(FX, FY)


def ring_surface(f0, rel_width, rms):
    amp = np.exp(-0.5 * ((F - f0) / (rel_width * f0)) ** 2)
    phase = np.exp(2j * np.pi * rng.random((n, n)))
    z = np.real(np.fft.ifft2(np.sqrt(amp) * phase))
    z *= rms / z.std()
    return z


def radial_psd(z):
    Z = np.fft.fft2(z) * dx * dx
    P2 = np.abs(Z) ** 2 / L ** 2          # m^4 -> C(q) in m^4 (2-D PSD per unit area)
    fb = np.logspace(4.3, 7.4, 90)
    idx = np.digitize(F.ravel(), fb)
    Pr = np.array([P2.ravel()[idx == i].mean() if np.any(idx == i) else np.nan for i in range(1, len(fb))])
    return np.sqrt(fb[1:] * fb[:-1]), Pr


targets = {  # (label, f0 = 1/P_opt (m^-1) from m6_gratings.txt, rel_width, rms)
    "PdD face in D2O, 785 nm (C1 cathode)": (1 / 562e-9, 0.08, 30e-9),
    "PdD face in vacuum, 785 nm (C3)": (1 / 761e-9, 0.05, 36e-9),
    "Au 20 nm over PdD in D2O, 785 nm": (1 / 566e-9, 0.01, 13e-9),
}
fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.5))
for i, (lab, (f0, w, rms)) in enumerate(targets.items()):
    z = ring_surface(f0, w, rms)
    fc, Pr = radial_psd(z)
    Ztot = np.sum(np.abs(np.fft.fft2(z) * dx * dx) ** 2 / L ** 2) / L ** 2 * 1.0
    out(f"{lab}: ring f0 = {f0:.3e} /m (period {1e9/f0:.0f} nm), sigma_f = {w*100:.0f} %, rms = {rms*1e9:.0f} nm; "
        f"realised rms = {z.std()*1e9:.2f} nm; Parseval rms from PSD = {np.sqrt(Ztot)*1e9:.2f} nm")
    ax[0].loglog(fc, Pr, color=f"C{i}", label=lab)
    if i == 0:
        im = ax[1].imshow(z[:256, :256] * 1e9, extent=[0, 10, 0, 10], cmap="viridis")
        plt.colorbar(im, ax=ax[1], label="height (nm)")
# generic etched Pd: ABC / k-correlation model C(f) = A / (1 + (B f)^2)^((1+H)) scaled to rms 50 nm
fb = np.logspace(4.3, 7.4, 200)
H, B = 0.8, 2e-6
C = 1 / (1 + (B * fb) ** 2) ** (1 + H)
norm = np.trapezoid(C * 2 * np.pi * fb, fb)
C *= (50e-9) ** 2 / norm
ax[0].loglog(fb, C, "k--", label="generic etched metal (ABC model, rms 50 nm, corr. 2 um) [BK]")
ax[0].axvspan(1e5, 1e7, color="orange", alpha=0.12, label="reconstructed ENEA band (f = 1e5-1e7 /m)")
for lam, nd in ((0.633, 1.378), (0.785, 1.363), (1.064, 1.346)):
    ax[0].axvline(nd / (lam * 1e-6), color="r", lw=0.7, ls=":")
ax[0].text(1.8e6, 1e-33, "SPP coupling f at 1064/785/633 nm (PdD/D2O)", rotation=90, fontsize=6, color="r")
ax[0].set(xlabel="spatial frequency f (1/m)", ylabel="radial 2-D PSD (m^4)", title="Target ring PSDs vs generic etched surface",
          ylim=(1e-36, 1e-24))
ax[0].legend(fontsize=6, loc="lower left")
ax[1].set(xlabel="um", ylabel="um", title="Realisation: ring PSD, P = 562 nm, rms 30 nm")
for a in ax[:1]:
    a.grid(alpha=0.3, which="both")
savefig(fig, "m6_psd.png")
frac = np.trapezoid((C * 2 * np.pi * fb)[(fb > 1.5e6) & (fb < 2.0e6)], fb[(fb > 1.5e6) & (fb < 2.0e6)]) / np.trapezoid(C * 2 * np.pi * fb, fb)
out(f"Generic etched surface (rms 50 nm): fraction of height variance within the 785-nm SPP ring (1.5-2.0e6 /m): {frac:.3f}")
out(f"  -> equivalent ring rms = {50*np.sqrt(frac):.1f} nm, i.e. ~1/3-1/2 of the ~30 nm needed for critical coupling.")
out("Rule of thumb (first-order perturbation, calibrated on the RCWA sinusoid): an isotropic ring PSD couples")
out("like a 1-D sinusoid of peak-to-valley h when its in-ring rms ~ h/2 (azimuthal average of cos^2 = 1/2).")
out.save()
