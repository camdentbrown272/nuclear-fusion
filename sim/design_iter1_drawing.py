"""Iteration-1 (rev C) drawings: DFM cell cross-section and membrane-type allocation.

Schematic, not a machining drawing. Vertical scale of thin layers is exaggerated
(labelled). Outputs docs/design/figs/iter1_dfm_cell.png and iter1_skins.png.
Geometry per docs/design/iteration-1.md rev C (ADR-005/006/007).
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon

OUT = os.environ.get("ITER1_OUT", os.path.join(os.path.dirname(__file__), "..", "docs", "design", "figs"))
os.makedirs(OUT, exist_ok=True)


def box(ax, x, y, w, h, fc, ec="k", lw=0.8, hatch=None, label=None, **kw):
    ax.add_patch(Rectangle((x, y), w, h, fc=fc, ec=ec, lw=lw, hatch=hatch, **kw))
    if label:
        ax.text(x + w / 2, y + h / 2, label, ha="center", va="center", fontsize=7)


def cell_section():
    fig, ax = plt.subplots(figsize=(10, 8.4))
    R = 10.0          # active radius (mm)
    wall = 4.0        # 316L wall incl. 1 mm PTFE liner
    # --- electrolyte cell (above membrane, y > 0) ---
    box(ax, -R - wall, 0.3, wall - 1, 22, "#c9ced6", label="316L\n+jacket")
    box(ax, R + 1, 0.3, wall - 1, 22, "#c9ced6", label="316L")
    box(ax, -R - 1, 0.3, 1, 22, "#f2efe6", lw=0.4)
    box(ax, R, 0.3, 1, 22, "#f2efe6", lw=0.4)
    ax.text(-R - 0.5, 12.0, "PTFE liner", rotation=90, fontsize=6, ha="center", va="center")
    box(ax, -R, 0.3, 2 * R, 14, "#cfe8ff")
    ax.text(0, 9.5, "1.0 M LiOD / D$_2$O, low-T, ≥99.9 % D  (~12 mL)", ha="center", fontsize=9)
    ax.plot([-R + 0.5, R - 0.5], [4.3, 4.3], color="#888", lw=3, ls=(0, (1, 0.6)))
    ax.text(R - 0.6, 5.0, "Pt mesh anode, g = 4.0 ± 0.2", ha="right", fontsize=8)
    ax.plot([-8.5, 8.5], [2.3, 2.3], ls="none", marker="o", ms=3, color="#888")
    ax.text(-R + 0.6, 2.9, "aux. Pt ring cathode (current steering)", fontsize=7)
    box(ax, -R, 14.3, 2 * R, 8, "#f4f4f4")
    ax.text(0, 16.4, "headspace ≤ 15 mL, D$_2$ 0.50 bar abs\n+1 % Kr tracer", ha="center", fontsize=8)
    box(ax, -4.5, 19.6, 9, 1.8, "#ffe0b3", label="recombiner (baffled)")
    box(ax, -R - wall - 3, 22.3, 2 * (R + wall + 3), 2.2, "#bfbfbf", label="316L lid, Au-wire seal")
    ax.plot([-6, -6], [24.5, 4.6], color="#b36b00", lw=1.2)
    ax.text(-6.3, 25.2, "anode feedthrough (+)", fontsize=7, ha="center")
    ax.plot([6, 6], [24.5, 27], color="k", lw=1)
    ax.text(6, 27.4, "P, T, burst disk 3.5 bar,\nHe-sample valve", fontsize=7, ha="center")

    # --- membrane (thickness exaggerated) ---
    t = 0.3
    box(ax, -R - 0.6, 0.0, 2 * R + 1.2, t, "#7f7f7f")
    # corrugated Pt annulus (non-hydriding), diffusion-bonded to the membrane rim
    import numpy as np
    for sgn in (-1, 1):
        xs = np.linspace(R + 0.3, R + 7.0, 60)
        ys = 0.15 + 0.5 * np.exp(-((xs - (R + 3.5)) / 0.9) ** 2)
        ax.plot(sgn * xs, ys, color="#8e8e8e", lw=2.2)
        box(ax, sgn * (R + 6.2) - 0.4, -0.35, 0.8, 0.3, "#d4a017", lw=0.3)
        box(ax, sgn * (R + 6.2) - 0.4, 0.65, 0.8, 0.3, "#d4a017", lw=0.3)
    ax.annotate("Pd membrane 12 ± 1 µm (drawn ×25)\nactive Ø20.0; entry face = cathode (ground)",
                xy=(R + 1.5, t / 2), xytext=(R + 4.5, 3.0), fontsize=8,
                arrowprops=dict(arrowstyle="->", lw=0.7))
    ax.text(R + 7.6, 0.9, "Pt annulus 0.15, two folds, press-bonded Ø23–25;\nPd-only gauge rim Ø20–23; Au-wire seals", fontsize=7, va="center")
    # exit-face skins (bottom surface), quadrants shown as two halves in section
    colors = {"L": "#d4a017", "M": "#6fa8dc", "F": "#b7b7b7", "X": "#93c47d"}
    box(ax, -R, -0.12, 2 * R, 0.12, colors["L"], lw=0.3)
    ax.annotate("exit finish (one per membrane):\nH-L Au 50 nm | H-M Ni 20 nm | FX bare/CaO", xy=(-5, -0.12),
                xytext=(-18, -2.8), fontsize=7, ha="right", arrowprops=dict(arrowstyle="->", lw=0.6))

    # --- grid, septum, telescope, front volume (below) ---
    y_grid = -1.2
    for xg in [x * 1.1 for x in range(-9, 10)]:
        box(ax, xg - 0.05, y_grid, 0.1, 0.25, "#5b3a29", lw=0)
    box(ax, -0.25, y_grid - 0.1, 0.5, 0.45, "#5b3a29", lw=0)
    ax.text(-18, y_grid + 0.1, "Mo catch grid 0.30 (proof 1.5 bar), 1.0 hex", fontsize=7, ha="right")
    box(ax, -0.15, -5.0, 0.3, 3.7, "#a0522d", lw=0.4)
    ax.text(0.5, -3.2, "cross septum\n4.0 × 0.3", fontsize=7)
    y_de = -5.3
    box(ax, -12.5, y_de, 12.2, 0.3, "#3c78d8", lw=0.4)
    box(ax, 0.3, y_de, 12.2, 0.3, "#3c78d8", lw=0.4)
    ax.text(-18, y_de + 0.15, "ΔE 25 µm, 600 mm², 4 quadrants\n(4 ± 1 below exit face)", fontsize=7, ha="right", va="center")
    box(ax, -14, y_de - 1.8, 28, 0.5, "#1c4587", lw=0.4)
    ax.text(-18, y_de - 1.55, "E 500 µm", fontsize=7, ha="right", va="center")
    box(ax, -R - wall - 3, -15.5, 2 * (R + wall + 3), 15.0, "none", ec="#444", lw=1.5, ls="--")
    ax.text(0, -12.6, "FRONT VOLUME ~30 cm³, D$_2$ 0.50 bar abs\nP$_{cell}$ − P$_{front}$ = +30 ± 10 mbar\n"
            "metal seals with pumped interspaces; no glass/epoxy/ion gauge\nforward relief +150 mbar; reverse disk 20 ± 5 mbar\nPd–Ag exhaust → mass-flow → bellows pump → Pd–Ag → headspace (D₂ recycle)",
            ha="center", fontsize=8)
    for xp, lab in [(-10, "Pd–Ag element 350 °C\n(only D$_2$ path)"), (0, "all-metal valve →\nHe manifold / HR-QMS"),
                    (10, "capacitance gauge")]:
        ax.plot([xp, xp], [-15.5, -18], color="k", lw=1)
        ax.text(xp, -19.3, lab, ha="center", fontsize=7)

    # proton tracks: exit-face and entry-face origin
    ax.annotate("", xy=(-7, y_de + 0.3), xytext=(-5, 0.3), arrowprops=dict(arrowstyle="->", color="r", lw=1))
    ax.text(-9.5, 1.0, "entry-face p: 3.02 → 2.23 MeV (normal)", color="r", fontsize=7)
    ax.annotate("", xy=(7.5, y_de + 0.3), xytext=(6, -0.1), arrowprops=dict(arrowstyle="->", color="m", lw=1))
    ax.text(8.0, -3.6, "exit-face p: 3.02 MeV", color="m", fontsize=7)

    ax.set_xlim(-34, 34)
    ax.set_ylim(-21, 29)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Iteration 1 rev C — DFM cell cross-section (mm; thin layers exaggerated)", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "iter1_dfm_cell.png"), dpi=160)
    plt.close(fig)


def skins():
    """Membrane types (single-regime, ADR-005) and two-stage allocation."""
    import numpy as np
    colors = {"H-L": "#d4a017", "H-M": "#6fa8dc", "F": "#d9d9d9", "X": "#93c47d"}
    stages = {"Stage 1": [("H-L", "a"), ("H-L", "b"), ("H-M", "a"), ("FX", "b")],
              "Stage 2": [("FX", "a"), ("H-M", "b"), ("H-L", "a", "Al"), ("H-L", "b")]}
    fig, axs = plt.subplots(2, 4, figsize=(12, 6.4))
    for row, (st, cells) in enumerate(stages.items()):
        for col, c in enumerate(cells):
            ax = axs[row, col]
            typ, lot = c[0], c[1]
            if typ == "FX":
                for k, (a0, a1) in enumerate([(0, 90), (90, 180), (180, 270), (270, 360)]):
                    th = np.linspace(np.radians(a0), np.radians(a1), 40)
                    pts = [(0, 0)] + [(10 * np.cos(x), 10 * np.sin(x)) for x in th]
                    lab = "F" if k % 2 == 0 else "X"
                    ax.add_patch(Polygon(pts, fc=colors[lab], ec="k", lw=0.6))
                    am = np.radians((a0 + a1) / 2)
                    ax.text(6 * np.cos(am), 6 * np.sin(am), lab, ha="center", va="center", fontsize=12, weight="bold")
            else:
                ax.add_patch(plt.Circle((0, 0), 10, fc=colors[typ], ec="k", lw=0.6))
                ax.text(0, 0, typ, ha="center", va="center", fontsize=14, weight="bold")
            ax.plot([-10, 10], [0, 0], color="k", lw=0.5, ls=":")
            ax.plot([0, 0], [-10, 10], color="k", lw=0.5, ls=":")
            extra = " + Al additive" if len(c) > 2 else ""
            ax.set_title(f"{st} A{col + 1}: lot {lot}{extra}", fontsize=9)
            ax.set_xlim(-11, 11); ax.set_ylim(-11, 11); ax.set_aspect("equal"); ax.axis("off")
    fig.suptitle("Single-regime membranes (exit face, viewed from telescope; dotted = telescope quadrants).\n"
                 "H-L = Au 50 nm (x≥0.9, J≈0) · H-M = Ni 20 nm (x≈0.91, small flux) · FX = bare F / Pd-CaO X "
                 "(x≈0.65, max flux). Controls each stage: T-Pt, T-H(L), T-H(FX), B", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "iter1_skins.png"), dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    cell_section()
    skins()
    print("wrote", OUT)
