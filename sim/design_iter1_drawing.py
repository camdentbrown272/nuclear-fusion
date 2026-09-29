"""Iteration-1 drawings: DFM cell cross-section and exit-face skin rotation.

Schematic, not a machining drawing. Vertical scale of thin layers is exaggerated
(labelled). Outputs docs/design/figs/iter1_dfm_cell.png and iter1_skins.png.
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
    wall = 4.0        # PCTFE wall
    # --- electrolyte cell (above membrane, y > 0) ---
    box(ax, -R - wall, 0.3, wall, 22, "#d9d9f3", label="PCTFE")
    box(ax, R, 0.3, wall, 22, "#d9d9f3", label="PCTFE")
    box(ax, -R, 0.3, 2 * R, 14, "#cfe8ff")
    ax.text(0, 9.5, "1.0 M LiOD / D$_2$O  (~12 mL)", ha="center", fontsize=9)
    ax.plot([-R + 0.5, R - 0.5], [4.3, 4.3], color="#888", lw=3, ls=(0, (1, 0.6)))
    ax.text(R - 0.6, 5.0, "Pt mesh anode, g = 4.0 ± 0.2", ha="right", fontsize=8)
    box(ax, -R, 14.3, 2 * R, 8, "#f4f4f4")
    ax.text(0, 16.4, "headspace ≤ 15 mL, D$_2$ 0.50 bar abs\n+1 % Kr tracer", ha="center", fontsize=8)
    box(ax, -4.5, 19.6, 9, 1.8, "#ffe0b3", label="recombiner (baffled)")
    box(ax, -R - wall - 3, 22.3, 2 * (R + wall + 3), 2.2, "#bfbfbf", label="316L lid, Cu-gasket seal")
    ax.plot([-6, -6], [24.5, 4.6], color="#b36b00", lw=1.2)
    ax.text(-6.3, 25.2, "anode feedthrough (+)", fontsize=7, ha="center")
    ax.plot([6, 6], [24.5, 27], color="k", lw=1)
    ax.text(6, 27.4, "P, T, burst disk 3.5 bar,\nHe-sample valve", fontsize=7, ha="center")

    # --- membrane (thickness exaggerated) ---
    t = 0.3
    box(ax, -R - 3, 0.0, 2 * R + 6, t, "#7f7f7f")
    ax.annotate("Pd membrane 15 ± 2 µm (drawn ×20)\nactive Ø20.0, blank Ø26; entry face = cathode (ground)",
                xy=(R + 1.5, t / 2), xytext=(R + 4.5, 3.0), fontsize=8,
                arrowprops=dict(arrowstyle="->", lw=0.7))
    box(ax, -R - wall, 0.3, wall, 0.9, "#333333")
    box(ax, R, 0.3, wall, 0.9, "#333333")
    ax.text(R + wall + 0.3, 0.9, "floating FFKM seal\n(≥1 mm radial travel)", fontsize=7, va="center")
    # exit-face skins (bottom surface), quadrants shown as two halves in section
    colors = {"L": "#d4a017", "M": "#6fa8dc", "F": "#b7b7b7", "X": "#93c47d"}
    box(ax, -R, -0.12, R - 0.5, 0.12, colors["F"], lw=0.3)
    box(ax, 0.5, -0.12, R - 0.5, 0.12, colors["L"], lw=0.3)
    ax.annotate("exit skin Q3: F (bare Pd)", xy=(-5, -0.12), xytext=(-18, -2.6), fontsize=7, ha="right",
                arrowprops=dict(arrowstyle="->", lw=0.6))
    ax.annotate("exit skin Q4: L (Au 20–50 nm)", xy=(8.5, -0.12), xytext=(18, -2.6), fontsize=7, ha="left",
                arrowprops=dict(arrowstyle="->", lw=0.6))

    # --- grid, septum, telescope, front volume (below) ---
    y_grid = -1.2
    for xg in [x * 1.1 for x in range(-9, 10)]:
        box(ax, xg - 0.05, y_grid, 0.1, 0.25, "#5b3a29", lw=0)
    box(ax, -0.25, y_grid - 0.1, 0.5, 0.45, "#5b3a29", lw=0)
    ax.text(-18, y_grid + 0.1, "Mo catch grid 0.10, 1.0 hex, open 0.80", fontsize=7, ha="right")
    box(ax, -0.15, -5.0, 0.3, 3.7, "#a0522d", lw=0.4)
    ax.text(0.5, -3.2, "cross septum\n4.0 × 0.3", fontsize=7)
    y_de = -5.3
    box(ax, -12.5, y_de, 12.2, 0.3, "#3c78d8", lw=0.4)
    box(ax, 0.3, y_de, 12.2, 0.3, "#3c78d8", lw=0.4)
    ax.text(-18, y_de + 0.15, "ΔE 25 µm, 600 mm², 4 quadrants\n(4 ± 1 below exit face)", fontsize=7, ha="right", va="center")
    box(ax, -14, y_de - 1.8, 28, 0.5, "#1c4587", lw=0.4)
    ax.text(-18, y_de - 1.55, "E 500 µm", fontsize=7, ha="right", va="center")
    box(ax, -R - wall - 3, -13.5, 2 * (R + wall + 3), 13.0, "none", ec="#444", lw=1.5, ls="--")
    ax.text(0, -10.5, "FRONT VOLUME ~30 cm³, D$_2$ 0.50 bar abs\nP$_{cell}$ − P$_{front}$ = +30 ± 20 mbar\n"
            "vacuum-fired 316L, ≤10 metal seals, no glass/epoxy/ion gauge",
            ha="center", fontsize=8)
    for xp, lab in [(-10, "Pd–Ag element 350 °C\n(only D$_2$ path)"), (0, "all-metal valve →\nHe manifold / HR-QMS"),
                    (10, "capacitance gauge")]:
        ax.plot([xp, xp], [-13.5, -16], color="k", lw=1)
        ax.text(xp, -17.3, lab, ha="center", fontsize=7)

    # proton tracks: exit-face and entry-face origin
    ax.annotate("", xy=(-7, y_de + 0.3), xytext=(-5, 0.3), arrowprops=dict(arrowstyle="->", color="r", lw=1))
    ax.text(-9.5, 1.0, "entry-face p: 3.02 → 2.05 MeV", color="r", fontsize=7)
    ax.annotate("", xy=(7.5, y_de + 0.3), xytext=(6, -0.1), arrowprops=dict(arrowstyle="->", color="m", lw=1))
    ax.text(8.0, -3.6, "exit-face p: 3.00 MeV", color="m", fontsize=7)

    ax.set_xlim(-34, 34)
    ax.set_ylim(-19, 29)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Iteration 1 — DFM cell cross-section (mm; thin layers exaggerated)", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "iter1_dfm_cell.png"), dpi=160)
    plt.close(fig)


def skins():
    layout = {"A1": ["L", "M", "F", "X"], "A2": ["F", "L", "X", "ED"],
              "A3": ["X", "ED", "L", "M"], "A4": ["M", "F", "ED", "L"]}
    colors = {"L": "#d4a017", "M": "#6fa8dc", "F": "#d9d9d9", "X": "#93c47d", "ED": "#e06666"}
    # quadrant order NE, NW, SW, SE -> angles
    ang = [(0, 90), (90, 180), (180, 270), (270, 360)]
    fig, axs = plt.subplots(1, 4, figsize=(12, 3.4))
    import numpy as np
    for ax, (name, q) in zip(axs, layout.items()):
        for (a0, a1), s in zip(ang, q):
            th = np.linspace(np.radians(a0), np.radians(a1), 40)
            pts = [(0, 0)] + [(10 * np.cos(x), 10 * np.sin(x)) for x in th]
            ax.add_patch(Polygon(pts, fc=colors[s], ec="k", lw=0.6))
            am = np.radians((a0 + a1) / 2)
            ax.text(6 * np.cos(am), 6 * np.sin(am), s, ha="center", va="center", fontsize=12, weight="bold")
        ax.plot([-10, 10], [0, 0], color="#d4a017", lw=3)
        ax.plot([0, 0], [-10, 10], color="#d4a017", lw=3)
        ax.add_patch(plt.Circle((0, 0), 0.8, fc="#444", ec="none"))
        ax.set_xlim(-11, 11); ax.set_ylim(-11, 11); ax.set_aspect("equal"); ax.axis("off")
        ax.set_title(name)
    fig.suptitle("Exit-face skins (viewed from telescope). L=Au cap, M=Ni/Au-mask, F=bare, X=Pd/CaO, ED=defect-rich Pd;"
                 " gold cross = 1 mm Au web, dot = $^{10}$B marker", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "iter1_skins.png"), dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    cell_section()
    skins()
    print("wrote", OUT)
