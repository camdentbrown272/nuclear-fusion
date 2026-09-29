"""M7 reference geometry of the dual-chamber permeation-foil target
(R7 section 7.1 / master-plan C3) plus gamma-detector layouts.

Coordinates in cm. The foil normal is +z: vacuum (detector) side z > 0,
electrolyte side z < -L. Everything is cylindrically symmetric about z except
the gamma detectors, which sit on the +-x (and optionally +-y) axes.

Assumptions (documented in M7 section 3):
  * Pd foil, 20 mm active diameter (brief: 10-25 mm), clamped out to r = 16 mm.
  * Vacuum-side clamp flange: 316L SS ring, r = 10-35 mm, 5 mm thick.
  * Electrolyte-side clamp: PTFE ring, r = 10-35 mm, 5 mm thick.
  * Electrolyte: D2O (+0.1 M LiOD, neglected), 26 mL: r < 10 mm for the
    5 mm next to the foil, then r < 15 mm for 35 mm. PTFE cell wall 5 mm.
    Pt anode (thin mesh) neglected.
  * UHV chamber: 316L tube ID 60 mm, wall 1.5 mm, top plate 3 mm at z = 60 mm.
  * Si Delta-E/E telescope on axis: 25 um + 1000 um Si, 24 mm diameter
    (450 mm2) at z = 20 mm and 22 mm, FR4 carrier 1.6 mm behind.
  * Enclosure: closed 5 cm HDPE box, inner half-size 16 cm, standing in for
    the neutron moderator / inner shield liner that surrounds everything.
  * Gamma detectors: cylinders along x, front face at |x| = xf (default
    3.7 cm, i.e. 2 mm clear of the flange), Al can 0.5 mm front / 1 mm side.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from m7_mc import Geometry, Tube, Box  # noqa: E402

UM = 1e-4

DETECTORS = {
    #  name      material  radius  length  can(front, side)  gap(can->crystal)
    'BGO3x3':  ('BGO',   3.81, 7.62, 0.05, 0.10, 0.05),
    'NAI3x3':  ('NAI',   3.81, 7.62, 0.05, 0.10, 0.05),
    'NAI5x5':  ('NAI',   6.35, 12.7, 0.05, 0.10, 0.05),
    'LABR3x3': ('LABR3', 3.81, 7.62, 0.05, 0.10, 0.05),
    'BGO4x4':  ('BGO',   5.08, 10.16, 0.05, 0.10, 0.05),
    # HPGe ~100 % relative: 80 mm dia x 80 mm crystal, 1 mm Al endcap, 5 mm vacuum gap
    'HPGE100': ('GE',    4.00, 8.00, 0.10, 0.10, 0.50),
}


def build(L_um=100.0, det='BGO3x3', xf=3.7, zc=-0.5, four=False, world=30.0,
          with_si=True, chamber_wall=0.15, chamber_mat='SS316', encl=16.0, encl_t=5.0,
          axes=None):
    """axes: optional list of (axis, sign, face_distance_cm) for the gamma detectors;
    default = +-x at xf (plus +-y if four). 'zpair' shortcut: +z above the top
    plate and -z below the cell bottom."""
    L = L_um * UM
    g = Geometry([world, world, world], 'AIR')
    # --- Si telescope (vacuum side)
    if with_si:
        g.add(Tube(2, [0, 0, 2.0 + 12.5 * UM], 0, 1.2, 12.5 * UM), 'SI', 'si_dE')
        g.add(Tube(2, [0, 0, 2.2 + 500 * UM], 0, 1.2, 500 * UM), 'SI', 'si_E')
        g.add(Tube(2, [0, 0, 2.35 + 0.08], 0, 2.0, 0.08), 'FR4', 'si_pcb')
    # --- foil
    g.add(Tube(2, [0, 0, -L / 2], 0, 1.6, L / 2), 'PD', 'foil')
    # --- vacuum volumes
    g.add(Tube(2, [0, 0, 0.25], 0, 1.0, 0.25), 'VAC', 'vac_ap')
    g.add(Tube(2, [0, 0, 0.5 + 2.75], 0, 3.0, 2.75), 'VAC', 'vac')
    # --- vacuum-side hardware
    g.add(Tube(2, [0, 0, 0.25], 1.0, 3.5, 0.25), 'SS316', 'flange')
    g.add(Tube(2, [0, 0, 0.5 + 2.75], 3.0, 3.0 + chamber_wall, 2.75), chamber_mat, 'chamber')
    g.add(Tube(2, [0, 0, 6.0 + 0.15], 0, 3.0 + chamber_wall, 0.15), chamber_mat, 'topplate')
    # --- electrolyte side
    g.add(Tube(2, [0, 0, -L - 0.25], 0, 1.0, 0.25), 'D2O', 'electrolyte')
    g.add(Tube(2, [0, 0, -L - 0.5 - 1.75], 0, 1.5, 1.75), 'D2O', 'electrolyte2')
    g.add(Tube(2, [0, 0, -L - 0.25], 1.0, 3.5, 0.25), 'PTFE', 'cellclamp')
    g.add(Tube(2, [0, 0, -L - 0.5 - 1.75], 1.5, 2.0, 1.75), 'PTFE', 'cellwall')
    g.add(Tube(2, [0, 0, -L - 4.0 - 0.25], 0, 2.0, 0.25), 'PTFE', 'cellbottom')
    # --- gamma detectors
    sens = ['si_dE', 'si_E'] if with_si else []
    if det is not None:
        mat, R, Lc, cf, cs, gap = DETECTORS[det]
        if axes is None:
            axes = [(0, 1, xf), (0, -1, xf)] + ([(1, 1, xf), (1, -1, xf)] if four else [])
        elif axes == 'zpair':
            axes = [(2, 1, 6.4), (2, -1, L + 4.6)]
        elif axes == 'xz':
            axes = [(0, 1, xf), (0, -1, xf), (2, 1, 6.4), (2, -1, L + 4.6)]
        for k, (ax, sg, xf) in enumerate(axes):
            nm = f'det{k}'
            cc = [0.0, 0.0, zc]
            cc[ax] = sg * (xf + cf + gap + Lc / 2)
            g.add(Tube(ax, cc, 0, R, Lc / 2), mat, nm)
            cn = [0.0, 0.0, zc]
            cn[ax] = sg * (xf + (cf + gap + Lc) / 2)
            hl = (cf + gap + Lc) / 2
            if gap < 0.3:      # scintillator: Al can + reflector, modelled as Al
                g.add(Tube(ax, cn, 0, R + cs, hl), 'AL', f'can{k}')
            else:              # HPGe: Al endcap (front disk + side shell) around vacuum
                cap = [0.0, 0.0, zc]
                cap[ax] = sg * (xf + cf / 2)
                g.add(Tube(ax, cap, 0, R + cs + 0.1, cf / 2), 'AL', f'cap{k}')
                g.add(Tube(ax, cn, R + cs, R + cs + 0.1, hl), 'AL', f'side{k}')
                g.add(Tube(ax, cn, 0, R + cs, hl), 'VAC', f'can{k}')
            sens.append(nm)
    # --- enclosure: 5 cm HDPE shell (neutron moderator / inner shield liner)
    if encl:
        ext = max([xf + DETECTORS[det][2] + 0.3 for (_, _, xf) in axes]) if det is not None else 0
        encl = max(encl, ext + 0.5)
        g.add(Box([0, 0, 0], [encl, encl, encl]), 'AIR', 'cavity')
        g.add(Box([0, 0, 0], [encl + encl_t] * 3), 'HDPE', 'enclosure')
    return g, sens


REGION_GROUPS = {
    'foil': ['foil'],
    'electrolyte': ['electrolyte', 'electrolyte2'],
    'PTFE cell': ['cellclamp', 'cellwall', 'cellbottom'],
    'SS flange/chamber': ['flange', 'chamber', 'topplate'],
    'Si telescope + PCB': ['si_dE', 'si_E', 'si_pcb'],
    'gamma detectors (crystal+can)': ['det0', 'det1', 'det2', 'det3', 'can0', 'can1', 'can2', 'can3',
                                      'cap0', 'cap1', 'side0', 'side1'],
    'air / outside': ['fill', 'cavity'],
    'enclosure (HDPE liner/moderator)': ['enclosure'],
}
