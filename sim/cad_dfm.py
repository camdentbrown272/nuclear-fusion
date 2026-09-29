"""Parametric 3D model of iteration-1 rev B (DFM cell + DFM-8 house).

Every dimension comes from docs/design/iteration-1.md (rev B). Units: mm. z = 0 is the
membrane exit face; the electrolytic cell is above (+z), the front volume/telescope below.

Outputs (docs/design/cad/):
  dfm_cell.step / .glb / .stl   detailed single cell (true 12 um membrane, hex catch grid)
  dfm_house.step / .glb         8-position array in its Cu / 3He / borated-HDPE house (simplified cells)

Requires: pip install cadquery
"""
import math
import os

import cadquery as cq

OUT = os.environ.get("CAD_OUT", os.path.join(os.path.dirname(__file__), "..", "docs", "design", "cad"))
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------- design parameters (rev B §3.2)
P = dict(
    membrane_t=0.012,        # 12 um Pd
    membrane_d=22.0,         # blank dia (bonded to Pt annulus over r 10.25..11)
    active_d=20.0,           # wetted/active
    annulus_id=20.5, annulus_od=34.0, annulus_t=0.10, corr_r=13.5, corr_h=0.6,
    grid_d=20.4, grid_t=0.10, hex_af=1.0, hex_web=0.10, rib_w=0.5,
    liner_t=1.0, bore_d=20.0, cell_wall=4.0, cell_h=30.0,
    anode_gap=4.0, anode_d=19.0, anode_t=0.10,
    electrolyte_h=14.0,
    flange_od=64.0, flange_t=6.0, bolt_pcd=56.0, bolt_d=5.5, n_bolts=8,
    front_id=48.0, front_od=56.0, front_depth=16.0,
    septum_h=4.0, septum_t=0.3, septum_l=26.0,
    de_gap=4.0, de_side=12.1, de_t=0.025, e_gap=1.5, e_side=27.0, e_t=0.5,
    board_d=46.0, board_t=1.0,
)

C = dict(  # RGBA colours
    steel=(0.72, 0.74, 0.77, 1.0), ptfe=(0.96, 0.95, 0.90, 1.0), electrolyte=(0.35, 0.62, 0.95, 0.35),
    pd=(0.50, 0.51, 0.54, 1.0), pt=(0.86, 0.86, 0.88, 1.0), au=(0.85, 0.65, 0.13, 1.0),
    mo=(0.45, 0.40, 0.36, 1.0), si=(0.20, 0.35, 0.80, 1.0), alumina=(0.93, 0.91, 0.85, 1.0),
    cu=(0.78, 0.45, 0.25, 1.0), hdpe=(0.95, 0.95, 0.93, 0.25), bhdpe=(0.75, 0.88, 0.75, 0.18),
    he3=(0.80, 0.82, 0.85, 1.0), veto=(0.55, 0.35, 0.75, 0.22), labr=(0.95, 0.90, 0.35, 1.0),
    recomb=(1.0, 0.85, 0.6, 1.0), cd=(0.6, 0.6, 0.65, 1.0), pdag=(0.62, 0.62, 0.66, 1.0),
)


def ring(od, id_, h, z0):
    return cq.Workplane("XY").workplane(offset=z0).circle(od / 2).circle(id_ / 2).extrude(h)


def disk(d, h, z0):
    return cq.Workplane("XY").workplane(offset=z0).circle(d / 2).extrude(h)


def bolt_holes(shape, z0, h):
    pts = [(P["bolt_pcd"] / 2 * math.cos(2 * math.pi * k / P["n_bolts"]),
            P["bolt_pcd"] / 2 * math.sin(2 * math.pi * k / P["n_bolts"])) for k in range(P["n_bolts"])]
    cut = cq.Workplane("XY").workplane(offset=z0 - 1).pushPoints(pts).circle(P["bolt_d"] / 2).extrude(h + 2)
    return shape.cut(cut)


# ---------------------------------------------------------------- parts
def membrane():
    return disk(P["membrane_d"], P["membrane_t"], 0.0)


def exit_finish(kind="H-L"):
    """50 nm Au cap drawn 0.02 mm thick for visibility (true 50 nm)."""
    return disk(P["active_d"], 0.02, -0.02) if kind == "H-L" else None


def pt_annulus():
    """0.10 mm Pt annulus with one convolution (radial compliance >= 1 mm), revolved profile."""
    ri, ro, t, rc, hc = P["annulus_id"] / 2, P["annulus_od"] / 2, P["annulus_t"], P["corr_r"], P["corr_h"]
    w = 1.2  # half-width of the convolution
    pts = [(ri, -t), (rc - w, -t), (rc - w / 2, hc - t), (rc + w / 2, hc - t), (rc + w, -t), (ro, -t),
           (ro, 0.0), (rc + w, 0.0), (rc + w / 2, hc), (rc - w / 2, hc), (rc - w, 0.0), (ri, 0.0)]
    return cq.Workplane("XZ").polyline(pts).close().revolve(360, (0, 0, 0), (0, 1, 0))


def au_wire(r, z):
    return cq.Workplane("XY").workplane(offset=z).center(r, 0).circle(0.25).revolve(360, (-r, 0, 0), (-r, 1, 0))


def catch_grid():
    """Photo-etched Mo grid: 1.0 mm A/F hex holes, 0.10 webs, 0.5 mm ribs along the septum cross."""
    g = disk(P["grid_d"], P["grid_t"], -P["grid_t"])
    pitch = P["hex_af"] + P["hex_web"]
    pts = []
    R = P["grid_d"] / 2 - 0.6
    ny = int(R / (pitch * math.sqrt(3) / 2)) + 1
    for j in range(-ny, ny + 1):
        y = j * pitch * math.sqrt(3) / 2
        off = (pitch / 2) if j % 2 else 0.0
        for i in range(-ny - 1, ny + 2):
            x = i * pitch + off
            if x * x + y * y > R * R:
                continue
            if abs(x) < P["rib_w"] / 2 + 0.6 or abs(y) < P["rib_w"] / 2 + 0.6:
                continue  # solid ribs under the cross septum
            pts.append((x, y))
    holes = (cq.Workplane("XY").workplane(offset=-P["grid_t"] - 0.1).pushPoints(pts)
             .polygon(6, P["hex_af"] / math.cos(math.pi / 6)).extrude(P["grid_t"] + 0.2))
    return g.cut(holes)


def septum():
    z0 = -P["grid_t"] - P["septum_h"]
    a = cq.Workplane("XY").workplane(offset=z0).rect(P["septum_l"], P["septum_t"]).extrude(P["septum_h"])
    b = cq.Workplane("XY").workplane(offset=z0).rect(P["septum_t"], P["septum_l"]).extrude(P["septum_h"])
    return a.union(b)


def telescope():
    """Quadrant dE (25 um) on an alumina board + E (500 um) on a second board."""
    zde = -P["de_gap"]
    s, gap = P["de_side"], 0.3
    de = None
    for sx in (-1, 1):
        for sy in (-1, 1):
            q = (cq.Workplane("XY").workplane(offset=zde - P["de_t"])
                 .center(sx * (s / 2 + gap / 2), sy * (s / 2 + gap / 2)).rect(s, s).extrude(P["de_t"]))
            de = q if de is None else de.union(q)
    board1 = (disk(P["board_d"], P["board_t"], zde - P["de_t"] - P["board_t"])
              .cut(cq.Workplane("XY").workplane(offset=zde - 3).rect(2 * s + gap - 0.6, 2 * s + gap - 0.6).extrude(5)))
    ze = zde - P["e_gap"]
    e = cq.Workplane("XY").workplane(offset=ze - P["e_t"]).rect(P["e_side"], P["e_side"]).extrude(P["e_t"])
    board2 = disk(P["board_d"], P["board_t"], ze - P["e_t"] - P["board_t"])
    return de, board1, e, board2


def cell_body():
    z0 = 0.35
    liner = ring(P["bore_d"] + 2 * P["liner_t"], P["bore_d"], P["cell_h"] - z0, z0)
    tube = ring(P["bore_d"] + 2 * P["liner_t"] + 2 * P["cell_wall"], P["bore_d"] + 2 * P["liner_t"],
                P["cell_h"] - z0, z0)
    flange = ring(P["flange_od"], P["bore_d"] + 2 * P["liner_t"], P["flange_t"], z0)
    flange = bolt_holes(flange, z0, P["flange_t"])
    top_fl = ring(44.0, P["bore_d"] + 2 * P["liner_t"], 4.0, P["cell_h"] - 4.0)
    body = tube.union(flange).union(top_fl)
    lid = disk(44.0, 5.0, P["cell_h"])
    for x, y, d, h in [(-6.0, 0.0, 4.0, 8.0), (6.0, 4.0, 6.0, 12.0), (6.0, -5.0, 6.0, 10.0)]:
        lid = lid.union(cq.Workplane("XY").workplane(offset=P["cell_h"] + 5.0).center(x, y).circle(d / 2).extrude(h))
    return liner, body, lid


def cell_internals():
    z0 = P["membrane_t"]
    electrolyte = disk(P["bore_d"] - 0.02, P["electrolyte_h"], z0)
    # Pt mesh anode (disk with square perforation pattern, >= 50 % open)
    anode = disk(P["anode_d"], P["anode_t"], z0 + P["anode_gap"])
    pts = [(x * 1.0, y * 1.0) for x in range(-9, 10) for y in range(-9, 10) if (x * x + y * y) < 81]
    anode = anode.cut(cq.Workplane("XY").workplane(offset=z0 + P["anode_gap"] - 0.1)
                      .pushPoints(pts).rect(0.72, 0.72).extrude(0.4))
    feed = cq.Workplane("XY").workplane(offset=z0 + P["anode_gap"]).center(-6.0, 0).circle(0.5).extrude(
        P["cell_h"] + 5.0 - z0 - P["anode_gap"])
    recomb = cq.Workplane("XY").workplane(offset=P["cell_h"] - 5.0).center(0, 0).rect(9.0, 5.0).extrude(2.0)
    return electrolyte, anode, feed, recomb


def front_body():
    top = -0.35
    fl = ring(P["flange_od"], P["front_id"], P["flange_t"], top - P["flange_t"])
    fl = bolt_holes(fl, top - P["flange_t"], P["flange_t"])
    tube = ring(P["front_od"], P["front_id"], P["front_depth"] - P["flange_t"], top - P["front_depth"])
    bottom = disk(P["flange_od"], 6.0, top - P["front_depth"] - 6.0)
    zb = top - P["front_depth"] - 6.0
    ports = None
    for x, y in [(-14.0, 0.0), (0.0, 0.0), (14.0, 0.0)]:
        p = cq.Workplane("XY").workplane(offset=zb - 18.0).center(x, y).circle(3.0).circle(2.2).extrude(18.0)
        ports = p if ports is None else ports.union(p)
    return fl.union(tube).union(bottom), ports, zb


def front_devices(zb):
    pdag = cq.Workplane("XY").workplane(offset=zb - 18.0 - 28.0).center(-14.0, 0).circle(7.0).extrude(28.0)
    valve = cq.Workplane("XY").workplane(offset=zb - 18.0 - 16.0).center(0.0, 0).polygon(6, 16.0).extrude(16.0)
    gauge = cq.Workplane("XY").workplane(offset=zb - 18.0 - 24.0).center(14.0, 0).circle(9.0).extrude(24.0)
    relief = (cq.Workplane("XY").workplane(offset=-12.0).center(0, P["front_od"] / 2 + 2.0).circle(1.5).extrude(18.0))
    return pdag, valve, gauge, relief


def build_cell(detailed=True, finish="H-L"):
    a = cq.Assembly(name="DFM_cell")
    liner, body, lid = cell_body()
    a.add(body, name="cell_body_316L", color=cq.Color(*C["steel"]))
    a.add(liner, name="PTFE_liner", color=cq.Color(*C["ptfe"]))
    a.add(lid, name="lid_316L", color=cq.Color(*C["steel"]))
    el, anode, feed, recomb = cell_internals()
    a.add(el, name="electrolyte_LiOD_D2O", color=cq.Color(*C["electrolyte"]))
    a.add(anode, name="Pt_mesh_anode", color=cq.Color(*C["pt"]))
    a.add(feed, name="anode_feedthrough", color=cq.Color(*C["cu"]))
    a.add(recomb, name="recombiner", color=cq.Color(*C["recomb"]))
    a.add(membrane(), name="Pd_membrane_12um", color=cq.Color(*C["pd"]))
    fin = exit_finish(finish)
    if fin is not None:
        a.add(fin, name="exit_finish_Au50nm", color=cq.Color(*C["au"]))
    a.add(pt_annulus(), name="Pt_corrugated_annulus", color=cq.Color(*C["pt"]))
    a.add(au_wire(16.5, 0.35), name="Au_wire_seal_upper", color=cq.Color(*C["au"]))
    a.add(au_wire(16.5, -0.35 - P["annulus_t"]), name="Au_wire_seal_lower", color=cq.Color(*C["au"]))
    a.add(catch_grid() if detailed else disk(P["grid_d"], P["grid_t"], -P["grid_t"]),
          name="Mo_catch_grid", color=cq.Color(*C["mo"]))
    a.add(septum(), name="cross_septum", color=cq.Color(*C["cu"]))
    de, b1, e, b2 = telescope()
    a.add(de, name="Si_dE_25um_quadrants", color=cq.Color(*C["si"]))
    a.add(b1, name="dE_alumina_board", color=cq.Color(*C["alumina"]))
    a.add(e, name="Si_E_500um", color=cq.Color(*C["si"]))
    a.add(b2, name="E_alumina_board", color=cq.Color(*C["alumina"]))
    fb, ports, zb = front_body()
    a.add(fb, name="front_volume_316L", color=cq.Color(*C["steel"]))
    a.add(ports, name="front_ports", color=cq.Color(*C["steel"]))
    pdag, valve, gauge, relief = front_devices(zb)
    a.add(pdag, name="PdAg_element_350C", color=cq.Color(*C["pdag"]))
    a.add(valve, name="all_metal_valve_He_manifold", color=cq.Color(*C["steel"]))
    a.add(gauge, name="capacitance_gauge", color=cq.Color(*C["steel"]))
    a.add(relief, name="reverse_relief_line", color=cq.Color(*C["steel"]))
    return a


# ---------------------------------------------------------------- house (rev B §4.2)
H = dict(cav=(450.0, 300.0, 300.0), cu_t=50.0, hdpe_t=105.0, bhdpe_t=200.0, pitch=90.0)


def box_shell(outer, inner, open_top=True):
    ox, oy, oz = outer
    ix, iy, iz = inner
    b = cq.Workplane("XY").box(ox, oy, oz)
    cut = cq.Workplane("XY").box(ix, iy, iz)
    if open_top:  # remove the lid so the model can be looked into
        cut = cut.union(cq.Workplane("XY").workplane(offset=iz / 2 - 1).rect(ix, iy).extrude(oz))
    return b.cut(cut)


def build_house():
    a = cq.Assembly(name="DFM8_house")
    cx, cy, cz = H["cav"]
    cu_o = (cx + 2 * H["cu_t"], cy + 2 * H["cu_t"], cz + 2 * H["cu_t"])
    a.add(box_shell(cu_o, H["cav"]), name="Cu_inner_shield_5cm", color=cq.Color(*C["cu"]))
    hd_o = (cu_o[0] + 2 * H["hdpe_t"], cu_o[1] + 2 * H["hdpe_t"], cu_o[2] + 2 * H["hdpe_t"])
    a.add(box_shell(hd_o, cu_o), name="HDPE_moderator_3He_bank", color=cq.Color(*C["hdpe"]))
    bh_o = (hd_o[0] + 2 * H["bhdpe_t"], hd_o[1] + 2 * H["bhdpe_t"], hd_o[2] + 2 * H["bhdpe_t"])
    a.add(box_shell(bh_o, (hd_o[0] + 2, hd_o[1] + 2, hd_o[2] + 2)), name="borated_HDPE_20cm_plus_Cd",
          color=cq.Color(*C["bhdpe"]))
    # 24 3He tubes, 1" x 400 mm, vertical, 45 mm outside the Cu (8 per long wall, 4 per short wall)
    tubes = None
    d_off = H["cu_t"] + 45.0
    for k in range(8):
        x = -cu_o[0] / 2 + 40 + k * (cu_o[0] - 80) / 7
        for y in (-(cy / 2 + d_off), cy / 2 + d_off):
            t = cq.Workplane("XY").workplane(offset=-200).center(x, y).circle(12.7).extrude(400)
            tubes = t if tubes is None else tubes.union(t)
    for k in range(4):
        y = -cu_o[1] / 2 + 40 + k * (cu_o[1] - 80) / 3
        for x in (-(cx / 2 + d_off), cx / 2 + d_off):
            tubes = tubes.union(cq.Workplane("XY").workplane(offset=-200).center(x, y).circle(12.7).extrude(400))
    a.add(tubes, name="He3_tubes_24x", color=cq.Color(*C["he3"]))
    # muon veto panels (50 mm plastic scintillator): top + two long sides
    vx, vy, vz = bh_o
    a.add(cq.Workplane("XY").box(vx, vy, 50).translate((0, 0, vz / 2 + 55)), name="muon_veto_top",
          color=cq.Color(*C["veto"]))
    for s_ in (-1, 1):
        a.add(cq.Workplane("XY").box(vx, 50, vz).translate((0, s_ * (vy / 2 + 55), 0)),
              name=f"muon_veto_side_{'N' if s_ > 0 else 'S'}", color=cq.Color(*C["veto"]))
    # cells: 2 rows x 4 at 90 mm pitch; labels per rev B §4.1
    labels = [["A1", "A2", "A3", "A4"], ["T-Pt", "T-H(L)", "T-H(FX)", "B"]]
    simple = build_cell(detailed=False)
    for r, row in enumerate(labels):
        for c, lab in enumerate(row):
            x = (c - 1.5) * H["pitch"]
            y = (0.5 - r) * H["pitch"]
            a.add(simple, name=f"cell_{lab}", loc=cq.Location(cq.Vector(x, y, -20.0)))
    # LaBr3 pair (2" x 2") back-to-back across the array, on a rail
    for s in (-1, 1):
        a.add(cq.Workplane("YZ").workplane(offset=s * 210).center(0, -20).circle(25.4).extrude(50.8 * s),
              name=f"LaBr3_{'E' if s > 0 else 'W'}", color=cq.Color(*C["labr"]))
    a.add(cq.Workplane("XY").workplane(offset=-cz / 2).rect(cx - 10, 20).extrude(10), name="LaBr3_rail",
          color=cq.Color(*C["steel"]))
    return a


def glb_to_json(stem):
    """Write <stem>.gltf.json: the GLB as self-contained glTF JSON (base64 buffer), servable as .json."""
    import base64
    import json
    import struct
    raw = open(os.path.join(OUT, stem + ".glb"), "rb").read()
    n_json = struct.unpack("<I", raw[12:16])[0]
    doc = json.loads(raw[20:20 + n_json])
    off = 20 + n_json
    n_bin = struct.unpack("<I", raw[off:off + 4])[0]
    bin_ = raw[off + 8:off + 8 + n_bin]
    doc["buffers"][0]["uri"] = "data:application/octet-stream;base64," + base64.b64encode(bin_).decode()
    with open(os.path.join(OUT, stem + ".gltf.json"), "w") as f:
        json.dump(doc, f, separators=(",", ":"))


def export(assy, stem, stl=False, tol=0.02):
    assy.export(os.path.join(OUT, stem + ".step"))
    assy.export(os.path.join(OUT, stem + ".glb"), tolerance=tol, angularTolerance=0.2)
    glb_to_json(stem)
    if stl:
        cq.exporters.export(assy.toCompound(), os.path.join(OUT, stem + ".stl"), tolerance=0.02, angularTolerance=0.2)


if __name__ == "__main__":
    cell = build_cell(detailed=True)
    export(cell, "dfm_cell", stl=True)
    export(build_house(), "dfm_house", tol=0.3)
    for f in sorted(os.listdir(OUT)):
        print(f, round(os.path.getsize(os.path.join(OUT, f)) / 1e6, 2), "MB")
