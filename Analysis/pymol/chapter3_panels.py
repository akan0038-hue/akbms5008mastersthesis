# =============================================================================
#  chapter3_panels.py   Chapter 3 structural panels
#
#  Renders each panel TWICE from an identical camera:
#     <tag>.png        clean geometry, no text at all      <- use this one
#     <tag>_guide.png  same view with labels               <- positioning aid
#
#  Workflow: place <tag>_guide.png as a bottom layer in Illustrator or Inkscape,
#  put <tag>.png on top, type the labels where the guide shows them, then delete
#  the guide layer. PyMOL has no label collision avoidance, so typesetting in
#  PyMOL will always pile text on top of itself.
#
#  A measurements.txt is written with the exact strings to type.
# =============================================================================

from pymol import cmd
import math, os, glob, re

ROOT = "/path/to/your/pymolstuff"
CIF  = os.path.join(ROOT, "cif")
OUT  = os.path.join(ROOT, "panels")
os.makedirs(OUT, exist_ok=True)

W, H, DPI = 2000, 1500, 300   # default; override per panel with size=(w,h)

# =============================================================================
#  PANELS
#    zoom      buffer in Angstrom around the measured atoms. Raise to pull back.
#    site_only True frames the active site and lets a distant nucleophile sit
#              near the edge; False frames the whole span.
# =============================================================================
PANELS = [
 dict(tag="fig3_2C", size=(2600, 1300), file="fig3_2C_chersinamycin_seed66799_sample0.cif",
      nuc="N130", carb="C2",   his=218, ala=91, zoom=6.0, site_only=False,
      expect=(22.92, 23.87, 108.7),
      caption="Synthetic construct on the chersinamycin TE domain"),

 dict(tag="fig3_3A", size=(2600, 1300), file="fig3_3A_synthpep1_chersinamycin_seed2039_sample2.cif",
      nuc="N1",   carb="C105", his=218, ala=91, zoom=6.0, site_only=False,
      expect=(26.04, 24.80, 142.1),
      caption="Full-length peptide on the chersinamycin TE domain"),

 dict(tag="fig3_3B", file="fig3_3B_synthpep6_chersinamycin_seed90416_sample4.cif",
      nuc="N1",   carb="C44",  his=218, ala=91, zoom=4.5, site_only=True,
      expect=(3.94, 4.85, 63.7),
      caption="Shortest peptide on the chersinamycin TE domain"),

 dict(tag="fig3_4B", file="fig3_4B_pris2_native_seed21853_sample1.cif",
      nuc="O3",   carb="C38",  his=224, ala=72, zoom=4.5, site_only=True,
      expect=(2.813, 3.889, 107.2),
      caption="Native pristinamycin IB substrate on Pristinamycin 2"),

 dict(tag="fig3_5A", file="fig3_5A_pris2_lser_seed3540_sample3.cif",
      nuc="O9",   carb="C20",  his=224, ala=72, zoom=4.5, site_only=True,
      expect=(2.651, 3.992, 76.1),
      caption="L-serine position-2 variant on Pristinamycin 2"),

 dict(tag="fig3_5B", size=(2400, 1400), file="fig3_5B_pris2_dap_seed23700_sample3.cif",
      nuc="N7",   carb="C20",  his=224, ala=72, zoom=5.5, site_only=False,
      expect=(10.85, 11.58, 76.7),
      caption="Dap position-2 variant on Pristinamycin 2"),
]

LIG = "chain L"

# =============================================================================
def sset(name, value, obj=None):
    """Set a PyMOL setting, ignoring names this build does not have."""
    try:
        cmd.set(name, value) if obj is None else cmd.set(name, value, obj)
    except Exception:
        print(f"   (skipped unknown setting: {name})")


def house_style():
    cmd.bg_color("white")
    sset("ray_opaque_background", 1)
    sset("ray_shadows", 0)
    sset("antialias", 2)
    sset("orthoscopic", 1)
    sset("depth_cue", 0)
    sset("specular", 0.15)
    sset("cartoon_transparency", 0.82)
    sset("stick_radius", 0.15)
    sset("sphere_scale", 0.32)
    sset("dash_gap", 0.35)
    sset("dash_length", 0.45)
    sset("dash_width", 4.0)
    sset("dash_radius", 0.05)
    sset("dash_round_ends", 0)
    sset("angle_size", 1.2)
    cmd.set_color("enzyme_c",   [0.88, 0.88, 0.90])
    cmd.set_color("ligand_c",   [0.13, 0.55, 0.55])
    cmd.set_color("triad_c",    [0.30, 0.30, 0.34])
    cmd.set_color("nuc_c",      [0.84, 0.11, 0.11])
    cmd.set_color("elec_c",     [0.10, 0.35, 0.80])
    cmd.set_color("met_dash",   [0.10, 0.35, 0.80])
    cmd.set_color("unmet_dash", [0.84, 0.11, 0.11])
    cmd.set_color("angle_c",    [0.55, 0.55, 0.58])

def label_style():
    sset("label_size", 26)
    sset("label_font_id", 7)
    sset("label_color", "black")
    sset("label_outline_color", "white")
    sset("label_distance_digits", 2)
    sset("label_angle_digits", 0)

def autocrop(path, margin=60):
    """Trim the white border so the molecule fills the frame. Framing a 26 A
    span and a 3 A span from one rule inevitably leaves whitespace; cropping
    after the fact is more reliable than tuning zoom per panel."""
    try:
        from PIL import Image, ImageChops
    except ImportError:
        return
    try:
        im = Image.open(path).convert("RGB")
        bg = Image.new("RGB", im.size, (255, 255, 255))
        diff = ImageChops.difference(im, bg).convert("L").point(lambda v: 255 if v > 6 else 0)
        box = diff.getbbox()
        if not box:
            return
        l, t, r, b = box
        l, t = max(0, l - margin), max(0, t - margin)
        r, b = min(im.width, r + margin), min(im.height, b + margin)
        im.crop((l, t, r, b)).save(path)
    except Exception as e:
        print(f"   (autocrop skipped: {e})")


def resolve(fname):
    stem = fname[:-4] if fname.lower().endswith(".cif") else fname
    for cand in (fname, fname + ".cif", stem, stem + ".cif.cif"):
        p = os.path.join(CIF, cand)
        if os.path.exists(p): return p
    hits = glob.glob(os.path.join(CIF, stem + "*"))
    if not hits:
        m = re.search(r"seed(\d+)_sample(\d+)", stem)
        if m:
            hits = [f for f in glob.glob(os.path.join(CIF, "**", "*.cif"), recursive=True)
                    if m.group(1) in f and ("sample-" + m.group(2)) in f]
    return hits[0] if hits else None

def frame_on_axis(a1, a2, a3, sel, buf):
    p1, p2, p3 = (cmd.get_atom_coords(a) for a in (a1, a2, a3))
    sub = lambda a, b: [a[i] - b[i] for i in range(3)]
    def norm(v):
        m = math.sqrt(sum(c*c for c in v)) or 1.0
        return [c/m for c in v]
    cross = lambda a, b: [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
    x = norm(sub(p2, p1)); z = norm(cross(x, norm(sub(p3, p1)))); y = cross(z, x)
    centre = [(p1[i]+p2[i])/2.0 for i in range(3)]
    v = cmd.get_view()
    cmd.set_view((x[0],y[0],z[0], x[1],y[1],z[1], x[2],y[2],z[2],
                  0.0,0.0,v[11]) + tuple(centre) + (v[15],v[16],v[17]))
    cmd.zoom(sel, buf)

def find_carbonyl_o(carb_sel, ala_cb_sel):
    cmd.select("_oxy", f"({LIG}) and elem O within 1.6 of ({carb_sel})")
    names = []
    cmd.iterate("_oxy", "names.append(name)", space={"names": names})
    if not names: return None
    if len(names) == 1: return names[0]
    return sorted((cmd.get_distance(f"{LIG} and name {n}", ala_cb_sel), n) for n in names)[-1][1]

# =============================================================================
report = []

def draw_panel(cfg):
    cmd.reinitialize(); house_style()
    path = resolve(cfg["file"])
    if path is None:
        print(f"{cfg['tag']:9s} SKIPPED, no file matching {cfg['file']}"); return
    cmd.load(path, "m")

    nuc  = f"m and {LIG} and name {cfg['nuc']}"
    carb = f"m and {LIG} and name {cfg['carb']}"
    ne2  = f"m and chain A and resi {cfg['his']} and name NE2"
    hisr = f"m and chain A and resi {cfg['his']}"
    alar = f"m and chain A and resi {cfg['ala']}"
    acb  = f"{alar} and name CB"
    oname = find_carbonyl_o(carb, acb)
    ocarb = f"m and {LIG} and name {oname}" if oname else None

    cmd.hide("everything"); cmd.remove("solvent or inorganic")
    cmd.show("cartoon", "m and polymer"); cmd.color("enzyme_c", "m and polymer")
    cmd.show("sticks", f"m and {LIG}"); cmd.color("ligand_c", f"m and {LIG} and elem C")
    cmd.show("sticks", f"({hisr} or {alar}) and sidechain")
    cmd.color("triad_c", f"({hisr} or {alar}) and elem C")
    cmd.util.cnc(f"m and ({LIG} or {hisr} or {alar})")
    cmd.show("spheres", f"{nuc} or {carb}")
    cmd.color("nuc_c", nuc); cmd.color("elec_c", carb)

    d_his  = cmd.get_distance(nuc, ne2)
    d_carb = cmd.get_distance(nuc, carb)
    cmd.distance("d_his",  nuc, ne2)
    cmd.distance("d_carb", nuc, carb)
    # Colour ON the object rather than with cmd.color. The angle object draws a
    # leg along exactly the same axis as d_carb and was winning the depth test,
    # which is why every nucleophile-to-carbon line came out grey.
    sset("dash_color", "met_dash" if d_his  < 3.5 else "unmet_dash", "d_his")
    sset("dash_color", "met_dash" if d_carb < 4.0 else "unmet_dash", "d_carb")
    sset("dash_radius", 0.055, "d_his")
    sset("dash_radius", 0.055, "d_carb")
    ang = None
    if ocarb:
        ang = cmd.get_angle(nuc, carb, ocarb)
        cmd.angle("bd", nuc, carb, ocarb)
        sset("angle_color", "angle_c", "bd")
        sset("angle_size", 1.5, "bd")       # large enough to read, small enough to stay attached
        sset("dash_radius", 0.018, "bd")    # legs recede to construction lines
        sset("dash_gap", 0.20, "bd")
        cmd.show("sticks", ocarb)

    sel = f"{nuc} or {carb} or {ne2} or {hisr} or {alar}" + (f" or {ocarb}" if ocarb else "")
    if not cfg.get("site_only", True):
        sel = f"{sel} or (m and {LIG})"
    cmd.select("measured", sel)
    frame_on_axis(nuc, ne2, carb, "measured", cfg.get("zoom", 5.0))
    cmd.deselect()

    # ---- clean render, no text ----
    cmd.hide("labels")
    sset("label_size", 0)
    w, h = cfg.get("size", (W, H))
    clean = os.path.join(OUT, cfg["tag"] + ".png")
    cmd.png(clean, width=w, height=h, dpi=DPI, ray=1)

    # ---- guide render, identical camera ----
    view = cmd.get_view()
    label_style()
    cmd.label(ne2, f'"His{cfg["his"]}"')
    cmd.label(acb, f'"Ala{cfg["ala"]}"')
    cmd.label(nuc,  '"NUC"')
    cmd.label(carb, '"ELEC"')
    if ocarb: cmd.label(ocarb, '"C=O"')
    # push each label a different way so the guide is actually usable
    sset("label_position", (0.0, -2.2, 0.0), nuc)
    sset("label_position", (2.4,  1.6, 0.0), carb)
    sset("label_position", (-2.4, 2.4, 0.0), ocarb) if ocarb else None
    sset("label_position", (2.6,  0.0, 0.0), ne2)
    sset("label_position", (2.4, -1.6, 0.0), acb)
    cmd.set_view(view)
    guide = os.path.join(OUT, cfg["tag"] + "_guide.png")
    cmd.png(guide, width=w, height=h, dpi=DPI, ray=1)
    # crop both to the SAME box so the guide still overlays the clean panel
    try:
        from PIL import Image, ImageChops
        a, b_ = Image.open(clean).convert("RGB"), Image.open(guide).convert("RGB")
        bg = Image.new("RGB", a.size, (255, 255, 255))
        d1 = ImageChops.difference(a,  bg).convert("L").point(lambda v: 255 if v > 6 else 0)
        d2 = ImageChops.difference(b_, bg).convert("L").point(lambda v: 255 if v > 6 else 0)
        box = ImageChops.lighter(d1, d2).getbbox()
        if box:
            m = 70
            l, t, r, bt = box
            box = (max(0, l-m), max(0, t-m), min(a.width, r+m), min(a.height, bt+m))
            a.crop(box).save(clean)
            b_.crop(box).save(guide)
    except ImportError:
        print("   (Pillow not available, panels left uncropped)")
    except Exception as e:
        print(f"   (crop skipped: {e})")

    e = cfg["expect"]
    ok = abs(d_his-e[0]) < 0.05 and abs(d_carb-e[1]) < 0.05
    a  = f"{ang:5.1f}" if ang is not None else "  n/a"
    print(f"{cfg['tag']:9s} NE2 {d_his:6.2f} (csv {e[0]:6.2f})   C {d_carb:6.2f} "
          f"(csv {e[1]:6.2f})   angle {a} (csv {e[2]:5.1f}){'' if ok else '   <-- CHECK'}")

    report.append((cfg, d_his, d_carb, ang, oname))

# =============================================================================
if not os.path.isdir(CIF):
    print(f"\nFolder does not exist: {CIF}\n")
else:
    print(f"\nFiles in {CIF}:")
    for f in sorted(os.listdir(CIF)): print("   ", f)

print("\npanel     measured vs CSV")
print("-" * 88)
for cfg in PANELS: draw_panel(cfg)
print("-" * 88)

# ---- text to type onto each panel ----
lines = ["LABELS TO ADD IN ILLUSTRATOR / INKSCAPE", "=" * 60, ""]
for cfg, dh, dc, ang, oname in report:
    met = lambda v, lim: "met" if v < lim else "not met"
    lines += [
        f"{cfg['tag']}   {cfg['caption']}",
        f"    red sphere    : ring-closing nucleophile ({cfg['nuc']})",
        f"    blue sphere   : electrophilic carbon ({cfg['carb']})",
        f"    grey sticks   : His{cfg['his']} and Ala{cfg['ala']}"
        f" (the catalytic serine, modelled as alanine)",
        f"    blue/red dash to His NE2 : {dh:.2f} A   [{met(dh,3.5)}, threshold 3.5]",
        f"    blue/red dash to carbon  : {dc:.2f} A   [{met(dc,4.0)}, threshold 4.0]",
    ]
    if ang is not None:
        v = "met" if 95 <= ang <= 115 else "not met"
        lines.append(f"    grey arc at the carbon   : {ang:.0f} deg  [{v}, window 95-115]"
                     f"  (measured to {oname})")
    lines.append("")
open(os.path.join(OUT, "measurements.txt"), "w").write("\n".join(lines))

print(f"\nClean panels : {OUT}\\<tag>.png")
print(f"Guides       : {OUT}\\<tag>_guide.png")
print(f"Label text   : {OUT}\\measurements.txt\n")
